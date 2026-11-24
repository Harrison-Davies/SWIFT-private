/*******************************************************************************
 * This file is part of SWIFT.
 *
 * DEBUG ONLY (Harrison): dump every neighbour of a particle that has asked for
 * a time-step below dt_min, together with everything that feeds the REMIX
 * geometric moment m2_bar, so that a near-singular m2 can be diagnosed.
 *
 * Called from get_part_timestep() in timestep.h just before the dt_min error.
 *
 * Output:
 *   - a summary via message()
 *   - one message() line per neighbour (core columns)
 *   - a whitespace-separated table with all columns, written to
 *       dtmin_neighbours_id<id>_step<step>_rank<nodeID>.txt
 *     Read it with pandas.read_csv(fname, sep=r"\s+", comment="#"): comment
 *     and summary lines start with '#', the column names are the first
 *     non-comment line.
 *
 * Caveats:
 *   - Brute-force search over all local (and, with MPI, foreign) particles.
 *     O(N) but only runs once, right before the code aborts.
 *   - Other threads may still be writing to particles while this runs; values
 *     are read without locks (fine for a post-mortem).
 *   - m2 is rebuilt here from the current positions / h / rho_evol. These can
 *     differ slightly from what the gradient loop saw this step (e.g. rho_evol
 *     has since been kicked), so expect small, not large, differences from the
 *     stored debug_m2_* values.
 *   - The gradient loop includes a pair if r < gamma*h_i OR r < gamma*h_j for
 *     symmetric interactions, but only r < gamma*h_i for non-symmetric ones
 *     (inactive or foreign neighbours). Both cases are included here and are
 *     flagged with the in_hi / in_hj columns.
 *   - Force-struct values of neighbours (A, vac_switch, debug_m2_q) are from
 *     their last prepare_force, and may be meaningless for foreign particles.
 ******************************************************************************/
#ifndef SWIFT_REMIX_HYDRO_DEBUG_NEIGHBOURS_H
#define SWIFT_REMIX_HYDRO_DEBUG_NEIGHBOURS_H

#include <float.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "active.h"
#include "engine.h"
#include "kernel_hydro.h"
#include "periodic.h"
#include "space.h"
#include "symmetric_matrix.h"

/* Running totals while scanning neighbours. */
struct debug_ngb_totals {
  int n_in_hi;      /* r < gamma h_i */
  int n_in_either;  /* r < gamma max(h_i, h_j) */
  int n_solid, n_fluid, n_other_phase;
  int n_other_mat;
  int n_inactive;
  int n_h_max;
  int n_foreign;
  int n_coincident; /* r/h_i < 1e-3 */
  float min_r_over_hi;
  float max_hj_over_hi, min_hj_over_hi;
  float max_vol_over_voli, min_vol_over_voli;
  double m0_bar;
  double m1_bar[3];
  double m2_bar[6]; /* xx yy zz xy xz yz, matching struct sym_matrix */
};

__attribute__((always_inline)) INLINE static const char *debug_ngb_phase_name(
    const int phase) {
  if (phase == mat_phase_solid) return "solid";
  if (phase == mat_phase_fluid) return "fluid";
  return "other";
}

/**
 * @brief Write the column header for the neighbour table.
 */
INLINE static void debug_ngb_write_header(FILE *f) {
  fprintf(f,
          "id foreign self x y z dx_hi dy_hi dz_hi r_hi r_hj in_hi in_hj "
          "h h_ratio mass rho rho_evol vol vol_ratio mat_id same_mat phase "
          "phase_full is_h_max m0 damage u vx vy vz time_bin active "
          "w_term_vol m2_trace_share A vac_switch ngb_debug_m2_q\n");
}

/**
 * @brief Process one candidate neighbour pj of pi: decide whether it is in
 * range, accumulate its m2 contribution and print / write its row.
 *
 * @param xpj The xpart of pj, or NULL if unavailable (foreign particles).
 */
INLINE static void debug_ngb_process(const struct part *pi,
                                     const struct part *pj,
                                     const struct xpart *xpj,
                                     const struct engine *e,
                                     const int is_foreign, FILE *f,
                                     struct debug_ngb_totals *t) {

  /* Skip particles that do not exist (any more). */
  if (pj->time_bin >= time_bin_inhibited) return;

  const struct space *s = e->s;
  const int is_self = (pj == pi);

  /* Separation, pi - pj, as in the interaction functions. */
  double dx[3];
  for (int k = 0; k < 3; k++) {
    dx[k] = pi->x[k] - pj->x[k];
    if (s->periodic) dx[k] = nearest(dx[k], s->dim[k]);
  }
  const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float r = (float)sqrt(r2);

  const float hi = pi->h;
  const float hj = pj->h;
  const int in_hi = (r < kernel_gamma * hi);
  const int in_hj = (r < kernel_gamma * hj);
  if (!is_self && !in_hi && !in_hj) return;

  /* Kernel weights exactly as in hydro_runner_iact_gradient_extra_kernel. */
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;
  float wi = 0.f, wj = 0.f;
  if (in_hi) kernel_eval(r * hi_inv, &wi);
  if (in_hj) kernel_eval(r * hj_inv, &wj);
  const float w_term =
      0.5f * (wi * pow_dimension(hi_inv) + wj * pow_dimension(hj_inv));

  const float vol_i = pi->mass / pi->rho_evol;
  const float vol_j = pj->mass / pj->rho_evol;
  const float w_term_vol = w_term * vol_j;

  /* Contribution to m2_bar / h_i^2 (trace), for ranking which neighbours
   * dominate the moment matrix. */
  const float m2_trace_contrib = (float)(r2 * w_term_vol) * hi_inv * hi_inv;

  if (!is_self) {
    /* Accumulate moments (the self term only enters m0_bar, added below). */
    t->m0_bar += w_term_vol;
    for (int k = 0; k < 3; k++) t->m1_bar[k] += dx[k] * w_term_vol;
    t->m2_bar[0] += dx[0] * dx[0] * w_term_vol;
    t->m2_bar[1] += dx[1] * dx[1] * w_term_vol;
    t->m2_bar[2] += dx[2] * dx[2] * w_term_vol;
    t->m2_bar[3] += dx[0] * dx[1] * w_term_vol;
    t->m2_bar[4] += dx[0] * dx[2] * w_term_vol;
    t->m2_bar[5] += dx[1] * dx[2] * w_term_vol;

    /* Counters. */
    if (in_hi) t->n_in_hi++;
    t->n_in_either++;
    if (pj->phase == mat_phase_solid)
      t->n_solid++;
    else if (pj->phase == mat_phase_fluid)
      t->n_fluid++;
    else
      t->n_other_phase++;
    if (pj->mat_id != pi->mat_id) t->n_other_mat++;
    if (!part_is_active(pj, e)) t->n_inactive++;
    if (pj->is_h_max) t->n_h_max++;
    if (is_foreign) t->n_foreign++;
    if (r * hi_inv < 1e-3f) t->n_coincident++;
    t->min_r_over_hi = fminf(t->min_r_over_hi, r * hi_inv);
    t->max_hj_over_hi = fmaxf(t->max_hj_over_hi, hj * hi_inv);
    t->min_hj_over_hi = fminf(t->min_hj_over_hi, hj * hi_inv);
    t->max_vol_over_voli = fmaxf(t->max_vol_over_voli, vol_j / vol_i);
    t->min_vol_over_voli = fminf(t->min_vol_over_voli, vol_j / vol_i);
  }

  const char *phase_full =
      (xpj != NULL) ? debug_ngb_phase_name(xpj->phase_full) : "NA";

  /* Full row to the file. */
  if (f != NULL) {
    fprintf(f,
            "%lld %d %d %.9e %.9e %.9e %.6e %.6e %.6e %.6e %.6e %d %d "
            "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %d %d %s "
            "%s %d %.6e %.6e %.6e %.6e %.6e %.6e %d %d "
            "%.6e %.6e %.6e %.6e %.6e\n",
            pj->id, is_foreign, is_self, pj->x[0], pj->x[1], pj->x[2],
            dx[0] * hi_inv, dx[1] * hi_inv, dx[2] * hi_inv, r * hi_inv,
            r * hj_inv, in_hi, in_hj, hj, hj * hi_inv, pj->mass, pj->rho,
            pj->rho_evol, vol_j, vol_j / vol_i, (int)pj->mat_id,
            pj->mat_id == pi->mat_id, debug_ngb_phase_name(pj->phase),
            phase_full,
            (int)pj->is_h_max, pj->m0, strength_get_damage(pj), pj->u,
            pj->v[0], pj->v[1], pj->v[2], (int)pj->time_bin,
            part_is_active(pj, e), w_term_vol, m2_trace_contrib, pj->force.A,
            pj->force.vac_switch, pj->force.debug_m2_q);
  }

  /* Compact row to stdout. */
  message(
      "ngb of %lld: id=%lld%s%s dx/h=(%+.3f,%+.3f,%+.3f) r/hi=%.3f "
      "r/hj=%.3f hj/hi=%.3f mat=%d %s%s vol/vol_i=%.3f m0=%.3f "
      "tbin=%d%s m2tr=%.3e",
      pi->id, pj->id, is_self ? " [SELF]" : "", is_foreign ? " [foreign]" : "",
      dx[0] * hi_inv, dx[1] * hi_inv, dx[2] * hi_inv, r * hi_inv, r * hj_inv,
      hj * hi_inv, (int)pj->mat_id, debug_ngb_phase_name(pj->phase),
      pj->is_h_max ? " h_max" : "", vol_j / vol_i, pj->m0, (int)pj->time_bin,
      part_is_active(pj, e) ? "" : " inactive", m2_trace_contrib);
}

/**
 * @brief Dump every neighbour of p and a reconstruction of its m2_bar.
 */
INLINE static void debug_dump_neighbours_remix(const struct part *p,
                                               const struct engine *e) {

  const struct space *s = e->s;

  /* Open the output file (failure is not fatal: we still print). */
  char fname[256];
  snprintf(fname, sizeof(fname), "dtmin_neighbours_id%lld_step%d_rank%d.txt",
           p->id, e->step, e->nodeID);
  FILE *f = fopen(fname, "w");
  if (f == NULL) message("Could not open %s, printing only", fname);

  if (f != NULL) {
    fprintf(f, "# Neighbours of part id=%lld at step=%d (rank %d)\n", p->id,
            e->step, e->nodeID);
    fprintf(f,
            "# Columns: dx_* = (x_i - x_j)/h_i (periodic-wrapped); in_hi/in_hj "
            "= r < gamma*h_i / gamma*h_j; phase_full (kick-time phase) is NA "
            "for foreign particles; w_term_vol = "
            "0.5(W_i+W_j)*m_j/rho_evol_j as used for m2_bar; m2_trace_share = "
            "r^2 * w_term_vol / h_i^2. The first data row (self=1) is the "
            "particle itself.\n");
    debug_ngb_write_header(f);
  }

  struct debug_ngb_totals t;
  memset(&t, 0, sizeof(t));
  t.min_r_over_hi = FLT_MAX;
  t.min_hj_over_hi = FLT_MAX;
  t.min_vol_over_voli = FLT_MAX;

  /* The particle itself first (pointer equality marks it as self). Find its
   * xpart via its index in the local array. */
  const ptrdiff_t ind = p - s->parts;
  const struct xpart *xp =
      (ind >= 0 && (size_t)ind < s->nr_parts) ? &s->xparts[ind] : NULL;
  debug_ngb_process(p, p, xp, e, /*is_foreign=*/0, f, &t);

  /* All local particles. */
  for (size_t k = 0; k < s->nr_parts; k++) {
    const struct part *pj = &s->parts[k];
    if (pj == p) continue;
    debug_ngb_process(p, pj, &s->xparts[k], e, /*is_foreign=*/0, f, &t);
  }

#ifdef WITH_MPI
  /* Foreign particles currently held on this rank. */
  for (size_t k = 0; k < s->nr_parts_foreign; k++) {
    const struct part *pj = &s->parts_foreign[k];
    debug_ngb_process(p, pj, NULL, e, /*is_foreign=*/1, f, &t);
  }
#endif

  /* Self term of m0_bar (hydro_end_gradient_extra_kernel). */
  const float hi = p->h;
  const float hi_inv = 1.f / hi;
  float w_root;
  kernel_eval(0.f, &w_root);
  t.m0_bar += (p->mass / p->rho_evol) * w_root * pow_dimension(hi_inv);

  /* Rebuilt m2_bar / h^2 and its shape. */
  struct sym_matrix m2;
  for (int k = 0; k < 6; k++)
    m2.elements[k] = (float)(t.m2_bar[k] * hi_inv * hi_inv);
  const float tr = m2.xx + m2.yy + m2.zz;
  const float det = m2.xx * (m2.yy * m2.zz - m2.yz * m2.yz) -
                    m2.xy * (m2.xy * m2.zz - m2.yz * m2.xz) +
                    m2.xz * (m2.xy * m2.yz - m2.yy * m2.xz);
  const float q = (tr > 0.f) ? 27.f * det / (tr * tr * tr) : 0.f;
  float eig[3];
  sym_matrix_compute_eigenvalues(eig, m2);
  float eig_min = fminf(eig[0], fminf(eig[1], eig[2]));
  float eig_max = fmaxf(eig[0], fmaxf(eig[1], eig[2]));

  /* One-sidedness of the neighbour cloud: |m1_bar| / (m0_bar h). ~0 in the
   * bulk, O(0.1-1) at a free surface or a hole. */
  const double m1_norm = sqrt(t.m1_bar[0] * t.m1_bar[0] +
                              t.m1_bar[1] * t.m1_bar[1] +
                              t.m1_bar[2] * t.m1_bar[2]);
  const float one_sided =
      (t.m0_bar > 0.) ? (float)(m1_norm / (t.m0_bar * hi)) : 0.f;

  char summary[8][512];
  snprintf(summary[0], 512,
           "part %lld: %d neighbours within gamma*max(hi,hj) (%d within "
           "gamma*hi), %d foreign, %d inactive",
           p->id, t.n_in_either, t.n_in_hi, t.n_foreign, t.n_inactive);
  snprintf(summary[1], 512,
           "part %lld: phases: %d solid, %d fluid, %d other; %d of other "
           "material; %d at h_max; %d coincident (r/hi<1e-3)",
           p->id, t.n_solid, t.n_fluid, t.n_other_phase, t.n_other_mat,
           t.n_h_max, t.n_coincident);
  snprintf(summary[2], 512,
           "part %lld: min r/hi=%.4f  hj/hi in [%.3f, %.3f]  vol_j/vol_i in "
           "[%.3e, %.3e]",
           p->id, t.min_r_over_hi, t.min_hj_over_hi, t.max_hj_over_hi,
           t.min_vol_over_voli, t.max_vol_over_voli);
  snprintf(summary[3], 512,
           "part %lld: rebuilt m2/h^2 = [[%.4e %.4e %.4e] [. %.4e %.4e] "
           "[. . %.4e]]",
           p->id, m2.xx, m2.xy, m2.xz, m2.yy, m2.yz, m2.zz);
  snprintf(summary[4], 512,
           "part %lld: rebuilt m2/h^2 eigenvalues = (%.4e, %.4e, %.4e) "
           "min/max=%.4e  q=%.4e tr=%.4e det=%.4e",
           p->id, eig[0], eig[1], eig[2],
           eig_max != 0.f ? eig_min / eig_max : 0.f, q, tr, det);
  snprintf(summary[5], 512,
           "part %lld: stored at prepare_force: q=%.4e tr=%.4e det=%.4e "
           "(compare with rebuilt values above)",
           p->id, p->force.debug_m2_q, p->force.debug_m2_tr,
           p->force.debug_m2_det);
  snprintf(summary[6], 512,
           "part %lld: m0_bar=%.4e (~1 when fully surrounded)  "
           "|m1_bar|/(m0_bar h)=%.4e (one-sidedness, ~0 in the bulk)",
           p->id, t.m0_bar, one_sided);
  snprintf(summary[7], 512, "part %lld: full neighbour table written to %s",
           p->id, f != NULL ? fname : "(file could not be opened)");

  for (int k = 0; k < 8; k++) {
    message("%s", summary[k]);
    if (f != NULL) fprintf(f, "# SUMMARY %s\n", summary[k]);
  }

  if (f != NULL) fclose(f);
}

#endif /* SWIFT_REMIX_HYDRO_DEBUG_NEIGHBOURS_H */
