import h5py
import numpy as np


def seed_flaws(num_part, volume, m_flaw, k_flaw, hardcoded_max_flaws, seed, chunk_size):
    """Assign Weibull-distributed activation thresholds to particles, adding flaws until each particle has at least one flaw."""
    rng = np.random.default_rng(seed)
    activation_thresholds = np.zeros((num_part, hardcoded_max_flaws))

    # Use chunks to avoid memory issues because of storing (num_draws, num_part) arrays
    chunks = []
    seen = np.zeros(num_part, dtype=bool)
    num_seen = 0

    # Loop over chunks until all particles have at least one flaw
    while num_seen < num_part:
        # Draw chunk worth of particle IDs
        chunk = rng.integers(0, num_part, size=chunk_size)

        # Find all the unique particle IDs that haven't been drawn before
        unique_in_chunk = np.unique(chunk)
        newly_seen = (~seen[unique_in_chunk]).sum()

        # If all particles have been drawn, finish up
        if num_seen + newly_seen >= num_part:
            # Mark which draws in this chunk are the first occurrence of each particle
            _, first_occurrences = np.unique(chunk, return_index=True)
            new_in_chunk = np.zeros(len(chunk), dtype=bool)
            new_in_chunk[first_occurrences] = True
            new_in_chunk &= ~seen[
                chunk
            ]  # exclude particles already seen before this chunk

            # Find the cutoff point where the final missing ID gets drawn for the first time
            cumulative_new = np.cumsum(new_in_chunk)
            cutoff = np.argmax(cumulative_new == (num_part - num_seen))

            # Add this chunk, up to the cutoff, to the list keeping track of IDs
            chunks.append(chunk[: cutoff + 1])
            break

        # Update quantities keeping track of the IDs that we've seen and move onto next chunk
        seen[unique_in_chunk] = True
        num_seen += newly_seen
        chunks.append(chunk)

    # Array of drawn indices
    sample_indices = np.concatenate(chunks)

    # Threshold for the nth flaw
    flaw_numbers = np.arange(1, len(sample_indices) + 1)
    thresholds = (flaw_numbers / (k_flaw * volume)) ** (1.0 / m_flaw)

    # For each draw, how many flaws did that particle already have
    order = np.argsort(sample_indices, kind="stable")
    prior_flaw_counts = np.empty(len(sample_indices), dtype=np.int32)
    prior_flaw_counts[order] = np.concatenate(
        [np.arange(count) for count in np.bincount(sample_indices, minlength=num_part)]
    )

    # Fill in the thresholds and number of flaws for each particle in the right orther
    activation_thresholds[sample_indices, prior_flaw_counts] = thresholds
    N_flaws = np.bincount(sample_indices, minlength=num_part)

    return N_flaws, activation_thresholds

def add_flaws_to_file(hdf5_file: str, nflaws, thresholds):
    with h5py.File(hdf5_file, "a") as f:
        grp: h5py.Group = f["PartType0"]  # pyright: ignore[reportAssignmentType]
        for name in ["NumFlaws", "ActivationThresholds"]:
            if name in grp:
                raise ValueError(f"{name} already exists in file")
        grp.create_dataset("NumFlaws", data=nflaws.reshape(-1, 1), dtype="i")
        grp.create_dataset("ActivationThresholds", data=thresholds, dtype="f")

if __name__ == "__main__":
    from tools.py_tools import unit_converter
    # Constants
    R_earth = 6.371e6   # m
    M_earth = 5.9724e24  # kg
    G = 6.67408e-11  # m^3 kg^-1 s^-2
    uc_earth = unit_converter(U_M=M_earth, U_L=R_earth)

    # Flaw parameters from B&A 1994
    m_flaw = 8
    k_flaw = 5e28 * 100**3  # cm^-3 => m^-3
    seed = 12345
    chunk_size = 100000  # to avoid memory issues with large arrays
    hardcoded_max_flaws = 40

    hdf5_file = "./demo_impact_n50.hdf5"
    with h5py.File(hdf5_file, "r") as f:
        num_part = f["Header"].attrs["NumPart_Total"][0] # pyright: ignore[reportIndexIssue]
        box_size = f["Header"].attrs["BoxSize"]
        masses = np.array(f["PartType0/Masses"]) * uc_earth.M
        densities = np.array(f["PartType0/Densities"]) * uc_earth.rho

    flaw_volume = np.mean(masses / densities)
    nflaws, thresholds = seed_flaws(
        num_part,
        flaw_volume,
        m_flaw,
        k_flaw,
        hardcoded_max_flaws,
        seed,
        chunk_size,
    )

    add_flaws_to_file(hdf5_file, nflaws, thresholds)