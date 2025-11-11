#!/bin/bash
PAGER=cat

make clean
# ./configure \
#   --with-hydro=remix \
#   --with-equation-of-state=planetary \
#   --with-kernel=wendland-C2 \
#   --enable-compiler-warnings=yes \
#   --with-gravity=basic \
#   --with-parmetis=/usr LDFLAGS="-L/usr/lib -lmetis" \
#   --enable-debug=yes \
#   --enable-optimization=no \
#   --enable-sanitizer=yes \
#   --enable-undefined-sanitizer=yes \
#   --enable-task-debugging=yes \
#   --enable-threadpool-debugging=yes \
#   --enable-debugging-checks=yes \
#   --enable-my-def\

./configure \
  --with-hydro=remix \
  --with-equation-of-state=planetary \
  --with-kernel=wendland-C2 \
  --enable-compiler-warnings=yes \
  --with-gravity=basic \
  --with-parmetis=/usr LDFLAGS="-L/usr/lib -lmetis" \
  --enable-debug=yes \
  --enable-optimization=no \
  --enable-my-def \

cd /mnt/c/SWIFT_Private || exit 1
make -j4 || { echo "Make failed"; exit 1; }
cd /mnt/c/SWIFT_Private/examples/Planetary/DemoImpact || exit 1

export ASAN_OPTIONS="halt_on_error=1:abort_on_error=1:detect_leaks=0:print_stacktrace=1:intercept_tls=0"
gdb ../../../swift \
    -ex "set args --hydro --self-gravity --threads=1 --steps=10 demo_impact_n50.yml" \
    -ex start