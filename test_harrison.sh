#!/bin/bash

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
  --enable-my-def \

cd /mnt/c/SWIFT_Private || exit 1
make -j4 || { echo "Make failed"; exit 1; }
cd /mnt/c/SWIFT_Private/examples/Planetary/DemoImpact || exit 1

../../../swift --hydro --self-gravity --threads=4 demo_impact_"n50".yml