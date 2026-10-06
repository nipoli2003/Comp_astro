#!/bin/bash
set -e

mkdir -p build
cd build
cmake ..
make -j4
cd ..

# Run from Ex1 root, NOT from inside build
./build/galaxy_grid

# Run the plotting script
python3 plot/plot_results.py
python3 plot/plot_sanity_check.py