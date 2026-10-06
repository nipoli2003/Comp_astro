mkdir -p build
cd build
cmake ..
make -j4
./galaxy_grid
cd ..
python3 plot_results.py