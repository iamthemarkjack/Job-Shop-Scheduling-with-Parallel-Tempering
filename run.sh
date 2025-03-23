#!/bin/bash
mkdir -p build
cd build
cmake ..
make
cd ..
./jssp_solver --processing-times data/instances/processing_times.txt --machines data/instances/machines.txt --iterations 10 --rho -10 --min-temp 0.01