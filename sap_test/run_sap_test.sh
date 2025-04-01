#!/bin/bash

# compile and execute adaptive_test.cpp
g++ -o  adaptive_test adaptive_test.cpp
./adaptive_test

# compile and execute unadaptive_test.cpp
g++ -o unadaptive_test unadaptive_test.cpp
./unadaptive_test

# run the python script to generate plots
python3 plot_adap.py
python3 plot_proximity.py