#!/bin/bash

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# And up a directory should be the main codebase.
pwd

###
# Compile
###
mkdir -p ../build
cd ../build

cmake -DCOSMO_N=50 -DCOSMO_NY=6 -DCOSMO_NZ=6 -DCOSMO_H_LEN_FRAC=1 .. && make -j16
./cosmo ../config/AwA/gauge_wave_test_r1.txt

cmake -DCOSMO_N=100 -DCOSMO_NY=6 -DCOSMO_NZ=6 -DCOSMO_H_LEN_FRAC=1 .. && make -j16
./cosmo ../config/AwA/gauge_wave_test_r2.txt

cmake -DCOSMO_N=200 -DCOSMO_NY=6 -DCOSMO_NZ=6 -DCOSMO_H_LEN_FRAC=1 .. && make -j16
./cosmo ../config/AwA/gauge_wave_test_r4.txt
