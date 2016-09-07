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

cmake -DCOSMO_N=50 -DCOSMO_NY=4 -DCOSMO_NZ=4 .. && make -j16
./cosmo ../config/AwA/stability_test_r1.txt

cmake -DCOSMO_N=100 -DCOSMO_NY=4 -DCOSMO_NZ=4 .. && make -j16
./cosmo ../config/AwA/stability_test_r2.txt

cmake -DCOSMO_N=200 -DCOSMO_NY=4 -DCOSMO_NZ=4 .. && make -j16
./cosmo ../config/AwA/stability_test_r4.txt
