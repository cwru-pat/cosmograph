#!/bin/bash

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# And up 2 directories should be the main codebase.
pwd

mkdir -p ../../build
cd ../../build

# Temporary file for many misc. tests
# Fewer steps here
TMP_CONFIG_FILE=gauge_wave_test_r4.text.tmp.test
cp ../config/AwA/gauge_wave_test_r4.txt $TMP_CONFIG_FILE

###
# Compile
###

# Use "K"-damping
cmake -DCOSMO_H_LEN_FRAC=1 -DCOSMO_N=200 -DCOSMO_NY=6 -DCOSMO_NZ=6 .. && make -j16
echo "Running with K Damping..."
sed -i -E "s/k_damping_amp = [0-9]+/k_damping_amp = 1/g" $TMP_CONFIG_FILE
./cosmo $TMP_CONFIG_FILE
sed -i -E "s/k_damping_amp = [0-9]+/k_damping_amp = 0/g" $TMP_CONFIG_FILE
# Use "A"-adjustment
echo "Running with A-Adjustment..."
sed -i -E "s/a_adj_amp = [0-9]+/a_adj_amp = 1/g" $TMP_CONFIG_FILE
./cosmo $TMP_CONFIG_FILE
sed -i -E "s/a_adj_amp = [0-9]+/a_adj_amp = 0/g" $TMP_CONFIG_FILE


# Gauge wave in y-direction (lowest res only)
cmake -DCOSMO_N=200 -DCOSMO_NX=6 -DCOSMO_NZ=6 -DCOSMO_H_LEN_FRAC=1 .. && make -j16
sed -i -E "s/gauge_wave_dir = [0-9]+/gauge_wave_dir = 2/g" $TMP_CONFIG_FILE
./cosmo $TMP_CONFIG_FILE

# Gauge wave in z-direction (lowest res only)
cmake -DCOSMO_N=200 -DCOSMO_NX=6 -DCOSMO_NY=6 -DCOSMO_H_LEN_FRAC=1 .. && make -j16
sed -i -E "s/gauge_wave_dir = [0-9]+/gauge_wave_dir = 3/g" $TMP_CONFIG_FILE
./cosmo $TMP_CONFIG_FILE

# Gauge wave in x-direction (vary resolutions)
cmake -DCOSMO_N=50 -DCOSMO_NY=6 -DCOSMO_NZ=6 -DCOSMO_H_LEN_FRAC=1 .. && make -j16


./cosmo ../config/AwA/gauge_wave_test_r1.txt
cmake -DCOSMO_N=100 -DCOSMO_NY=6 -DCOSMO_NZ=6 -DCOSMO_H_LEN_FRAC=1 .. && make -j16
./cosmo ../config/AwA/gauge_wave_test_r2.txt
cmake -DCOSMO_N=200 -DCOSMO_NY=6 -DCOSMO_NZ=6 -DCOSMO_H_LEN_FRAC=1 .. && make -j16
./cosmo ../config/AwA/gauge_wave_test_r4.txt

rm $TMP_CONFIG_FILE
