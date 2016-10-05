#!/bin/bash

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# And up two directories should be the main codebase.
pwd

###
# Compile
###
mkdir -p ../../build
cd ../../build
TMP_CONFIG_FILE=stability_test_r1.txt.test
cp ../config/AwA/stability_test_r1.txt $TMP_CONFIG_FILE

# KO Dissipation run
echo "Running with KO Dissipation..."
cmake -DCOSMO_N=50 -DCOSMO_NY=6 -DCOSMO_NZ=6 .. && make -j16
sed -i -E "s/KO_damping_coefficient = [0-9]+/KO_damping_coefficient = 1/g" $TMP_CONFIG_FILE
./cosmo $TMP_CONFIG_FILE
sed -i -E "s/KO_damping_coefficient = [0-9]+/KO_damping_coefficient = 0/g" $TMP_CONFIG_FILE
# Use "K"-damping
echo "Running with K Damping..."
cmake -DCOSMO_N=50 -DCOSMO_NY=4 -DCOSMO_NZ=4 .. && make -j16
sed -i -E "s/k_damping_amp = [0-9]+/k_damping_amp = 1/g" $TMP_CONFIG_FILE
./cosmo $TMP_CONFIG_FILE
sed -i -E "s/k_damping_amp = [0-9]+/k_damping_amp = 0/g" $TMP_CONFIG_FILE
# Use "A"-adjustment
echo "Running with A-Adjustment..."
sed -i -E "s/a_adj_amp = [0-9]+/a_adj_amp = 1/g" $TMP_CONFIG_FILE
./cosmo $TMP_CONFIG_FILE
sed -i -E "s/a_adj_amp = [0-9]+/a_adj_amp = 0/g" $TMP_CONFIG_FILE
# No metric norming
echo "Running with No metric norming..."
sed -i -E "s/normalize_metric = [0-9]+/normalize_metric = 0/g" $TMP_CONFIG_FILE
./cosmo $TMP_CONFIG_FILE
sed -i -E "s/normalize_metric = [0-9]+/normalize_metric = 1/g" $TMP_CONFIG_FILE
# No non-linear behavior (small amplitude)
sed -i -E "s/noise_amp = 1.0e-10/noise_amp = 1.0e-20/g" $TMP_CONFIG_FILE
./cosmo $TMP_CONFIG_FILE
sed -i -E "s/noise_amp = 1.0e-20/noise_amp = 1.0e-10/g" $TMP_CONFIG_FILE

# Use Z4c terms
echo "Using Z4c terms..."
cmake -DCOSMO_N=50 -DCOSMO_NY=6 -DCOSMO_NZ=6 -DCOSMO_USE_Z4c_DAMPING=1 .. && make -j16
./cosmo $TMP_CONFIG_FILE

# Vary resolutions
cmake -DCOSMO_N=50 -DCOSMO_NY=4 -DCOSMO_NZ=4 .. && make -j16
./cosmo ../config/AwA/stability_test_r1.txt
cmake -DCOSMO_N=100 -DCOSMO_NY=4 -DCOSMO_NZ=4 .. && make -j16
./cosmo ../config/AwA/stability_test_r2.txt
cmake -DCOSMO_N=200 -DCOSMO_NY=4 -DCOSMO_NZ=4 .. && make -j16
./cosmo ../config/AwA/stability_test_r4.txt

rm $TMP_CONFIG_FILE
