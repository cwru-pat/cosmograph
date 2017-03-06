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

# Temporary file for many misc. tests
# Fewer steps here
TMP_CONFIG_FILE=linear_wave_test_r4.text.tmp.test
cp ../config/AwA/linear_wave_test_r4.txt $TMP_CONFIG_FILE

# KO Dissipation run
cmake -DCOSMO_H_LEN_FRAC=1 -DCOSMO_N=200 -DCOSMO_NY=6 -DCOSMO_NZ=6 .. && make -j16

echo "Running with KO Dissipation..."
sed -i -E "s/KO_damping_coefficient = [0-9]+/KO_damping_coefficient = 1/g" $TMP_CONFIG_FILE
./cosmo $TMP_CONFIG_FILE
sed -i -E "s/KO_damping_coefficient = [0-9]+/KO_damping_coefficient = 0/g" $TMP_CONFIG_FILE
# Use "K"-damping
cmake -DCOSMO_H_LEN_FRAC=1 -DCOSMO_N=200 -DCOSMO_NY=4 -DCOSMO_NZ=4 .. && make -j16
echo "Running with K Damping..."
sed -i -E "s/k_damping_amp = [0-9]+/k_damping_amp = 1/g" $TMP_CONFIG_FILE
./cosmo $TMP_CONFIG_FILE
sed -i -E "s/k_damping_amp = [0-9]+/k_damping_amp = 0/g" $TMP_CONFIG_FILE
# Use "A"-adjustment
echo "Running with A-Adjustment..."
sed -i -E "s/a_adj_amp = [0-9]+/a_adj_amp = 1/g" $TMP_CONFIG_FILE
./cosmo $TMP_CONFIG_FILE
sed -i -E "s/a_adj_amp = [0-9]+/a_adj_amp = 0/g" $TMP_CONFIG_FILE

# Full, "normal" runs
cmake -DCOSMO_H_LEN_FRAC=1 -DCOSMO_N=50 -DCOSMO_NY=4 -DCOSMO_NZ=4 .. && make -j16
./cosmo ../config/AwA/linear_wave_test_r1.txt
cmake -DCOSMO_H_LEN_FRAC=1 -DCOSMO_N=100 -DCOSMO_NY=4 -DCOSMO_NZ=4 .. && make -j16
./cosmo ../config/AwA/linear_wave_test_r2.txt
cmake -DCOSMO_H_LEN_FRAC=1 -DCOSMO_N=200 -DCOSMO_NY=4 -DCOSMO_NZ=4 .. && make -j16
./cosmo ../config/AwA/linear_wave_test_r4.txt

rm $TMP_CONFIG_FILE
