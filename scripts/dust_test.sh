#!/bin/bash

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# And up two directories should be the main codebase.
pwd

###
# Compile
###
mkdir -p ../build
cd ../build
TMP_CONFIG_FILE=long_dust_test.txt.test
cp ../config/tests/long_dust_test.txt $TMP_CONFIG_FILE

do_runs()
{
  # "normal" run
  ./cosmo $TMP_CONFIG_FILE
  # KO Dissipation run
  echo "Running with KO Dissipation..."
  sed -i -E "s/KO_damping_coefficient = [0-9]+/KO_damping_coefficient = 0.1/g" $TMP_CONFIG_FILE
  ./cosmo $TMP_CONFIG_FILE
  sed -i -E "s/KO_damping_coefficient = [0-9\.]+/KO_damping_coefficient = 0/g" $TMP_CONFIG_FILE
  # Use "K"-damping
  echo "Running with K Damping..."
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
}

cmake -DCOSMO_N=16 .. && make -j32
sed -i -E "s/steps = [0-9]+/steps = 5000/g" $TMP_CONFIG_FILE
sed -i -E "s/output_dir = long_dust_test/output_dir = long_dust_test_r16/g" $TMP_CONFIG_FILE
do_runs

cmake -DCOSMO_N=32 .. && make -j32
sed -i -E "s/steps = [0-9]+/steps = 10000/g" $TMP_CONFIG_FILE
sed -i -E "s/output_dir = long_dust_test_r16/output_dir = long_dust_test_r32/g" $TMP_CONFIG_FILE
do_runs

cmake -DCOSMO_N=64 .. && make -j32
sed -i -E "s/steps = [0-9]+/steps = 20000/g" $TMP_CONFIG_FILE
sed -i -E "s/output_dir = long_dust_test_r32/output_dir = long_dust_test_r64/g" $TMP_CONFIG_FILE
do_runs

rm $TMP_CONFIG_FILE
