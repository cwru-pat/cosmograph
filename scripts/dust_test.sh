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
}

cmake -DCOSMO_N=20 .. && make -j32
sed -i -E "s/steps = [0-9]+/steps = 6000/g" $TMP_CONFIG_FILE
sed -i -E "s/output_dir = long_dust_test/output_dir = long_dust_test_r20/g" $TMP_CONFIG_FILE
do_runs

cmake -DCOSMO_N=30 .. && make -j32
sed -i -E "s/steps = [0-9]+/steps = 9000/g" $TMP_CONFIG_FILE
sed -i -E "s/output_dir = long_dust_test_r20/output_dir = long_dust_test_r30/g" $TMP_CONFIG_FILE
do_runs

cmake -DCOSMO_N=40 .. && make -j32
sed -i -E "s/steps = [0-9]+/steps = 12000/g" $TMP_CONFIG_FILE
sed -i -E "s/output_dir = long_dust_test_r30/output_dir = long_dust_test_r4/g" $TMP_CONFIG_FILE
do_runs

rm $TMP_CONFIG_FILE
