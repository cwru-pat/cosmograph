#!/bin/bash

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# up a directory should be the main codebase.

cd ..

mkdir -p build
cd build

mkdir -p sinusoid_gauge_test
cd sinusoid_gauge_test

cp ../../config/particles_sinusoid.txt config.txt


do_runs () {
  sed -i.bak -r "s/k_driver_coeff = [\.0-9]+/k_driver_coeff = 0.00/" config.txt
  ./cosmo config.txt

  sed -i.bak -r "s/k_driver_coeff = [\.0-9]+/k_driver_coeff = 0.01/" config.txt
  ./cosmo config.txt

  sed -i.bak -r "s/k_driver_coeff = [\.0-9]+/k_driver_coeff = 0.1/" config.txt
  ./cosmo config.txt
}

# "Standard" run

sed -i.bak -r "s/output_dir = [_A-Za-z0-9]+/output_dir = particle_sinusoid_test_Hfrac1/" config.txt
sed -i.bak -r "s/steps = [0-9]+/steps = 10000/" config.txt
sed -i.bak -r "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = 10/" config.txt
sed -i.bak -r "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = 10/" config.txt
sed -i.bak -r "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = 10/" config.txt
cmake -DCOSMO_N=8 -DCOSMO_STENCIL_ORDER=8 -DCOSMO_H_LEN_FRAC=1 ../.. && make -j16
do_runs

sed -i.bak -r "s/output_dir = [_A-Za-z0-9]+/output_dir = particle_sinusoid_test_Hfrac2/" config.txt
sed -i.bak -r "s/steps = [0-9]+/steps = 20000/" config.txt
sed -i.bak -r "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = 20/" config.txt
sed -i.bak -r "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = 20/" config.txt
sed -i.bak -r "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = 20/" config.txt
cmake -DCOSMO_N=8 -DCOSMO_STENCIL_ORDER=8 -DCOSMO_H_LEN_FRAC=0.5 ../.. && make -j16
do_runs

sed -i.bak -r "s/output_dir = [_A-Za-z0-9]+/output_dir = particle_sinusoid_test_Hfrac10/" config.txt
sed -i.bak -r "s/steps = [0-9]+/steps = 100000/" config.txt
sed -i.bak -r "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = 100/" config.txt
sed -i.bak -r "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = 100/" config.txt
sed -i.bak -r "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = 100/" config.txt
cmake -DCOSMO_N=8 -DCOSMO_STENCIL_ORDER=8 -DCOSMO_H_LEN_FRAC=0.1 ../.. && make -j16
do_runs

sed -i.bak -r "s/output_dir = [_A-Za-z0-9]+/output_dir = particle_sinusoid_test_Hfrac20/" config.txt
sed -i.bak -r "s/steps = [0-9]+/steps = 200000/" config.txt
sed -i.bak -r "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = 200/" config.txt
sed -i.bak -r "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = 200/" config.txt
sed -i.bak -r "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = 200/" config.txt
cmake -DCOSMO_N=8 -DCOSMO_STENCIL_ORDER=8 -DCOSMO_H_LEN_FRAC=0.05 ../.. && make -j16
do_runs

sed -i.bak -r "s/output_dir = [_A-Za-z0-9]+/output_dir = particle_sinusoid_test_Hfrac100/" config.txt
sed -i.bak -r "s/steps = [0-9]+/steps = 1000000/" config.txt
sed -i.bak -r "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = 1000/" config.txt
sed -i.bak -r "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = 1000/" config.txt
sed -i.bak -r "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = 1000/" config.txt
cmake -DCOSMO_N=8 -DCOSMO_STENCIL_ORDER=8 -DCOSMO_H_LEN_FRAC=0.01 ../.. && make -j16
do_runs

sed -i.bak -r "s/output_dir = [_A-Za-z0-9]+/output_dir = particle_sinusoid_test_Hfrac100/" config.txt
sed -i.bak -r "s/steps = [0-9]+/steps = 1000000/" config.txt
sed -i.bak -r "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = 1000/" config.txt
sed -i.bak -r "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = 1000/" config.txt
sed -i.bak -r "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = 1000/" config.txt
cmake -DCOSMO_N=8 -DCOSMO_STENCIL_ORDER=8 -DCOSMO_H_LEN_FRAC=0.01 ../.. && make -j16
do_runs

