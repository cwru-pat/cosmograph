#!/bin/bash

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# up a directory should be the main codebase.

cd ..

mkdir -p build
cd build

mkdir -p SVT_test
cd SVT_test

cp ../../config/dust_SVT_test.txt config.txt


# "Standard" run

sed -i.bak -r "s/output_dir = [_A-Za-z0-9]+/output_dir = SVT_test_Hfrac001/" config.txt
sed -i.bak -r "s/steps = [0-9]+/steps = 60000/" config.txt
sed -i.bak -r "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = 60/" config.txt
sed -i.bak -r "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = 60/" config.txt
sed -i.bak -r "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = 60/" config.txt
cmake -DCOSMO_N=48 -DCOSMO_NY=4 -DCOSMO_NZ=4 -DCOSMO_STENCIL_ORDER=8 -DCOSMO_H_LEN_FRAC=1 ../.. && make -j16
./cosmo config.txt

sed -i.bak -r "s/output_dir = [_A-Za-z0-9]+/output_dir = SVT_test_Hfrac002/" config.txt
sed -i.bak -r "s/steps = [0-9]+/steps = 120000/" config.txt
sed -i.bak -r "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = 120/" config.txt
sed -i.bak -r "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = 120/" config.txt
sed -i.bak -r "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = 120/" config.txt
cmake -DCOSMO_N=48 -DCOSMO_NY=4 -DCOSMO_NZ=4 -DCOSMO_STENCIL_ORDER=8 -DCOSMO_H_LEN_FRAC=0.5 ../.. && make -j16
./cosmo config.txt

sed -i.bak -r "s/output_dir = [_A-Za-z0-9]+/output_dir = SVT_test_Hfrac010/" config.txt
sed -i.bak -r "s/steps = [0-9]+/steps = 600000/" config.txt
sed -i.bak -r "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = 600/" config.txt
sed -i.bak -r "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = 600/" config.txt
sed -i.bak -r "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = 600/" config.txt
cmake -DCOSMO_N=48 -DCOSMO_NY=4 -DCOSMO_NZ=4 -DCOSMO_STENCIL_ORDER=8 -DCOSMO_H_LEN_FRAC=0.1 ../.. && make -j16
./cosmo config.txt

sed -i.bak -r "s/output_dir = [_A-Za-z0-9]+/output_dir = SVT_test_Hfrac020/" config.txt
sed -i.bak -r "s/steps = [0-9]+/steps = 1200000/" config.txt
sed -i.bak -r "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = 1200/" config.txt
sed -i.bak -r "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = 1200/" config.txt
sed -i.bak -r "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = 1200/" config.txt
cmake -DCOSMO_N=48 -DCOSMO_NY=4 -DCOSMO_NZ=4 -DCOSMO_STENCIL_ORDER=8 -DCOSMO_H_LEN_FRAC=0.05 ../.. && make -j16
./cosmo config.txt

sed -i.bak -r "s/output_dir = [_A-Za-z0-9]+/output_dir = SVT_test_Hfrac100/" config.txt
sed -i.bak -r "s/steps = [0-9]+/steps = 6000000/" config.txt
sed -i.bak -r "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = 6000/" config.txt
sed -i.bak -r "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = 6000/" config.txt
sed -i.bak -r "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = 6000/" config.txt
cmake -DCOSMO_N=48 -DCOSMO_NY=4 -DCOSMO_NZ=4 -DCOSMO_STENCIL_ORDER=8 -DCOSMO_H_LEN_FRAC=0.01 ../.. && make -j16
./cosmo config.txt
