#!/bin/bash

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# And up a directory should be the main codebase.
cd ..

cd build

## OpenMP Scaling

## Strong scaling runs
cmake -DCOSMO_N=64 .. && make -j4
sed -i 's/omp_num_threads = 4/omp_num_threads = 1/g' ../config/stability_test.txt
./cosmo ../config/stability_test.txt
sed -i 's/omp_num_threads = 1/omp_num_threads = 2/g' ../config/stability_test.txt
./cosmo ../config/stability_test.txt
sed -i 's/omp_num_threads = 2/omp_num_threads = 4/g' ../config/stability_test.txt
./cosmo ../config/stability_test.txt
sed -i 's/omp_num_threads = 4/omp_num_threads = 8/g' ../config/stability_test.txt
./cosmo ../config/stability_test.txt
sed -i 's/omp_num_threads = 8/omp_num_threads = 16/g' ../config/stability_test.txt
./cosmo ../config/stability_test.txt
sed -i 's/omp_num_threads = 16/omp_num_threads = 32/g' ../config/stability_test.txt
./cosmo ../config/stability_test.txt
sed -i 's/omp_num_threads = 32/omp_num_threads = 64/g' ../config/stability_test.txt
./cosmo ../config/stability_test.txt
sed -i 's/omp_num_threads = 64/omp_num_threads = 128/g' ../config/stability_test.txt
./cosmo ../config/stability_test.txt
sed -i 's/omp_num_threads = 128/omp_num_threads = 4/g' ../config/stability_test.txt

## Weak scaling runs
cmake -DCOSMO_N=16 .. && make -j4
sed -i 's/omp_num_threads = 4/omp_num_threads = 1/g' ../config/stability_test.txt
./cosmo ../config/stability_test.txt
cmake -DCOSMO_N=16 -DCOSMO_NX=32 .. && make -j4
sed -i 's/omp_num_threads = 1/omp_num_threads = 2/g' ../config/stability_test.txt
./cosmo ../config/stability_test.txt
cmake -DCOSMO_N=16 -DCOSMO_NX=32 -DCOSMO_NY=32 .. && make -j4
sed -i 's/omp_num_threads = 2/omp_num_threads = 4/g' ../config/stability_test.txt
./cosmo ../config/stability_test.txt
cmake -DCOSMO_N=32 .. && make -j4
sed -i 's/omp_num_threads = 4/omp_num_threads = 8/g' ../config/stability_test.txt
./cosmo ../config/stability_test.txt
cmake -DCOSMO_N=32 -DCOSMO_NX=64 .. && make -j4
sed -i 's/omp_num_threads = 8/omp_num_threads = 16/g' ../config/stability_test.txt
./cosmo ../config/stability_test.txt
cmake -DCOSMO_N=32 -DCOSMO_NX=64 -DCOSMO_NY=64 .. && make -j4
sed -i 's/omp_num_threads = 16/omp_num_threads = 32/g' ../config/stability_test.txt
./cosmo ../config/stability_test.txt
cmake -DCOSMO_N=64 .. && make -j4
sed -i 's/omp_num_threads = 32/omp_num_threads = 64/g' ../config/stability_test.txt
./cosmo ../config/stability_test.txt
cmake -DCOSMO_N=64 -DCOSMO_NX=128 .. && make -j4
sed -i 's/omp_num_threads = 64/omp_num_threads = 128/g' ../config/stability_test.txt
./cosmo ../config/stability_test.txt
sed -i 's/omp_num_threads = 128/omp_num_threads = 4/g' ../config/stability_test.txt

