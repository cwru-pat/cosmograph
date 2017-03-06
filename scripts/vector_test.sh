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
TMP_CONFIG_FILE=particles_vector.txt.test
cp ../config/particles_vector.txt $TMP_CONFIG_FILE

cmake -DCOSMO_N=16 -DCOSMO_NX=1 -DCOSMO_NZ=1 -DCOSMO_STENCIL_ORDER=2 .. && make -j32
sed -i -E "s/steps = [0-9]+/steps = 50/g" $TMP_CONFIG_FILE
sed -i -E "s/particle_vectorpert_test/particle_vectorpert_test_r0016/g" $TMP_CONFIG_FILE
./cosmo $TMP_CONFIG_FILE

cmake -DCOSMO_N=32 -DCOSMO_NX=1 -DCOSMO_NZ=1 -DCOSMO_STENCIL_ORDER=2 .. && make -j32
sed -i -E "s/steps = [0-9]+/steps = 100/g" $TMP_CONFIG_FILE
sed -i -E "s/particle_vectorpert_test_r0016/particle_vectorpert_test_r0032/g" $TMP_CONFIG_FILE
./cosmo $TMP_CONFIG_FILE

cmake -DCOSMO_N=64 -DCOSMO_NX=1 -DCOSMO_NZ=1 -DCOSMO_STENCIL_ORDER=2 .. && make -j32
sed -i -E "s/steps = [0-9]+/steps = 200/g" $TMP_CONFIG_FILE
sed -i -E "s/particle_vectorpert_test_r0032/particle_vectorpert_test_r0064/g" $TMP_CONFIG_FILE
./cosmo $TMP_CONFIG_FILE

cmake -DCOSMO_N=128 -DCOSMO_NX=1 -DCOSMO_NZ=1 -DCOSMO_STENCIL_ORDER=2 .. && make -j32
sed -i -E "s/steps = [0-9]+/steps = 400/g" $TMP_CONFIG_FILE
sed -i -E "s/particle_vectorpert_test_r0064/particle_vectorpert_test_r0128/g" $TMP_CONFIG_FILE
./cosmo $TMP_CONFIG_FILE

# cmake -DCOSMO_N=1024 -DCOSMO_NX=1 -DCOSMO_NZ=1 -DCOSMO_STENCIL_ORDER=2 .. && make -j32
# sed -i -E "s/steps = [0-9]+/steps = 3200/g" $TMP_CONFIG_FILE
# sed -i -E "s/particle_vectorpert_test_r128/particle_vectorpert_test_r1024/g" $TMP_CONFIG_FILE
# ./cosmo $TMP_CONFIG_FILE

rm $TMP_CONFIG_FILE
