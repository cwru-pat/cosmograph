#!/bin/bash

echo "Running some tests."

###
# Tests in tests directory
###
cd tests
g++ --std=c++11 reference_frw.cc && ./a.out
if [ $? -ne 0 ]; then
    echo "Error: reference FRW integrator check failed!"
    exit 1
fi
g++ --std=c++11 timer.cc ../utils/Timer.cc -lrt -O0 && ./a.out
if [ $? -ne 0 ]; then
    echo "Error: timer class check failed!"
    exit 1
fi
g++ --std=c++11 array.cc -O0 && ./a.out
if [ $? -ne 0 ]; then
    echo "Error: array class check failed!"
    exit 1
fi
g++ --std=c++11 rk4.cc -O0 && ./a.out
if [ $? -ne 0 ]; then
    echo "Error: RK4 class check failed!"
    exit 1
fi
rm a.out

###
# Compile
###
mkdir -p ../build
cd ../build
# Run a 16^3 test
cmake -DCMAKE_CXX_COMPILER=g++ -DCOSMO_N=16 -DCOSMO_USE_REFERENCE_FRW=1 ..
if [ $? -ne 0 ]; then
    echo "Error: cmake failed!"
    exit 1
fi
make
if [ $? -ne 0 ]; then
    echo "Error: make failed!"
    exit 1
fi
# Test runs
./cosmo ../config/dust_test.txt
if [ $? -ne 0 ]; then
    echo "Error: dust test run failed!"
    exit 1
fi
./cosmo ../config/lambda_test.txt
if [ $? -ne 0 ]; then
    echo "Error: dust+lambda test run failed!"
    exit 1
fi
./cosmo ../config/scalar_test.txt
if [ $? -ne 0 ]; then
    echo "Error: run failed!"
    exit 1
fi
