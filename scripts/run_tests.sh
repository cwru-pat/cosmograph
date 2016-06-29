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
sed -i.bak "s/define N [0-9]\+/define N 16/" ../cosmo_macros.h
cmake ..
if [ $? -ne 0 ]; then
    echo "Error: cmake failed!"
    exit 1
fi
make
if [ $? -ne 0 ]; then
    echo "Error: make failed!"
    exit 1
fi
# Test run
./cosmo ../config/dust_test.txt
if [ $? -ne 0 ]; then
    echo "Error: run failed!"
    exit 1
fi
