#!/bin/bash

echo "Running some tests."

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# And up a directory should be the main codebase.
cd ..

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
cmake -DCMAKE_CXX_COMPILER=g++ ..
if [ $? -ne 0 ]; then
    echo "Error: cmake failed!"
    exit 1
fi
make -j16
if [ $? -ne 0 ]; then
    echo "Error: make failed!"
    exit 1
fi

# Test runs

echo ""
echo "Running stability test"
echo "----------------------"
./cosmo ../config/AwA/stability_test_r1.txt
if [ $? -ne 0 ]; then
    echo "Error: stability run failed!"
    exit 1
fi

echo ""
echo "Running dust test"
echo "-----------------"
./cosmo ../config/dust_test.txt
if [ $? -ne 0 ]; then
    echo "Error: dust test run failed!"
    exit 1
fi

echo ""
echo "Running dust+lambda test"
echo "------------------------"
./cosmo ../config/lambda_test.txt
if [ $? -ne 0 ]; then
    echo "Error: dust+lambda test run failed!"
    exit 1
fi

echo ""
echo "Running particles test"
echo "----------------------"
./cosmo ../config/particles_test.txt
if [ $? -ne 0 ]; then
    echo "Error: particle run failed!"
    exit 1
fi

echo ""
echo "Running scalar test"
echo "-------------------"
./cosmo ../config/scalar_test.txt
if [ $? -ne 0 ]; then
    echo "Error: scalar run failed!"
    exit 1
fi
