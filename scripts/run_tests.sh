#!/bin/bash

echo "Running some tests."

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# And up a directory should be the main codebase.
cd ..

if [ -z "$CXX" ]; then
    echo "C compiler used is g++"
    CXX=g++
else
    echo "C compiler used is ${CXX}"
fi
$CXX --version

###
# Tests in tests directory
###
cd tests
$CXX --std=c++11 reference_frw.cc && ./a.out
if [ $? -ne 0 ]; then
    echo "Error: reference FRW integrator check failed!"
    exit 1
fi
$CXX --std=c++11 timer.cc ../utils/Timer.cc -lrt -O0 && ./a.out
if [ $? -ne 0 ]; then
    echo "Error: timer class check failed!"
    exit 1
fi
$CXX --std=c++11 array.cc -O0 && ./a.out
if [ $? -ne 0 ]; then
    echo "Error: array class check failed!"
    exit 1
fi
$CXX --std=c++11 rk4.cc -O0 && ./a.out
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
cmake ..
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
./cosmo ../config/tests/vacuum_test.txt
if [ $? -ne 0 ]; then
    echo "Error: stability run failed!"
    exit 1
fi

echo ""
echo "Running dust test"
echo "-----------------"
./cosmo ../config/tests/dust_ray_test.txt
if [ $? -ne 0 ]; then
    echo "Error: dust test run failed!"
    exit 1
fi

echo ""
echo "Running dust+lambda test"
echo "------------------------"
./cosmo ../config/tests/lambda_test.txt
if [ $? -ne 0 ]; then
    echo "Error: dust+lambda test run failed!"
    exit 1
fi

echo ""
echo "Running particles test"
echo "----------------------"
./cosmo ../config/tests/particles_test.txt
if [ $? -ne 0 ]; then
    echo "Error: particle run failed!"
    exit 1
fi

echo ""
echo "Running scalar test"
echo "-------------------"
./cosmo ../config/tests/scalar_multigrid_test.txt
if [ $? -ne 0 ]; then
    echo "Error: scalar run failed!"
    exit 1
fi



echo ""
echo "Compiling with shift & gamma driver support"
echo "-------------------------------------------"
cmake -DCOSMO_USE_GAMMA_DRIVER=1 -DCOSMO_USE_BSSN_SHIFT=1 ..
if [ $? -ne 0 ]; then
    echo "Error: cmake failed!"
    exit 1
fi
make -j16
if [ $? -ne 0 ]; then
    echo "Error: make failed!"
    exit 1
fi

echo ""
echo "Running harmonic gauge & gamma driver test"
echo "------------------------------------------"
./cosmo ../config/tests/gauges/scalar_harmonic_gammadriver_test.txt
if [ $? -ne 0 ]; then
    echo "Error: scalar run failed!"
    exit 1
fi

echo ""
echo "Running 1+log gauge & gamma driver test"
echo "---------------------------------------"
./cosmo ../config/tests/gauges/scalar_1pluslog_gammadriver_test.txt
if [ $? -ne 0 ]; then
    echo "Error: scalar run failed!"
    exit 1
fi

echo ""
echo "Running damped wave gauge test"
echo "------------------------------"
./cosmo ../config/tests/gauges/scalar_dampedwave_test.txt
if [ $? -ne 0 ]; then
    echo "Error: scalar run failed!"
    exit 1
fi

