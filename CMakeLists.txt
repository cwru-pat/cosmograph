
project(COSMO)
cmake_minimum_required(VERSION 2.6.0)

set(OPT_LEVEL        "-O3 -g -ffast-math")
set(WARNINGS         "-pedantic -Wall -Wno-unused-variable -Wno-unused-function")
set(CC_OPTS          "-march=native")
set(CMAKE_CXX_FLAGS  "${CC_OPTS} ${OPT_LEVEL} ${WARNINGS}")

# Check for OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Need these standards for random support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# FFT libraries for ICs
find_library(FFTW_LIBRARY
     NAMES fftw3 libfftw3 fftw libfftw
     PATHS /usr/lib)
set(FFTW_LIBRARIES "${FFTW_LIBRARY}")
message(STATUS " FFTW_LIBRARY: ${FFTW_LIBRARY}")

file(GLOB COSMO_SOURCES *.cc utils/*.cc)

add_executable(cosmo ${COSMO_SOURCES})
target_link_libraries(cosmo m rt z hdf5 ${FFTW_LIBRARIES})
