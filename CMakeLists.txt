
project(COSMO)
cmake_minimum_required(VERSION 2.6.0)

# add -g for valgrind
set(OPT_LEVEL         "-O3 -ffast-math")
set(WARNINGS          "-pedantic -Wall -Wno-unused-variable -Wno-unused-function")
set(CC_OPTS           "-march=native -g")
set(CMAKE_CXX_FLAGS   "${CC_OPTS} ${OPT_LEVEL} ${WARNINGS}")

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
     PATHS /home/jbm120/fftw/lib /usr/lib /usr/local/fftw/3.3.4/lib)
set(FFTW_LIBRARIES "${FFTW_LIBRARY}")
message(STATUS " FFTW_LIBRARY: ${FFTW_LIBRARY}")
include_directories(/home/jbm120/fftw/include)

# HDF5 libraries
find_library(HDF5_LIBRARY
     NAMES hdf5 libhdf5
     PATHS /home/jbm120/hdf5/build/lib /usr/lib)
set(HDF5_LIBRARIES "${HDF5_LIBRARY}")
message(STATUS " HDF5_LIBRARY: ${HDF5_LIBRARY}")
include_directories(/home/jbm120/hdf5/build/include)


file(GLOB COSMO_SOURCES *.cc utils/*.cc)

add_executable(cosmo ${COSMO_SOURCES})
target_link_libraries(cosmo m rt z ${HDF5_LIBRARY} ${FFTW_LIBRARIES})
