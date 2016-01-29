
project(COSMO)
cmake_minimum_required(VERSION 2.6.0)

# add -g for valgrind
set(OPT_LEVEL         " -O4 -fast")
set(CC_OPTS           "-std=c++11 -acc -Minfo=accel")
set(CMAKE_CXX_FLAGS   "${CC_OPTS} ${OPT_LEVEL} ${WARNINGS}")

# FFT libraries for ICs
find_library(FFTW_LIBRARY
     NAMES fftw3 libfftw3 fftw libfftw
     PATHS /usr/lib)
set(FFTW_LIBRARIES "${FFTW_LIBRARY}")
message(STATUS " FFTW_LIBRARY: ${FFTW_LIBRARY}")

file(GLOB COSMO_SOURCES *.cc utils/*.cc)

SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
add_executable(cosmo ${COSMO_SOURCES})
target_link_libraries(cosmo m rt z hdf5 ${FFTW_LIBRARIES})
