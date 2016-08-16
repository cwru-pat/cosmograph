# FFT libraries for ICs

find_library(FFTW_LIBRARY
     NAMES fftw3 libfftw3 fftw libfftw)

if(NOT FFTW_LIBRARY)
  message(FATAL_ERROR "${Red}The FFTW library was not found. Please make sure you have FFTW installed/loaded, and that the library and include directories can be found in, eg, CPLUS_INCLUDE_PATH and LD_LIBRARY_PATH.${ColorReset}")
else()
  set(FFTW_LIBRARIES "${FFTW_LIBRARY}")
  message(STATUS "FFTW_LIBRARY: ${FFTW_LIBRARY}")
endif()
