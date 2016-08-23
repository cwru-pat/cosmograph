# FFT libraries for ICs

find_library(FFTW_LIBRARY
  NAMES fftw3 libfftw3 fftw libfftw
  HINTS ENV LD_LIBRARY_PATH)

find_path(FFTW_INCLUDE_DIRS
  NAMES fftw.h
  HINTS ENV CPLUS_INCLUDE_PATH)

if(NOT FFTW_LIBRARY)
  message(FATAL_ERROR "${Red}The FFTW library and include locations were not found. Please make sure you have FFTW installed or loaded, and that the library and include directories can be found in CPLUS_INCLUDE_PATH and LD_LIBRARY_PATH environment variables.${ColorReset}")
else()
  set(FFTW_LIBRARIES "${FFTW_LIBRARY}")
  if(NOT FFTW_INCLUDE_DIRS)
    message(STATUS "FFTW partially found. An include directory was not found for FFTW in the CPLUS_INCLUDE_PATH environment variable.")
    message(STATUS "If compilation or linking fails, you may need to add the include path to this variable.")
    message(STATUS " FFTW_LIBRARY: ${FFTW_LIBRARY}")
  else()
    include_directories("${FFTW_INCLUDE_DIRS}")
    message(STATUS "Found FFTW:")
    message(STATUS " FFTW_LIBRARY: ${FFTW_LIBRARY}")
    message(STATUS " FFTW_INCLUDE_DIRS: ${FFTW_INCLUDE_DIRS}")
  endif()
endif()
