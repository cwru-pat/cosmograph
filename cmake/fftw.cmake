# FFT libraries for ICs

find_library(FFTW_LIBRARY
     NAMES fftw3 libfftw3 fftw libfftw
     PATHS /home/jbm120/fftw/lib /usr/lib /usr/local/fftw/3.3.4/lib)

set(FFTW_LIBRARIES "${FFTW_LIBRARY}")
message(STATUS " FFTW_LIBRARY: ${FFTW_LIBRARY}")

# Custom fftw dir for cluster
if(EXISTS /home/jbm120/fftw/include)
  include_directories(/home/jbm120/fftw/include)
endif()
