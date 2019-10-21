
project(COSMO)
cmake_minimum_required(VERSION 2.6.0)

# coloring for messages
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColorReset  "${Esc}[m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Cyan        "${Esc}[36m")
endif()

# optimization and debugging
include(cmake/optimization.cmake)

# check for OpenMP
include(cmake/openmp.cmake)

# check compiler support
include(cmake/compiler.cmake)

# compilation options
include(cmake/options.cmake)

# check for raytracing, multigrid modules
include(cmake/cosmotrace.cmake)
include(cmake/multigrid.cmake)

# FFT libraries
include(cmake/fftw.cmake)

# HDF5 libraries
include(cmake/hdf5.cmake)

unset(COSMO_SOURCES CACHE)
file(GLOB COSMO_SOURCES *.cc IO/*.cc utils/*.cc
  components/static/*.cc components/particles/*.cc
  components/scalar/*.cc components/bssn/*.cc components/phase_space_sheet/*.cc
  components/Lambda/*.cc components/dust_fluid/*.cc sims/*.cc ICs/*.cc)

add_executable(cosmo ${COSMO_SOURCES} ${MULTIGRID_SOURCES})
target_link_libraries(cosmo m rt z ${HDF5_LINK_LIBRARY} ${FFTW_LIBRARIES} ${FFTWL_LIBRARIES})
