SET(COSMO_DEBUG FALSE CACHE STRING "Turn on debug mode (turn off optimizations, turn on valgrind")
SET(COSMO_PROFILE FALSE CACHE STRING "Turn on profiling (gprof flag)")
SET(COSMO_STATIC FALSE CACHE STRING "Static build (-static flag)")


if(COSMO_DEBUG)
  # add -g for valgrind
  message(STATUS "${Cyan} Debug build enabled (-g flag).${ColorReset}")
  set(PROFILING     "-g")
  set(OPT_LEVEL     "-O1")
  set(CC_OPTS       "")
else()
  set(PROFILING     "")
  set(OPT_LEVEL     "-O3 -ffast-math -flto")
  set(CC_OPTS       "-march=native")
endif()

# Statically linked build?
if(COSMO_STATIC)
  message(STATUS "${Cyan} Static build enabled (-static flag).${ColorReset}")
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

# add -pg for gprof
if(COSMO_PROFILE)
  message(STATUS "${Cyan} Profiling enabled (-pg flag).${ColorReset}")
  set(PROFILING "${PROFILING} -pg")
endif()


set(WARNINGS          "-pedantic -Wall -Wno-unused-variable -Wno-unused-function")
set(CMAKE_CXX_FLAGS   "${CC_OPTS} ${OPT_LEVEL} ${WARNINGS} ${PROFILING}")
set(CMAKE_EXE_LINKER_FLAGS  "${PROFILING}")

unset(COSMO_DEBUG CACHE)
unset(COSMO_STATIC CACHE)
unset(COSMO_PROFILE CACHE)
