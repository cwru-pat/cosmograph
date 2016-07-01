
# Make sure compiler supports c++11 standards

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  message(WARNING "${Yellow}The compiler ${CMAKE_CXX_COMPILER} has no C++11 support, but has c++0x support. Attempting to build anyways; compile may fail.${ColorReset}")
else()
  message(FATAL_ERROR "${Red}The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different compiler.${ColorReset}")
endif()
