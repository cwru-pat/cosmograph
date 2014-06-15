
project(COSMO)
cmake_minimum_required(VERSION 2.6.0)

set(OPT_LEVEL "-O3 -ffast-math")
set(WARNINGS "-pedantic -Wall -Wno-unused-variable -Wno-unused-function")
set(CC_OPTS "-march=native")
set(CMAKE_CXX_FLAGS "${CC_OPTS} ${OPT_LEVEL} ${WARNINGS}")

file(GLOB COSMO_SOURCES *.cc utils/*.cc)

add_executable(cosmo ${COSMO_SOURCES})
target_link_libraries(cosmo m rt z)
