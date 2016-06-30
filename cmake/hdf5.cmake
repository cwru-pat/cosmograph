# HDF5 libraries
find_library(HDF5_LIBRARY
     NAMES hdf5 libhdf5
     PATHS /home/jbm120/hdf5/build/lib /usr/lib)
set(HDF5_LIBRARIES "${HDF5_LIBRARY}")
message(STATUS " HDF5_LIBRARY: ${HDF5_LIBRARY}")

if(EXISTS /home/jbm120/hdf5/build/include)
  include_directories(/home/jbm120/hdf5/build/include)
endif()
