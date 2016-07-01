# cosmograph

[![Build Status](https://travis-ci.com/jbcm627/cosmograph.svg?token=j5zJrjKFZL3UXL3HwPp6&branch=master)](https://travis-ci.com/jbcm627/cosmograph)

## Cosmological General Relativity And (Perfect fluid | Particle) Hydrodynamics

*But we'll see about the "particle" part.*

### Software dependencies

 - git
 - hdf5
 - cmake
 - fftw3
 - libz
 - a compiler with c++11 support

To install these on on linux (Ubuntu, Mint), run a command like:

```{r, engine='bash', compile}
sudo apt-get install git cmake libhdf5-dev fftw3-dev libzip-dev g++
```

### Setting up the code
 
 - 1) Clone the code: `git clone https://github.com/jbcm627/cosmograph`
 - 2) Change into the cloned repository directory: `cd cosmograph`
 - 3) Initialize submodules: `git submodule update --init --recursive`
 - 4) Create a build directory `mkdir build` and chage into it: `cd build`
 - 5) Compile code using `cmake ..` then `make`

To the extent different compile-time options are supported, they should be
enumerated in the corresponding [cmake](https://github.com/jbcm627/cosmograph/blob/master/cmake/options.cmake)
file. You can adjust these by supplying an appropriate argument to cmake. For
example, to change the resolution to `N=32`, the compile command will be

```{r, engine='bash', compile}
cmake -DCOSMO_N=32 ..
```

### Running the code

The executable should accept a configuration file as a parameter. The
configuration parameter contains various runtime options that can be 
updated without needing to recompile. For example, in order to run a
'test' simulation of a pure-dust universe, the corresponding command
could be

```{r, engine='bash', compile}
./cosmo ../config/dust_test.txt
```

### Support

CosmoGRaPH is under heavy development, and as such, all functionality is
subject to change. You are nevertheless welcome to use and modify the code
in its present form. Some documentation may be available in the `docs`
directory, although this is not guaranteed to be up-to-date.
