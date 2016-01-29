cosmograph
==========

Cosmological General Relativity And Particle Hydrodynamics

But we'll see about the "particle" part.

Software needed:
 - git
 - hdf5
 - cmake
 - fftw
 - libz
 
To install on on linux (Ubuntu or linux mint): `sudo apt-get install git cmake libhdf5-dev fftw-dev libzip-dev`

How to set up:
 - 1) Clone code: `git clone https://github.com/jbcm627/cosmograph`
 - 2) CD into directory: `cd cosmograph`
 - 3) Update submodules: `git submodule update --init --recursive`
 - 4) Make build directory `mkdir build` and chage into it: `cd build`
 - 5) Compile code `cmake ..` then `make`
