#ifndef COSMO_GLOBALS_H
#define COSMO_GLOBALS_H

//#include "cosmo_types.h"
#include "cosmo_macros.h"

#ifndef dt
  extern double dt;
#endif
#ifndef dx
  extern double dx;
#endif

// define a global variable for box size
// will only be usded in sheets and particles 1D
// simulation to reduce it to 1D
// dx will always be set equally for all 3 directions
#ifndef Lx
  extern double Lx;
#endif

#ifndef Ly
  extern double Ly;
#endif

#ifndef Lz
  extern double Lz;
#endif

#ifndef Nx
  extern int Nx;
#endif

#ifndef Ny
  extern int Ny;
#endif

#ifndef Nz
  extern int Nz;
#endif



#include "utils/Timer.h"
#include "utils/ConfigParser.h"

extern cosmo::TimerManager _timer;
extern cosmo::ConfigParser _config;

#endif
