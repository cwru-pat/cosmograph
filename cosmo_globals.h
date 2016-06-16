#ifndef COSMO_GLOBALS_H
#define COSMO_GLOBALS_H

#include "cosmo_types.h"
#include "cosmo_macros.h"

#ifndef dt
  extern cosmo::real_t dt;
#endif
#ifndef dx
  extern cosmo::real_t dx;
#endif

#include "utils/Timer.h"
#include "utils/ConfigParser.h"

extern cosmo::TimerManager _timer;
extern cosmo::ConfigParser _config;

#endif
