#ifndef COSMO_WAVE_H
#define COSMO_WAVE_H

#include "cosmo.h"

namespace cosmo
{

class Wave
{
public:
  Wave();
  ~Wave();

  void step();
  //void step_boundary();

private:
  /* wave equation fields */
  real_t *RESTRICT phi, *RESTRICT phi_p;
  real_t *RESTRICT www, *RESTRICT www_p;
};

}

#endif
