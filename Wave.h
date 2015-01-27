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

  void init();
  void step();

private:
  /* wave equation fields */
  real_t *RESTRICT phi, *RESTRICT phi_p;
  real_t *RESTRICT www, *RESTRICT www_p;

  /* easy way to access/iterate over fields?
   * (will eventually have many many of these...) */
  std::map <std::string, real_t *> fields;

};

}

#endif
