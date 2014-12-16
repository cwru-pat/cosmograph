#ifndef COSMO_ICS
#define COSMO_ICS

#include "cosmo.h"
#include "globals.h"

#include "ICs_data.h"

namespace cosmo
{

real_t cosmo_power_spectrum(real_t k, ICsData *icd);
void set_gaussian_random_field(void (*Pk)(real_t, real_t, real_t),
  real_t *field, ICsData *icd);


} // namespace cosmo

#endif
