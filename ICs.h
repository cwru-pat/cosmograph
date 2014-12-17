#ifndef COSMO_ICS
#define COSMO_ICS

#include "cosmo.h"
#include "globals.h"

#include "ICs_data.h"

namespace cosmo
{

real_t cosmo_power_spectrum(real_t k, ICsData *icd);

void set_gaussian_random_field(real_t *field, Fourier *fourier, ICsData *icd);

void set_physical_metric_and_density(
  std::map <std::string, real_t *> & bssn_fields,
  std::map <std::string, real_t *> & hydro_fields,
  Fourier *fourier);

} // namespace cosmo

#endif
