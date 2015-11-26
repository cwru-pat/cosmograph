#ifndef COSMO_ICS
#define COSMO_ICS

#include "cosmo.h"
#include "globals.h"

#include "ICs_data.h"
#include "io_data.h"
#include "utils/reference_frw.h"

namespace cosmo
{

ICsData cosmo_get_ICsData();

real_t cosmo_power_spectrum(real_t k, ICsData *icd);

void set_gaussian_random_field(real_t *field, Fourier *fourier, ICsData *icd);

void set_conformal_ICs(
  std::map <std::string, arr_t *> & bssn_fields,
  Fourier *fourier, IOData *iod, FRW<real_t> *frw);

void set_stability_test_ICs(
  std::map <std::string, arr_t *> & bssn_fields);

void set_linear_wave_ICs(
  std::map <std::string, arr_t *> & bssn_fields);

} // namespace cosmo

#endif
