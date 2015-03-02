#ifndef COSMO_ICS
#define COSMO_ICS

#include "cosmo.h"
#include "globals.h"

#include "ICs_data.h"

namespace cosmo
{

real_t cosmo_power_spectrum(real_t k, ICsData *icd);

void set_gaussian_random_field(real_t *field, Fourier *fourier, ICsData *icd);

void set_conformal_ICs(
  std::map <std::string, real_t *> & bssn_fields,
  std::map <std::string, real_t *> & hydro_fields,
  Fourier *fourier,
  ICsData *icd,
  real_t rho_K_matter,
  real_t rho_K_lambda);

void set_flat_dynamic_ICs(
  std::map <std::string, real_t *> & bssn_fields,
  std::map <std::string, real_t *> & hydro_fields,
  Fourier *fourier,
  ICsData *icd,
  real_t rho_K_lambda,
  real_t rho_K_matter);

void set_flat_static_ICs(
  std::map <std::string, real_t *> & bssn_fields,
  std::map <std::string, real_t *> & hydro_fields,
  Fourier *fourier,
  ICsData *icd,
  real_t rho_K_lambda,
  real_t rho_K_matter);

} // namespace cosmo

#endif
