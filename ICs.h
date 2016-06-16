#ifndef COSMO_ICS
#define COSMO_ICS

#include "cosmo_includes.h"
#include "cosmo_types.h"
#include "globals.h"

#include "ICs_data.h"

#include "utils/Fourier.h"

#include "cosmotrace/raytrace.h"

namespace cosmo
{

ICsData cosmo_get_ICsData();

real_t cosmo_power_spectrum(real_t k, ICsData *icd);
void set_gaussian_random_field(real_t *field, Fourier *fourier, ICsData *icd);

void init_ray_vector(std::vector<RayTrace<real_t, idx_t> *> * rays,
  idx_t n_rays);

void set_stability_test_ICs( map_t & bssn_fields, map_t & static_field);
void set_linear_wave_ICs(map_t & bssn_fields);

} // namespace cosmo

#endif
