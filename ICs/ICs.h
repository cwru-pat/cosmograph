#ifndef COSMO_ICS
#define COSMO_ICS

#include "../cosmo_includes.h"
#include "../cosmo_types.h"
#include "../cosmo_globals.h"

#include "ICs_data.h"

#include "../utils/Fourier.h"

#include "../cosmotrace/raytrace.h"

namespace cosmo
{

ICsData cosmo_get_ICsData();

real_t cosmo_power_spectrum(real_t k, ICsData *icd);
void set_gaussian_random_field(arr_t & field, Fourier *fourier, ICsData *icd);

void init_ray_vector(std::vector<RayTrace<real_t, idx_t> *> * rays,
  idx_t n_rays);

} // namespace cosmo

#endif
