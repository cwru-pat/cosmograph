#ifndef COSMO_ICS
#define COSMO_ICS

#include "../cosmo_includes.h"
#include "../cosmo_types.h"
#include "../cosmo_globals.h"

#include "ICs_data.h"

#include "../utils/Fourier.h"

#if USE_COSMOTRACE
#include "../components/cosmotrace/raytrace.h"
#endif

namespace cosmo
{

ICsData cosmo_get_ICsData();

real_t cosmo_power_spectrum(real_t k, ICsData *icd);
void set_gaussian_random_field(arr_t & field, Fourier *fourier, ICsData *icd);

#if USE_COSMOTRACE
void init_ray_vector(std::vector<RayTrace<real_t, idx_t> *> * rays,
  idx_t n_rays);
void init_healpix_ray_vectors(std::vector<RayTrace<real_t, idx_t> *> * rays);
void init_random_ray_vectors(std::vector<RayTrace<real_t, idx_t> *> * rays, idx_t n_rays);
#endif

} // namespace cosmo

#endif
