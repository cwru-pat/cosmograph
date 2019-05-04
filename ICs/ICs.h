#ifndef COSMO_ICS
#define COSMO_ICS

#include "../cosmo_includes.h"
#include "../cosmo_types.h"
#include "../cosmo_globals.h"

#include "../utils/Fourier.h"

#if USE_COSMOTRACE
#include "../components/cosmotrace/raytrace.h"
#endif

namespace cosmo
{

real_t cosmo_power_spectrum(real_t k, real_t A, real_t k0);
void set_gaussian_random_Phi_N(arr_t & field, Fourier *fourier,
  real_t A, real_t k0, real_t p_cut);

#if USE_COSMOTRACE
void init_ray_vector(std::vector<RayTrace<real_t, idx_t> *> * rays);
void init_healpix_ray_vectors(std::vector<RayTrace<real_t, idx_t> *> * rays);
void init_random_ray_vectors(std::vector<RayTrace<real_t, idx_t> *> * rays);
#endif

} // namespace cosmo

#endif
