#ifndef COSMO_ICS
#define COSMO_ICS

#include "cosmo_includes.h"
#include "cosmo_types.h"
#include "globals.h"

#include "utils/Fourier.h"
#include "cosmotrace/raytrace.h"

#include "ICs_data.h"
#include "IO/io.h"
#include "particles/particles.h"
#include "scalar/scalar.h"

namespace cosmo
{

ICsData cosmo_get_ICsData();

real_t cosmo_power_spectrum(real_t k, ICsData *icd);
void set_gaussian_random_field(real_t *field, Fourier *fourier, ICsData *icd);

void ICs_set_dust(map_t & bssn_fields, map_t & static_field,
  Fourier *fourier, IOData *iod, FRW<real_t> *frw);

void ICs_set_particle(Particles * particles, map_t & bssn_fields,
  Fourier *fourier, IOData *iod);

void ICs_set_vacuum(map_t & bssn_fields, IOData *iod);

void init_ray_vector(std::vector<RayTrace<real_t, idx_t> *> * rays,
  idx_t n_rays);

void ICs_set_scalar_wave(Scalar * scalar);
void ICs_set_scalar_inflation(map_t & bssn_fields, map_t & scalar_fields);

void set_stability_test_ICs( map_t & bssn_fields, map_t & static_field);
void set_linear_wave_ICs(map_t & bssn_fields);

} // namespace cosmo

#endif
