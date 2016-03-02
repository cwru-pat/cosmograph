#ifndef COSMO_MAIN_H
#define COSMO_MAIN_H

#include "cosmo.h"
#include "globals.h"

namespace cosmo
{

// main evolution functions

void init_ray_vector(std::vector<RayTrace<real_t, idx_t> *> * rays, idx_t n_rays);

void call_io_routines(BSSN * bssnSim, Static * staticSim,
                             IOData *iodata, idx_t step, idx_t steps,
                             Fourier *fourier, FRW<real_t> *frw,
                             std::vector<RayTrace<real_t, idx_t> *> * rays);

} /* namespace cosmo */

#endif
