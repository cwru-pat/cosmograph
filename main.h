#ifndef COSMO_MAIN_H
#define COSMO_MAIN_H

#include "cosmo_includes.h"
#include "cosmo_types.h"
#include "globals.h"

#include "utils/Fourier.h"
#include "utils/reference_frw.h"

#include "ICs.h"
#include "io.h"
#include "bssn.h"
#include "static.h"

namespace cosmo
{

void call_io_routines(BSSN * bssnSim, Static * staticSim,
                      IOData *iodata, idx_t step, idx_t steps,
                      Fourier *fourier, FRW<real_t> *frw,
                      std::vector<RayTrace<real_t, idx_t> *> * rays);

} /* namespace cosmo */

#endif
