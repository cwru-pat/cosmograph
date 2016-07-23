/** @file static_ic.h
 * @brief Functions to set initial conditions for the Static (dust) class.
 * Functions should be made callable via a config setting in the DustSim
 * class.
 */

#ifndef COSMO_STATIC_ICS
#define COSMO_STATIC_ICS

#include "../bssn/bssn.h"
#include "../../IO/IOData.h"
#include "../../utils/Fourier.h"
#include "static.h"

namespace cosmo
{

void dust_ic_set_random(BSSN * bssn, Static * dust, Fourier * fourier,
  IOData * iodata);

}

#endif
