/** @file bssn_ic.h
 * @brief Functions to set initial conditions for the BSSN class.
 * Functions are called via a config setting in VacuumSim
 * class.
 */

#ifndef COSMO_SCALAR_ICS
#define COSMO_SCALAR_ICS

#include "bssn.h"

namespace cosmo
{

void bssn_ic_awa_stability(BSSN * bssn);
void bssn_ic_awa_linear_wave(BSSN * bssn);
void bssn_ic_awa_linear_wave_desitter(BSSN * bssn);

}

#endif
