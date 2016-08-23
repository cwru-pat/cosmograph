/** @file scalar_ic.h
 * @brief Functions to set initial conditions for the Scalar class.
 * Functions should be made callable via a config setting in ScalarSim
 * class.
 */

#ifndef COSMO_SCALAR_ICS
#define COSMO_SCALAR_ICS

#include "../bssn/bssn.h"
#include "scalar.h"
#include "../../IO/io.h"

namespace cosmo
{

void scalar_ic_set_wave(BSSN * bssn, Scalar * scalar);
void scalar_ic_set_Lambda(BSSN * bssn, Scalar * scalar);
void scalar_ic_set_semianalytic_test(BSSN * bssn, Scalar * scalar, IOData * iodata);

#if USE_MULTIGRID
void scalar_ic_set_multigrid(BSSN * bssn, Scalar * scalar, IOData * iodata);
#endif

}

#endif
