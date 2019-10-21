/** @file static_ic.h
 * @brief Functions to set initial conditions for the Static class.
 * Functions should be made callable via a config setting in the StaticSim
 * class.
 */

#ifndef COSMO_STATIC_ICS
#define COSMO_STATIC_ICS

#include "../bssn/bssn.h"
#include "../Lambda/lambda.h"
#include "../../IO/IOData.h"
#include "../../utils/Fourier.h"
#include "static.h"

namespace cosmo
{

void static_ic_set_random( BSSN * bssn, Static * stat, Lambda * lambda,
  Fourier * fourier, IOData * iodata );

void static_ic_set_sinusoid_3d( BSSN * bssn, Static * stat, Lambda * lambda,
  Fourier * fourier, IOData * iodata );
 
void static_ic_set_sinusoid( BSSN * bssn, Static * stat, Lambda * lambda,
  Fourier * fourier, IOData * iodata );

void static_ic_set_sphere( BSSN * bssn, Static * stat, IOData * iodata );

void static_ic_set_semianalytic( BSSN * bssn, Static * stat, Lambda * lambda,
  Fourier * fourier, IOData * iodata );
  
}

#endif
