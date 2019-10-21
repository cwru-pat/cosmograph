#ifndef COSMO_DUST_ICS
#define COSMO_DUST_ICS

#include "../bssn/bssn.h"
#include "dust.h"
#include "../Lambda/lambda.h"
#include "../../utils/Fourier.h"
#include "../../IO/IOData.h"

namespace cosmo
{

void dust_ic_set_random( BSSN * bssn, Dust * dust, Lambda * lambda,
  Fourier * fourier, IOData * iodata );

}

#endif
