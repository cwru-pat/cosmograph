#ifndef COSMO_ICS_DATA
#define COSMO_ICS_DATA

namespace cosmo
{

typedef struct {

  // spectrum data
  real_t peak_k;
  real_t peak_amplitude;

  // plan for taking FFTs
  fftw_plan plan;

} ICsData;

} /* namespace cosmo */

#endif
