#ifndef COSMO_ICS_DATA
#define COSMO_ICS_DATA

namespace cosmo
{

typedef struct {

  // spectrum data
  real_t peak_k;
  real_t peak_amplitude;
  idx_t ic_spec_cut;

} ICsData;

} /* namespace cosmo */

#endif
