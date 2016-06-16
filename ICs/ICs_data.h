#ifndef COSMO_ICS_DATA
#define COSMO_ICS_DATA

namespace cosmo
{

typedef struct {

  // spectrum data
  real_t peak_k;
  real_t peak_amplitude;
  idx_t ic_spec_cut;

  // matter params
  real_t rho_K_matter;
  real_t rho_K_lambda;

  // constraint violating term amplitude
  real_t viol_amp;

} ICsData;

} /* namespace cosmo */

#endif
