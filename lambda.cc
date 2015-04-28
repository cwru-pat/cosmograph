#include "lambda.h"

namespace cosmo
{

Lambda::Lambda()
{
  ICsData icd = cosmo_get_ICsData();
  lambda = icd.rho_K_lambda;
}

Lambda::~Lambda()
{
  /* xxx */
}

void Lambda::addBSSNSrc(std::map <std::string, real_t *> & bssn_fields)
{
  idx_t i, j, k;

  real_t * const r_a = bssn_fields["r_a"];
  real_t * const S_a = bssn_fields["S_a"];
  #pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    r_a[idx] += lambda;
    S_a[idx] += -3.0*lambda;

    // no S_i or STF_ij contributions.
  }

}


} /* namespace cosmo */
