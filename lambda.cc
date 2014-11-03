#include "lambda.h"

namespace cosmo
{

Lambda::Lambda(real_t L)
{
  lambda = L;
}

Lambda::~Lambda()
{
  /* xxx */
}

void Lambda::addBSSNSrc(std::map <std::string, real_t *> & bssn_fields)
{

  LOOP3(i,j,k)
  {
    idx_t idx = INDEX(i,j,k);

    bssn_fields["r_a"][idx] += lambda;
    bssn_fields["S_a"][idx] += -3.0*lambda;

    // no S_i or STF_ij contributions.

  }

}


} /* namespace cosmo */
