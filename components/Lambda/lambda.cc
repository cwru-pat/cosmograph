#include "lambda.h"

namespace cosmo
{

Lambda::Lambda()
{
  lambda = 0.0;
}

Lambda::~Lambda()
{
  /* xxx */
}

void Lambda::setLambda(real_t lambda_in)
{
  lambda = lambda_in;
}

void Lambda::addBSSNSource(BSSN *bssn)
{
  arr_t & DIFFr_a = *bssn->fields["DIFFr_a"];
  arr_t & DIFFS_a = *bssn->fields["DIFFS_a"];

  idx_t i, j, k;
# pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    DIFFr_a[idx] += lambda;
    DIFFS_a[idx] += -3.0*lambda;

    // no S_i or STF_ij contributions.
  }

}


} /* namespace cosmo */
