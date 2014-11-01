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

    real_t g11 = bssn_fields["gamma11_f"][idx];
    real_t g12 = bssn_fields["gamma12_f"][idx];
    real_t g13 = bssn_fields["gamma13_f"][idx];
    real_t g22 = bssn_fields["gamma22_f"][idx];
    real_t g23 = bssn_fields["gamma23_f"][idx];
    real_t g33 = bssn_fields["gamma33_f"][idx];

    bssn_fields["r_a"][idx] += lambda;

    bssn_fields["S11_a"][idx] += -g11*lambda;
    bssn_fields["S12_a"][idx] += -g12*lambda;
    bssn_fields["S13_a"][idx] += -g13*lambda;
    bssn_fields["S22_a"][idx] += -g22*lambda;
    bssn_fields["S23_a"][idx] += -g23*lambda;
    bssn_fields["S33_a"][idx] += -g33*lambda;

    bssn_fields["S_a"][idx] += -3.0*lambda;
  }

}


} /* namespace cosmo */
