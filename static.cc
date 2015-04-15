#include "static.h"

namespace cosmo
{

Static::Static()
{
  GEN1_ARRAY_ALLOC(D);
  GEN1_ARRAY_ADDMAP(D);
}

Static::~Static()
{
  GEN1_ARRAY_DELETE(D);
}

void Static::addBSSNSrc(std::map <std::string, real_t *> & bssn_fields)
{
  idx_t i, j, k;
  #pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    // \rho = \gamma^{-1/2}*D
    //      = exp(-6*\phi)*D
    bssn_fields["r_a"][idx] += exp(-6.0*bssn_fields["phi_a"][idx])*D_a[idx];
  }

}

void Static::init()
{
  // initialize values
  idx_t i, j, k;
  LOOP3(i,j,k)
  {
    idx_t idx = INDEX(i,j,k);

    // empty universe
    D_a[idx] = 0.0;
  }

}

} /* namespace cosmo */
