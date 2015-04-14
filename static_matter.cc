#include "static_matter.h"

namespace cosmo
{

Static::Static()
{
  GEN1_ARRAY_ALLOC(r);
  GEN1_ARRAY_ADDMAP(r);
}

Static::~Static()
{
  GEN1_ARRAY_DELETE(r);
}

void Static::addBSSNSrc(std::map <std::string, real_t *> & bssn_fields)
{
  idx_t i, j, k;
  #pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = INDEX(i,j,k);
    bssn_fields["r_a"][idx] += r_a[idx];
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
    r_a[idx] = 0.0;
  }

}

} /* namespace cosmo */
