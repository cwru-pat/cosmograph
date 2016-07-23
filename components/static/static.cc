#include "static.h"
#include "../../cosmo_includes.h"
#include "../../cosmo_globals.h"

namespace cosmo
{

Static::Static()
{
  GEN1_ARRAY_ALLOC(DIFFD);
  GEN1_ARRAY_ADDMAP(DIFFD);
}

Static::~Static()
{
  GEN1_ARRAY_DELETE(DIFFD);
}

void Static::addBSSNSrc(map_t & bssn_fields, FRW<real_t> *frw)
{
  idx_t i=0, j=0, k=0;
  arr_t & DIFFalpha_a = *bssn_fields["DIFFalpha_a"];
  arr_t & DIFFr_a = *bssn_fields["DIFFr_a"];
  arr_t & DIFFphi_a = *bssn_fields["DIFFphi_a"];

  real_t phi_FRW = frw->get_phi();
  real_t rho_FRW = frw->get_rho();

  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    // \Delta \rho = \rho - \rho_FRW = ...
    DIFFr_a[idx] += exp(-6.0*(DIFFphi_a[idx] + phi_FRW))*( DIFFD_a[idx] / (1.0 + DIFFalpha_a[idx]) )
      + rho_FRW*expm1(-6.0*DIFFphi_a[idx]);
  }

}

void Static::init()
{
  // initialize values
  idx_t i, j, k;
  LOOP3(i,j,k)
  {
    idx_t idx = INDEX(i,j,k);
    DIFFD_a[idx] = 0.0;
  }

}

} /* namespace cosmo */
