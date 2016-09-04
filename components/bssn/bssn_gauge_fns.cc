#include "bssn_gauge_fns.h"
#include "../../cosmo_types.h"
#include "../../cosmo_globals.h"
#include <map>
#include <iostream>
#include <cmath>

namespace cosmo
{

/**
 * @brief Get a function pointer for a lapse function
 * associated with a string
 * 
 * @param gauge string describing lapse gauge function
 * @return pointer to gauge function
 */
bssn_func_t bssn_gauge_get_lapse_fn(std::string gauge)
{
  if(gauge == "") return bssn_gauge_static;

  // list of possible lapse functions
  std::map<std::string, bssn_func_t> lapse_functions;
  lapse_functions["static"] = bssn_gauge_static;
  lapse_functions["harmonic"] = bssn_gauge_alpha_harmonic;
  lapse_functions["1pluslog"] = bssn_gauge_alpha_1pluslog;
  lapse_functions["dampedwave"] = bssn_gauge_alpha_dampedwave;
  // TODO: gauges need more work:
  lapse_functions["anharmonic"] = bssn_gauge_alpha_anharmonic;
  lapse_functions["conformalsync"] = bssn_gauge_alpha_conformalsync;
  lapse_functions["AwA_gaugewave"] = bssn_gauge_alpha_AwA_gaugewave;

  if( lapse_functions.find(gauge) == lapse_functions.end() )
  {
    // gauge doesn't exist
    std::cerr << "Lapse gauge `" << gauge << "` not found!";
    throw -1;
  }

  return lapse_functions[gauge];
}

/**
 * @brief Get a function pointer for a shift function
 * associated with a string
 * 
 * @param gauge string describing shift gauge function
 * @return pointer to gauge function
 */
bssn_func_t bssn_gauge_get_shift_fn(std::string gauge, int i)
{
  if(gauge == "") return bssn_gauge_static;
  // No map here; only a couple possible shift functions for now


# if USE_BSSN_SHIFT

  if(gauge == "gammadriver")
  {
    if(!USE_GAMMA_DRIVER)
    {
      std::cerr << "Code must be compiled with gamma driver enabled to use the gamma driver gauge.";
      throw -1;
    }
#   if USE_GAMMA_DRIVER
      if(i == 1) return bssn_gauge_beta1_gammadriver;
      if(i == 2) return bssn_gauge_beta2_gammadriver;
      if(i == 3) return bssn_gauge_beta3_gammadriver;
#   endif
  }

  if(gauge == "dampedwave")
  {
    if(i == 1) return bssn_gauge_beta1_dampedwave;
    if(i == 2) return bssn_gauge_beta2_dampedwave;
    if(i == 3) return bssn_gauge_beta3_dampedwave;
  }

  // no gauge found?
  std::cerr << "Shift gauge `" << gauge << "` not found!";
  throw -1;

# else // from if USE_BSSN_SHIFT
  std::cerr << "Unable to use shift, code was not compiled with shift option.";
  throw -1;
# endif

}

/**
 * @brief Don't evolve anything
 * @return 0
 */
real_t bssn_gauge_static(BSSNData *bd)
{
  return 0.0;
}

/** Lapse gauge functions **/

real_t bssn_gauge_alpha_harmonic(BSSNData *bd)
{
  // TODO: Generalize K0
  return -1.0*pw2(bd->alpha)*( bd->K - bd->K0 );
}

real_t bssn_gauge_alpha_anharmonic(BSSNData *bd)
{
  //<<< TODO: generalize K "offset" in harmonic gauge.
  // Ref. showing presence of offset:
  // http://relativity.livingreviews.org/open?pubNo=lrr-2012-9&amp;page=articlesu7.html
  // for FRW (+ perturbation) sims, having no offset leads to lapse blowing up?
  real_t K_FRW_0 = -3.0;
  return 1.0*pw2(bd->alpha)*( bd->K - K_FRW_0 );
}

real_t bssn_gauge_alpha_1pluslog(BSSNData *bd)
{
  return -2.0*bd->alpha*( bd->K  )*GD_C
      + bd->beta1*bd->d1a + bd->beta2*bd->d2a + bd->beta3*bd->d3a;
}

real_t bssn_gauge_alpha_dampedwave(BSSNData *bd)
{
  return pw2(bd->alpha) * (DW_MU_L * (12.0 * bd->phi * DW_P - std::log(bd->alpha)) - bd->K)
    + bd->beta1 * bd->d1a + bd->beta2 * bd->d2a + bd->beta3 * bd->d3a;
}

real_t bssn_gauge_alpha_conformalsync(BSSNData *bd)
{
  return -1.0/3.0*bd->alpha*bd->K_FRW;
}

real_t bssn_gauge_alpha_AwA_gaugewave(BSSNData *bd)
{
  return -1.0*pw2(bd->alpha)*bd->DIFFK;
}

/** Shift gauge functions **/

#if USE_BSSN_SHIFT
#if USE_GAMMA_DRIVER
real_t bssn_gauge_beta1_gammadriver(BSSNData *bd)
{
  return bd->auxB1;
}
real_t bssn_gauge_beta2_gammadriver(BSSNData *bd)
{
  return bd->auxB2;
}
real_t bssn_gauge_beta3_gammadriver(BSSNData *bd)
{
  return bd->auxB3;
}
#endif

real_t bssn_gauge_beta1_dampedwave(BSSNData *bd)
{
  return bd->beta1*bd->d1beta1 + bd->beta2*bd->d2beta1 + bd->beta3*bd->d3beta1
    - DW_MU_S*bd->alpha*bd->beta1
    + bd->alpha * ( -DW_MU_L*(12.0*bd->phi*DW_P - std::log(bd->alpha))*bd->beta1
        + std::exp(-4.0*bd->phi)*(
          - bd->gammai11*bd->d1a - bd->gammai12*bd->d2a - bd->gammai13*bd->d3a
          + bd->alpha*(bd->Gamma1
              -2.0*(bd->gammai11*bd->d1phi + bd->gammai12*bd->d2phi + bd->gammai13*bd->d3phi)
            )
        )
      );
}
real_t bssn_gauge_beta2_dampedwave(BSSNData *bd)
{
  return bd->beta1*bd->d1beta2 + bd->beta2*bd->d2beta2 + bd->beta3*bd->d3beta2
    - DW_MU_S*bd->alpha*bd->beta2
    + bd->alpha * ( -DW_MU_L*(12.0*bd->phi*DW_P - std::log(bd->alpha))*bd->beta2
        + std::exp(-4.0*bd->phi)*(
          - bd->gammai21*bd->d1a - bd->gammai22*bd->d2a - bd->gammai23*bd->d3a
          + bd->alpha*(bd->Gamma3
              -2.0*(bd->gammai21*bd->d1phi + bd->gammai22*bd->d2phi + bd->gammai23*bd->d3phi)
            )
        )
      );
}
real_t bssn_gauge_beta3_dampedwave(BSSNData *bd)
{
  return bd->beta1*bd->d1beta3 + bd->beta2*bd->d2beta3 + bd->beta3*bd->d3beta3
    - DW_MU_S*bd->alpha*bd->beta3
    + bd->alpha * ( -DW_MU_L*(12.0*bd->phi*DW_P - std::log(bd->alpha))*bd->beta3
        + std::exp(-4.0*bd->phi)*(
          - bd->gammai31*bd->d1a - bd->gammai32*bd->d2a - bd->gammai33*bd->d3a
          + bd->alpha*(bd->Gamma3
              -2.0*(bd->gammai31*bd->d1phi + bd->gammai32*bd->d2phi + bd->gammai33*bd->d3phi)
            )
        )
      );
}
#endif

} // namespace cosmo
