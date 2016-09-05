/** @file bssn_gauge_fns.h
 * @brief Functions to determine gauge evolution for the BSSN class.
 * Functions are determined via a config setting in the CosmoSim class.
 */

#ifndef COSMO_BSSN_GAUGE_FNS
#define COSMO_BSSN_GAUGE_FNS

#include "bssn_data.h"
#include <string>

namespace cosmo
{

typedef real_t (*bssn_func_t)(BSSNData *bd); /**< (gauge) function pointer type */

bssn_func_t bssn_gauge_get_lapse_fn(std::string gauge);
bssn_func_t bssn_gauge_get_shift_fn(std::string gauge, int i);

real_t bssn_gauge_static(BSSNData *bd);

real_t bssn_gauge_alpha_harmonic(BSSNData *bd);
real_t bssn_gauge_alpha_anharmonic(BSSNData *bd);
real_t bssn_gauge_alpha_1pluslog(BSSNData *bd);
real_t bssn_gauge_alpha_dampedwave(BSSNData *bd);
real_t bssn_gauge_alpha_conformalsync(BSSNData *bd);

real_t bssn_gauge_alpha_AwA_gaugewave(BSSNData *bd);

#if USE_BSSN_SHIFT
#if USE_GAMMA_DRIVER
real_t bssn_gauge_beta1_gammadriver(BSSNData *bd);
real_t bssn_gauge_beta2_gammadriver(BSSNData *bd);
real_t bssn_gauge_beta3_gammadriver(BSSNData *bd);
#endif
real_t bssn_gauge_beta1_dampedwave(BSSNData *bd);
real_t bssn_gauge_beta2_dampedwave(BSSNData *bd);
real_t bssn_gauge_beta3_dampedwave(BSSNData *bd);
#endif

}

#endif
