#ifndef COSMO_BSSN
#define COSMO_BSSN

#include "cosmo.h"
#include "globals.h"

#include "bssn_data.h"
#include "bssn_macros.h"

namespace cosmo
{

/** BSSN class **/
class BSSN
{
  /* arrays for storing fields */
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_CREATE)
  BSSN_APPLY_TO_SOURCES(GEN1_ARRAY_CREATE)
  // additional arrays for output of calculated quantities
  GEN1_ARRAY_CREATE(ricci);
  GEN1_ARRAY_CREATE(AijAij);

  GEN1_ARRAY_CREATE(gammai11);
  GEN1_ARRAY_CREATE(gammai12);
  GEN1_ARRAY_CREATE(gammai13);
  GEN1_ARRAY_CREATE(gammai22);
  GEN1_ARRAY_CREATE(gammai23);
  GEN1_ARRAY_CREATE(gammai33);

public:
  std::map <std::string, real_t *> fields;

  BSSN();
  ~BSSN();

  void regSwap_c_a();

  void stepInit();
  void K1Calc(BSSNData *paq);
  void K1CalcPt(idx_t i, idx_t j, idx_t k, BSSNData *paq);
  void K2Calc(BSSNData *paq);
  void K2CalcPt(idx_t i, idx_t j, idx_t k, BSSNData *paq);
  void K3Calc(BSSNData *paq);
  void K3CalcPt(idx_t i, idx_t j, idx_t k, BSSNData *paq);
  void K4Calc(BSSNData *paq);
  void K4CalcPt(idx_t i, idx_t j, idx_t k, BSSNData *paq);
  void stepTerm();

// DEBUGGING FUNCTIONS
void set_gammai_values(idx_t i, idx_t j, idx_t k);
void calculate_dK(BSSNData *paq);
void set_detgamma(idx_t i, idx_t j, idx_t k);

  /* set current local field values */
  void set_source_vals(BSSNData *paq);
  void set_local_vals(BSSNData *paq);
  void calculate_Acont(BSSNData *paq);

  /* Calculate metric derivatives */
  void calculate_dgamma(BSSNData *paq);
  void calculate_ddgamma(BSSNData *paq);

  /* Calculate metric derivatives */
  void calculate_dgammai(BSSNData *paq);
  void calculate_christoffels(BSSNData *paq);
  void calculate_dphi(BSSNData *paq);

  /* Calculate trace-free ricci tensor components */
  void calculateRicciTF(BSSNData *paq);
  void calculateDDphi(BSSNData *paq);

  /* Evolution functions */
  real_t ev_gamma11(BSSNData *paq);
  real_t ev_gamma12(BSSNData *paq);
  real_t ev_gamma13(BSSNData *paq);
  real_t ev_gamma22(BSSNData *paq);
  real_t ev_gamma23(BSSNData *paq);
  real_t ev_gamma33(BSSNData *paq);
  real_t ev_A11(BSSNData *paq);
  real_t ev_A12(BSSNData *paq);
  real_t ev_A13(BSSNData *paq);
  real_t ev_A22(BSSNData *paq);
  real_t ev_A23(BSSNData *paq);
  real_t ev_A33(BSSNData *paq);
  real_t ev_K(BSSNData *paq);
  real_t ev_phi(BSSNData *paq);
  real_t ev_Gamma1(BSSNData *paq);
  real_t ev_Gamma2(BSSNData *paq);
  real_t ev_Gamma3(BSSNData *paq);

  // calculate needed quantities (need the inverse metric set everywhere first)
  void set_paq_values(idx_t i, idx_t j, idx_t k, BSSNData *paq);
  void set_full_metric(BSSNData *paq);
  void set_full_metric_der(BSSNData *paq);

  real_t hamiltonianConstraintCalc(BSSNData *paq);
  real_t hamiltonianConstraintScale(BSSNData *paq);
  real_t hamiltonianConstraintMagMean();
  real_t hamiltonianConstraintMagStDev(real_t mean);
  real_t hamiltonianConstraintMagMax();

  real_t momentumConstraintCalc(BSSNData *paq, idx_t d);
  real_t momentumConstraintScale(BSSNData *paq, idx_t d);
  real_t momentumConstraintMagMean();
  real_t momentumConstraintMagStDev(real_t mean);
  real_t momentumConstraintMagMax();

  real_t metricConstraintTotalMag();

  void clearSrc();

  void step(BSSNData *paq);
  void init();

};

}

#endif
