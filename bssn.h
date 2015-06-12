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

  real_t der(real_t field_adj[3][3][3], int d);
  real_t dder(real_t field_adj[3][3][3], int d1, int d2);
  real_t der_ext(real_t field_adj[3][3][3], real_t field_adj_ext[3][2], int d);
  real_t dder_ext(real_t field_adj[3][3][3], real_t field_adj_ext[3][2], int d1, int d2);
  
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
  real_t hamiltonianConstraintMag(idx_t i, idx_t j, idx_t k);
  real_t momentumConstraintCalc(BSSNData *paq, idx_t i);
  real_t momentumConstraintMag(BSSNData *paq, idx_t i);

  void clearSrc();

  void step(BSSNData *paq);
  void init();

};

}

#endif
