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

public:
  std::map <std::string, real_t *> fields;

  /* local values */
  cosmo::BSSNData paq;

  BSSN();
  ~BSSN();

  void regSwap_c_a();

  void stepInit();
  void K1Calc();
  void K1CalcPt(idx_t i, idx_t j, idx_t k);
  void K2Calc();
  void K2CalcPt(idx_t i, idx_t j, idx_t k);
  void K3Calc();
  void K3CalcPt(idx_t i, idx_t j, idx_t k);
  void K4Calc();
  void K4CalcPt(idx_t i, idx_t j, idx_t k);
  void stepTerm();

  real_t der(real_t field_adj[3][3][3], int d);
  real_t dder(real_t field_adj[3][3][3], int d1, int d2);
  
  /* set current local field values */
  void set_local_vals(cosmo::BSSNData *paq);
  void calculate_Acont(cosmo::BSSNData *paq);

  /* Calculate metric derivatives */
  void calculate_dgamma(cosmo::BSSNData *paq);
  void calculate_ddgamma(cosmo::BSSNData *paq);

  /* Calculate metric derivatives */
  void calculate_dgammai(cosmo::BSSNData *paq);
  void calculate_christoffels(cosmo::BSSNData *paq);
  void calculate_dalpha_dphi(cosmo::BSSNData *paq);
  void calculate_dbeta(cosmo::BSSNData *paq);

  /* Calculate trace-free ricci tensor components */
  void calculateRicciTF(cosmo::BSSNData *paq);
  void calculateDDphi(cosmo::BSSNData *paq);
  void calculateDDalphaTF(cosmo::BSSNData *paq);

  /* Evolution functions */
  real_t ev_gamma11(cosmo::BSSNData *paq);
  real_t ev_gamma12(cosmo::BSSNData *paq);
  real_t ev_gamma13(cosmo::BSSNData *paq);
  real_t ev_gamma22(cosmo::BSSNData *paq);
  real_t ev_gamma23(cosmo::BSSNData *paq);
  real_t ev_gamma33(cosmo::BSSNData *paq);
  real_t ev_A11(cosmo::BSSNData *paq);
  real_t ev_A12(cosmo::BSSNData *paq);
  real_t ev_A13(cosmo::BSSNData *paq);
  real_t ev_A22(cosmo::BSSNData *paq);
  real_t ev_A23(cosmo::BSSNData *paq);
  real_t ev_A33(cosmo::BSSNData *paq);
  real_t ev_K(cosmo::BSSNData *paq);
  real_t ev_phi(cosmo::BSSNData *paq);
  real_t ev_Gamma1(cosmo::BSSNData *paq);
  real_t ev_Gamma2(cosmo::BSSNData *paq);
  real_t ev_Gamma3(cosmo::BSSNData *paq);
  // static gauge for now:
  real_t ev_alpha(cosmo::BSSNData *paq) { return 0; }
  real_t ev_beta1(cosmo::BSSNData *paq) { return 0; }
  real_t ev_beta2(cosmo::BSSNData *paq) { return 0; }
  real_t ev_beta3(cosmo::BSSNData *paq) { return 0; }

  // calculate needed quantities (need the inverse metric set everywhere first)
  void set_paq_values(idx_t i, idx_t j, idx_t k, cosmo::BSSNData *paq);

  void step();
  void init();

};

}

#endif
