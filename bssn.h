#ifndef COSMO_BSSN
#define COSMO_BSSN

#include "cosmo.h"
#include "globals.h"

#include "utils/reference_frw.h"
#include "bssn_data.h"
#include "bssn_macros.h"

namespace cosmo
{

/** BSSN class **/
class BSSN
{
  /* arrays for storing fields */
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_CREATE)
  BSSN_APPLY_TO_GEN1_EXTRAS(GEN1_ARRAY_CREATE)

public:
  std::map <std::string, arr_t *> fields;

  // Standard FRW spacetime integrator - for a
  // reference metric
  FRW<real_t> * frw;

  BSSN();
  ~BSSN();

  void init();

  /* RK integrator functions */
    void regSwap_c_a();

    void stepInit();
    void WedgeStep();
    void WedgeK1Calc(idx_t i_K1);
    void WedgeK2Calc(idx_t i_K2);
    void WedgeK3Calc(idx_t i_K3);
    void WedgeTailCalc(idx_t i_tail);

    void DrawWedgeSlice(idx_t i_p, idx_t i_K1, idx_t i_K2, idx_t i_K3, idx_t i_tail);

  /* calculating quantities during an RK step */
    void set_paq_values(idx_t i, idx_t j, idx_t k, BSSNData *paq);

    /* set current local field values */
      void set_local_vals(BSSNData *paq);
      void set_gammai_values(idx_t i, idx_t j, idx_t k, BSSNData *paq);
      void set_DIFFgamma_Aij_norm();
      void set_DIFFgamma_Aij_norm_pt(idx_t i, idx_t j, idx_t k);

    /* Calculate quantities only dependent on FRW soln in paq*/
      void calculate_Acont(BSSNData *paq);
      void calculate_dgamma(BSSNData *paq);
      void calculate_ddgamma(BSSNData *paq);
      void calculate_dalpha_dphi(BSSNData *paq);
      void calculate_dK(BSSNData *paq);
      #if USE_Z4c_DAMPING
        void calculate_dtheta(BSSNData *paq);
      #endif
      #if USE_BSSN_SHIFT
        void calculate_dbeta(BSSNData *paq);
      #endif

    /* Calculate "dependent" quantities (depend on previously calc'd vals) */
      void calculate_conformal_christoffels(BSSNData *paq);

    /* Calculate doubly-"dependent" quantities (depend on previously calc'd vals) */
      void calculateDDphi(BSSNData *paq);
      void calculateRicciTF(BSSNData *paq);
      void calculateDDalphaTF(BSSNData *paq);

    /* (optional) Calculations of additional quantities */
      void set_KillingDelta(idx_t i, idx_t j, idx_t k, BSSNData *paq);
      void set_full_metric(BSSNData *paq);
      void set_full_metric_der(BSSNData *paq);

  /* Evolution functions */
    real_t ev_DIFFgamma11(BSSNData *paq);
    real_t ev_DIFFgamma12(BSSNData *paq);
    real_t ev_DIFFgamma13(BSSNData *paq);
    real_t ev_DIFFgamma22(BSSNData *paq);
    real_t ev_DIFFgamma23(BSSNData *paq);
    real_t ev_DIFFgamma33(BSSNData *paq);
    real_t ev_A11(BSSNData *paq);
    real_t ev_A12(BSSNData *paq);
    real_t ev_A13(BSSNData *paq);
    real_t ev_A22(BSSNData *paq);
    real_t ev_A23(BSSNData *paq);
    real_t ev_A33(BSSNData *paq);
    real_t ev_DIFFK(BSSNData *paq);
    real_t ev_DIFFphi(BSSNData *paq);
    real_t ev_Gamma1(BSSNData *paq);
    real_t ev_Gamma2(BSSNData *paq);
    real_t ev_Gamma3(BSSNData *paq);

    real_t ev_DIFFalpha(BSSNData *paq);
    real_t ev_DIFFdustrho(BSSNData *paq);
    real_t ev_eta(BSSNData *paq);

    #if USE_Z4c_DAMPING
      real_t ev_theta(BSSNData *paq);
    #endif

    #if USE_BSSN_SHIFT
      real_t ev_beta1(BSSNData *paq);
      real_t ev_beta2(BSSNData *paq);
      real_t ev_beta3(BSSNData *paq);
    #endif

  /* constraint violation calculations */
    void setHamiltonianConstraintCalcs(real_t H_values[7], bool reset_paq);
    real_t hamiltonianConstraintCalc(idx_t idx);
    real_t hamiltonianConstraintCalc(BSSNData *paq);
    real_t hamiltonianConstraintScale(idx_t idx);

    void setMomentumConstraintCalcs(real_t M_values[7]);
    real_t momentumConstraintCalc(BSSNData *paq, idx_t d);
    real_t momentumConstraintScale(BSSNData *paq, idx_t d);

    real_t metricConstraintTotalMag();

};

}

#endif
