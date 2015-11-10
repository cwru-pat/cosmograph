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
  BSSN_APPLY_TO_GEN1_EXTRAS(GEN1_ARRAY_CREATE)

public:
  std::map <std::string, real_t *> fields;

  BSSN();
  ~BSSN();

  void init();

  /* RK integrator functions */
    void clearSrc();
    void regSwap_c_a();

    /* functions to outline / perform RK integration */
    void step(BSSNData *paq);
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

  /* functions to pre-calculate quantities before RK steps */
    void set_gammai_values(idx_t i, idx_t j, idx_t k);
    void set_detgamma(idx_t i, idx_t j, idx_t k);

  /* calculating quantities during an RK step */
    void set_paq_values(idx_t i, idx_t j, idx_t k, BSSNData *paq);

    /* set current local field values */
      void set_local_vals(BSSNData *paq);

    /* Calculate independent quantities */
      void calculate_Acont(BSSNData *paq);
      void calculate_dgamma(BSSNData *paq);
      void calculate_ddgamma(BSSNData *paq);
      void calculate_dgammai(BSSNData *paq);
      void calculate_dalpha_dphi(BSSNData *paq);
      void calculate_dK(BSSNData *paq);
      #if Z4c_DAMPING > 0
        void calculate_dtheta(BSSNData *paq);
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

    real_t ev_alpha(BSSNData *paq);

    #if Z4c_DAMPING > 0
      real_t ev_Z1(BSSNData *paq);
      real_t ev_Z2(BSSNData *paq);
      real_t ev_Z3(BSSNData *paq);
      real_t ev_theta(BSSNData *paq);
    #endif

  /* constraint violation calculations */
    real_t hamiltonianConstraintCalc(BSSNData *paq);
    real_t hamiltonianConstraintScale(BSSNData *paq);
    real_t hamiltonianConstraintMagMean();
    real_t hamiltonianConstraintMagStDev(real_t mean);
    real_t hamiltonianConstraintMagMax();

    real_t hamiltonianConstraintMax();

    real_t momentumConstraintCalc(BSSNData *paq, idx_t d);
    real_t momentumConstraintScale(BSSNData *paq, idx_t d);
    real_t momentumConstraintMagMean();
    real_t momentumConstraintMagStDev(real_t mean);
    real_t momentumConstraintMagMax();

    real_t metricConstraintTotalMag();

};

}

#endif
