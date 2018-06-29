#ifndef COSMO_BSSN
#define COSMO_BSSN

#include "../../cosmo_includes.h"
#include "../../cosmo_types.h"
#include "bssn_data.h"
#include "bssn_macros.h"
#include "BSSNGaugeHandler.h"
#include "../../utils/Array.h"
#include "../../utils/FRW.h"
#include "../../utils/ConfigParser.h"

#if USE_COSMOTRACE
#include "../cosmotrace/raytrace.h"
#endif

namespace cosmo
{

/**
 * @brief BSSN Class: evolves BSSN metric fields, computes derived quantities
 */
class BSSN
{
  /* arrays for storing fields */
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_CREATE)
  BSSN_APPLY_TO_SOURCES(GEN1_ARRAY_CREATE)
  BSSN_APPLY_TO_GEN1_EXTRAS(GEN1_ARRAY_CREATE)

  BSSNGaugeHandler * gaugeHandler;

  real_t KO_damping_coefficient; ///< KO_dissipation coefficient amplitude (default 0)
  real_t k_damping_amp; ///< Additional "K"-damping (default 0)
  real_t a_adj_amp; ///< A-adjusted term amplitude (default 0)
  real_t gd_eta; ///< Gamma driver "eta" parameter
  int normalize_metric; ///< Normalize A_ij and \gamma_ij? Default: 1 (true)

public:
  map_t fields; ///< Public map from names to internal arrays

  // Standard FRW spacetime integrator - for a reference metric
  FRW<real_t> * frw; ///< FRW reference metric instance
  // Average K, rho for reference
  real_t K_avg, rho_avg;
  real_t K_min;

  real_t avg_vol;

  BSSN(ConfigParser * config);
  ~BSSN();

  void init();

  void setDt(real_t dt);
  void setKODampingCoefficient(real_t KO_damping_coefficient_in);

  /* RK integrator functions */
    void stepInit();
    void RKEvolve();
    void RKEvolvePt(idx_t i, idx_t j, idx_t k, BSSNData * bd);
    void K1Finalize();
    void K2Finalize();
    void K3Finalize();
    void K4Finalize();
    void clearSrc();
    void step();

  /* calculating quantities during an RK step */
    void set_bd_values(idx_t i, idx_t j, idx_t k, BSSNData *bd);

    /* set current local field values */
      void set_local_vals(BSSNData *bd);
      void set_gammai_values(idx_t i, idx_t j, idx_t k, BSSNData *bd);
      void set_DIFFgamma_Aij_norm();

    /* Calculate quantities only dependent on FRW soln in bd*/
      void calculate_Acont(BSSNData *bd);
      void calculate_dgamma(BSSNData *bd);
      void calculate_ddgamma(BSSNData *bd);
      void calculate_dalpha_dphi(BSSNData *bd);
      void calculate_dK(BSSNData *bd);
#     if USE_Z4c_DAMPING
        void calculate_dtheta(BSSNData *bd);
#     endif
#     if USE_BSSN_SHIFT
        void calculate_dbeta(BSSNData *bd);
        void calculate_dexpN(BSSNData *bd);
#     endif

    /* Calculate "dependent" quantities (depend on previously calc'd vals) */
      void calculate_conformal_christoffels(BSSNData *bd);

    /* Calculate doubly-"dependent" quantities (depend on previously calc'd vals) */
      void calculateDDphi(BSSNData *bd);
      void calculateRicciTF(BSSNData *bd);
      void calculateDDalphaTF(BSSNData *bd);

      void enforceTFSIJ(BSSNData *bd);

    /* (optional) Calculations of additional quantities */
      void set_full_metric(BSSNData *bd);
      void set_full_metric_der(BSSNData *bd);

  /* Evolution functions */
    real_t ev_DIFFgamma11(BSSNData *bd);
    real_t ev_DIFFgamma12(BSSNData *bd);
    real_t ev_DIFFgamma13(BSSNData *bd);
    real_t ev_DIFFgamma22(BSSNData *bd);
    real_t ev_DIFFgamma23(BSSNData *bd);
    real_t ev_DIFFgamma33(BSSNData *bd);
    real_t ev_A11(BSSNData *bd);
    real_t ev_A12(BSSNData *bd);
    real_t ev_A13(BSSNData *bd);
    real_t ev_A22(BSSNData *bd);
    real_t ev_A23(BSSNData *bd);
    real_t ev_A33(BSSNData *bd);
    real_t ev_DIFFK(BSSNData *bd);
    real_t ev_DIFFphi(BSSNData *bd);
    real_t ev_Gamma1(BSSNData *bd);
    real_t ev_Gamma2(BSSNData *bd);
    real_t ev_Gamma3(BSSNData *bd);

    real_t ev_DIFFalpha(BSSNData *bd);

#   if USE_Z4c_DAMPING
      real_t ev_theta(BSSNData *bd);
#   endif

#   if USE_BSSN_SHIFT
      real_t ev_beta1(BSSNData *bd);
      real_t ev_beta2(BSSNData *bd);
      real_t ev_beta3(BSSNData *bd);
      real_t ev_expN(BSSNData *bd);
#   endif

#   if USE_GAMMA_DRIVER
      real_t ev_auxB1(BSSNData *bd);
      real_t ev_auxB2(BSSNData *bd);
      real_t ev_auxB3(BSSNData *bd);
#   endif

  /* constraint violation calculations */
    void set1DConstraintOutput(
      real_t H_values[], real_t M_values[], int axis, idx_t n1, idx_t n2);

    void setConstraintCalcs(real_t H_values[7], real_t M_values[7],
      real_t G_values[7], real_t A_values[7], real_t S_values[7]);

    real_t hamiltonianConstraintCalc(BSSNData *bd);
    real_t hamiltonianConstraintScale(BSSNData *bd);

    real_t momentumConstraintCalc(BSSNData *bd, idx_t d);
    real_t momentumConstraintScale(BSSNData *bd, idx_t d);

    real_t christoffelConstraintCalc(BSSNData *bd, idx_t d);
    real_t christoffelConstraintScale(BSSNData *bd, idx_t d);

    real_t AijTFConstraintCalc(BSSNData *bd);
    real_t AijTFConstraintScale(BSSNData *bd);

    real_t unitDetConstraintCalc(BSSNData *bd);
    real_t unitDetConstraintScale(BSSNData *bd);

# if USE_COSMOTRACE
  /* Raytracing functionality */
    RaytracePrimitives<real_t> getRaytraceData(BSSNData *bd);
    void setRaytraceCornerPrimitives(RayTrace<real_t, idx_t> *rt);
    void setRaytracePrimitives(RayTrace<real_t, idx_t> *rt);
# endif

};

}

#endif
