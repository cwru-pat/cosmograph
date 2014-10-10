#ifndef COSMO_BSSN
#define COSMO_BSSN

#include "cosmo.h"
#include "globals.h"

#include "bssn_macros.h"

namespace cosmo
{

typedef struct {

  idx_t i, j, k, idx;

  // generic var for misc. expressions
  real_t trace, expression;

  // trace-free ricci tensor components
  real_t ricciTF11, ricciTF12, ricciTF13, ricciTF22, ricciTF23, ricciTF33;

  // derivatives of \alpha
    // covariant double-derivatives 
    real_t D1D1aTF, D1D2aTF, D1D3aTF, D2D2aTF, D2D3aTF, D3D3aTF;
    real_t DDaTR;
    // normal derivatives of
    real_t d1a, d2a, d3a;

  // derivatives of \beta
    real_t d1beta1, d1beta2, d1beta3,
           d2beta1, d2beta2, d2beta3,
           d3beta1, d3beta2, d3beta3;

  // derivatives of phi
    // covariant double-derivatives 
    real_t D1D1phi, D1D2phi, D1D3phi, D2D2phi, D2D3phi, D3D3phi;
    // normal derivatives of
    real_t d1phi, d2phi, d3phi;

  // Contravariant (upstairs index) ext. curvature
  real_t Acont11, Acont12, Acont13, Acont22, Acont23, Acont33;

  // inverse metric
  real_t gammai11, gammai12, gammai13, gammai22, gammai23, gammai33;
  real_t gammai11_adj[3][3][3], gammai12_adj[3][3][3], gammai13_adj[3][3][3],
         gammai22_adj[3][3][3], gammai23_adj[3][3][3], gammai33_adj[3][3][3];

  // Christoffel symbols
  real_t G111, G112, G113, G122, G123, G133,
         G211, G212, G213, G222, G223, G233,
         G311, G312, G313, G322, G323, G333;

  // derivatives of the metric, d_i g_jk
  real_t d1g11, d1g12, d1g13, d1g22, d1g23, d1g33,
         d2g11, d2g12, d2g13, d2g22, d2g23, d2g33,
         d3g11, d3g12, d3g13, d3g22, d3g23, d3g33;

  // derivatives of the inverse metric d_i g^jk
  real_t d1gi11, d1gi12, d1gi13, d1gi22, d1gi23, d1gi33,
         d2gi11, d2gi12, d2gi13, d2gi22, d2gi23, d2gi33,
         d3gi11, d3gi12, d3gi13, d3gi22, d3gi23, d3gi33;

  // second derivatives of the metric d_i d_j g_kl
  real_t d1d1g11, d1d1g12, d1d1g13, d1d1g22, d1d1g23, d1d1g33,
         d1d2g11, d1d2g12, d1d2g13, d1d2g22, d1d2g23, d1d2g33,
         d1d3g11, d1d3g12, d1d3g13, d1d3g22, d1d3g23, d1d3g33,
         d2d2g11, d2d2g12, d2d2g13, d2d2g22, d2d2g23, d2d2g33,
         d2d3g11, d2d3g12, d2d3g13, d2d3g22, d2d3g23, d2d3g33,
         d3d3g11, d3d3g12, d3d3g13, d3d3g22, d3d3g23, d3d3g33;

  // local copies of current field values
  BSSN_APPLY_TO_FIELDS(DECLARE_REAL_T)

  // local copies of adjacent current field values for fast derivatives
  BSSN_APPLY_TO_FIELDS(DECLARE_ADJACENT_REAL_T)

  // Source terms
  real_t rho;
  real_t S11, S12, S13, S23, S33;

} BSSNData;


/** BSSN class **/
class BSSN
{
  /* arrays for storing fields */
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_CREATE)
  /* local values */
  BSSNData paq;

public:
  std::map <std::string, real_t *> fields;

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
  void set_local_vals(BSSNData *paq);
  void calculate_Acont(BSSNData *paq);

  /* Calculate metric derivatives */
  void calculate_dgamma(BSSNData *paq);
  void calculate_ddgamma(BSSNData *paq);

  /* Calculate metric derivatives */
  void calculate_dgammai(BSSNData *paq);
  void calculate_christoffels(BSSNData *paq);
  void calculate_dalpha_dphi(BSSNData *paq);
  void calculate_dbeta(BSSNData *paq);

  /* Calculate trace-free ricci tensor components */
  void calculateRicciTF(BSSNData *paq);
  void calculateDDphi(BSSNData *paq);
  void calculateDDalphaTF(BSSNData *paq);

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
  // static gauge for now:
  real_t ev_alpha(BSSNData *paq) { return 0; }
  real_t ev_beta1(BSSNData *paq) { return 0; }
  real_t ev_beta2(BSSNData *paq) { return 0; }
  real_t ev_beta3(BSSNData *paq) { return 0; }

  // calculate needed quantities (need the inverse metric set everywhere first)
  void set_paq_values(idx_t i, idx_t j, idx_t k, BSSNData *paq);

  void step();
  void init();

};

}

#endif
