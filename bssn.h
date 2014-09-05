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

  // source?
  real_t rho; // energy density

} PointData;


/** BSSN class **/
class BSSN
{
  /* arrays for storing fields */
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_CREATE)
  /* local values */
  PointData paq;

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

  real_t der(real_t field_adj[3][3][3], int d)
  {
   // return 0;
    
    switch (d) {
      case 1:
        return (field_adj[2][1][1] - field_adj[0][1][1])/dx/2.0;
        break;
      case 2:
        return (field_adj[1][2][1] - field_adj[1][0][1])/dx/2.0;
        break;
      case 3:
        return (field_adj[1][1][2] - field_adj[1][1][0])/dx/2.0;
        break;
    }

    /* XXX */
    return 0;
  }

  real_t dder(real_t field_adj[3][3][3], int d1, int d2)
  {
    //return 0;

    switch (d1) {
      case 1:
        switch (d2) {
          case 1:
            return (field_adj[0][1][1] + field_adj[2][1][1] - 2.0*field_adj[1][1][1])/dx/dx;
            break;
          case 2:
            return (field_adj[0][0][1] + field_adj[2][2][1] - field_adj[0][2][1] - field_adj[2][0][1])/dx/dx/4.0;
            break;
          case 3:
            return (field_adj[0][1][0] + field_adj[2][1][2] - field_adj[0][1][2] - field_adj[2][1][0])/dx/dx/4.0;
            break;
        }
        break;
      case 2:
        switch (d2) {
          case 1:
            return (field_adj[0][0][1] + field_adj[2][2][1] - field_adj[0][2][1] - field_adj[2][0][1])/dx/dx/4.0;
            break;
          case 2:
            return (field_adj[1][0][1] + field_adj[1][2][1] - 2.0*field_adj[1][1][1])/dx/dx;
            break;
          case 3:
            return (field_adj[1][0][0] + field_adj[1][2][2] - field_adj[1][0][2] - field_adj[1][2][0])/dx/dx/4.0;
            break;
        }
        break;
      case 3:
        switch (d2) {
          case 1:
            return (field_adj[0][1][0] + field_adj[2][1][2] - field_adj[0][1][2] - field_adj[2][1][0])/dx/dx/4.0;
            break;
          case 2:
            return (field_adj[1][0][0] + field_adj[1][2][2] - field_adj[1][0][2] - field_adj[1][2][0])/dx/dx/4.0;
            break;
          case 3:
            return (field_adj[1][1][0] + field_adj[1][1][2] - 2.0*field_adj[1][1][1])/dx/dx;
            break;
        }
        break;
    }

    /* XXX */
    return 0.0;
  }
  

  /* set current local field values */
  void set_local_vals(PointData *paq)
  {
    SET_LOCAL_INDEXES;

    // Need full mixed 2nd derivatives for all of these
    SET_LOCAL_VALUES_PFE(gamma11);
    SET_LOCAL_VALUES_PFE(gamma12);
    SET_LOCAL_VALUES_PFE(gamma13);
    SET_LOCAL_VALUES_PFE(gamma22);
    SET_LOCAL_VALUES_PFE(gamma23);
    SET_LOCAL_VALUES_PFE(gamma33);
    SET_LOCAL_VALUES_PFE(phi);
    SET_LOCAL_VALUES_PFE(K);

    SET_LOCAL_VALUES_PFE(beta1);
    SET_LOCAL_VALUES_PFE(beta2);
    SET_LOCAL_VALUES_PFE(beta3);
    
    SET_LOCAL_VALUES_PFE(alpha);

    // Seems faster to compute all of these at each step, rather than pulling them out of memory.
    BSSN_COMPUTE_LOCAL_GAMMAI_PF(11, 22, 33, 23, 23);
    BSSN_COMPUTE_LOCAL_GAMMAI_PF(12, 13, 23, 12, 33);
    BSSN_COMPUTE_LOCAL_GAMMAI_PF(13, 12, 23, 13, 22);
    BSSN_COMPUTE_LOCAL_GAMMAI_PF(22, 11, 33, 13, 13);
    BSSN_COMPUTE_LOCAL_GAMMAI_PF(23, 12, 13, 23, 11);
    BSSN_COMPUTE_LOCAL_GAMMAI_PF(33, 11, 22, 12, 12);

    // adjacent values only
    SET_LOCAL_VALUES_PF(Gamma1);
    SET_LOCAL_VALUES_PF(Gamma2);
    SET_LOCAL_VALUES_PF(Gamma3);
    SET_LOCAL_VALUES_PF(A11);
    SET_LOCAL_VALUES_PF(A12);
    SET_LOCAL_VALUES_PF(A13);
    SET_LOCAL_VALUES_PF(A22);
    SET_LOCAL_VALUES_PF(A23);
    SET_LOCAL_VALUES_PF(A33);
  }

  void calculate_Acont(PointData *paq)
  {
    BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_ACONT)
  }

  /* Calculate metric derivatives */
  void calculate_dgamma(PointData *paq)
  {
    BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_DGAMMA)
  }

  void calculate_ddgamma(PointData *paq)
  {
    BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_DIDJGAMMA_PERMS)
  }

  /* Calculate metric derivatives */
  void calculate_dgammai(PointData *paq)
  {
    BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_DGAMMAI);
  }

  void calculate_christoffels(PointData *paq)
  {
    // christoffel symbols: \Gamma^i_{jk} = Gijk
    BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_CHRISTOFFEL)
  }

  void calculate_dalpha_dphi(PointData *paq)
  {
    // normal derivatives of phi
    paq->d1phi = der(paq->phi_adj, 1);
    paq->d2phi = der(paq->phi_adj, 2);
    paq->d3phi = der(paq->phi_adj, 3);

    // normal derivatives of alpha
    paq->d1a = der(paq->alpha_adj, 1);
    paq->d2a = der(paq->alpha_adj, 2);
    paq->d3a = der(paq->alpha_adj, 3);
  }

  void calculate_dbeta(PointData *paq)
  {
    paq->d1beta1 = der(paq->beta1_adj, 1);
    paq->d1beta2 = der(paq->beta2_adj, 1);
    paq->d1beta3 = der(paq->beta3_adj, 1);
    paq->d2beta1 = der(paq->beta1_adj, 2);
    paq->d2beta2 = der(paq->beta2_adj, 2);
    paq->d2beta3 = der(paq->beta3_adj, 2);
    paq->d3beta1 = der(paq->beta1_adj, 3);
    paq->d3beta2 = der(paq->beta2_adj, 3);
    paq->d3beta3 = der(paq->beta3_adj, 3);
  }

  /* Calculate trace-free ricci tensor components */
  void calculateRicciTF(PointData *paq)
  {
    // unitary pieces
    BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_RICCITF_UNITARY)

    real_t expression = (
      paq->gammai11*(paq->D1D1phi - 2.0*paq->d1phi*paq->d1phi)
      + paq->gammai22*(paq->D2D2phi - 2.0*paq->d2phi*paq->d2phi)
      + paq->gammai33*(paq->D3D3phi - 2.0*paq->d3phi*paq->d3phi)
      + 2.0*(
        paq->gammai12*(paq->D1D2phi - 2.0*paq->d1phi*paq->d2phi)
        + paq->gammai13*(paq->D1D3phi - 2.0*paq->d1phi*paq->d3phi)
        + paq->gammai23*(paq->D2D3phi - 2.0*paq->d2phi*paq->d3phi)
      )
    );

    /* phi-piece */
    paq->ricciTF11 += -2.0*( paq->D1D1phi - 2.0*paq->d1phi*paq->d1phi + paq->gamma11*(expression) );
    paq->ricciTF12 += -2.0*( paq->D1D2phi - 2.0*paq->d1phi*paq->d2phi + paq->gamma12*(expression) );
    paq->ricciTF13 += -2.0*( paq->D1D3phi - 2.0*paq->d1phi*paq->d3phi + paq->gamma13*(expression) );
    paq->ricciTF22 += -2.0*( paq->D2D2phi - 2.0*paq->d2phi*paq->d2phi + paq->gamma22*(expression) );
    paq->ricciTF23 += -2.0*( paq->D2D3phi - 2.0*paq->d2phi*paq->d3phi + paq->gamma23*(expression) );
    paq->ricciTF33 += -2.0*( paq->D3D3phi - 2.0*paq->d3phi*paq->d3phi + paq->gamma33*(expression) );

    /* remove trace... */ 
    paq->trace = paq->gammai11*paq->ricciTF11 + paq->gammai22*paq->ricciTF22 + paq->gammai33*paq->ricciTF33
        + 2.0*(paq->gammai12*paq->ricciTF12 + paq->gammai13*paq->ricciTF13 + paq->gammai23*paq->ricciTF23);
  
    paq->ricciTF11 -= (1.0/3.0)*paq->gamma11*paq->trace;
    paq->ricciTF12 -= (1.0/3.0)*paq->gamma12*paq->trace;
    paq->ricciTF13 -= (1.0/3.0)*paq->gamma13*paq->trace;
    paq->ricciTF22 -= (1.0/3.0)*paq->gamma22*paq->trace;
    paq->ricciTF23 -= (1.0/3.0)*paq->gamma23*paq->trace;
    paq->ricciTF33 -= (1.0/3.0)*paq->gamma33*paq->trace;
  }

  void calculateDDphi(PointData *paq)
  {
    // double covariant derivatives, using normal metric
    paq->D1D1phi = dder(paq->phi_adj, 1, 1) - (paq->G111*paq->d1phi + paq->G211*paq->d2phi + paq->G311*paq->d3phi);
    paq->D1D2phi = dder(paq->phi_adj, 1, 2) - (paq->G112*paq->d1phi + paq->G212*paq->d2phi + paq->G312*paq->d3phi);
    paq->D1D3phi = dder(paq->phi_adj, 1, 3) - (paq->G113*paq->d1phi + paq->G213*paq->d2phi + paq->G313*paq->d3phi);
    paq->D2D2phi = dder(paq->phi_adj, 2, 2) - (paq->G122*paq->d1phi + paq->G222*paq->d2phi + paq->G322*paq->d3phi);
    paq->D2D3phi = dder(paq->phi_adj, 2, 3) - (paq->G123*paq->d1phi + paq->G223*paq->d2phi + paq->G323*paq->d3phi);
    paq->D3D3phi = dder(paq->phi_adj, 3, 3) - (paq->G133*paq->d1phi + paq->G233*paq->d2phi + paq->G333*paq->d3phi);
  }

  void calculateDDalphaTF(PointData *paq)
  {
    // double covariant derivatives - use non-unitary metric - extra pieces that depend on phi!
    // the gamma*ldlphi are needed for the BSSN_CALCULATE_DIDJALPHA macro
    real_t gamma1ldlphi = paq->gammai11*paq->d1phi + paq->gammai12*paq->d2phi + paq->gammai13*paq->d3phi;
    real_t gamma2ldlphi = paq->gammai21*paq->d1phi + paq->gammai22*paq->d2phi + paq->gammai23*paq->d3phi;
    real_t gamma3ldlphi = paq->gammai31*paq->d1phi + paq->gammai32*paq->d2phi + paq->gammai33*paq->d3phi;
    BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_DIDJALPHA)

    // subtract trace
    paq->DDaTR = paq->gammai11*paq->D1D1aTF + paq->gammai22*paq->D2D2aTF + paq->gammai33*paq->D3D3aTF
        + 2.0*(paq->gammai12*paq->D1D2aTF + paq->gammai13*paq->D1D3aTF + paq->gammai23*paq->D2D3aTF);
    
    paq->D1D1aTF -= (1/3.0)*paq->gamma11*paq->DDaTR;
    paq->D1D2aTF -= (1/3.0)*paq->gamma12*paq->DDaTR;
    paq->D1D3aTF -= (1/3.0)*paq->gamma13*paq->DDaTR;
    paq->D2D2aTF -= (1/3.0)*paq->gamma22*paq->DDaTR;
    paq->D2D3aTF -= (1/3.0)*paq->gamma23*paq->DDaTR;
    paq->D3D3aTF -= (1/3.0)*paq->gamma33*paq->DDaTR;
  }

  real_t ev_gamma11(PointData *paq) { return BSSN_DT_GAMMAIJ(1, 1); }
  real_t ev_gamma12(PointData *paq) { return BSSN_DT_GAMMAIJ(1, 2); }
  real_t ev_gamma13(PointData *paq) { return BSSN_DT_GAMMAIJ(1, 3); }
  real_t ev_gamma22(PointData *paq) { return BSSN_DT_GAMMAIJ(2, 2); }
  real_t ev_gamma23(PointData *paq) { return BSSN_DT_GAMMAIJ(2, 3); }
  real_t ev_gamma33(PointData *paq) { return BSSN_DT_GAMMAIJ(3, 3); }

  real_t ev_A11(PointData *paq) { return BSSN_DT_AIJ(1, 1); }
  real_t ev_A12(PointData *paq) { return BSSN_DT_AIJ(1, 2); }
  real_t ev_A13(PointData *paq) { return BSSN_DT_AIJ(1, 3); }
  real_t ev_A22(PointData *paq) { return BSSN_DT_AIJ(2, 2); }
  real_t ev_A23(PointData *paq) { return BSSN_DT_AIJ(2, 3); }
  real_t ev_A33(PointData *paq) { return BSSN_DT_AIJ(3, 3); }

  real_t ev_K(PointData *paq)
  {
    return (
      - paq->DDaTR
      + paq->alpha*(
          paq->A11*paq->A11 + paq->A22*paq->A22 + paq->A33*paq->A33
          + 2.0*(paq->A12*paq->A12 + paq->A13*paq->A13 + paq->A23*paq->A23)
          + (1.0/3.0)*paq->K*paq->K
        )
      + paq->beta1*der(paq->K_adj, 1) + paq->beta2*der(paq->K_adj, 2) + paq->beta3*der(paq->K_adj, 3)
      + 4.0*PI*paq->alpha*paq->rho
    );
  }

  real_t ev_phi(PointData *paq)
  {
    return (
      1.0/6.0*(der(paq->beta1_adj, 1) + der(paq->beta2_adj, 2) + der(paq->beta3_adj, 3) - paq->alpha*paq->K)
      + paq->beta1*paq->d1phi + paq->beta2*paq->d2phi + paq->beta3*paq->d3phi
    );
  }

  real_t ev_Gamma1(PointData *paq) { return BSSN_DT_GAMMAI(1); }
  real_t ev_Gamma2(PointData *paq) { return BSSN_DT_GAMMAI(2); }
  real_t ev_Gamma3(PointData *paq) { return BSSN_DT_GAMMAI(3); }

  // static gauge for now:
  real_t ev_alpha(PointData *paq) { return 0; }
  real_t ev_beta1(PointData *paq) { return 0; }
  real_t ev_beta2(PointData *paq) { return 0; }
  real_t ev_beta3(PointData *paq) { return 0; }

  // calculate needed quantities (need the inverse metric set everywhere first)
  void set_paq_values(idx_t i, idx_t j, idx_t k, PointData *paq)
  {
    paq->i = i;
    paq->j = j;
    paq->k = k;
    paq->idx = INDEX(i,j,k);

    // source terms?
    paq->rho = -3.0/PI/8.0;

    // draw data from cache
    set_local_vals(paq);

    // pre-compute re-used quantities
    // gammas & derivs first
    calculate_Acont(paq);
    calculate_dgamma(paq);
    calculate_ddgamma(paq);
    calculate_dgammai(paq);
    calculate_dbeta(paq);
    calculate_dalpha_dphi(paq);
    // Christoffels depend on metric & derivs.
    calculate_christoffels(paq);
    // Ricci, DDa, DDw depend on christoffels, metric, and derivs
    calculateRicciTF(paq);
    calculateDDphi(paq);
    calculateDDalphaTF(paq);
  }

  void step();
  void init();

};



}

#endif
