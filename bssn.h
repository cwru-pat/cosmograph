#ifndef COSMO_BSSN
#define COSMO_BSSN

#include "cosmo.h"

namespace cosmo
{

#include "BSSN_macros.h"

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

  // derivatives of phi
    // covariant double-derivatives 
    real_t D1D1phi, D1D2phi, D1D3phi, D2D2phi, D2D3phi, D3D3phi;
    // normal derivatives of
    real_t d1phi, d2phi, d3phi;

  // Contravariant (upstairs index) ext. curvature
  real_t Acont11, Acont12, Acont13, Acont22, Acont23, Acont33;

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

  BSSN()
  {
    BSSN_APPLY_TO_FIELDS(RK4_ARRAY_ALLOC)
    BSSN_APPLY_TO_FIELDS(RK4_ARRAY_ADDMAP)
  }

  ~BSSN()
  {
    BSSN_APPLY_TO_FIELDS(RK4_ARRAY_DELETE)
  }

  void step() {

    // Init arr_a with _p values
    BSSN_COPY_ARRAYS(_p, _a);

    // First RK step: add in k_1 coeff
    LOOP3(i, j, k)
    {
      set_paq_values(i, j, k, &paq);

      // evolve fields: arr_c = arr_p + dt/2*k_1
        // arr_c[idx] = arr_p[idx] + dt/2.0*evfn(arr_a);
        BSSN_COMPUTE_RK_STEP(0.5);
        // algebraically calculate gammai values in _c register
        BSSN_COMPUTE_GAMMAI(_c);

      // add computation to _f array: arr_f = arr_p + dt/2*k_1
        // arr_f += arr_c
        BSSN_ADD_C_TO_F(1.0);
    }
    // swap _c <-> _a registers
    BSSN_SWAP_ARRAYS(_c, _a);

    // Calculate k_2 coeff (using k_1 values now in _a register):
    LOOP3(i, j, k)
    {
      set_paq_values(i, j, k, &paq);

      // evolve fields: arr_c = arr_p + dt/2*k_2
        // arr_c[idx] = arr_p[idx] + dt/2.0*evfn(arr_a);
        BSSN_COMPUTE_RK_STEP(0.5);
        // algebraically calculate gammai values in _c register
        BSSN_COMPUTE_GAMMAI(_c);

      // add computation to _f array: arr_f = 3*arr_p + dt/2*(k_1 + 2*k_2)
        // arr_f += 2.0*arr_c
        BSSN_ADD_C_TO_F(2.0);
    }
    // swap _c <-> _a
    BSSN_SWAP_ARRAYS(_c, _a);

    // Calculate k_3 coeff (using k_2 values now in _a register):
    LOOP3(i, j, k)
    {
      set_paq_values(i, j, k, &paq);

      // evolve fields: arr_c = arr_p + dt*k_3
        // arr_c[idx] = arr_p[idx] + dt*evfn(arr_a);
        BSSN_COMPUTE_RK_STEP(1.0);
        // algebraically calculate gammai values in _c register
        BSSN_COMPUTE_GAMMAI(_c);

      // add computation to _f array: arr_f = 4*arr_p + dt/2*(k_1 + 2*k_2 + 2*k_3)
        // arr_f += arr_c
        BSSN_ADD_C_TO_F(1.0);
    }
    // swap _c <-> _a
    BSSN_SWAP_ARRAYS(_c, _a);

    // Add in k_4 contribution to _f register and "weight" correctly:
    LOOP3(i, j, k)
    {
      set_paq_values(i, j, k, &paq);

      // evolve fields and add to _f register:
      // arr_f = arr_p + dt/6*(k_1 + 2*k_2 + 2*k_3 + k_4)
        // arr_f = (1.0/3.0)*(arr_f - arr_p) + (1.0/6.0)*evfn(arr_a)
      
      // calculate gamma_i in _f register
      BSSN_COMPUTE_GAMMAI(_f);
    }
    // arr_f register now holds "final" calculation; move back to _p register:
    // swap _f <-> _p
    BSSN_SWAP_ARRAYS(_f, _p);

    // done!
  }

  inline real_t der(real_t *field, int d, PointData *paq)
  {
    return derivative_stencil(paq->i, paq->j, paq->k, d, field);
  }

  inline real_t dder(real_t *field, int d1, int d2, PointData *paq)
  {
    return double_derivative(paq->i, paq->j, paq->k, d1, d2, field);
  }

  /* set current local field values */
  inline void set_local_vals(PointData *paq)
  {
    BSSN_APPLY_TO_FIELDS(SET_LOCAL_VALUES)
  }

  inline void calculate_Acont(PointData *paq)
  {
    BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_ACONT)
  }

  /* Calculate metric derivatives */
  inline void calculate_dgamma(PointData *paq)
  {
    BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_DGAMMA)
  }

  inline void calculate_ddgamma(PointData *paq)
  {
    BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_DIDJGAMMA_PERMS)
  }

  /* Calculate metric derivatives */
  inline void calculate_dgammai(PointData *paq)
  {
    BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_DGAMMAI);
  }

  inline void calculate_christoffels(PointData *paq)
  {
    // christoffel symbols: \Gamma^i_{jk} = Gijk
    BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_CHRISTOFFEL)
  }

  inline void calculate_dalpha_dphi(PointData *paq)
  {
    // normal derivatives of phi
    paq->d1phi = der(phi_a, 1, paq);
    paq->d2phi = der(phi_a, 2, paq);
    paq->d3phi = der(phi_a, 3, paq);

    // normal derivatives of alpha
    paq->d1a = der(alpha_a, 1, paq);
    paq->d2a = der(alpha_a, 2, paq);
    paq->d3a = der(alpha_a, 3, paq);
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
    paq->ricciTF11 += -2.0*( paq->D1D1phi -2.0*paq->d1phi*paq->d1phi + paq->gamma11*(expression) );
    paq->ricciTF12 += -2.0*( paq->D1D2phi -2.0*paq->d1phi*paq->d2phi + paq->gamma12*(expression) );
    paq->ricciTF13 += -2.0*( paq->D1D3phi -2.0*paq->d1phi*paq->d3phi + paq->gamma13*(expression) );
    paq->ricciTF22 += -2.0*( paq->D2D2phi -2.0*paq->d2phi*paq->d2phi + paq->gamma22*(expression) );
    paq->ricciTF23 += -2.0*( paq->D2D3phi -2.0*paq->d2phi*paq->d3phi + paq->gamma23*(expression) );
    paq->ricciTF33 += -2.0*( paq->D3D3phi -2.0*paq->d3phi*paq->d3phi + paq->gamma33*(expression) );

    /* remove trace... */ 
    paq->trace = paq->gammai11*paq->ricciTF11 + paq->gammai22*paq->ricciTF22 + paq->gammai33*paq->ricciTF33
        + 2.0*(paq->gammai12*paq->ricciTF12 + paq->gammai13*paq->ricciTF13 + paq->gammai23*paq->ricciTF23);
    
    paq->ricciTF11 -= (1/3.0)*paq->gamma11*paq->trace;
    paq->ricciTF12 -= (1/3.0)*paq->gamma12*paq->trace;
    paq->ricciTF13 -= (1/3.0)*paq->gamma13*paq->trace;
    paq->ricciTF22 -= (1/3.0)*paq->gamma22*paq->trace;
    paq->ricciTF23 -= (1/3.0)*paq->gamma23*paq->trace;
    paq->ricciTF33 -= (1/3.0)*paq->gamma33*paq->trace;
  }

  void calculateDDphi(PointData *paq)
  {
    // double covariant derivatives, using normal metric
    paq->D1D1phi = dder(phi_a, 1, 1, paq) - (paq->G111*paq->d1phi + paq->G211*paq->d2phi + paq->G311*paq->d3phi);
    paq->D1D2phi = dder(phi_a, 1, 2, paq) - (paq->G112*paq->d1phi + paq->G212*paq->d2phi + paq->G312*paq->d3phi);
    paq->D1D3phi = dder(phi_a, 1, 3, paq) - (paq->G113*paq->d1phi + paq->G213*paq->d2phi + paq->G313*paq->d3phi);
    paq->D2D2phi = dder(phi_a, 2, 2, paq) - (paq->G122*paq->d1phi + paq->G222*paq->d2phi + paq->G322*paq->d3phi);
    paq->D2D3phi = dder(phi_a, 2, 3, paq) - (paq->G123*paq->d1phi + paq->G223*paq->d2phi + paq->G323*paq->d3phi);
    paq->D3D3phi = dder(phi_a, 3, 3, paq) - (paq->G133*paq->d1phi + paq->G233*paq->d2phi + paq->G333*paq->d3phi);
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
      + paq->beta1*der(K_a, 1, paq) + paq->beta2*der(K_a, 2, paq) + paq->beta3*der(K_a, 3, paq)
    );
  }

  real_t ev_phi(PointData *paq)
  {
    return (
      1.0/6.0*(der(beta1_a, 1, paq) + der(beta2_a, 2, paq) + der(beta3_a, 3, paq) - paq->alpha*paq->K)
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

  /* calculate needed quantities (need the inverse metric set everywhere first) */
  inline void set_paq_values(idx_t i, idx_t j, idx_t k, PointData *paq)
  {
    set_paq_index(i, j, k, paq);
    set_local_vals(paq);
    // gammas & derivs first
    calculate_Acont(paq);
    calculate_dgamma(paq);
    calculate_ddgamma(paq);
    calculate_dgammai(paq);
    calculate_dalpha_dphi(paq);
    // Christoffels depend on metric & derivs.
    calculate_christoffels(paq);
    // Ricci, DDa, DDw depend on christoffels, metric, and derivs
    calculateRicciTF(paq);
    calculateDDphi(paq);
    calculateDDalphaTF(paq);
  }

  inline void set_paq_index(idx_t i, idx_t j, idx_t k, PointData *paq)
  {
    paq->i = i;
    paq->i = j;
    paq->i = k;
    paq->idx = INDEX(i,j,k);
  }

};

}

#endif
