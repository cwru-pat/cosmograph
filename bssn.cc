#include "bssn.h"

namespace cosmo
{

BSSN::BSSN()
{
  // BSSN fields
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_ALLOC)
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_ADDMAP)
  // BSSN source fields
  BSSN_APPLY_TO_SOURCES(GEN1_ARRAY_ALLOC)
  BSSN_APPLY_TO_SOURCES(GEN1_ARRAY_ADDMAP)

  // any additional arrays for calcuated quantities
  GEN1_ARRAY_ALLOC(dk0_slice_phi);
  GEN1_ARRAY_ADDMAP(dk0_slice_phi);

  GEN1_ARRAY_ALLOC(KDx);
  GEN1_ARRAY_ADDMAP(KDx);

  GEN1_ARRAY_ALLOC(KDy);
  GEN1_ARRAY_ADDMAP(KDy);

  GEN1_ARRAY_ALLOC(KDz);
  GEN1_ARRAY_ADDMAP(KDz);

  GEN1_ARRAY_ALLOC(ricci);
  GEN1_ARRAY_ADDMAP(ricci);

  GEN1_ARRAY_ALLOC(AijAij);
  GEN1_ARRAY_ADDMAP(AijAij);

  GEN1_ARRAY_ALLOC(gammai11);
  GEN1_ARRAY_ADDMAP(gammai11);

  GEN1_ARRAY_ALLOC(gammai12);
  GEN1_ARRAY_ADDMAP(gammai12);

  GEN1_ARRAY_ALLOC(gammai13);
  GEN1_ARRAY_ADDMAP(gammai13);

  GEN1_ARRAY_ALLOC(gammai22);
  GEN1_ARRAY_ADDMAP(gammai22);

  GEN1_ARRAY_ALLOC(gammai23);
  GEN1_ARRAY_ADDMAP(gammai23);

  GEN1_ARRAY_ALLOC(gammai33);
  GEN1_ARRAY_ADDMAP(gammai33);
}

BSSN::~BSSN()
{
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_DELETE)
  BSSN_APPLY_TO_SOURCES(GEN1_ARRAY_DELETE)
  GEN1_ARRAY_DELETE(dk0_slice_phi);
  GEN1_ARRAY_DELETE(KDx);
  GEN1_ARRAY_DELETE(KDy);
  GEN1_ARRAY_DELETE(KDz);
  GEN1_ARRAY_DELETE(ricci);
  GEN1_ARRAY_DELETE(AijAij);
  GEN1_ARRAY_DELETE(gammai11);
  GEN1_ARRAY_DELETE(gammai12);
  GEN1_ARRAY_DELETE(gammai13);
  GEN1_ARRAY_DELETE(gammai22);
  GEN1_ARRAY_DELETE(gammai23);
  GEN1_ARRAY_DELETE(gammai33);
}

void BSSN::set_gammai_values(idx_t i, idx_t j, idx_t k)
{
  // Compute the inverse metric algebraically at each point
  idx_t idx = NP_INDEX(i,j,k);
  gammai11_a[idx] = gamma22_a[idx]*gamma33_a[idx] - gamma23_a[idx]*gamma23_a[idx];
  gammai12_a[idx] = gamma13_a[idx]*gamma23_a[idx] - gamma12_a[idx]*gamma33_a[idx];
  gammai13_a[idx] = gamma12_a[idx]*gamma23_a[idx] - gamma13_a[idx]*gamma22_a[idx];
  gammai22_a[idx] = gamma11_a[idx]*gamma33_a[idx] - gamma13_a[idx]*gamma13_a[idx];
  gammai23_a[idx] = gamma12_a[idx]*gamma13_a[idx] - gamma23_a[idx]*gamma11_a[idx];
  gammai33_a[idx] = gamma11_a[idx]*gamma22_a[idx] - gamma12_a[idx]*gamma12_a[idx];
}

void BSSN::set_paq_values(idx_t i, idx_t j, idx_t k, BSSNData *paq)
{
  paq->i = i;
  paq->j = j;
  paq->k = k;
  paq->idx = NP_INDEX(i,j,k);

  // draw data from cache
  set_local_vals(paq);

  // source values
  set_source_vals(paq);

  // pre-compute re-used quantities
  // gammas & derivs first
  calculate_Acont(paq);
  calculate_dgamma(paq);
  calculate_ddgamma(paq);
  calculate_dgammai(paq);
  calculate_dphi(paq);
  calculate_dK(paq);
  #if Z4c_DAMPING > 0
    calculate_dtheta(paq);
  #endif

  // Christoffels depend on metric & derivs.
  calculate_conformal_christoffels(paq);
  // DDw depend on christoffels, metric, and derivs
  calculateDDphi(paq);
  // Ricci depends on DDphi
  calculateRicciTF(paq);

  paq->H = hamiltonianConstraintCalc(paq);
}

void BSSN::set_detgamma(idx_t i, idx_t j, idx_t k)
{
  idx_t idx = NP_INDEX(i,j,k);

  real_t det_p = pow(
    - gamma13_p[idx] * gamma13_p[idx] * gamma22_p[idx]
    + 2.0*gamma12_p[idx] * gamma13_p[idx] * gamma23_p[idx]
    - gamma11_p[idx] * gamma23_p[idx] * gamma23_p[idx]
    - gamma12_p[idx] * gamma12_p[idx] * gamma33_p[idx]
    + gamma11_p[idx] * gamma22_p[idx] * gamma33_p[idx]
  , 1.0/3.0);

  gamma11_p[idx] /= det_p; gamma12_p[idx] /= det_p; gamma13_p[idx] /= det_p; gamma22_p[idx] /= det_p; gamma23_p[idx] /= det_p; gamma33_p[idx] /= det_p;
  A11_p[idx] /= det_p; A12_p[idx] /= det_p; A13_p[idx] /= det_p; A22_p[idx] /= det_p; A23_p[idx] /= det_p; A33_p[idx] /= det_p;

  real_t det_a = pow(
    - gamma13_a[idx] * gamma13_a[idx] * gamma22_a[idx]
    + 2.0*gamma12_a[idx] * gamma13_a[idx] * gamma23_a[idx]
    - gamma11_a[idx] * gamma23_a[idx] * gamma23_a[idx]
    - gamma12_a[idx] * gamma12_a[idx] * gamma33_a[idx]
    + gamma11_a[idx] * gamma22_a[idx] * gamma33_a[idx]
  , 1.0/3.0);

  gamma11_a[idx] /= det_a; gamma12_a[idx] /= det_a; gamma13_a[idx] /= det_a; gamma22_a[idx] /= det_a; gamma23_a[idx] /= det_a; gamma33_a[idx] /= det_a;
  A11_a[idx] /= det_a; A12_a[idx] /= det_a; A13_a[idx] /= det_a; A22_a[idx] /= det_a; A23_a[idx] /= det_a; A33_a[idx] /= det_a;

}

// Full RK step (More useful when not evolving the source simultaneously)
void BSSN::step(BSSNData *paq)
{
  _timer["RK Init"].start();
  stepInit();
  _timer["RK Init"].stop();

  _timer["RK K1 Calc"].start();
  K1Calc(paq);
  _timer["RK K1 Calc"].stop();

  _timer["RK K2 Calc"].start();
  K2Calc(paq);
  _timer["RK K2 Calc"].stop();

  _timer["RK K3 Calc"].start();
  K3Calc(paq);
  _timer["RK K3 Calc"].stop();

  _timer["RK K4 Calc"].start();
  K4Calc(paq);
  _timer["RK K4 Calc"].stop();

  stepTerm();

  // done!
}

void BSSN::regSwap_c_a()
{
  BSSN_SWAP_ARRAYS(_c, _a);
}

// Init _a register with _p values, _f with 0
void BSSN::stepInit()
{
  BSSN_COPY_ARRAYS(_p, _a);
  idx_t i, j, k;
  #pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i, j, k)
  {
    idx_t idx = NP_INDEX(i, j, k);
    BSSN_ZERO_ARRAYS(_f, idx);
  }
}


// First RK step: calculate and add k_1 coeff to _f array
void BSSN::K1Calc(BSSNData *paq)
{
  idx_t i, j, k;
  LOOP3(i, j, k)
  {
    K1CalcPt(i, j, k, paq);
  }
  // swap _c <-> _a registers
  BSSN_SWAP_ARRAYS(_c, _a);
}
void BSSN::K1CalcPt(idx_t i, idx_t j, idx_t k, BSSNData *paq)
{
  set_paq_values(i, j, k, paq);

  // evolve fields: arr_c = arr_p + dt/2*k_1
    // arr_c[idx] = arr_p[idx] + dt/2.0*evfn(arr_a);
    BSSN_COMPUTE_RK_STEP(0.5);

  // add computation to _f array: arr_f = arr_p + dt/2*k_1
    // arr_f += arr_c
    BSSN_ADD_C_TO_F(1.0);
}


// Calculate k_2 coeff (using k_1 values now in _a register)
// and add to _f array
void BSSN::K2Calc(BSSNData *paq)
{
  idx_t i, j, k;
  LOOP3(i, j, k)
  {
    K2CalcPt(i, j, k, paq);
  }
  // swap _c <-> _a
  BSSN_SWAP_ARRAYS(_c, _a);
}
void BSSN::K2CalcPt(idx_t i, idx_t j, idx_t k, BSSNData *paq)
{
  set_paq_values(i, j, k, paq);

  // evolve fields: arr_c = arr_p + dt/2*k_2
    // arr_c[idx] = arr_p[idx] + dt/2.0*evfn(arr_a);
    BSSN_COMPUTE_RK_STEP(0.5);

  // add computation to _f array: arr_f = 3*arr_p + dt/2*(k_1 + 2*k_2)
    // arr_f += 2.0*arr_c
    BSSN_ADD_C_TO_F(2.0);
}


// Calculate k_3 coeff (using k_2 values now in _a register)
// and add to _f array
void BSSN::K3Calc(BSSNData *paq)
{
  idx_t i, j, k;
  LOOP3(i, j, k)
  {
    K3CalcPt(i, j, k, paq);
  }
  // swap _c <-> _a
  BSSN_SWAP_ARRAYS(_c, _a);
}
void BSSN::K3CalcPt(idx_t i, idx_t j, idx_t k, BSSNData *paq)
{
  set_paq_values(i, j, k, paq);

  // evolve fields: arr_c = arr_p + dt*k_3
    // arr_c[idx] = arr_p[idx] + dt*evfn(arr_a);
    BSSN_COMPUTE_RK_STEP(1.0);

  // add computation to _f array: arr_f = 4*arr_p + dt/2*(k_1 + 2*k_2 + 2*k_3)
    // arr_f += arr_c
    BSSN_ADD_C_TO_F(1.0);
}


// Add in k_4 contribution to _f register,
// and "weight" the final calculation correctly:
void BSSN::K4Calc(BSSNData *paq)
{
  idx_t i, j, k;
  LOOP3(i, j, k)
  {
    K4CalcPt(i, j, k, paq);
  }
}
void BSSN::K4CalcPt(idx_t i, idx_t j, idx_t k, BSSNData *paq)
{
  set_paq_values(i, j, k, paq);

  // evolve fields and add to _f register:
  // arr_f = arr_p + dt/6*(k_1 + 2*k_2 + 2*k_3 + k_4)
  //       = (1.0/3.0)*(arr_f - arr_p) + (dt/6.0)*evfn(arr_a)
  BSSN_FINAL_RK4_STEP();
}


// arr_f register now holds "final" calculation; move back to _p register:
// swap _f <-> _p
void BSSN::stepTerm()
{
  BSSN_SWAP_ARRAYS(_f, _p);
}

void BSSN::clearSrc()
{
  idx_t i, j, k;
  #pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i, j, k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    r_a[idx]   = 0.0;
    S_a[idx]   = 0.0;
    S1_a[idx]  = 0.0;
    S2_a[idx]  = 0.0;
    S3_a[idx]  = 0.0;
    STF11_a[idx] = 0.0;
    STF12_a[idx] = 0.0;
    STF13_a[idx] = 0.0;
    STF22_a[idx] = 0.0;
    STF23_a[idx] = 0.0;
    STF33_a[idx] = 0.0;
  }
}

void BSSN::init()
{
  idx_t idx;

  idx_t i, j, k;
  LOOP3(i, j, k)
  {
    idx = NP_INDEX(i,j,k);

    // default flat static vacuum spacetime.
    gamma11_p[idx] = gamma11_a[idx] = gamma11_c[idx] = gamma11_f[idx]  = 1.0;
    gamma12_p[idx] = gamma12_a[idx] = gamma12_c[idx] = gamma12_f[idx]  = 0.0;
    gamma13_p[idx] = gamma13_a[idx] = gamma13_c[idx] = gamma13_f[idx]  = 0.0;
    gamma22_p[idx] = gamma22_a[idx] = gamma22_c[idx] = gamma22_f[idx]  = 1.0;
    gamma23_p[idx] = gamma23_a[idx] = gamma23_c[idx] = gamma23_f[idx]  = 0.0;
    gamma33_p[idx] = gamma33_a[idx] = gamma33_c[idx] = gamma33_f[idx]  = 1.0;
    
    phi_p[idx]     = phi_a[idx]     = phi_c[idx]     = phi_f[idx]      = 0.0;
    
    A11_p[idx]     = A11_a[idx]     = A11_c[idx]     = A11_f[idx]      = 0.0;
    A12_p[idx]     = A12_a[idx]     = A12_c[idx]     = A12_f[idx]      = 0.0;
    A13_p[idx]     = A13_a[idx]     = A13_c[idx]     = A13_f[idx]      = 0.0;
    A22_p[idx]     = A22_a[idx]     = A22_c[idx]     = A22_f[idx]      = 0.0;
    A23_p[idx]     = A23_a[idx]     = A23_c[idx]     = A23_f[idx]      = 0.0;
    A33_p[idx]     = A33_a[idx]     = A33_c[idx]     = A33_f[idx]      = 0.0;

    K_p[idx]       = K_a[idx]       = K_c[idx]       = K_f[idx]        = 0.0;

    Gamma1_p[idx]  = Gamma1_a[idx]  = Gamma1_c[idx]  = Gamma1_f[idx]   = 0.0;
    Gamma2_p[idx]  = Gamma2_a[idx]  = Gamma2_c[idx]  = Gamma2_f[idx]   = 0.0;
    Gamma3_p[idx]  = Gamma3_a[idx]  = Gamma3_c[idx]  = Gamma3_f[idx]   = 0.0;

    #if Z4c_DAMPING > 0
      theta_p[idx]  = theta_a[idx]  = theta_c[idx]  = theta_f[idx]   = 0.0;
      Z1_p[idx]  = Z1_a[idx]  = Z1_c[idx]  = Z1_f[idx]   = 0.0;
      Z2_p[idx]  = Z2_a[idx]  = Z2_c[idx]  = Z2_f[idx]   = 0.0;
      Z3_p[idx]  = Z3_a[idx]  = Z3_c[idx]  = Z3_f[idx]   = 0.0;
    #endif

    dk0_slice_phi_a[idx] = 0.0;
    KDx_a[idx] = 0.0;
    KDy_a[idx] = 0.0;
    KDz_a[idx] = 0.0;

    r_a[idx]        = 0.0;
    S_a[idx]        = 0.0;
    S1_a[idx]       = 0.0;
    S2_a[idx]       = 0.0;
    S3_a[idx]       = 0.0;
    STF11_a[idx]      = 0.0;
    STF12_a[idx]      = 0.0;
    STF13_a[idx]      = 0.0;
    STF22_a[idx]      = 0.0;
    STF23_a[idx]      = 0.0;
    STF33_a[idx]      = 0.0;
  }
}

/* Set local source term values */
void BSSN::set_source_vals(BSSNData *paq)
{
  idx_t idx = paq->idx;

  paq->rho = r_a[idx];
  paq->S = S_a[idx];

  paq->S1 = S1_a[idx];
  paq->S2 = S2_a[idx];
  paq->S3 = S3_a[idx];

  paq->STF11 = STF11_a[idx];
  paq->STF12 = STF12_a[idx];
  paq->STF13 = STF13_a[idx];
  paq->STF22 = STF22_a[idx];
  paq->STF23 = STF23_a[idx];
  paq->STF33 = STF33_a[idx];
}


/* set current local field values */
void BSSN::set_local_vals(BSSNData *paq)
{
  SET_LOCAL_INDEXES;

  // Pull out & store the metric at a point and adjacent points
  SET_LOCAL_VALUES_P(gamma11);
  SET_LOCAL_VALUES_P(gamma12);
  SET_LOCAL_VALUES_P(gamma13);
  SET_LOCAL_VALUES_P(gamma22);
  SET_LOCAL_VALUES_P(gamma23);
  SET_LOCAL_VALUES_P(gamma33);

  SET_LOCAL_VALUES_P(gammai11);
  SET_LOCAL_VALUES_P(gammai12);
  SET_LOCAL_VALUES_P(gammai13);
  SET_LOCAL_VALUES_P(gammai22);
  SET_LOCAL_VALUES_P(gammai23);
  SET_LOCAL_VALUES_P(gammai33);

  // Pull out values of quantities at a single point
  SET_LOCAL_VALUES_P(phi);
  SET_LOCAL_VALUES_P(K);
  SET_LOCAL_VALUES_P(Gamma1);
  SET_LOCAL_VALUES_P(Gamma2);
  SET_LOCAL_VALUES_P(Gamma3);
  SET_LOCAL_VALUES_P(A11);
  SET_LOCAL_VALUES_P(A12);
  SET_LOCAL_VALUES_P(A13);
  SET_LOCAL_VALUES_P(A22);
  SET_LOCAL_VALUES_P(A23);
  SET_LOCAL_VALUES_P(A33);
}



/*
******************************************************************************
Calculate independent quantities for later use (minimize # times calc'd)
******************************************************************************
*/

void BSSN::calculate_Acont(BSSNData *paq)
{
  // A^ij is calculated from A_ij by raising wrt. the conformal metric
  BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_ACONT)

  // calculate A_ij A^ij term
  AijAij_a[paq->idx] = paq->Acont11*paq->A11 + paq->Acont22*paq->A22 + paq->Acont33*paq->A33
      + 2.0*(paq->Acont12*paq->A12 + paq->Acont13*paq->A13 + paq->Acont23*paq->A23);
  paq->AijAij = AijAij_a[paq->idx];
}

/* Calculate metric derivatives */
void BSSN::calculate_dgamma(BSSNData *paq)
{
  BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_DGAMMA)
}

void BSSN::calculate_ddgamma(BSSNData *paq)
{
  BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_DIDJGAMMA_PERMS)
}

/* Calculate metric derivatives */
void BSSN::calculate_dgammai(BSSNData *paq)
{
  BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_DGAMMAI);
}

void BSSN::calculate_dphi(BSSNData *paq)
{
  // normal derivatives of phi
  paq->d1phi = derivative(paq->i, paq->j, paq->k, 1, phi_a);
  paq->d2phi = derivative(paq->i, paq->j, paq->k, 2, phi_a);
  paq->d3phi = derivative(paq->i, paq->j, paq->k, 3, phi_a);
}

void BSSN::calculate_dK(BSSNData *paq)
{
  // normal derivatives of phi
  paq->d1K = derivative(paq->i, paq->j, paq->k, 1, K_a);
  paq->d2K = derivative(paq->i, paq->j, paq->k, 2, K_a);
  paq->d3K = derivative(paq->i, paq->j, paq->k, 3, K_a);
}

#if Z4c_DAMPING > 0
void BSSN::calculate_dtheta(BSSNData *paq)
{
  // normal derivatives of phi
  paq->d1theta = derivative(paq->i, paq->j, paq->k, 1, theta_a);
  paq->d2theta = derivative(paq->i, paq->j, paq->k, 2, theta_a);
  paq->d3theta = derivative(paq->i, paq->j, paq->k, 3, theta_a);
}
#endif



/*
******************************************************************************
"dependent" quantities (depend on previously calc'd vals)
******************************************************************************
*/

void BSSN::calculate_conformal_christoffels(BSSNData *paq)
{
  // christoffel symbols: \Gamma^i_{jk} = Gijk
  BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_CHRISTOFFEL)
  // "lowered" christoffel symbols: \Gamma_{ijk} = GLijk
  BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_CHRISTOFFEL_LOWER)
}

void BSSN::calculateDDphi(BSSNData *paq)
{
  idx_t i = paq->i;
  idx_t j = paq->j;
  idx_t k = paq->k;

  // double covariant derivatives, using unitary metric
  paq->D1D1phi = double_derivative(i, j, k, 1, 1, phi_a) - (paq->G111*paq->d1phi + paq->G211*paq->d2phi + paq->G311*paq->d3phi);
  paq->D2D2phi = double_derivative(i, j, k, 2, 2, phi_a) - (paq->G122*paq->d1phi + paq->G222*paq->d2phi + paq->G322*paq->d3phi);
  paq->D3D3phi = double_derivative(i, j, k, 3, 3, phi_a) - (paq->G133*paq->d1phi + paq->G233*paq->d2phi + paq->G333*paq->d3phi);

  paq->D1D2phi = double_derivative(i, j, k, 1, 2, phi_a) - (paq->G112*paq->d1phi + paq->G212*paq->d2phi + paq->G312*paq->d3phi);
  paq->D1D3phi = double_derivative(i, j, k, 1, 3, phi_a) - (paq->G113*paq->d1phi + paq->G213*paq->d2phi + paq->G313*paq->d3phi);
  paq->D2D3phi = double_derivative(i, j, k, 2, 3, phi_a) - (paq->G123*paq->d1phi + paq->G223*paq->d2phi + paq->G323*paq->d3phi);  
}

/* Calculate trace-free ricci tensor components */
void BSSN::calculateRicciTF(BSSNData *paq)
{
  // unitary pieces
  /*BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_RICCITF_UNITARY)*/
  // may be more accurate but computationally intensive:
   BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_RICCITF_UNITARY_ALT) 

  paq->Uricci11 = paq->ricciTF11;
  paq->Uricci12 = paq->ricciTF12;
  paq->Uricci13 = paq->ricciTF13;
  paq->Uricci22 = paq->ricciTF22;
  paq->Uricci23 = paq->ricciTF23;
  paq->Uricci33 = paq->ricciTF33;

  /* calculate unitary Ricci scalar at this point; ricciTF isn't actually TF yet. */
  paq->unitRicci = paq->ricciTF11*paq->gammai11 + paq->ricciTF22*paq->gammai22 + paq->ricciTF33*paq->gammai33
            + 2.0*(paq->ricciTF12*paq->gammai12 + paq->ricciTF13*paq->gammai13 + paq->ricciTF23*paq->gammai23);

  real_t expression = (
    paq->gammai11*(paq->D1D1phi + 2.0*paq->d1phi*paq->d1phi)
    + paq->gammai22*(paq->D2D2phi + 2.0*paq->d2phi*paq->d2phi)
    + paq->gammai33*(paq->D3D3phi + 2.0*paq->d3phi*paq->d3phi)
    + 2.0*(
      paq->gammai12*(paq->D1D2phi + 2.0*paq->d1phi*paq->d2phi)
      + paq->gammai13*(paq->D1D3phi + 2.0*paq->d1phi*paq->d3phi)
      + paq->gammai23*(paq->D2D3phi + 2.0*paq->d2phi*paq->d3phi)
    )
  );

  /* phi-piece */
  paq->ricciTF11 += -2.0*( paq->D1D1phi - 2.0*paq->d1phi*paq->d1phi + paq->gamma11*(expression) );
  paq->ricciTF12 += -2.0*( paq->D1D2phi - 2.0*paq->d1phi*paq->d2phi + paq->gamma12*(expression) );
  paq->ricciTF13 += -2.0*( paq->D1D3phi - 2.0*paq->d1phi*paq->d3phi + paq->gamma13*(expression) );
  paq->ricciTF22 += -2.0*( paq->D2D2phi - 2.0*paq->d2phi*paq->d2phi + paq->gamma22*(expression) );
  paq->ricciTF23 += -2.0*( paq->D2D3phi - 2.0*paq->d2phi*paq->d3phi + paq->gamma23*(expression) );
  paq->ricciTF33 += -2.0*( paq->D3D3phi - 2.0*paq->d3phi*paq->d3phi + paq->gamma33*(expression) );

  /* calculate Ricci scalar at this point; ricciTF isn't TF at this point */
  paq->ricci = paq->ricciTF11*paq->gammai11 + paq->ricciTF22*paq->gammai22 + paq->ricciTF33*paq->gammai33
            + 2.0*(paq->ricciTF12*paq->gammai12 + paq->ricciTF13*paq->gammai13 + paq->ricciTF23*paq->gammai23);
  paq->ricci *= exp(-4.0*paq->phi);
  /* store ricci scalar here too. */
  ricci_a[paq->idx] = paq->ricci;

  /* remove trace. Note that \bar{gamma}_{ij}*\bar{gamma}^{kl}R_{kl} = (unbarred gammas). */
  paq->trace = paq->gammai11*paq->ricciTF11 + paq->gammai22*paq->ricciTF22 + paq->gammai33*paq->ricciTF33
      + 2.0*(paq->gammai12*paq->ricciTF12 + paq->gammai13*paq->ricciTF13 + paq->gammai23*paq->ricciTF23);
  paq->ricciTF11 -= (1.0/3.0)*paq->gamma11*paq->trace;
  paq->ricciTF12 -= (1.0/3.0)*paq->gamma12*paq->trace;
  paq->ricciTF13 -= (1.0/3.0)*paq->gamma13*paq->trace;
  paq->ricciTF22 -= (1.0/3.0)*paq->gamma22*paq->trace;
  paq->ricciTF23 -= (1.0/3.0)*paq->gamma23*paq->trace;
  paq->ricciTF33 -= (1.0/3.0)*paq->gamma33*paq->trace;
}



/*
******************************************************************************
(optional) Calculations of additional quantities
******************************************************************************
*/

void BSSN::set_KillingDelta(idx_t i, idx_t j, idx_t k, BSSNData *paq)
{
  idx_t idx = NP_INDEX(i,j,k);
  // Delta_mu^mu = 2 Grad_mu K^mu_i   (in i direction)
  //             = 2 _^(3)Gamma^j_ij
  //             = 2 ( \bar{Gamma}^j_ij + 6 d_i \phi )
  KDx_a[idx] = 2.0*( paq->G111 + paq->G221 + paq->G331 + 6.0*paq->d1phi );
  KDy_a[idx] = 2.0*( paq->G112 + paq->G222 + paq->G332 + 6.0*paq->d2phi );
  KDz_a[idx] = 2.0*( paq->G113 + paq->G223 + paq->G333 + 6.0*paq->d3phi );
}

void BSSN::set_full_metric(BSSNData *paq)
{
  SET_M00();

  SET_M0I(1); SET_M0I(2); SET_M0I(3);

  SET_MIJ(1, 1); SET_MIJ(2, 2); SET_MIJ(3, 3);
  SET_MIJ(1, 2); SET_MIJ(1, 3); SET_MIJ(2, 3);


  SET_Mi00();

  SET_Mi0I(1); SET_Mi0I(2); SET_Mi0I(3);

  SET_MiIJ(1, 1); SET_MiIJ(2, 2); SET_MiIJ(3, 3);
  SET_MiIJ(1, 2); SET_MiIJ(1, 3); SET_MiIJ(2, 3);
}

void BSSN::set_full_metric_der(BSSNData *paq)
{
  SET_DKM0I(1, 1); SET_DKM0I(2, 1); SET_DKM0I(3, 1);
  SET_DKM0I(1, 2); SET_DKM0I(2, 2); SET_DKM0I(3, 2);
  SET_DKM0I(1, 3); SET_DKM0I(2, 3); SET_DKM0I(3, 3);

  SET_DKMIJ(1, 1, 1); SET_DKMIJ(1, 2, 2); SET_DKMIJ(1, 3, 3);
  SET_DKMIJ(1, 1, 2); SET_DKMIJ(1, 1, 3); SET_DKMIJ(1, 2, 3);
  SET_DKMIJ(2, 1, 1); SET_DKMIJ(2, 2, 2); SET_DKMIJ(2, 3, 3);
  SET_DKMIJ(2, 1, 2); SET_DKMIJ(2, 1, 3); SET_DKMIJ(2, 2, 3);
  SET_DKMIJ(3, 1, 1); SET_DKMIJ(3, 2, 2); SET_DKMIJ(3, 3, 3);
  SET_DKMIJ(3, 1, 2); SET_DKMIJ(3, 1, 3); SET_DKMIJ(3, 2, 3);

  // macro depends on SET_DKMIJ first
  SET_DKM00(1); SET_DKM00(2); SET_DKM00(3);
}



/*
******************************************************************************
Evolution equation calculations
******************************************************************************
*/

real_t BSSN::ev_gamma11(BSSNData *paq) { return BSSN_DT_GAMMAIJ(1, 1) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->gamma11 - KO_dissipation_Q(paq->i, paq->j, paq->k, gamma11_a); }
real_t BSSN::ev_gamma12(BSSNData *paq) { return BSSN_DT_GAMMAIJ(1, 2) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->gamma12 - KO_dissipation_Q(paq->i, paq->j, paq->k, gamma12_a); }
real_t BSSN::ev_gamma13(BSSNData *paq) { return BSSN_DT_GAMMAIJ(1, 3) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->gamma13 - KO_dissipation_Q(paq->i, paq->j, paq->k, gamma13_a); }
real_t BSSN::ev_gamma22(BSSNData *paq) { return BSSN_DT_GAMMAIJ(2, 2) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->gamma22 - KO_dissipation_Q(paq->i, paq->j, paq->k, gamma22_a); }
real_t BSSN::ev_gamma23(BSSNData *paq) { return BSSN_DT_GAMMAIJ(2, 3) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->gamma23 - KO_dissipation_Q(paq->i, paq->j, paq->k, gamma23_a); }
real_t BSSN::ev_gamma33(BSSNData *paq) { return BSSN_DT_GAMMAIJ(3, 3) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->gamma33 - KO_dissipation_Q(paq->i, paq->j, paq->k, gamma33_a); }

real_t BSSN::ev_A11(BSSNData *paq) { return BSSN_DT_AIJ(1, 1) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A11*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, A11_a); }
real_t BSSN::ev_A12(BSSNData *paq) { return BSSN_DT_AIJ(1, 2) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A12*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, A12_a); }
real_t BSSN::ev_A13(BSSNData *paq) { return BSSN_DT_AIJ(1, 3) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A13*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, A13_a); }
real_t BSSN::ev_A22(BSSNData *paq) { return BSSN_DT_AIJ(2, 2) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A22*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, A22_a); }
real_t BSSN::ev_A23(BSSNData *paq) { return BSSN_DT_AIJ(2, 3) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A23*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, A23_a); }
real_t BSSN::ev_A33(BSSNData *paq) { return BSSN_DT_AIJ(3, 3) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A33*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, A33_a); }

real_t BSSN::ev_Gamma1(BSSNData *paq) { return BSSN_DT_GAMMAI(1) - KO_dissipation_Q(paq->i, paq->j, paq->k, Gamma1_a); }
real_t BSSN::ev_Gamma2(BSSNData *paq) { return BSSN_DT_GAMMAI(2) - KO_dissipation_Q(paq->i, paq->j, paq->k, Gamma2_a); }
real_t BSSN::ev_Gamma3(BSSNData *paq) { return BSSN_DT_GAMMAI(3) - KO_dissipation_Q(paq->i, paq->j, paq->k, Gamma3_a); }

real_t BSSN::ev_K(BSSNData *paq)
{
  return (
    pw2(paq->K + 2.0*paq->theta)/3.0 + paq->AijAij
    + 4.0*PI*(paq->rho + paq->S)
    - 1.0*JM_K_DAMPING_AMPLITUDE*paq->H*exp(-5.0*paq->phi)
    + Z4c_K1_DAMPING_AMPLITUDE*(1.0 - Z4c_K2_DAMPING_AMPLITUDE)*paq->theta
    - KO_dissipation_Q(paq->i, paq->j, paq->k, K_a)
  );
}

real_t BSSN::ev_phi(BSSNData *paq)
{
  return (
    0.1*BS_H_DAMPING_AMPLITUDE*dt*paq->H
    -1.0/6.0*(
      paq->K + 2.0*paq->theta
    )
    - KO_dissipation_Q(paq->i, paq->j, paq->k, phi_a)
  );
}


#if Z4c_DAMPING > 0
real_t BSSN::ev_theta(BSSNData *paq)
{
  return (
    0.5*paq->H
    - Z4c_K1_DAMPING_AMPLITUDE*(2.0 + Z4c_K2_DAMPING_AMPLITUDE)*paq->theta
  ) - KO_dissipation_Q(paq->i, paq->j, paq->k, theta_a);
}

real_t BSSN::ev_Z1(BSSNData *paq)
{
  return (
    BSSN_MI(1)
    + derivative(paq->i, paq->j, paq->k, 1, theta_a)
    - Z4c_K1_DAMPING_AMPLITUDE*paq->Z1
  ) - KO_dissipation_Q(paq->i, paq->j, paq->k, Z1_a);
}

real_t BSSN::ev_Z2(BSSNData *paq)
{
  return (
    BSSN_MI(2)
    + derivative(paq->i, paq->j, paq->k, 2, theta_a)
    - Z4c_K1_DAMPING_AMPLITUDE*paq->Z2
  ) - KO_dissipation_Q(paq->i, paq->j, paq->k, Z2_a);
}

real_t BSSN::ev_Z3(BSSNData *paq)
{
  return (
    BSSN_MI(3)
    + derivative(paq->i, paq->j, paq->k, 3, theta_a)
    - Z4c_K1_DAMPING_AMPLITUDE*paq->Z3
  ) - KO_dissipation_Q(paq->i, paq->j, paq->k, Z3_a);
}
#endif


/*
******************************************************************************
Constraint violtion calculations
******************************************************************************
*/

real_t BSSN::hamiltonianConstraintMagMean()
{
  idx_t i, j, k;
  real_t mean_hamiltonian_constraint = 0.0;

  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:mean_hamiltonian_constraint)
  LOOP3(i,j,k)
  {
    BSSNData b_paq = {0};
    set_paq_values(i, j, k, &b_paq);
    real_t ham = hamiltonianConstraintCalc(&b_paq);
    real_t ham2 = hamiltonianConstraintScale(&b_paq);
    real_t violation_fraction = ham/ham2;

    mean_hamiltonian_constraint += violation_fraction;
  }

  return mean_hamiltonian_constraint/POINTS;
}

real_t BSSN::hamiltonianConstraintMagStDev(real_t mean)
{
  idx_t i, j, k;
  real_t stdev_hamiltonian_constraint = 0.0;

  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:stdev_hamiltonian_constraint)
  LOOP3(i,j,k)
  {
    BSSNData b_paq = {0};
    set_paq_values(i, j, k, &b_paq);
    real_t ham = hamiltonianConstraintCalc(&b_paq);
    real_t ham2 = hamiltonianConstraintScale(&b_paq);
    real_t violation_fraction = ham/ham2;
    stdev_hamiltonian_constraint += pw2(violation_fraction - mean);
  }

  return sqrt(stdev_hamiltonian_constraint/(POINTS-1.0));
}

real_t BSSN::hamiltonianConstraintMagMax()
{
  idx_t i, j, k;
  real_t max_hamiltonian_constraint = 0.0;

  #pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    BSSNData b_paq = {0};
    set_paq_values(i, j, k, &b_paq);
    real_t ham = hamiltonianConstraintCalc(&b_paq);
    real_t ham2 = hamiltonianConstraintScale(&b_paq);
    real_t violation_fraction = ham/ham2;

    #pragma omp critical
    {
      if(fabs(violation_fraction) > max_hamiltonian_constraint)
      {
        max_hamiltonian_constraint = fabs(violation_fraction);
      }
    }
  }

  return max_hamiltonian_constraint;
}

real_t BSSN::hamiltonianConstraintCalc(BSSNData *paq)
{
  return -exp(5.0*paq->phi)/8.0*(paq->ricci + 2.0/3.0*pw2(paq->K + 2.0*paq->theta) - paq->AijAij - 16.0*PI*paq->rho);
}

real_t BSSN::hamiltonianConstraintScale(BSSNData *paq)
{
  // sqrt sum of sq. of terms for appx. mag / scale
  return (exp(5.0*paq->phi)/8.0)*sqrt( pw2(paq->ricci) + pw2(paq->AijAij) + pw2(2.0/3.0*pw2(paq->K + 2.0*paq->theta)) + pw2(16.0*PI*paq->rho) );
}

real_t BSSN::metricConstraintTotalMag()
{
  idx_t i, j, k;
  real_t constraint_mag = 0.0;

  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:constraint_mag)
  LOOP3(i,j,k)
  {
    BSSNData b_paq = {0};
    set_paq_values(i, j, k, &b_paq);

    constraint_mag += fabs(
      -1.0 +
      b_paq.gamma11*b_paq.gamma22*b_paq.gamma33 + 2.0*b_paq.gamma12*b_paq.gamma13*b_paq.gamma23
      - b_paq.gamma11*pw2(b_paq.gamma23) - b_paq.gamma22*pw2(b_paq.gamma13) - b_paq.gamma33*pw2(b_paq.gamma12)
    );
  }

  return constraint_mag;
}

real_t BSSN::momentumConstraintMagMean()
{
  idx_t i, j, k;
  real_t momentum_constraint_mag = 0.0;

  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:momentum_constraint_mag)
  LOOP3(i,j,k)
  {
    BSSNData b_paq = {0};
    set_paq_values(i, j, k, &b_paq);

    real_t Mi2 = sqrt(
      pw2(momentumConstraintCalc(&b_paq, 1))
      + pw2(momentumConstraintCalc(&b_paq, 2))
      + pw2(momentumConstraintCalc(&b_paq, 3))
    );
    real_t MiScale2 = sqrt(
      pw2(momentumConstraintScale(&b_paq, 1))
      + pw2(momentumConstraintScale(&b_paq, 2))
      + pw2(momentumConstraintScale(&b_paq, 3))
    );
    momentum_constraint_mag += Mi2/MiScale2;
  }

  return momentum_constraint_mag/POINTS;
}

real_t BSSN::momentumConstraintMagStDev(real_t mean)
{
  idx_t i, j, k;
  real_t stdev_momentum_constraint = 0.0;

  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:stdev_momentum_constraint)
  LOOP3(i,j,k)
  {
    BSSNData b_paq = {0};
    set_paq_values(i, j, k, &b_paq);

    real_t momentum_constraint_mag = sqrt(
      pw2(momentumConstraintCalc(&b_paq, 1))
      + pw2(momentumConstraintCalc(&b_paq, 2))
      + pw2(momentumConstraintCalc(&b_paq, 3))
    ) / sqrt(
      pw2(momentumConstraintScale(&b_paq, 1))
      + pw2(momentumConstraintScale(&b_paq, 2))
      + pw2(momentumConstraintScale(&b_paq, 3))
    );

    stdev_momentum_constraint += pw2(momentum_constraint_mag - mean);
  }

  return sqrt(stdev_momentum_constraint/(POINTS-1.0));
}

real_t BSSN::momentumConstraintMagMax()
{
  idx_t i, j, k;
  real_t momentum_constraint_max = 0.0;

  #pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    BSSNData b_paq = {0};
    set_paq_values(i, j, k, &b_paq);

    real_t Mi2 = sqrt(
      pw2(momentumConstraintCalc(&b_paq, 1))
      + pw2(momentumConstraintCalc(&b_paq, 2))
      + pw2(momentumConstraintCalc(&b_paq, 3))
    );
    real_t MiScale2 = sqrt(
      pw2(momentumConstraintScale(&b_paq, 1))
      + pw2(momentumConstraintScale(&b_paq, 2))
      + pw2(momentumConstraintScale(&b_paq, 3))
    );
    
    #pragma omp critical
    {
      if(Mi2/MiScale2 > momentum_constraint_max) {
        momentum_constraint_max = Mi2/MiScale2;
      }
    }
  }

  return momentum_constraint_max;
}

real_t BSSN::momentumConstraintCalc(BSSNData *paq, idx_t d)
{
  // needs paq vals and aijaij calc'd first
  switch(d)
  {
    case 1:
      return BSSN_MI(1);
    case 2:
      return BSSN_MI(2);
    case 3:
      return BSSN_MI(3);
  }

  /* xxx */
  throw -1;
  return 0;
}

real_t BSSN::momentumConstraintScale(BSSNData *paq, idx_t d)
{
  // needs paq vals and aijaij calc'd first
  switch(d)
  {
    case 1:
      return BSSN_MI_SCALE(1);
    case 2:
      return BSSN_MI_SCALE(2);
    case 3:
      return BSSN_MI_SCALE(3);
  }

  /* xxx */
  throw -1;
  return 0;
}

}
