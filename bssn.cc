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
  BSSN_APPLY_TO_GEN1_EXTRAS(GEN1_ARRAY_ALLOC)
  BSSN_APPLY_TO_GEN1_EXTRAS(GEN1_ARRAY_ADDMAP)

  // FRW reference integrator
  frw = new FRW<real_t> (0.0, 0.0);
}

BSSN::~BSSN()
{
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_DELETE)
  BSSN_APPLY_TO_SOURCES(GEN1_ARRAY_DELETE)
  BSSN_APPLY_TO_GEN1_EXTRAS(GEN1_ARRAY_DELETE)
}

void BSSN::set_gammai_values(idx_t i, idx_t j, idx_t k, BSSNData *paq)
{
  // Compute the inverse metric algebraically at each point
  // assumes det(gamma) = 1
  idx_t idx = NP_INDEX(i,j,k);
  paq->gammai11 = 1.0 + DIFFgamma22_a[idx] + DIFFgamma33_a[idx] - pw2(DIFFgamma23_a[idx]) + DIFFgamma22_a[idx]*DIFFgamma33_a[idx];
  paq->gammai22 = 1.0 + DIFFgamma11_a[idx] + DIFFgamma33_a[idx] - pw2(DIFFgamma13_a[idx]) + DIFFgamma11_a[idx]*DIFFgamma33_a[idx];
  paq->gammai33 = 1.0 + DIFFgamma11_a[idx] + DIFFgamma22_a[idx] - pw2(DIFFgamma12_a[idx]) + DIFFgamma11_a[idx]*DIFFgamma22_a[idx];
  paq->gammai12 = DIFFgamma13_a[idx]*DIFFgamma23_a[idx] - DIFFgamma12_a[idx]*(1.0 + DIFFgamma33_a[idx]);
  paq->gammai13 = DIFFgamma12_a[idx]*DIFFgamma23_a[idx] - DIFFgamma13_a[idx]*(1.0 + DIFFgamma22_a[idx]);
  paq->gammai23 = DIFFgamma12_a[idx]*DIFFgamma13_a[idx] - DIFFgamma23_a[idx]*(1.0 + DIFFgamma11_a[idx]);
}

void BSSN::set_DIFFgamma_Aij_norm()
{
  idx_t i, j, k;

  /* This potentially breaks conservation of trace:
   * need to come up with something else to preserve both
   * trace + determinant constraints.
   */
  #pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i, j, k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    // 1 - det(1 + DiffGamma)
    real_t one_minus_det_gamma = -1.0*(
      DIFFgamma11_a[idx] + DIFFgamma22_a[idx] + DIFFgamma33_a[idx]
      - pw2(DIFFgamma12_a[idx]) - pw2(DIFFgamma13_a[idx]) - pw2(DIFFgamma23_a[idx])
      + DIFFgamma11_a[idx]*DIFFgamma22_a[idx] + DIFFgamma11_a[idx]*DIFFgamma33_a[idx] + DIFFgamma22_a[idx]*DIFFgamma33_a[idx]
      - pw2(DIFFgamma23_a[idx])*DIFFgamma11_a[idx] - pw2(DIFFgamma13_a[idx])*DIFFgamma22_a[idx] - pw2(DIFFgamma12_a[idx])*DIFFgamma33_a[idx]
      + 2.0*DIFFgamma12_a[idx]*DIFFgamma13_a[idx]*DIFFgamma23_a[idx] + DIFFgamma11_a[idx]*DIFFgamma22_a[idx]*DIFFgamma33_a[idx]
    );

    // accurately compute 1 - det(g)^(1/3), without roundoff error
    // = -( det(g)^(1/3) - 1 )
    // = -( exp{log[det(g)^(1/3)]} - 1 )
    // = -( expm1{log[det(g)]/3} )
    // = -expm1{log1p[-one_minus_det_gamma]/3.0}
    real_t one_minus_det_gamma_thirdpow = -1.0*expm1(log1p(-1.0*one_minus_det_gamma)/3.0);

    // Perform the equivalent of re-scaling the conformal metric so det(gamma) = 1
    // gamma -> gamma / det(gamma)^(1/3)
    // DIFFgamma -> (delta + DiffGamma) / det(gamma)^(1/3) - delta
    //            = ( DiffGamma + delta*[1 - det(gamma)^(1/3)] ) / ( 1 - [1 - det(1 + DiffGamma)^1/3] )
    DIFFgamma11_a[idx] = (DIFFgamma11_a[idx] + one_minus_det_gamma_thirdpow) / (1.0 - one_minus_det_gamma_thirdpow);
    DIFFgamma22_a[idx] = (DIFFgamma22_a[idx] + one_minus_det_gamma_thirdpow) / (1.0 - one_minus_det_gamma_thirdpow);
    DIFFgamma33_a[idx] = (DIFFgamma33_a[idx] + one_minus_det_gamma_thirdpow) / (1.0 - one_minus_det_gamma_thirdpow);
    DIFFgamma12_a[idx] = (DIFFgamma12_a[idx]) / (1.0 - one_minus_det_gamma_thirdpow);
    DIFFgamma13_a[idx] = (DIFFgamma13_a[idx]) / (1.0 - one_minus_det_gamma_thirdpow);
    DIFFgamma23_a[idx] = (DIFFgamma23_a[idx]) / (1.0 - one_minus_det_gamma_thirdpow);

    // re-scale A_ij / ensure it is trace-free
    // need inverse gamma for finding Tr(A)
    real_t gammai11 = 1.0 + DIFFgamma22_a[idx] + DIFFgamma33_a[idx] - pw2(DIFFgamma23_a[idx]) + DIFFgamma22_a[idx]*DIFFgamma33_a[idx];
    real_t gammai22 = 1.0 + DIFFgamma11_a[idx] + DIFFgamma33_a[idx] - pw2(DIFFgamma13_a[idx]) + DIFFgamma11_a[idx]*DIFFgamma33_a[idx];
    real_t gammai33 = 1.0 + DIFFgamma11_a[idx] + DIFFgamma22_a[idx] - pw2(DIFFgamma12_a[idx]) + DIFFgamma11_a[idx]*DIFFgamma22_a[idx];
    real_t gammai12 = DIFFgamma13_a[idx]*DIFFgamma23_a[idx] - DIFFgamma12_a[idx]*(1.0 + DIFFgamma33_a[idx]);
    real_t gammai13 = DIFFgamma12_a[idx]*DIFFgamma23_a[idx] - DIFFgamma13_a[idx]*(1.0 + DIFFgamma22_a[idx]);
    real_t gammai23 = DIFFgamma12_a[idx]*DIFFgamma13_a[idx] - DIFFgamma23_a[idx]*(1.0 + DIFFgamma11_a[idx]);
    real_t trA = gammai11*A11_a[idx] + gammai22*A22_a[idx] + gammai33*A33_a[idx]
      + 2.0*(gammai12*A12_a[idx] + gammai13*A13_a[idx] + gammai23*A23_a[idx]);
    // A_ij -> ( A_ij - 1/3 gamma_ij A )
    A11_a[idx] = ( A11_a[idx] - 1.0/3.0*(1.0 + DIFFgamma11_a[idx])*trA ) / (1.0 - one_minus_det_gamma_thirdpow);
    A22_a[idx] = ( A22_a[idx] - 1.0/3.0*(1.0 + DIFFgamma22_a[idx])*trA ) / (1.0 - one_minus_det_gamma_thirdpow);
    A33_a[idx] = ( A33_a[idx] - 1.0/3.0*(1.0 + DIFFgamma33_a[idx])*trA ) / (1.0 - one_minus_det_gamma_thirdpow);
    A12_a[idx] = ( A12_a[idx] - 1.0/3.0*DIFFgamma12_a[idx]*trA ) / (1.0 - one_minus_det_gamma_thirdpow);
    A13_a[idx] = ( A13_a[idx] - 1.0/3.0*DIFFgamma13_a[idx]*trA ) / (1.0 - one_minus_det_gamma_thirdpow);
    A23_a[idx] = ( A23_a[idx] - 1.0/3.0*DIFFgamma23_a[idx]*trA ) / (1.0 - one_minus_det_gamma_thirdpow);
  }
}

void BSSN::set_paq_values(idx_t i, idx_t j, idx_t k, BSSNData *paq)
{
  paq->i = i;
  paq->j = j;
  paq->k = k;
  paq->idx = NP_INDEX(i,j,k);

  // need to set FRW quantities first
  paq->phi_FRW = frw->get_phi();
  paq->K_FRW = frw->get_K();
  paq->rho_FRW = frw->get_rho();
  paq->S_FRW = frw->get_S();

  // draw data from cache
  set_local_vals(paq);
  set_gammai_values(i, j, k, paq);

  // non-DIFF quantities
  paq->phi      =   paq->DIFFphi + paq->phi_FRW;
  paq->K        =   paq->DIFFK + paq->K_FRW;
  paq->gamma11  =   paq->DIFFgamma11 + 1.0;
  paq->gamma12  =   paq->DIFFgamma12;
  paq->gamma13  =   paq->DIFFgamma13;
  paq->gamma22  =   paq->DIFFgamma22 + 1.0;
  paq->gamma23  =   paq->DIFFgamma23;
  paq->gamma33  =   paq->DIFFgamma33 + 1.0;
  paq->r        =   paq->DIFFr + paq->rho_FRW;
  paq->S        =   paq->DIFFS + paq->S_FRW;
  paq->alpha    =   paq->DIFFalpha + 1.0;

  // pre-compute re-used quantities
  // gammas & derivs first
  calculate_Acont(paq);
  calculate_dgamma(paq);
  calculate_ddgamma(paq);
  calculate_dalpha_dphi(paq);
  calculate_dK(paq);
  #if USE_Z4c_DAMPING
    calculate_dtheta(paq);
  #endif
  #if USE_BSSN_SHIFT
    calculate_dbeta(paq);
  #endif

  // Christoffels depend on metric & derivs.
  calculate_conformal_christoffels(paq);
  // DDw depend on christoffels, metric, and derivs
  calculateDDphi(paq);
  calculateDDalphaTF(paq);
  // Ricci depends on DDphi
  calculateRicciTF(paq);

  // H depends on AijAij and ricci arrays
  paq->H = hamiltonianConstraintCalc(paq->idx);
}

// Full RK step (More useful when not evolving the source simultaneously)
void BSSN::step(BSSNData *paq)
{
  K1Calc();
  K2Calc();
  K3Calc();
  K4Calc();
  stepTerm();
  // done!
}

void BSSN::regSwap_c_a()
{
  BSSN_SWAP_ARRAYS(_c, _a);
  #if NORMALIZE_GAMMAIJ_AIJ
    set_DIFFgamma_Aij_norm(); // norms _a register
  #endif
}

// Init _a register with _p values, _f with 0
void BSSN::stepInit()
{
  BSSN_COPY_ARRAYS(_p, _a);
  #if NORMALIZE_GAMMAIJ_AIJ
    set_DIFFgamma_Aij_norm(); // norms _a register
    BSSN_COPY_ARRAYS(_a, _p);
  #endif
  idx_t i, j, k;
  #pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i, j, k)
  {
    idx_t idx = NP_INDEX(i, j, k);
    BSSN_ZERO_ARRAYS(_f, idx);
  }
}


// First RK step: calculate and add k_1 coeff to _f array
void BSSN::K1Calc()
{
  BSSN_RK_PERFORM_KN_CALC(1);
  BSSN_SWAP_ARRAYS(_c, _a);
  frw->P1_step(dt);
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
void BSSN::K2Calc()
{
  BSSN_RK_PERFORM_KN_CALC(2);
  BSSN_SWAP_ARRAYS(_c, _a);
  frw->P2_step(dt);
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
void BSSN::K3Calc()
{
  BSSN_RK_PERFORM_KN_CALC(3);
  BSSN_SWAP_ARRAYS(_c, _a);
  frw->P3_step(dt);
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
void BSSN::K4Calc()
{
  BSSN_RK_PERFORM_KN_CALC(4);
  frw->RK_total_step(dt);
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
    BSSN_ZERO_SOURCES()
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
    BSSN_ZERO_ARRAYS(_p, idx)
    BSSN_ZERO_ARRAYS(_a, idx)
    BSSN_ZERO_ARRAYS(_c, idx)
    BSSN_ZERO_ARRAYS(_f, idx)

    BSSN_ZERO_GEN1_EXTRAS()
    BSSN_ZERO_SOURCES()
  }
}

/* set current local field values */
void BSSN::set_local_vals(BSSNData *paq)
{
  SET_LOCAL_INDEXES;

  // Pull out values of quantities at a single point
  BSSN_APPLY_TO_FIELDS(SET_LOCAL_VALUES_P);
  BSSN_APPLY_TO_GEN1_EXTRAS(SET_LOCAL_VALUES_P);
  BSSN_APPLY_TO_SOURCES(SET_LOCAL_VALUES_P);
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

void BSSN::calculate_dalpha_dphi(BSSNData *paq)
{
  // normal derivatives of phi
  paq->d1phi = derivative(paq->i, paq->j, paq->k, 1, DIFFphi_a);
  paq->d2phi = derivative(paq->i, paq->j, paq->k, 2, DIFFphi_a);
  paq->d3phi = derivative(paq->i, paq->j, paq->k, 3, DIFFphi_a);

  paq->d1d1phi = double_derivative(paq->i, paq->j, paq->k, 1, 1, DIFFphi_a);
  paq->d2d2phi = double_derivative(paq->i, paq->j, paq->k, 2, 2, DIFFphi_a);
  paq->d3d3phi = double_derivative(paq->i, paq->j, paq->k, 3, 3, DIFFphi_a);
  paq->d1d2phi = double_derivative(paq->i, paq->j, paq->k, 1, 2, DIFFphi_a);
  paq->d1d3phi = double_derivative(paq->i, paq->j, paq->k, 1, 3, DIFFphi_a);
  paq->d2d3phi = double_derivative(paq->i, paq->j, paq->k, 2, 3, DIFFphi_a);

  // normal derivatives of alpha
  paq->d1a = derivative(paq->i, paq->j, paq->k, 1, DIFFalpha_a);
  paq->d2a = derivative(paq->i, paq->j, paq->k, 2, DIFFalpha_a);
  paq->d3a = derivative(paq->i, paq->j, paq->k, 3, DIFFalpha_a);
}

void BSSN::calculate_dK(BSSNData *paq)
{
  // normal derivatives of phi
  paq->d1K = derivative(paq->i, paq->j, paq->k, 1, DIFFK_a);
  paq->d2K = derivative(paq->i, paq->j, paq->k, 2, DIFFK_a);
  paq->d3K = derivative(paq->i, paq->j, paq->k, 3, DIFFK_a);
}

#if USE_Z4c_DAMPING
void BSSN::calculate_dtheta(BSSNData *paq)
{
  // normal derivatives of phi
  paq->d1theta = derivative(paq->i, paq->j, paq->k, 1, theta_a);
  paq->d2theta = derivative(paq->i, paq->j, paq->k, 2, theta_a);
  paq->d3theta = derivative(paq->i, paq->j, paq->k, 3, theta_a);
}
#endif

#if USE_BSSN_SHIFT
void BSSN::calculate_dbeta(BSSNData *paq)
{
  paq->d1beta1 = derivative(paq->i, paq->j, paq->k, 1, beta1_a);
  paq->d1beta2 = derivative(paq->i, paq->j, paq->k, 1, beta2_a);
  paq->d1beta3 = derivative(paq->i, paq->j, paq->k, 1, beta3_a);
  paq->d2beta1 = derivative(paq->i, paq->j, paq->k, 2, beta1_a);
  paq->d2beta2 = derivative(paq->i, paq->j, paq->k, 2, beta2_a);
  paq->d2beta3 = derivative(paq->i, paq->j, paq->k, 2, beta3_a);
  paq->d3beta1 = derivative(paq->i, paq->j, paq->k, 3, beta1_a);
  paq->d3beta2 = derivative(paq->i, paq->j, paq->k, 3, beta2_a);
  paq->d3beta3 = derivative(paq->i, paq->j, paq->k, 3, beta3_a);
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

  paq->Gammad1 = paq->G111*paq->gammai11 + paq->G122*paq->gammai22 + paq->G133*paq->gammai33
    + 2.0*(paq->G112*paq->gammai12 + paq->G113*paq->gammai13 + paq->G123*paq->gammai23);
  paq->Gammad2 = paq->G211*paq->gammai11 + paq->G222*paq->gammai22 + paq->G233*paq->gammai33
    + 2.0*(paq->G212*paq->gammai12 + paq->G213*paq->gammai13 + paq->G223*paq->gammai23);
  paq->Gammad3 = paq->G311*paq->gammai11 + paq->G322*paq->gammai22 + paq->G333*paq->gammai33
    + 2.0*(paq->G312*paq->gammai12 + paq->G313*paq->gammai13 + paq->G323*paq->gammai23);
}

void BSSN::calculateDDphi(BSSNData *paq)
{
  idx_t i = paq->i;
  idx_t j = paq->j;
  idx_t k = paq->k;

  // double covariant derivatives, using unitary metric
  paq->D1D1phi = paq->d1d1phi - (paq->G111*paq->d1phi + paq->G211*paq->d2phi + paq->G311*paq->d3phi);
  paq->D2D2phi = paq->d2d2phi - (paq->G122*paq->d1phi + paq->G222*paq->d2phi + paq->G322*paq->d3phi);
  paq->D3D3phi = paq->d3d3phi - (paq->G133*paq->d1phi + paq->G233*paq->d2phi + paq->G333*paq->d3phi);

  paq->D1D2phi = paq->d1d2phi - (paq->G112*paq->d1phi + paq->G212*paq->d2phi + paq->G312*paq->d3phi);
  paq->D1D3phi = paq->d1d3phi - (paq->G113*paq->d1phi + paq->G213*paq->d2phi + paq->G313*paq->d3phi);
  paq->D2D3phi = paq->d2d3phi - (paq->G123*paq->d1phi + paq->G223*paq->d2phi + paq->G323*paq->d3phi);  
}

void BSSN::calculateDDalphaTF(BSSNData *paq)
{
  // double covariant derivatives - use non-unitary metric - extra pieces that depend on phi!
  // the gammaIldlphi are needed for the BSSN_CALCULATE_DIDJALPHA macro
  real_t gammai1ldlphi = paq->gammai11*paq->d1phi + paq->gammai12*paq->d2phi + paq->gammai13*paq->d3phi;
  real_t gammai2ldlphi = paq->gammai21*paq->d1phi + paq->gammai22*paq->d2phi + paq->gammai23*paq->d3phi;
  real_t gammai3ldlphi = paq->gammai31*paq->d1phi + paq->gammai32*paq->d2phi + paq->gammai33*paq->d3phi;
  // Calculates full (not trace-free) piece:
  BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_DIDJALPHA)

  // subtract trace (traced with full spatial metric but subtracted later)
  paq->DDaTR = paq->gammai11*paq->D1D1aTF + paq->gammai22*paq->D2D2aTF + paq->gammai33*paq->D3D3aTF
      + 2.0*(paq->gammai12*paq->D1D2aTF + paq->gammai13*paq->D1D3aTF + paq->gammai23*paq->D2D3aTF);
  paq->D1D1aTF -= (1.0/3.0)*paq->gamma11*paq->DDaTR;
  paq->D1D2aTF -= (1.0/3.0)*paq->gamma12*paq->DDaTR;
  paq->D1D3aTF -= (1.0/3.0)*paq->gamma13*paq->DDaTR;
  paq->D2D2aTF -= (1.0/3.0)*paq->gamma22*paq->DDaTR;
  paq->D2D3aTF -= (1.0/3.0)*paq->gamma23*paq->DDaTR;
  paq->D3D3aTF -= (1.0/3.0)*paq->gamma33*paq->DDaTR;

  // scale trace back (=> contracted with "real" metric)
  paq->DDaTR *= exp(-4.0*paq->phi);
}

/* Calculate trace-free ricci tensor components */
void BSSN::calculateRicciTF(BSSNData *paq)
{
  // unitary pieces
  BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_RICCI_UNITARY)

  paq->Uricci11 = paq->ricci11;
  paq->Uricci12 = paq->ricci12;
  paq->Uricci13 = paq->ricci13;
  paq->Uricci22 = paq->ricci22;
  paq->Uricci23 = paq->ricci23;
  paq->Uricci33 = paq->ricci33;

  /* calculate unitary Ricci scalar at this point. */
  paq->unitRicci = paq->ricci11*paq->gammai11 + paq->ricci22*paq->gammai22 + paq->ricci33*paq->gammai33
            + 2.0*(paq->ricci12*paq->gammai12 + paq->ricci13*paq->gammai13 + paq->ricci23*paq->gammai23);

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
  paq->ricci11 += -2.0*( paq->D1D1phi - 2.0*paq->d1phi*paq->d1phi + paq->gamma11*(expression) );
  paq->ricci12 += -2.0*( paq->D1D2phi - 2.0*paq->d1phi*paq->d2phi + paq->gamma12*(expression) );
  paq->ricci13 += -2.0*( paq->D1D3phi - 2.0*paq->d1phi*paq->d3phi + paq->gamma13*(expression) );
  paq->ricci22 += -2.0*( paq->D2D2phi - 2.0*paq->d2phi*paq->d2phi + paq->gamma22*(expression) );
  paq->ricci23 += -2.0*( paq->D2D3phi - 2.0*paq->d2phi*paq->d3phi + paq->gamma23*(expression) );
  paq->ricci33 += -2.0*( paq->D3D3phi - 2.0*paq->d3phi*paq->d3phi + paq->gamma33*(expression) );

  /* calculate full Ricci scalar at this point */
  paq->ricci = paq->ricci11*paq->gammai11 + paq->ricci22*paq->gammai22 + paq->ricci33*paq->gammai33
            + 2.0*(paq->ricci12*paq->gammai12 + paq->ricci13*paq->gammai13 + paq->ricci23*paq->gammai23);
  paq->ricci *= exp(-4.0*paq->phi);
  /* store ricci scalar here too. */
  ricci_a[paq->idx] = paq->ricci;

  /* remove trace. Note that \bar{gamma}_{ij}*\bar{gamma}^{kl}R_{kl} = (unbarred gammas). */
  paq->trace = paq->gammai11*paq->ricci11 + paq->gammai22*paq->ricci22 + paq->gammai33*paq->ricci33
      + 2.0*(paq->gammai12*paq->ricci12 + paq->gammai13*paq->ricci13 + paq->gammai23*paq->ricci23);
  paq->ricciTF11 = paq->ricci11 - (1.0/3.0)*paq->gamma11*paq->trace;
  paq->ricciTF12 = paq->ricci12 - (1.0/3.0)*paq->gamma12*paq->trace;
  paq->ricciTF13 = paq->ricci13 - (1.0/3.0)*paq->gamma13*paq->trace;
  paq->ricciTF22 = paq->ricci22 - (1.0/3.0)*paq->gamma22*paq->trace;
  paq->ricciTF23 = paq->ricci23 - (1.0/3.0)*paq->gamma23*paq->trace;
  paq->ricciTF33 = paq->ricci33 - (1.0/3.0)*paq->gamma33*paq->trace;
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

real_t BSSN::ev_DIFFgamma11(BSSNData *paq) { return BSSN_DT_DIFFGAMMAIJ(1, 1) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->DIFFgamma11 - KO_dissipation_Q(paq->i, paq->j, paq->k, DIFFgamma11_a); }
real_t BSSN::ev_DIFFgamma12(BSSNData *paq) { return BSSN_DT_DIFFGAMMAIJ(1, 2) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->DIFFgamma12 - KO_dissipation_Q(paq->i, paq->j, paq->k, DIFFgamma12_a); }
real_t BSSN::ev_DIFFgamma13(BSSNData *paq) { return BSSN_DT_DIFFGAMMAIJ(1, 3) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->DIFFgamma13 - KO_dissipation_Q(paq->i, paq->j, paq->k, DIFFgamma13_a); }
real_t BSSN::ev_DIFFgamma22(BSSNData *paq) { return BSSN_DT_DIFFGAMMAIJ(2, 2) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->DIFFgamma22 - KO_dissipation_Q(paq->i, paq->j, paq->k, DIFFgamma22_a); }
real_t BSSN::ev_DIFFgamma23(BSSNData *paq) { return BSSN_DT_DIFFGAMMAIJ(2, 3) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->DIFFgamma23 - KO_dissipation_Q(paq->i, paq->j, paq->k, DIFFgamma23_a); }
real_t BSSN::ev_DIFFgamma33(BSSNData *paq) { return BSSN_DT_DIFFGAMMAIJ(3, 3) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->DIFFgamma33 - KO_dissipation_Q(paq->i, paq->j, paq->k, DIFFgamma33_a); }

real_t BSSN::ev_A11(BSSNData *paq) { return BSSN_DT_AIJ(1, 1) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A11*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, A11_a); }
real_t BSSN::ev_A12(BSSNData *paq) { return BSSN_DT_AIJ(1, 2) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A12*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, A12_a); }
real_t BSSN::ev_A13(BSSNData *paq) { return BSSN_DT_AIJ(1, 3) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A13*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, A13_a); }
real_t BSSN::ev_A22(BSSNData *paq) { return BSSN_DT_AIJ(2, 2) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A22*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, A22_a); }
real_t BSSN::ev_A23(BSSNData *paq) { return BSSN_DT_AIJ(2, 3) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A23*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, A23_a); }
real_t BSSN::ev_A33(BSSNData *paq) { return BSSN_DT_AIJ(3, 3) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A33*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, A33_a); }

real_t BSSN::ev_Gamma1(BSSNData *paq) { return BSSN_DT_GAMMAI(1) - KO_dissipation_Q(paq->i, paq->j, paq->k, Gamma1_a); }
real_t BSSN::ev_Gamma2(BSSNData *paq) { return BSSN_DT_GAMMAI(2) - KO_dissipation_Q(paq->i, paq->j, paq->k, Gamma2_a); }
real_t BSSN::ev_Gamma3(BSSNData *paq) { return BSSN_DT_GAMMAI(3) - KO_dissipation_Q(paq->i, paq->j, paq->k, Gamma3_a); }

real_t BSSN::ev_DIFFK(BSSNData *paq)
{
  return (
    - paq->DDaTR
    + paq->alpha*(
        paq->AijAij
        + 1.0/3.0*(paq->DIFFK + 2.0*paq->theta)*(paq->DIFFK + 2.0*paq->theta + 2.0*paq->K_FRW)
    )
    + 4.0*PI*paq->alpha*(paq->DIFFr + paq->DIFFS)
    - paq->DIFFalpha*(
        1.0/3.0*pw2(paq->K_FRW)
        + 4.0*PI*(paq->rho_FRW + paq->S_FRW)
      )
    + paq->beta1*paq->d1K + paq->beta2*paq->d2K + paq->beta3*paq->d3K
    - 1.0*JM_K_DAMPING_AMPLITUDE*paq->H*exp(-5.0*paq->phi)
    + Z4c_K1_DAMPING_AMPLITUDE*(1.0 - Z4c_K2_DAMPING_AMPLITUDE)*paq->theta
    - KO_dissipation_Q(paq->i, paq->j, paq->k, DIFFK_a)
  );
}

real_t BSSN::ev_DIFFphi(BSSNData *paq)
{
  return (
    0.1*BS_H_DAMPING_AMPLITUDE*dt*paq->H
    -1.0/6.0*(
      paq->alpha*(paq->DIFFK + 2.0*paq->theta)
      - paq->DIFFalpha*paq->K_FRW
      - ( paq->d1beta1 + paq->d2beta2 + paq->d3beta3 )
    )
    + paq->beta1*paq->d1phi + paq->beta2*paq->d2phi + paq->beta3*paq->d3phi
    - KO_dissipation_Q(paq->i, paq->j, paq->k, DIFFphi_a)
  );
}

real_t BSSN::ev_DIFFalpha(BSSNData *paq)
{
  #if USE_HARMONIC_ALPHA
    return -1.0*pw2(paq->alpha)*paq->K - KO_dissipation_Q(paq->i, paq->j, paq->k, DIFFalpha_a);
  #endif

  #if USE_CONFORMAL_SYNC_ALPHA
    return -1.0/3.0*paq->alpha*paq->K_FRW;
  #endif

  return 0.0;
}

#if USE_Z4c_DAMPING
real_t BSSN::ev_theta(BSSNData *paq)
{
  return (
    0.5*paq->alpha*(
      paq->ricci + 2.0/3.0*pw2(paq->K + 2.0*paq->theta) - paq->AijAij - 16.0*PI*( paq->r )
    )
    - paq->alpha*Z4c_K1_DAMPING_AMPLITUDE*(2.0 + Z4c_K2_DAMPING_AMPLITUDE)*paq->theta
  ) - KO_dissipation_Q(paq->i, paq->j, paq->k, theta_a);
}
#endif

#if USE_BSSN_SHIFT
real_t BSSN::ev_beta1(BSSNData *paq)
{
  return 0.0;
}

real_t BSSN::ev_beta2(BSSNData *paq)
{
  return 0.0;
}

real_t BSSN::ev_beta3(BSSNData *paq)
{
  return 0.0;
}
#endif


/*
******************************************************************************
Constraint violtion calculations
******************************************************************************
*/

void BSSN::setHamiltonianConstraintCalcs(real_t H_values[7], bool reset_paq)
{
  idx_t i, j, k;
  // unscaled quantities
  real_t mean_H = 0.0;
  real_t stdev_H = 0.0;
  real_t max_H = 0.0;
  // the "scale"
  real_t mean_H_scale = 0.0;
  //scaled quantities
  real_t mean_H_scaled = 0.0;
  real_t stdev_H_scaled = 0.0;
  real_t max_H_scaled = 0.0;

  if(reset_paq)
  {
    #pragma omp parallel for default(shared) private(i, j, k)
    LOOP3(i,j,k)
    {
      BSSNData b_paq = {0};
      set_paq_values(i, j, k, &b_paq); // sets AijAij and Ricci too
    }
  }

  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:mean_H,mean_H_scale,mean_H_scaled)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    real_t H = hamiltonianConstraintCalc(idx);
    real_t H_scale = hamiltonianConstraintScale(idx);
    real_t H_scaled = H/H_scale;

    mean_H += H;
    mean_H_scale += H_scale;
    mean_H_scaled += H_scaled;

    #pragma omp critical
    {
      if(fabs(H) > max_H)
      {
        max_H = fabs(H);
      }
      if(fabs(H_scaled) > max_H_scaled)
      {
        max_H_scaled = fabs(H_scaled);
      }
    }
  }
  mean_H /= POINTS;
  mean_H_scaled /= POINTS;

  // stdev relies on mean calcs
  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:stdev_H,stdev_H_scaled)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    real_t H = hamiltonianConstraintCalc(idx);
    real_t H_scale = hamiltonianConstraintScale(idx);
    real_t H_scaled = H/H_scale;

    stdev_H += pw2(H - mean_H);
    stdev_H_scaled += pw2(H_scaled - mean_H_scaled);
  }

  stdev_H = sqrt(stdev_H/(POINTS-1.0));
  stdev_H_scaled = sqrt(stdev_H_scaled/(POINTS-1.0));

  H_values[0] = mean_H;
  H_values[1] = stdev_H;
  H_values[2] = max_H;
  H_values[3] = mean_H_scale;
  H_values[4] = mean_H_scaled;
  H_values[5] = stdev_H_scaled;
  H_values[6] = max_H_scaled;

  return;
}

real_t BSSN::hamiltonianConstraintCalc(idx_t idx)
{
  #if USE_Z4c_DAMPING
    real_t theta = theta_a[idx];
  #else
    real_t theta = 0.0;
  #endif

  #if USE_REFERENCE_FRW
    real_t K_FRW = frw->get_K();
    real_t phi_FRW = frw->get_phi();
    real_t rho_FRW = frw->get_rho();
    return -exp(5.0*(DIFFphi_a[idx] + phi_FRW))/8.0*(
      ricci_a[idx] + 2.0/3.0*pw2(K_FRW + DIFFK_a[idx] + 2.0*theta) - AijAij_a[idx] - 16.0*PI*(DIFFr_a[idx] + rho_FRW)
    );
  #else
    return -exp(5.0*DIFFphi_a[idx])/8.0*(
      ricci_a[idx] + 2.0/3.0*pw2(DIFFK_a[idx] + 2.0*theta) - AijAij_a[idx] - 16.0*PI*DIFFr_a[idx]
    );
  #endif
}

real_t BSSN::hamiltonianConstraintScale(idx_t idx)
{
  real_t K_FRW = frw->get_K();
  real_t phi_FRW = frw->get_phi();
  real_t rho_FRW = frw->get_rho();

  #if USE_Z4c_DAMPING
    real_t theta = theta_a[idx];
  #else
    real_t theta = 0.0;
  #endif

  // sqrt sum of sq. of terms for appx. mag / scale
  return (exp(5.0*(DIFFphi_a[idx] + phi_FRW))/8.0)*
    sqrt( pw2(ricci_a[idx]) + pw2(AijAij_a[idx]) + pw2(2.0/3.0*pw2(K_FRW + DIFFK_a[idx] + 2.0*theta)) + pw2(16.0*PI*(DIFFr_a[idx] + rho_FRW))
  );
}


void BSSN::setMomentumConstraintCalcs(real_t M_values[7])
{
  idx_t i, j, k;
  // unscaled quantities
  real_t mean_M = 0.0;
  real_t stdev_M = 0.0;
  real_t max_M = 0.0;
  // the "scale"
  real_t mean_M_scale = 0.0;
  //scaled quantities
  real_t mean_M_scaled = 0.0;
  real_t stdev_M_scaled = 0.0;
  real_t max_M_scaled = 0.0;

  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:mean_M,mean_M_scale,mean_M_scaled)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    BSSNData b_paq = {0};
    set_paq_values(i, j, k, &b_paq); // sets AijAij and Ricci too

    real_t M = sqrt(
      pw2(momentumConstraintCalc(&b_paq, 1))
      + pw2(momentumConstraintCalc(&b_paq, 2))
      + pw2(momentumConstraintCalc(&b_paq, 3))
    );
    real_t M_scale = sqrt(
      pw2(momentumConstraintScale(&b_paq, 1))
      + pw2(momentumConstraintScale(&b_paq, 2))
      + pw2(momentumConstraintScale(&b_paq, 3))
    );
    real_t M_scaled = M/M_scale;

    mean_M += M;
    mean_M_scale += M_scale;
    mean_M_scaled += M_scaled;

    #pragma omp critical
    {
      if(fabs(M) > max_M)
      {
        max_M = fabs(M);
      }
      if(fabs(M_scaled) > max_M_scaled)
      {
        max_M_scaled = fabs(M_scaled);
      }
    }
  }
  mean_M /= POINTS;
  mean_M_scaled /= POINTS;

  // stdev relies on mean calcs
  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:stdev_M,stdev_M_scaled)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    BSSNData b_paq = {0};
    set_paq_values(i, j, k, &b_paq); // sets AijAij and Ricci too

    real_t M = sqrt(
      pw2(momentumConstraintCalc(&b_paq, 1))
      + pw2(momentumConstraintCalc(&b_paq, 2))
      + pw2(momentumConstraintCalc(&b_paq, 3))
    );
    real_t M_scale = sqrt(
      pw2(momentumConstraintScale(&b_paq, 1))
      + pw2(momentumConstraintScale(&b_paq, 2))
      + pw2(momentumConstraintScale(&b_paq, 3))
    );
    real_t M_scaled = M/M_scale;

    stdev_M += pw2(M - mean_M);
    stdev_M_scaled += pw2(M_scaled - mean_M_scaled);
  }

  stdev_M = sqrt(stdev_M/(POINTS-1.0));
  stdev_M_scaled = sqrt(stdev_M_scaled/(POINTS-1.0));

  M_values[0] = mean_M;
  M_values[1] = stdev_M;
  M_values[2] = max_M;
  M_values[3] = mean_M_scale;
  M_values[4] = mean_M_scaled;
  M_values[5] = stdev_M_scaled;
  M_values[6] = max_M_scaled;

  return;
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


/*
******************************************************************************
Populate a RaytracePrimitives struct with values from a BSSN struct
 (plus derivatives on the A_ij field)
******************************************************************************
*/

void BSSN::setRaytracePrimitives(RayTrace<real_t, idx_t> *rt)
{
  setRaytraceCornerPrimitives(rt);
  rt->interpolatePrimitives();
}

void BSSN::setRaytraceCornerPrimitives(RayTrace<real_t, idx_t> *rt)
{
  BSSNData b_paq = {0};

  struct RaytracePrimitives<real_t> corner_rp[2][2][2];

  idx_t x_idx = rt->getRayIDX(1, dx, NX);
  idx_t y_idx = rt->getRayIDX(2, dx, NY);
  idx_t z_idx = rt->getRayIDX(3, dx, NZ);

  set_paq_values(x_idx + 0, y_idx + 0, z_idx + 0, &b_paq);
  corner_rp[0][0][0] = getRaytraceData(&b_paq);
  set_paq_values(x_idx + 0, y_idx + 0, z_idx + 1, &b_paq);
  corner_rp[0][0][1] = getRaytraceData(&b_paq);
  set_paq_values(x_idx + 0, y_idx + 1, z_idx + 0, &b_paq);
  corner_rp[0][1][0] = getRaytraceData(&b_paq);
  set_paq_values(x_idx + 0, y_idx + 1, z_idx + 1, &b_paq);
  corner_rp[0][1][1] = getRaytraceData(&b_paq);
  set_paq_values(x_idx + 1, y_idx + 0, z_idx + 0, &b_paq);
  corner_rp[1][0][0] = getRaytraceData(&b_paq);
  set_paq_values(x_idx + 1, y_idx + 0, z_idx + 1, &b_paq);
  corner_rp[1][0][1] = getRaytraceData(&b_paq);
  set_paq_values(x_idx + 1, y_idx + 1, z_idx + 0, &b_paq);
  corner_rp[1][1][0] = getRaytraceData(&b_paq);
  set_paq_values(x_idx + 1, y_idx + 1, z_idx + 1, &b_paq);
  corner_rp[1][1][1] = getRaytraceData(&b_paq);

  rt->copyInCornerPrimitives(corner_rp);
}

RaytracePrimitives<real_t> BSSN::getRaytraceData(BSSNData *paq)
{
  RaytracePrimitives<real_t> rp = {0};

  // normalization factor
  real_t P = exp(4.0*paq->phi);

  // metric
  rp.g[0] = P*paq->gamma11; rp.g[1] = P*paq->gamma12; rp.g[2] = P*paq->gamma13;
  rp.g[3] = P*paq->gamma22; rp.g[4] = P*paq->gamma23; rp.g[5] = P*paq->gamma33;
  // inverse metric
  rp.gi[0] = paq->gammai11/P; rp.gi[1] = paq->gammai12/P; rp.gi[2] = paq->gammai13/P;
  rp.gi[3] = paq->gammai22/P; rp.gi[4] = paq->gammai23/P; rp.gi[5] = paq->gammai33/P;
  // derivatives of metric
  rp.dg[0][0] = BSSN_RP_DG(1,1,1); rp.dg[0][1] = BSSN_RP_DG(1,2,1); rp.dg[0][2] = BSSN_RP_DG(1,3,1);
  rp.dg[0][3] = BSSN_RP_DG(2,2,1); rp.dg[0][4] = BSSN_RP_DG(2,3,1); rp.dg[0][5] = BSSN_RP_DG(3,3,1);
  rp.dg[1][0] = BSSN_RP_DG(1,1,2); rp.dg[1][1] = BSSN_RP_DG(1,2,2); rp.dg[1][2] = BSSN_RP_DG(1,3,2);
  rp.dg[1][3] = BSSN_RP_DG(2,2,2); rp.dg[1][4] = BSSN_RP_DG(2,3,2); rp.dg[1][5] = BSSN_RP_DG(3,3,2);
  rp.dg[2][0] = BSSN_RP_DG(1,1,3); rp.dg[2][1] = BSSN_RP_DG(1,2,3); rp.dg[2][2] = BSSN_RP_DG(1,3,3);
  rp.dg[2][3] = BSSN_RP_DG(2,2,3); rp.dg[2][4] = BSSN_RP_DG(2,3,3); rp.dg[2][5] = BSSN_RP_DG(3,3,3);
  // second derivatives of metric
  rp.ddg[0][0] = BSSN_RP_DDG(1,1,1,1); rp.ddg[0][1] = BSSN_RP_DDG(1,2,1,1); rp.ddg[0][2] = BSSN_RP_DDG(1,3,1,1);
  rp.ddg[0][3] = BSSN_RP_DDG(2,2,1,1); rp.ddg[0][4] = BSSN_RP_DDG(2,3,1,1); rp.ddg[0][5] = BSSN_RP_DDG(3,3,1,1);
  rp.ddg[1][0] = BSSN_RP_DDG(1,1,1,2); rp.ddg[1][1] = BSSN_RP_DDG(1,2,1,2); rp.ddg[1][2] = BSSN_RP_DDG(1,3,1,2);
  rp.ddg[1][3] = BSSN_RP_DDG(2,2,1,2); rp.ddg[1][4] = BSSN_RP_DDG(2,3,1,2); rp.ddg[1][5] = BSSN_RP_DDG(3,3,1,2);
  rp.ddg[2][0] = BSSN_RP_DDG(1,1,1,3); rp.ddg[2][1] = BSSN_RP_DDG(1,2,1,3); rp.ddg[2][2] = BSSN_RP_DDG(1,3,1,3);
  rp.ddg[2][3] = BSSN_RP_DDG(2,2,1,3); rp.ddg[2][4] = BSSN_RP_DDG(2,3,1,3); rp.ddg[2][5] = BSSN_RP_DDG(3,3,1,3);
  rp.ddg[3][0] = BSSN_RP_DDG(1,1,2,2); rp.ddg[3][1] = BSSN_RP_DDG(1,2,2,2); rp.ddg[3][2] = BSSN_RP_DDG(1,3,2,2);
  rp.ddg[3][3] = BSSN_RP_DDG(2,2,2,2); rp.ddg[3][4] = BSSN_RP_DDG(2,3,2,2); rp.ddg[3][5] = BSSN_RP_DDG(3,3,2,2);
  rp.ddg[4][0] = BSSN_RP_DDG(1,1,2,3); rp.ddg[4][1] = BSSN_RP_DDG(1,2,2,3); rp.ddg[4][2] = BSSN_RP_DDG(1,3,2,3);
  rp.ddg[4][3] = BSSN_RP_DDG(2,2,2,3); rp.ddg[4][4] = BSSN_RP_DDG(2,3,2,3); rp.ddg[4][5] = BSSN_RP_DDG(3,3,2,3);
  rp.ddg[5][0] = BSSN_RP_DDG(1,1,3,3); rp.ddg[5][1] = BSSN_RP_DDG(1,2,3,3); rp.ddg[5][2] = BSSN_RP_DDG(1,3,3,3);
  rp.ddg[5][3] = BSSN_RP_DDG(2,2,3,3); rp.ddg[5][4] = BSSN_RP_DDG(2,3,3,3); rp.ddg[5][5] = BSSN_RP_DDG(3,3,3,3);

  // extrinsic curvature:
  rp.K[0] = BSSN_RP_K(1,1); rp.K[1] = BSSN_RP_K(1,2); rp.K[2] = BSSN_RP_K(1,3);
  rp.K[3] = BSSN_RP_K(2,2); rp.K[4] = BSSN_RP_K(2,3); rp.K[5] = BSSN_RP_K(3,3);
  // derivatives of extrinsic curvature
  rp.dK[0][0] = BSSN_RP_DK(1,1,1); rp.dK[0][1] = BSSN_RP_DK(1,2,1); rp.dK[0][2] = BSSN_RP_DK(1,3,1);
  rp.dK[0][3] = BSSN_RP_DK(2,2,1); rp.dK[0][4] = BSSN_RP_DK(2,3,1); rp.dK[0][5] = BSSN_RP_DK(3,3,1);
  rp.dK[1][0] = BSSN_RP_DK(1,1,2); rp.dK[1][1] = BSSN_RP_DK(1,2,2); rp.dK[1][2] = BSSN_RP_DK(1,3,2);
  rp.dK[1][3] = BSSN_RP_DK(2,2,2); rp.dK[1][4] = BSSN_RP_DK(2,3,2); rp.dK[1][5] = BSSN_RP_DK(3,3,2);
  rp.dK[2][0] = BSSN_RP_DK(1,1,3); rp.dK[2][1] = BSSN_RP_DK(1,2,3); rp.dK[2][2] = BSSN_RP_DK(1,3,3);
  rp.dK[2][3] = BSSN_RP_DK(2,2,3); rp.dK[2][4] = BSSN_RP_DK(2,3,3); rp.dK[2][5] = BSSN_RP_DK(3,3,3);

  // 3-Christoffel symbols
  // raised first index
  rp.G[0][0] = BSSN_RP_GAMMA(1,1,1); rp.G[0][1] = BSSN_RP_GAMMA(1,2,1); rp.G[0][2] = BSSN_RP_GAMMA(1,3,1);
  rp.G[0][3] = BSSN_RP_GAMMA(2,2,1); rp.G[0][4] = BSSN_RP_GAMMA(2,3,1); rp.G[0][5] = BSSN_RP_GAMMA(3,3,1);
  rp.G[1][0] = BSSN_RP_GAMMA(1,1,2); rp.G[1][1] = BSSN_RP_GAMMA(1,2,2); rp.G[1][2] = BSSN_RP_GAMMA(1,3,2);
  rp.G[1][3] = BSSN_RP_GAMMA(2,2,2); rp.G[1][4] = BSSN_RP_GAMMA(2,3,2); rp.G[1][5] = BSSN_RP_GAMMA(3,3,2);
  rp.G[2][0] = BSSN_RP_GAMMA(1,1,3); rp.G[2][1] = BSSN_RP_GAMMA(1,2,3); rp.G[2][2] = BSSN_RP_GAMMA(1,3,3);
  rp.G[2][3] = BSSN_RP_GAMMA(2,2,3); rp.G[2][4] = BSSN_RP_GAMMA(2,3,3); rp.G[2][5] = BSSN_RP_GAMMA(3,3,3);
  // lower first index
  rp.GL[0][0] = BSSN_RP_GAMMAL(1,1,1); rp.GL[0][1] = BSSN_RP_GAMMAL(1,2,1); rp.GL[0][2] = BSSN_RP_GAMMAL(1,3,1);
  rp.GL[0][3] = BSSN_RP_GAMMAL(2,2,1); rp.GL[0][4] = BSSN_RP_GAMMAL(2,3,1); rp.GL[0][5] = BSSN_RP_GAMMAL(3,3,1);
  rp.GL[1][0] = BSSN_RP_GAMMAL(1,1,2); rp.GL[1][1] = BSSN_RP_GAMMAL(1,2,2); rp.GL[1][2] = BSSN_RP_GAMMAL(1,3,2);
  rp.GL[1][3] = BSSN_RP_GAMMAL(2,2,2); rp.GL[1][4] = BSSN_RP_GAMMAL(2,3,2); rp.GL[1][5] = BSSN_RP_GAMMAL(3,3,2);
  rp.GL[2][0] = BSSN_RP_GAMMAL(1,1,3); rp.GL[2][1] = BSSN_RP_GAMMAL(1,2,3); rp.GL[2][2] = BSSN_RP_GAMMAL(1,3,3);
  rp.GL[2][3] = BSSN_RP_GAMMAL(2,2,3); rp.GL[2][4] = BSSN_RP_GAMMAL(2,3,3); rp.GL[2][5] = BSSN_RP_GAMMAL(3,3,3);

  // 3-Ricci tensor
  rp.Ricci[0] = paq->ricci11; rp.Ricci[1] = paq->ricci12;
  rp.Ricci[2] = paq->ricci13; rp.Ricci[3] = paq->ricci22;
  rp.Ricci[4] = paq->ricci23; rp.Ricci[5] = paq->ricci33;

  rp.rho = paq->r;
  rp.trK = paq->K;

  return rp;
}

}