
#include "bssn.h"

namespace cosmo
{

BSSN::BSSN()
{
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_ALLOC)
  BSSN_APPLY_TO_SOURCES(GEN1_ARRAY_ALLOC)

  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_ADDMAP)
  BSSN_APPLY_TO_SOURCES(GEN1_ARRAY_ADDMAP)
}

BSSN::~BSSN()
{
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_DELETE)
  BSSN_APPLY_TO_SOURCES(GEN1_ARRAY_DELETE)
}


void BSSN::set_paq_values(idx_t i, idx_t j, idx_t k, BSSNData *paq)
{
  paq->i = i;
  paq->j = j;
  paq->k = k;
  paq->idx = INDEX(i,j,k);

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

  // source values
  set_source_vals(paq);  
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
  LOOP3(i, j, k)
  {
    idx_t idx = INDEX(i, j, k);
    BSSN_ZERO_ARRAY(_f, idx);
  }
}


// First RK step: calculate and add k_1 coeff to _f array
void BSSN::K1Calc(BSSNData *paq)
{
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
  LOOP3(i, j, k)
  {
    idx_t idx = INDEX(i,j,k);

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

  LOOP3(i, j, k)
  {
    idx = INDEX(i,j,k);

    gamma11_p[idx]  = gamma11_f[idx]  = 1.0;
    gamma12_p[idx]  = gamma12_f[idx]  = 0.0;
    gamma13_p[idx]  = gamma13_f[idx]  = 0.0;
    gamma22_p[idx]  = gamma22_f[idx]  = 1.0;
    gamma23_p[idx]  = gamma23_f[idx]  = 0.0;
    gamma33_p[idx]  = gamma33_f[idx]  = 1.0;
    
    phi_p[idx]      = phi_f[idx]      = 0.0;
    
    A11_p[idx]      = A11_f[idx]      = 0.0;
    A12_p[idx]      = A12_f[idx]      = 0.0;
    A13_p[idx]      = A13_f[idx]      = 0.0;
    A22_p[idx]      = A22_f[idx]      = 0.0;
    A23_p[idx]      = A23_f[idx]      = 0.0;
    A33_p[idx]      = A33_f[idx]      = 0.0;

    K_p[idx]        = K_f[idx]        = -sqrt(24.0*PI*(0.01));

    Gamma1_p[idx]   = Gamma1_f[idx]   = 0.0;
    Gamma2_p[idx]   = Gamma2_f[idx]   = 0.0;
    Gamma3_p[idx]   = Gamma3_f[idx]   = 0.0;

    beta1_p[idx]    = beta1_f[idx]    = 0.0;
    beta2_p[idx]    = beta2_f[idx]    = 0.0;
    beta3_p[idx]    = beta3_f[idx]    = 0.0;

    alpha_p[idx]    = alpha_f[idx]    = 1.0;

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


real_t BSSN::der(real_t field_adj[3][3][3], int d)
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

real_t BSSN::dder(real_t field_adj[3][3][3], int d1, int d2)
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

void BSSN::calculate_Acont(BSSNData *paq)
{
  BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_ACONT)
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

void BSSN::calculate_christoffels(BSSNData *paq)
{
  // christoffel symbols: \Gamma^i_{jk} = Gijk
  BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_CHRISTOFFEL)
}

void BSSN::calculate_dalpha_dphi(BSSNData *paq)
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

void BSSN::calculate_dbeta(BSSNData *paq)
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
void BSSN::calculateRicciTF(BSSNData *paq)
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

void BSSN::calculateDDphi(BSSNData *paq)
{
  // double covariant derivatives, using normal metric
  paq->D1D1phi = dder(paq->phi_adj, 1, 1) - (paq->G111*paq->d1phi + paq->G211*paq->d2phi + paq->G311*paq->d3phi);
  paq->D1D2phi = dder(paq->phi_adj, 1, 2) - (paq->G112*paq->d1phi + paq->G212*paq->d2phi + paq->G312*paq->d3phi);
  paq->D1D3phi = dder(paq->phi_adj, 1, 3) - (paq->G113*paq->d1phi + paq->G213*paq->d2phi + paq->G313*paq->d3phi);
  paq->D2D2phi = dder(paq->phi_adj, 2, 2) - (paq->G122*paq->d1phi + paq->G222*paq->d2phi + paq->G322*paq->d3phi);
  paq->D2D3phi = dder(paq->phi_adj, 2, 3) - (paq->G123*paq->d1phi + paq->G223*paq->d2phi + paq->G323*paq->d3phi);
  paq->D3D3phi = dder(paq->phi_adj, 3, 3) - (paq->G133*paq->d1phi + paq->G233*paq->d2phi + paq->G333*paq->d3phi);
}

void BSSN::calculateDDalphaTF(BSSNData *paq)
{
  // double covariant derivatives - use non-unitary metric - extra pieces that depend on phi!
  // the gamma*ldlphi are needed for the BSSN_CALCULATE_DIDJALPHA macro
  real_t gamma1ldlphi = paq->gammai11*paq->d1phi + paq->gammai12*paq->d2phi + paq->gammai13*paq->d3phi;
  real_t gamma2ldlphi = paq->gammai21*paq->d1phi + paq->gammai22*paq->d2phi + paq->gammai23*paq->d3phi;
  real_t gamma3ldlphi = paq->gammai31*paq->d1phi + paq->gammai32*paq->d2phi + paq->gammai33*paq->d3phi;
  // Calculates full (not trace-free) piece:
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

real_t BSSN::ev_gamma11(BSSNData *paq) { return BSSN_DT_GAMMAIJ(1, 1); }
real_t BSSN::ev_gamma12(BSSNData *paq) { return BSSN_DT_GAMMAIJ(1, 2); }
real_t BSSN::ev_gamma13(BSSNData *paq) { return BSSN_DT_GAMMAIJ(1, 3); }
real_t BSSN::ev_gamma22(BSSNData *paq) { return BSSN_DT_GAMMAIJ(2, 2); }
real_t BSSN::ev_gamma23(BSSNData *paq) { return BSSN_DT_GAMMAIJ(2, 3); }
real_t BSSN::ev_gamma33(BSSNData *paq) { return BSSN_DT_GAMMAIJ(3, 3); }

real_t BSSN::ev_A11(BSSNData *paq) { return BSSN_DT_AIJ(1, 1); }
real_t BSSN::ev_A12(BSSNData *paq) { return BSSN_DT_AIJ(1, 2); }
real_t BSSN::ev_A13(BSSNData *paq) { return BSSN_DT_AIJ(1, 3); }
real_t BSSN::ev_A22(BSSNData *paq) { return BSSN_DT_AIJ(2, 2); }
real_t BSSN::ev_A23(BSSNData *paq) { return BSSN_DT_AIJ(2, 3); }
real_t BSSN::ev_A33(BSSNData *paq) { return BSSN_DT_AIJ(3, 3); }

real_t BSSN::ev_K(BSSNData *paq)
{
  return (
    - paq->DDaTR
    + paq->alpha*(
        paq->A11*paq->A11 + paq->A22*paq->A22 + paq->A33*paq->A33
        + 2.0*(paq->A12*paq->A12 + paq->A13*paq->A13 + paq->A23*paq->A23)
        + (1.0/3.0)*paq->K*paq->K
      )
    + paq->beta1*der(paq->K_adj, 1) + paq->beta2*der(paq->K_adj, 2) + paq->beta3*der(paq->K_adj, 3)
    + 4.0*PI*paq->alpha*(paq->rho + paq->S)
  );
}

real_t BSSN::ev_phi(BSSNData *paq)
{
  return (
    1.0/6.0*(der(paq->beta1_adj, 1) + der(paq->beta2_adj, 2) + der(paq->beta3_adj, 3) - paq->alpha*paq->K)
    + paq->beta1*paq->d1phi + paq->beta2*paq->d2phi + paq->beta3*paq->d3phi
  );
}

real_t BSSN::ev_Gamma1(BSSNData *paq) { return BSSN_DT_GAMMAI(1); }
real_t BSSN::ev_Gamma2(BSSNData *paq) { return BSSN_DT_GAMMAI(2); }
real_t BSSN::ev_Gamma3(BSSNData *paq) { return BSSN_DT_GAMMAI(3); }

// static gauge for now:
// real_t BSSN::ev_alpha(BSSNData *paq) { return 0; }
// real_t BSSN::ev_beta1(BSSNData *paq) { return 0; }
// real_t BSSN::ev_beta2(BSSNData *paq) { return 0; }
// real_t BSSN::ev_beta3(BSSNData *paq) { return 0; }


}
