#include "bssn.h"

#include <chrono>
#include <thread>
#include <iomanip>

namespace cosmo
{

BSSN::BSSN()
{
  // BSSN fields
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_ALLOC)
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_ADDMAP)

  // any additional arrays for calcuated quantities
  BSSN_APPLY_TO_GEN1_EXTRAS(GEN1_ARRAY_ALLOC)
  BSSN_APPLY_TO_GEN1_EXTRAS(GEN1_ARRAY_ADDMAP)

  // FRW reference integrator
  frw = new FRW<real_t> (0.0, 0.0);
}

BSSN::~BSSN()
{
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_DELETE)
  BSSN_APPLY_TO_GEN1_EXTRAS(GEN1_ARRAY_DELETE)
}

/*
******************************************************************************
Functionality for special "wedge" RK4 integrator using 1 register

Integration procedure outline:
- For i = 0, i < NX + 2*Length(afterimage array), i++ :
   - Data starts in p register, _a references _p
   - Fill "rightmost" value in _K1 register (area loop)
     - point _a to _K1
   - Fill "rightmost" value in _K2 register
     - point _a to _K2
   - Fill "rightmost" value in _K3 register
     - point _a to _K3
   - Fill in "rightmost" point in tail
     - if i < WEDGE_SLICE_1_LEN + WEDGE_TAIL_LEN
         "burn-in" phase, do nothing
       Tail and wedge should be fully populated after this
     - if i == WEDGE_SLICE_1_LEN + WEDGE_TAIL_LEN
         populate afterimage
     - if i > WEDGE_SLICE_1_LEN + WEDGE_TAIL_LEN
         "leftmost" point in tail goes in _p register

******************************************************************************
*/

#define DISPP(v) (v==0.0?0.00:(-1.0*log10(fabs(v))))
#define BETWEEN(l, n1, n2) ( l>n1 && l<=n2 )
void BSSN::DrawWedgeSlice(idx_t i_p, idx_t i_K1, idx_t i_K2,
  idx_t i_K3, idx_t i_tail)
{
  return;
  int l;
  l = system("clear");

  std::cout.precision(3);
  std::cout << std::fixed << "\n\n";
  for(l=0; l<NX; l++) std::cout << DISPP(A22_p(l,0,0)) << " ";
  std::cout << "\n";
  for(l=0; l<NX; l++) {
    if(BETWEEN(l, i_K1-WEDGE_SLICE_1_LEN, i_K1)) { std::cout << DISPP(A22_K1(l,0,0)) << " "; }
    else if(BETWEEN(l+NX, i_K1-WEDGE_SLICE_1_LEN, i_K1)) { std::cout << DISPP(A22_K1(l+NX,0,0)) << " "; }
    else if(BETWEEN(l-NX, i_K1-WEDGE_SLICE_1_LEN, i_K1)) { std::cout << DISPP(A22_K1(l-NX,0,0)) << " "; }
    else { std::cout << "      "; }
  }
  std::cout << "\n";
  for(l=0; l<NX; l++) {
    if(BETWEEN(l, i_K2-WEDGE_SLICE_2_LEN, i_K2)) { std::cout << DISPP(A22_K2(l,0,0)) << " "; }
    else if(BETWEEN(l+NX, i_K2-WEDGE_SLICE_2_LEN, i_K2)) { std::cout << DISPP(A22_K2(l+NX,0,0)) << " "; }
    else if(BETWEEN(l-NX, i_K2-WEDGE_SLICE_2_LEN, i_K2)) { std::cout << DISPP(A22_K2(l-NX,0,0)) << " "; }
    else { std::cout << "      "; }
  }
  std::cout << "\n";
  for(l=0; l<NX; l++) {
    if(BETWEEN(l, i_K3-WEDGE_SLICE_3_LEN, i_K3)) { std::cout << DISPP(A22_K3(l,0,0)) << " "; }
    else if(BETWEEN(l+NX, i_K3-WEDGE_SLICE_3_LEN, i_K3)) { std::cout << DISPP(A22_K3(l+NX,0,0)) << " "; }
    else if(BETWEEN(l-NX, i_K3-WEDGE_SLICE_3_LEN, i_K3)) { std::cout << DISPP(A22_K3(l-NX,0,0)) << " "; }
    else { std::cout << "      "; }
  }
  std::cout << "\n";
  for(l=0; l<NX; l++) {
    if(BETWEEN(l, i_tail-WEDGE_TAIL_LEN, i_tail)) { std::cout << DISPP(A22_t(l,0,0)) << " "; }
    else if(BETWEEN(l+NX, i_tail-WEDGE_TAIL_LEN, i_tail)) { std::cout << DISPP(A22_t(l+NX,0,0)) << " "; }
    else if(BETWEEN(l-NX, i_tail-WEDGE_TAIL_LEN, i_tail)) { std::cout << DISPP(A22_t(l-NX,0,0)) << " "; }
    else { std::cout << "      "; }
  }
  std::cout << "\n";
  for(l=0; l<NX; l++) if(l>=WEDGE_AFTER_START_IDX && l<=WEDGE_AFTER_END_IDX) { std::cout << DISPP(A22_r(l,0,0)) << " "; } else { std::cout << "      "; }
  std::cout << "\n";
  
  std::this_thread::sleep_for(std::chrono::milliseconds(200));

  return;
}

void BSSN::WedgeStep()
{
  idx_t i, j, k, l;
  for(i=0; i<=NX + WEDGE_AFTER_END_IDX - WEDGE_SLICE_LEN_DIFF/2; ++i)
  {
    // "i" is location of calculation for wedge "peak"
    idx_t i_p = i;
    idx_t i_K1 = i;
    idx_t i_K2 = i - WEDGE_SLICE_LEN_DIFF/2;
    idx_t i_K3 = i - 2*WEDGE_SLICE_LEN_DIFF/2;
    idx_t i_tail = i - 3*WEDGE_SLICE_LEN_DIFF/2;

    frw->P0_step();
    WedgeK1Calc(i_K1); // uses _p
    frw->P1_step(dt);
    WedgeK2Calc(i_K2); // uses "_p + K1"
    frw->P2_step(dt);
    WedgeK3Calc(i_K3); // uses "_p + K2"
    frw->P3_step(dt);
    WedgeTailCalc(i_tail); // uses "_p + K3"

// calculate phi and rho on a constant-eta slice
real_t eta_val = 1.0;
PARALLEL_AREA_LOOP(j,k)
{
  real_t eta_ref = 2.0;
  // eta increases monatonically
  // check to see when simulation has crossed a reference eta
  real_t eta2 = eta_t(i_tail,j,k); // next step is in _t array
  real_t eta1 = eta_p(i_tail,j,k); // previous step is in _p array

  if( eta_ref >= eta1 && eta_ref <= eta2 )
  {
    // linear interpolation to determine "time" of K_ref
    real_t ref_time = (eta_ref - eta1) / (eta2 - eta1);
    
    // get phi value at this "time"
    real_t phi2 = DIFFphi_t(i_tail,j,k);
    real_t phi1 = DIFFphi_p(i_tail,j,k);
    real_t phi_at_ref = (phi2 - phi1)*ref_time + phi1;
    // store phi value on the slice
    eta_phi_slice_a(i_tail,j,k) = phi_at_ref;

    // get rho value at this "time"
    real_t rho2 = DIFFdustrho_t(i_tail,j,k);
    real_t rho1 = DIFFdustrho_p(i_tail,j,k);
    real_t rho_at_ref = (rho2 - rho1)*ref_time + rho1;
    // store rho value on the slice
    eta_rho_slice_a(i_tail,j,k) = rho_at_ref;
  }
}


    // Store snapshot to "afterimage"
    if(i_tail == WEDGE_AFTER_END_IDX)
    {
      // populate afterimage
      for(l = WEDGE_AFTER_START_IDX; l <= WEDGE_AFTER_END_IDX; ++l)
      {
        PARALLEL_AREA_LOOP(j,k)
        {
          BSSN_STORE_AFTER();
        }
      }
    }

    // store points from tail back in 
    if( i_tail >= WEDGE_AFTER_END_IDX + WEDGE_TAIL_LEN )
    {
      // "leftmost" point in tail goes in _p register
      PARALLEL_AREA_LOOP(j,k)
      {
        BSSN_STORE_TAILEND();
      }
    }
  }

  // restore afterimage and tail
  for(i = WEDGE_AFTER_START_IDX; i <= WEDGE_AFTER_END_IDX; ++i)
  {
    PARALLEL_AREA_LOOP(j,k)
    {
      BSSN_RE_STORE_AFTER();
    }
  }
  for(i = NX + WEDGE_AFTER_START_IDX - WEDGE_TAIL_LEN; i < NX + WEDGE_AFTER_START_IDX; ++i)
  {
    PARALLEL_AREA_LOOP(j,k)
    {
      BSSN_RE_STORE_TAIL();
    }
  }
  // advance FRW integrator
  frw->RK_total_step(dt);
}

// y_{K1} = y_n + h/2*f[y_n]
void BSSN::WedgeK1Calc(idx_t i_K1)
{
  idx_t j, k;
  BSSN_ASSIGN_TO_A_REGISTER(_p);
  PARALLEL_AREA_LOOP(j,k)
  {
    BSSNData paq = {0};
    set_paq_values(i_K1, j, k, &paq);
    BSSN_COMPUTE_WEDGE_STEP_K1();
  }
}

// y_{K2} = y_n + h/2*f[y_{K1}]
void BSSN::WedgeK2Calc(idx_t i_K2)
{
  idx_t j, k;
  BSSN_ASSIGN_TO_A_REGISTER(_K1);
  PARALLEL_AREA_LOOP(j,k)
  {
    BSSNData paq = {0};
    set_paq_values(i_K2, j, k, &paq);
    BSSN_COMPUTE_WEDGE_STEP_K2();
  }
}

// y_{K3} = y_n +   h*f[y_{K2}]
void BSSN::WedgeK3Calc(idx_t i_K3)
{
  idx_t j, k;
  BSSN_ASSIGN_TO_A_REGISTER(_K2);
  PARALLEL_AREA_LOOP(j,k)
  {
    BSSNData paq = {0};
    set_paq_values(i_K3, j, k, &paq);
    BSSN_COMPUTE_WEDGE_STEP_K3();
  }
}

// y_{n+1} = ( -y_n + y_{K1} + 2 y_{K2} + y_{K3} ) / 3 + h/6*f[y_{K3}]
void BSSN::WedgeTailCalc(idx_t i_tail)
{
  idx_t j, k;
  BSSN_ASSIGN_TO_A_REGISTER(_K3);
  PARALLEL_AREA_LOOP(j,k)
  {
    BSSNData paq = {0};
    set_paq_values(i_tail, j, k, &paq);
    BSSN_COMPUTE_WEDGE_STEP_TAIL();
  }
  BSSN_ASSIGN_TO_A_REGISTER(_p);
}


/*
******************************************************************************
Functionality for setting/calculating data values
******************************************************************************
*/

void BSSN::set_gammai_values(idx_t i, idx_t j, idx_t k, BSSNData *paq)
{
  // Compute the inverse metric algebraically at each point
  // assumes det(gamma) = 1
  idx_t idx = NP_INDEX(i,j,k);
  paq->gammai11 = 1.0 + paq->DIFFgamma22 + paq->DIFFgamma33 - pw2(paq->DIFFgamma23) + paq->DIFFgamma22*paq->DIFFgamma33;
  paq->gammai22 = 1.0 + paq->DIFFgamma11 + paq->DIFFgamma33 - pw2(paq->DIFFgamma13) + paq->DIFFgamma11*paq->DIFFgamma33;
  paq->gammai33 = 1.0 + paq->DIFFgamma11 + paq->DIFFgamma22 - pw2(paq->DIFFgamma12) + paq->DIFFgamma11*paq->DIFFgamma22;
  paq->gammai12 = paq->DIFFgamma13*paq->DIFFgamma23 - paq->DIFFgamma12*(1.0 + paq->DIFFgamma33);
  paq->gammai13 = paq->DIFFgamma12*paq->DIFFgamma23 - paq->DIFFgamma13*(1.0 + paq->DIFFgamma22);
  paq->gammai23 = paq->DIFFgamma12*paq->DIFFgamma13 - paq->DIFFgamma23*(1.0 + paq->DIFFgamma11);
}

void BSSN::set_DIFFgamma_Aij_norm()
{
  idx_t i, j, k;

  /* This potentially breaks conservation of trace:
   * need to come up with something else to preserve both
   * trace + determinant constraints.
   */
  PARALLEL_LOOP3(i, j, k)
  {
    set_DIFFgamma_Aij_norm_pt(i,j,k);
  }
}

void BSSN::set_DIFFgamma_Aij_norm_pt(idx_t i, idx_t j, idx_t k)
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

void BSSN::set_paq_values(idx_t i, idx_t j, idx_t k, BSSNData *paq)
{
  paq->i = i;
  paq->j = j;
  paq->k = k;

  // this might already be set
  if(paq->idx == 0)
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
  paq->r        =   paq->DIFFdustrho + paq->rho_FRW;
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
  paq->H = hamiltonianConstraintCalc(paq);
}

void BSSN::init()
{
  idx_t idx;

  idx_t i, j, k;
  PARALLEL_LOOP3(i, j, k)
  {
    idx = NP_INDEX(i,j,k);
    BSSN_ZERO_ARRAYS(_p, idx)
    BSSN_ZERO_GEN1_EXTRAS()
  }
}

/* set current local field values */
void BSSN::set_local_vals(BSSNData *paq)
{
  // Pull out values of quantities at a single point
  BSSN_APPLY_TO_FIELDS(SET_LOCAL_VALUES);
  BSSN_APPLY_TO_GEN1_EXTRAS(SET_LOCAL_VALUES);
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
  AijAij_a(paq->i, paq->j, paq->k)
    = paq->Acont11*paq->A11 + paq->Acont22*paq->A22 + paq->Acont33*paq->A33
      + 2.0*(paq->Acont12*paq->A12 + paq->Acont13*paq->A13 + paq->Acont23*paq->A23);
  paq->AijAij = AijAij_a(paq->i, paq->j, paq->k);
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
  paq->d1phi = derivative(paq->i, paq->j, paq->k, 1, &DIFFphi_a);
  paq->d2phi = derivative(paq->i, paq->j, paq->k, 2, &DIFFphi_a);
  paq->d3phi = derivative(paq->i, paq->j, paq->k, 3, &DIFFphi_a);

  // normal derivatives of alpha
  paq->d1a = derivative(paq->i, paq->j, paq->k, 1, &DIFFalpha_a);
  paq->d2a = derivative(paq->i, paq->j, paq->k, 2, &DIFFalpha_a);
  paq->d3a = derivative(paq->i, paq->j, paq->k, 3, &DIFFalpha_a);
}

void BSSN::calculate_dK(BSSNData *paq)
{
  // normal derivatives of phi
  paq->d1K = derivative(paq->i, paq->j, paq->k, 1, &DIFFK_a);
  paq->d2K = derivative(paq->i, paq->j, paq->k, 2, &DIFFK_a);
  paq->d3K = derivative(paq->i, paq->j, paq->k, 3, &DIFFK_a);
}

#if USE_Z4c_DAMPING
void BSSN::calculate_dtheta(BSSNData *paq)
{
  // normal derivatives of phi
  paq->d1theta = derivative(paq->i, paq->j, paq->k, 1, &theta_a);
  paq->d2theta = derivative(paq->i, paq->j, paq->k, 2, &theta_a);
  paq->d3theta = derivative(paq->i, paq->j, paq->k, 3, &theta_a);
}
#endif

#if USE_BSSN_SHIFT
void BSSN::calculate_dbeta(BSSNData *paq)
{
  paq->d1beta1 = derivative(paq->i, paq->j, paq->k, 1, &beta1_a);
  paq->d1beta2 = derivative(paq->i, paq->j, paq->k, 1, &beta2_a);
  paq->d1beta3 = derivative(paq->i, paq->j, paq->k, 1, &beta3_a);
  paq->d2beta1 = derivative(paq->i, paq->j, paq->k, 2, &beta1_a);
  paq->d2beta2 = derivative(paq->i, paq->j, paq->k, 2, &beta2_a);
  paq->d2beta3 = derivative(paq->i, paq->j, paq->k, 2, &beta3_a);
  paq->d3beta1 = derivative(paq->i, paq->j, paq->k, 3, &beta1_a);
  paq->d3beta2 = derivative(paq->i, paq->j, paq->k, 3, &beta2_a);
  paq->d3beta3 = derivative(paq->i, paq->j, paq->k, 3, &beta3_a);
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
  paq->D1D1phi = double_derivative(i, j, k, 1, 1, &DIFFphi_a) - (paq->G111*paq->d1phi + paq->G211*paq->d2phi + paq->G311*paq->d3phi);
  paq->D2D2phi = double_derivative(i, j, k, 2, 2, &DIFFphi_a) - (paq->G122*paq->d1phi + paq->G222*paq->d2phi + paq->G322*paq->d3phi);
  paq->D3D3phi = double_derivative(i, j, k, 3, 3, &DIFFphi_a) - (paq->G133*paq->d1phi + paq->G233*paq->d2phi + paq->G333*paq->d3phi);

  paq->D1D2phi = double_derivative(i, j, k, 1, 2, &DIFFphi_a) - (paq->G112*paq->d1phi + paq->G212*paq->d2phi + paq->G312*paq->d3phi);
  paq->D1D3phi = double_derivative(i, j, k, 1, 3, &DIFFphi_a) - (paq->G113*paq->d1phi + paq->G213*paq->d2phi + paq->G313*paq->d3phi);
  paq->D2D3phi = double_derivative(i, j, k, 2, 3, &DIFFphi_a) - (paq->G123*paq->d1phi + paq->G223*paq->d2phi + paq->G323*paq->d3phi);  
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
  BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_RICCITF_UNITARY)

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
  ricci_a(paq->i, paq->j, paq->k) = paq->ricci;

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

real_t BSSN::ev_DIFFgamma11(BSSNData *paq) { return BSSN_DT_DIFFGAMMAIJ(1, 1) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->DIFFgamma11 - KO_dissipation_Q(paq->i, paq->j, paq->k, &DIFFgamma11_a); }
real_t BSSN::ev_DIFFgamma12(BSSNData *paq) { return BSSN_DT_DIFFGAMMAIJ(1, 2) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->DIFFgamma12 - KO_dissipation_Q(paq->i, paq->j, paq->k, &DIFFgamma12_a); }
real_t BSSN::ev_DIFFgamma13(BSSNData *paq) { return BSSN_DT_DIFFGAMMAIJ(1, 3) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->DIFFgamma13 - KO_dissipation_Q(paq->i, paq->j, paq->k, &DIFFgamma13_a); }
real_t BSSN::ev_DIFFgamma22(BSSNData *paq) { return BSSN_DT_DIFFGAMMAIJ(2, 2) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->DIFFgamma22 - KO_dissipation_Q(paq->i, paq->j, paq->k, &DIFFgamma22_a); }
real_t BSSN::ev_DIFFgamma23(BSSNData *paq) { return BSSN_DT_DIFFGAMMAIJ(2, 3) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->DIFFgamma23 - KO_dissipation_Q(paq->i, paq->j, paq->k, &DIFFgamma23_a); }
real_t BSSN::ev_DIFFgamma33(BSSNData *paq) { return BSSN_DT_DIFFGAMMAIJ(3, 3) + 0.5*BS_H_DAMPING_AMPLITUDE*dt*paq->H*paq->DIFFgamma33 - KO_dissipation_Q(paq->i, paq->j, paq->k, &DIFFgamma33_a); }

real_t BSSN::ev_A11(BSSNData *paq) { return BSSN_DT_AIJ(1, 1) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A11*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, &A11_a); }
real_t BSSN::ev_A12(BSSNData *paq) { return BSSN_DT_AIJ(1, 2) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A12*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, &A12_a); }
real_t BSSN::ev_A13(BSSNData *paq) { return BSSN_DT_AIJ(1, 3) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A13*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, &A13_a); }
real_t BSSN::ev_A22(BSSNData *paq) { return BSSN_DT_AIJ(2, 2) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A22*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, &A22_a); }
real_t BSSN::ev_A23(BSSNData *paq) { return BSSN_DT_AIJ(2, 3) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A23*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, &A23_a); }
real_t BSSN::ev_A33(BSSNData *paq) { return BSSN_DT_AIJ(3, 3) - 1.0*BS_H_DAMPING_AMPLITUDE*dt*paq->A33*paq->H - KO_dissipation_Q(paq->i, paq->j, paq->k, &A33_a); }

real_t BSSN::ev_Gamma1(BSSNData *paq) { return BSSN_DT_GAMMAI(1) - KO_dissipation_Q(paq->i, paq->j, paq->k, &Gamma1_a); }
real_t BSSN::ev_Gamma2(BSSNData *paq) { return BSSN_DT_GAMMAI(2) - KO_dissipation_Q(paq->i, paq->j, paq->k, &Gamma2_a); }
real_t BSSN::ev_Gamma3(BSSNData *paq) { return BSSN_DT_GAMMAI(3) - KO_dissipation_Q(paq->i, paq->j, paq->k, &Gamma3_a); }

real_t BSSN::ev_DIFFK(BSSNData *paq)
{
  return (
    - paq->DDaTR
    + paq->alpha*(
        paq->AijAij
        + 1.0/3.0*(paq->DIFFK + 2.0*paq->theta)*(paq->DIFFK + 2.0*paq->theta + 2.0*paq->K_FRW)
    )
    + 4.0*PI*paq->alpha*(paq->DIFFdustrho)
    - paq->DIFFalpha*(
        1.0/3.0*pw2(paq->K_FRW)
        + 4.0*PI*(paq->rho_FRW + paq->S_FRW)
      )
    + paq->beta1*paq->d1K + paq->beta2*paq->d2K + paq->beta3*paq->d3K
    - 1.0*JM_K_DAMPING_AMPLITUDE*paq->H*exp(-5.0*paq->phi)
    + Z4c_K1_DAMPING_AMPLITUDE*(1.0 - Z4c_K2_DAMPING_AMPLITUDE)*paq->theta
    - KO_dissipation_Q(paq->i, paq->j, paq->k, &DIFFK_a)
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
    - KO_dissipation_Q(paq->i, paq->j, paq->k, &DIFFphi_a)
  );
}

real_t BSSN::ev_DIFFalpha(BSSNData *paq)
{
  #if USE_HARMONIC_ALPHA
    return -1.0*pw2(paq->alpha)*paq->K - KO_dissipation_Q(paq->i, paq->j, paq->k, &DIFFalpha_a);
  #endif

  #if USE_CONFORMAL_SYNC_ALPHA
    return -1.0/3.0*paq->alpha*paq->K_FRW;
  #endif

  return 0.0;
}

real_t BSSN::ev_DIFFdustrho(BSSNData *paq)
{
  // only good for alpha = 1
  return paq->K_FRW*paq->DIFFdustrho + paq->DIFFK*paq->rho_FRW + paq->DIFFK*paq->DIFFdustrho;
}

real_t BSSN::ev_eta(BSSNData *paq)
{
  // Appx. newtonian gauge time elapsed
  // a(eta) = a_0 * eta^(2/3)
  // d\eta / dt = eta^(4/3) * exp( 2 phi )
  // eta_0 = 
  return pow(paq->eta, 4.0/3.0)*exp(2.0*(paq->phi));
}

#if USE_Z4c_DAMPING
real_t BSSN::ev_theta(BSSNData *paq)
{
  return (
    0.5*paq->alpha*(
      paq->ricci + 2.0/3.0*pw2(paq->K + 2.0*paq->theta) - paq->AijAij - 16.0*PI*( paq->DIFFdustrho + paq->rho_FRW )
    )
    - paq->alpha*Z4c_K1_DAMPING_AMPLITUDE*(2.0 + Z4c_K2_DAMPING_AMPLITUDE)*paq->theta
  ) - KO_dissipation_Q(paq->i, paq->j, paq->k, &theta_a);
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
    PARALLEL_LOOP3(i,j,k)
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
      ricci_a[idx] + 2.0/3.0*pw2(K_FRW + DIFFK_a[idx] + 2.0*theta) - AijAij_a[idx] - 16.0*PI*(DIFFdustrho_a[idx] + rho_FRW)
    );
  #else
    return -exp(5.0*DIFFphi_a[idx])/8.0*(
      ricci_a[idx] + 2.0/3.0*pw2(DIFFK_a[idx] + 2.0*theta) - AijAij_a[idx] - 16.0*PI*DIFFdustrho_a[idx]
    );
  #endif
}

real_t BSSN::hamiltonianConstraintCalc(BSSNData *paq)
{
  #if USE_Z4c_DAMPING
    real_t theta = paq->theta;
  #else
    real_t theta = 0.0;
  #endif

  #if USE_REFERENCE_FRW
    real_t K_FRW = frw->get_K();
    real_t phi_FRW = frw->get_phi();
    real_t rho_FRW = frw->get_rho();
    return -exp(5.0*(paq->DIFFphi + phi_FRW))/8.0*(
      paq->ricci + 2.0/3.0*pw2(K_FRW + paq->DIFFK + 2.0*theta) - paq->AijAij - 16.0*PI*(paq->DIFFdustrho + rho_FRW)
    );
  #else
    return -exp(5.0*paq->DIFFphi)/8.0*(
      paq->ricci + 2.0/3.0*pw2(paq->DIFFK + 2.0*theta) - paq->AijAij - 16.0*PI*paq->DIFFdustrho
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
    sqrt( pw2(ricci_a[idx]) + pw2(AijAij_a[idx]) + pw2(2.0/3.0*pw2(K_FRW + DIFFK_a[idx] + 2.0*theta)) + pw2(16.0*PI*(DIFFdustrho_a[idx] + rho_FRW))
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

}
