
#include "bssn.h"

namespace cosmo
{

BSSN::BSSN()
{
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_ALLOC)
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_ADDMAP)
}

BSSN::~BSSN()
{
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_DELETE)
}

// Full RK step (More useful when not evolving the source simultaneously)
void BSSN::step()
{
  _timer["RK Init"].start();
  stepInit();
  _timer["RK Init"].stop();

  _timer["RK K1 Calc"].start();
  K1Calc();
  _timer["RK K1 Calc"].stop();

  _timer["RK K2 Calc"].start();
  K2Calc();
  _timer["RK K2 Calc"].stop();

  _timer["RK K3 Calc"].start();
  K3Calc();
  _timer["RK K3 Calc"].stop();

  _timer["RK K4 Calc"].start();
  K4Calc();
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
void BSSN::K1Calc()
{
  LOOP3(i, j, k)
  {
    K1CalcPt(i, j, k);
  }
  // swap _c <-> _a registers
  BSSN_SWAP_ARRAYS(_c, _a);
}
void BSSN::K1CalcPt(idx_t i, idx_t j, idx_t k)
{
  set_paq_values(i, j, k, &paq);

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
  LOOP3(i, j, k)
  {
    K2CalcPt(i, j, k);
  }
  // swap _c <-> _a
  BSSN_SWAP_ARRAYS(_c, _a);
}
void BSSN::K2CalcPt(idx_t i, idx_t j, idx_t k)
{
  set_paq_values(i, j, k, &paq);

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
  LOOP3(i, j, k)
  {
    K3CalcPt(i, j, k);
  }
  // swap _c <-> _a
  BSSN_SWAP_ARRAYS(_c, _a);
}
void BSSN::K3CalcPt(idx_t i, idx_t j, idx_t k)
{
  set_paq_values(i, j, k, &paq);

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
  LOOP3(i, j, k)
  {
    K4CalcPt(i, j, k);
  }
}
void BSSN::K4CalcPt(idx_t i, idx_t j, idx_t k)
{
  set_paq_values(i, j, k, &paq);

  // evolve fields and add to _f register:
  // arr_f = arr_p + dt/6*(k_1 + 2*k_2 + 2*k_3 + k_4)
  //       = (1.0/3.0)*(arr_f - arr_p) + (1.0/6.0)*evfn(arr_a)
  BSSN_FINAL_RK4_STEP();
}


// arr_f register now holds "final" calculation; move back to _p register:
// swap _f <-> _p
void BSSN::stepTerm()
{
  BSSN_SWAP_ARRAYS(_f, _p);
}


void BSSN::init()
{
  real_t eps = 0.01;

  LOOP3(i, j, k)
  {
    paq.idx = INDEX(i,j,k);

    gamma11_p[paq.idx]  = 1.0;
    gamma12_p[paq.idx]  = 0.0;
    gamma13_p[paq.idx]  = 0.0;
    gamma22_p[paq.idx]  = 1.0;
    gamma23_p[paq.idx]  = 0.0;
    gamma33_p[paq.idx]  = 1.0;
    
    phi_p[paq.idx]      = 0.0;
    
    A11_p[paq.idx]      = 0.0;
    A12_p[paq.idx]      = 0.0;
    A13_p[paq.idx]      = 0.0;
    A22_p[paq.idx]      = 0.0;
    A23_p[paq.idx]      = 0.0;
    A33_p[paq.idx]      = 0.0;

    K_p[paq.idx]        = -3.0;

    Gamma1_p[paq.idx]   = 0.0;
    Gamma2_p[paq.idx]   = 0.0;
    Gamma3_p[paq.idx]   = 0.0;

    beta1_p[paq.idx]    = 0.0;
    beta2_p[paq.idx]    = 0.0;
    beta3_p[paq.idx]    = 0.0;

    alpha_p[paq.idx]    = 1.0;
  }
}

}
