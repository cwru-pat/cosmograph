#ifndef BSSN_MACROS
#define BSSN_MACROS

/*
 * applying functions to lots of vars
 */

#if USE_Z4c_DAMPING
  #define Z4c_APPLY_TO_FIELDS(function)  \
    function(theta);
  #define Z4c_APPLY_TO_FIELDS_ARGS(function, ...) \
    function(theta, __VA_ARGS__);
#else
  #define Z4c_APPLY_TO_FIELDS(function)
  #define Z4c_APPLY_TO_FIELDS_ARGS(function, ...)
#endif

#if USE_BSSN_SHIFT
  #define BSSN_APPLY_TO_SHIFT(function)  \
    function(beta1); \
    function(beta2); \
    function(beta3);
  #define BSSN_APPLY_TO_SHIFT_ARGS(function, ...) \
    function(beta1, __VA_ARGS__); \
    function(beta2, __VA_ARGS__); \
    function(beta3, __VA_ARGS__);
#else
  #define BSSN_APPLY_TO_SHIFT(function)
  #define BSSN_APPLY_TO_SHIFT_ARGS(function, ...)
#endif

#define BSSN_APPLY_TO_FIELDS_ARGS(function, ...)   \
  function(DIFFgamma11, __VA_ARGS__);              \
  function(DIFFgamma12, __VA_ARGS__);              \
  function(DIFFgamma13, __VA_ARGS__);              \
  function(DIFFgamma22, __VA_ARGS__);              \
  function(DIFFgamma23, __VA_ARGS__);              \
  function(DIFFgamma33, __VA_ARGS__);              \
  function(DIFFphi, __VA_ARGS__);                  \
  function(A11, __VA_ARGS__);                      \
  function(A12, __VA_ARGS__);                      \
  function(A13, __VA_ARGS__);                      \
  function(A22, __VA_ARGS__);                      \
  function(A23, __VA_ARGS__);                      \
  function(A33, __VA_ARGS__);                      \
  function(DIFFK, __VA_ARGS__);                    \
  function(Gamma1, __VA_ARGS__);                   \
  function(Gamma2, __VA_ARGS__);                   \
  function(Gamma3, __VA_ARGS__);                   \
  function(DIFFalpha, __VA_ARGS__);                \
  Z4c_APPLY_TO_FIELDS_ARGS(function, __VA_ARGS__)  \
  BSSN_APPLY_TO_SHIFT_ARGS(function, __VA_ARGS__)

#define BSSN_APPLY_TO_FIELDS(function) \
  function(DIFFgamma11);               \
  function(DIFFgamma12);               \
  function(DIFFgamma13);               \
  function(DIFFgamma22);               \
  function(DIFFgamma23);               \
  function(DIFFgamma33);               \
  function(DIFFphi);                   \
  function(A11);                       \
  function(A12);                       \
  function(A13);                       \
  function(A22);                       \
  function(A23);                       \
  function(A33);                       \
  function(DIFFK);                     \
  function(Gamma1);                    \
  function(Gamma2);                    \
  function(Gamma3);                    \
  function(DIFFalpha);                 \
  Z4c_APPLY_TO_FIELDS(function)        \
  BSSN_APPLY_TO_SHIFT(function)

#define BSSN_APPLY_TO_SOURCES(function) \
  function(DIFFr);                      \
  function(DIFFS);                      \
  function(S1);                         \
  function(S2);                         \
  function(S3);                         \
  function(STF11);                      \
  function(STF12);                      \
  function(STF13);                      \
  function(STF22);                      \
  function(STF23);                      \
  function(STF33);

#define BSSN_APPLY_TO_GEN1_EXTRAS(function) \
  function(KDx);                            \
  function(KDy);                            \
  function(KDz);                            \
  function(ricci);                          \
  function(AijAij);

#define BSSN_APPLY_TO_IJ_PERMS(function) \
  function(1, 1);                        \
  function(1, 2);                        \
  function(1, 3);                        \
  function(2, 2);                        \
  function(2, 3);                        \
  function(3, 3);

// apply when the "I" index is special, e.g., d_i g_{jk}
#define BSSN_APPLY_TO_IJK_PERMS(function)   \
  function(1, 1, 1);                  \
  function(1, 1, 2);                  \
  function(1, 1, 3);                  \
  function(1, 2, 2);                  \
  function(1, 2, 3);                  \
  function(1, 3, 3);                  \
  function(2, 1, 1);                  \
  function(2, 1, 2);                  \
  function(2, 1, 3);                  \
  function(2, 2, 2);                  \
  function(2, 2, 3);                  \
  function(2, 3, 3);                  \
  function(3, 1, 1);                  \
  function(3, 1, 2);                  \
  function(3, 1, 3);                  \
  function(3, 2, 2);                  \
  function(3, 2, 3);                  \
  function(3, 3, 3);


#define BSSN_SWAP_REGISTERS(field, reg_prefix_1, reg_prefix_2) \
  cosmoArraySwap(field##reg_prefix_1, field##reg_prefix_2);

#define BSSN_SWAP_ARRAYS(reg_prefix_1, reg_prefix_2) \
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_SWAP_REGISTERS, reg_prefix_1, reg_prefix_2)



#define BSSN_COPY_FIELD(field, reg_prefix_from, reg_prefix_to) \
  field##reg_prefix_to = field##reg_prefix_from;

#define BSSN_COPY_ARRAYS(reg_prefix_from, reg_prefix_to) \
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_COPY_FIELD, reg_prefix_from, reg_prefix_to)



#define BSSN_ZERO_FIELD(field, reg, idx) \
  field##reg[idx] = 0.0;

#define BSSN_ZERO_ARRAYS(reg, idx) \
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_ZERO_FIELD, reg, idx)



// macro requires idx to be set
#define BSSN_ZERO_GEN1_FIELD(field) \
  field##_a[idx] = 0.0;

// macro requires idx to be set
#define BSSN_ZERO_GEN1_EXTRAS() \
  BSSN_APPLY_TO_GEN1_EXTRAS(BSSN_ZERO_GEN1_FIELD)

// macro requires idx to be set
#define BSSN_ZERO_SOURCES() \
  BSSN_APPLY_TO_SOURCES(BSSN_ZERO_GEN1_FIELD)



// arr_c[idx] = arr_p[idx] + dt*mult*evfn(arr_a);
#define BSSN_COMPUTE_RK_STEP_FIELD(field, mult) \
  field##_c[paq->idx] = field##_p[paq->idx] + dt*mult*ev_##field(paq);

#define BSSN_COMPUTE_RK_STEP(mult) \
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_COMPUTE_RK_STEP_FIELD, mult)



#define BSSN_ADD_C_TO_F_FIELD(field, mult) \
  field##_f[paq->idx] += mult*field##_c[paq->idx];

#define BSSN_ADD_C_TO_F(mult) \
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_ADD_C_TO_F_FIELD, mult)


// arr_f = (1.0/3.0)*(arr_f - arr_p) + (1.0/6.0)*evfn(arr_a)
#define BSSN_FINAL_RK4_STEP() \
  BSSN_APPLY_TO_FIELDS(BSSN_FINAL_RK4_STEP_FIELD)

#define BSSN_FINAL_RK4_STEP_FIELD(field) \
  field##_f[paq->idx] = (1.0/3.0)*(field##_f[paq->idx] - field##_p[paq->idx]) + (1.0/6.0)*dt*ev_##field(paq);



// common definition for RK step
#define BSSN_RK_PERFORM_KN_CALC(N) \
  idx_t i, j, k; \
  _Pragma("omp parallel for default(shared) private(i, j, k)") \
  LOOP3(i, j, k) \
  { \
    BSSNData b_paq = {0}; \
    K##N##CalcPt(i, j, k, &b_paq); \
  }



/*
 * Aux. variable calculations
 */

#define BSSN_CALCULATE_CHRISTOFFEL(I, J, K) paq->G##I##J##K = 0.5*( \
    paq->gammai##I##1 * (paq->d##J##g##K##1 + paq->d##K##g##J##1 - paq->d1g##J##K) + \
    paq->gammai##I##2 * (paq->d##J##g##K##2 + paq->d##K##g##J##2 - paq->d2g##J##K) + \
    paq->gammai##I##3 * (paq->d##J##g##K##3 + paq->d##K##g##J##3 - paq->d3g##J##K) \
  )

#define BSSN_CALCULATE_CHRISTOFFEL_LOWER(I, J, K) paq->GL##I##J##K = 0.5*( \
    paq->d##J##g##K##I + paq->d##K##g##J##I - paq->d##I##g##J##K \
  )

#define BSSN_CALCULATE_DGAMMA(I, J, K) paq->d##I##g##J##K = derivative(paq->i, paq->j, paq->k, I, DIFFgamma##J##K##_a);

#define BSSN_CALCULATE_ACONT(I, J) paq->Acont##I##J = ( \
    paq->gammai##I##1*paq->gammai##J##1*paq->A11 + paq->gammai##I##2*paq->gammai##J##1*paq->A21 + paq->gammai##I##3*paq->gammai##J##1*paq->A31 \
    + paq->gammai##I##1*paq->gammai##J##2*paq->A12 + paq->gammai##I##2*paq->gammai##J##2*paq->A22 + paq->gammai##I##3*paq->gammai##J##2*paq->A32 \
    + paq->gammai##I##1*paq->gammai##J##3*paq->A13 + paq->gammai##I##2*paq->gammai##J##3*paq->A23 + paq->gammai##I##3*paq->gammai##J##3*paq->A33 \
  );

// needs the gamma*ldlphi vars defined:
// not actually trace free yet!
#define BSSN_CALCULATE_DIDJALPHA(I, J) paq->D##I##D##J##aTF = double_derivative(paq->i, paq->j, paq->k, I, J, DIFFalpha_a) - ( \
    (paq->G1##I##J + 2.0*( (1==I)*paq->d##J##phi + (1==J)*paq->d##I##phi - paq->gamma##I##J*gammai1ldlphi))*paq->d1a + \
    (paq->G2##I##J + 2.0*( (2==I)*paq->d##J##phi + (2==J)*paq->d##I##phi - paq->gamma##I##J*gammai2ldlphi))*paq->d2a + \
    (paq->G3##I##J + 2.0*( (3==I)*paq->d##J##phi + (3==J)*paq->d##I##phi - paq->gamma##I##J*gammai3ldlphi))*paq->d3a \
  );

#define BSSN_CALCULATE_RICCI_UNITARY(I, J) paq->ricci##I##J = ( \
    - 0.5*( \
      paq->gammai11*paq->d1d1g##I##J + paq->gammai22*paq->d2d2g##I##J + paq->gammai33*paq->d3d3g##I##J \
      + 2.0*(paq->gammai12*paq->d1d2g##I##J + paq->gammai13*paq->d1d3g##I##J + paq->gammai23*paq->d2d3g##I##J) \
    ) \
    + 0.5*( \
      paq->gamma1##I*derivative(paq->i, paq->j, paq->k, J, Gamma1_a) + paq->gamma2##I*derivative(paq->i, paq->j, paq->k, J, Gamma2_a) + paq->gamma3##I*derivative(paq->i, paq->j, paq->k, J, Gamma3_a) + \
      paq->gamma1##J*derivative(paq->i, paq->j, paq->k, I, Gamma1_a) + paq->gamma2##J*derivative(paq->i, paq->j, paq->k, I, Gamma2_a) + paq->gamma3##J*derivative(paq->i, paq->j, paq->k, I, Gamma3_a) \
    ) \
    + 0.5*( \
      paq->Gammad1*paq->GL##I##J##1 + paq->Gammad2*paq->GL##I##J##2 + paq->Gammad3*paq->GL##I##J##3 \
      + paq->Gammad1*paq->GL##J##I##1 + paq->Gammad2*paq->GL##J##I##2 + paq->Gammad3*paq->GL##J##I##3 \
    ) \
    + paq->gammai11*( \
        paq->G11##I*paq->GL##J##11 + paq->G21##I*paq->GL##J##21 + paq->G31##I*paq->GL##J##31 \
        + paq->G11##J*paq->GL##I##11 + paq->G21##J*paq->GL##I##21 + paq->G31##J*paq->GL##I##31 \
        + paq->G1##I##1*paq->GL11##J + paq->G2##I##1*paq->GL21##J + paq->G3##I##1*paq->GL31##J \
      ) \
    + paq->gammai12*( \
        paq->G11##I*paq->GL##J##12 + paq->G21##I*paq->GL##J##22 + paq->G31##I*paq->GL##J##32 \
        + paq->G11##J*paq->GL##I##12 + paq->G21##J*paq->GL##I##22 + paq->G31##J*paq->GL##I##32 \
        + paq->G1##I##2*paq->GL11##J + paq->G2##I##2*paq->GL21##J + paq->G3##I##2*paq->GL31##J \
      ) \
    + paq->gammai13*( \
        paq->G11##I*paq->GL##J##13 + paq->G21##I*paq->GL##J##23 + paq->G31##I*paq->GL##J##33 \
        + paq->G11##J*paq->GL##I##13 + paq->G21##J*paq->GL##I##23 + paq->G31##J*paq->GL##I##33 \
        + paq->G1##I##3*paq->GL11##J + paq->G2##I##3*paq->GL21##J + paq->G3##I##3*paq->GL31##J \
      ) \
    \
    + paq->gammai21*( \
        paq->G12##I*paq->GL##J##11 + paq->G22##I*paq->GL##J##21 + paq->G32##I*paq->GL##J##31 \
        + paq->G12##J*paq->GL##I##11 + paq->G22##J*paq->GL##I##21 + paq->G32##J*paq->GL##I##31 \
        + paq->G1##I##1*paq->GL12##J + paq->G2##I##1*paq->GL22##J + paq->G3##I##1*paq->GL32##J \
      ) \
    + paq->gammai22*( \
        paq->G12##I*paq->GL##J##12 + paq->G22##I*paq->GL##J##22 + paq->G32##I*paq->GL##J##32 \
        + paq->G12##J*paq->GL##I##12 + paq->G22##J*paq->GL##I##22 + paq->G32##J*paq->GL##I##32 \
        + paq->G1##I##2*paq->GL12##J + paq->G2##I##2*paq->GL22##J + paq->G3##I##2*paq->GL32##J \
      ) \
    + paq->gammai23*( \
        paq->G12##I*paq->GL##J##13 + paq->G22##I*paq->GL##J##23 + paq->G32##I*paq->GL##J##33 \
        + paq->G12##J*paq->GL##I##13 + paq->G22##J*paq->GL##I##23 + paq->G32##J*paq->GL##I##33 \
        + paq->G1##I##3*paq->GL12##J + paq->G2##I##3*paq->GL22##J + paq->G3##I##3*paq->GL32##J \
      ) \
    \
    + paq->gammai31*( \
        paq->G13##I*paq->GL##J##11 + paq->G23##I*paq->GL##J##21 + paq->G33##I*paq->GL##J##31 \
        + paq->G13##J*paq->GL##I##11 + paq->G23##J*paq->GL##I##21 + paq->G33##J*paq->GL##I##31 \
        + paq->G1##I##1*paq->GL13##J + paq->G2##I##1*paq->GL23##J + paq->G3##I##1*paq->GL33##J \
      ) \
    + paq->gammai32*( \
        paq->G13##I*paq->GL##J##12 + paq->G23##I*paq->GL##J##22 + paq->G33##I*paq->GL##J##32 \
        + paq->G13##J*paq->GL##I##12 + paq->G23##J*paq->GL##I##22 + paq->G33##J*paq->GL##I##32 \
        + paq->G1##I##2*paq->GL13##J + paq->G2##I##2*paq->GL23##J + paq->G3##I##2*paq->GL33##J \
      ) \
    + paq->gammai33*( \
        paq->G13##I*paq->GL##J##13 + paq->G23##I*paq->GL##J##23 + paq->G33##I*paq->GL##J##33 \
        + paq->G13##J*paq->GL##I##13 + paq->G23##J*paq->GL##I##23 + paq->G33##J*paq->GL##I##33 \
        + paq->G1##I##3*paq->GL13##J + paq->G2##I##3*paq->GL23##J + paq->G3##I##3*paq->GL33##J \
      ) \
  );


#define BSSN_CALCULATE_DIDJGAMMA_PERMS(I, J)           \
  paq->d##I##d##J##g11 = double_derivative(paq->i, paq->j, paq->k, I, J, DIFFgamma11_a); \
  paq->d##I##d##J##g12 = double_derivative(paq->i, paq->j, paq->k, I, J, DIFFgamma12_a); \
  paq->d##I##d##J##g13 = double_derivative(paq->i, paq->j, paq->k, I, J, DIFFgamma13_a); \
  paq->d##I##d##J##g22 = double_derivative(paq->i, paq->j, paq->k, I, J, DIFFgamma22_a); \
  paq->d##I##d##J##g23 = double_derivative(paq->i, paq->j, paq->k, I, J, DIFFgamma23_a); \
  paq->d##I##d##J##g33 = double_derivative(paq->i, paq->j, paq->k, I, J, DIFFgamma33_a)


/*
 * Evolution equations for indexed components
 */

#define BSSN_DT_DIFFGAMMAIJ(I, J) ( \
    - 2.0*paq->alpha*paq->A##I##J \
    + paq->gamma##I##1*paq->d##J##beta1 + paq->gamma##I##2*paq->d##J##beta2 + paq->gamma##I##3*paq->d##J##beta3 \
    + paq->gamma##J##1*paq->d##I##beta1 + paq->gamma##J##2*paq->d##I##beta2 + paq->gamma##J##3*paq->d##I##beta3 \
    - (2.0/3.0)*paq->gamma##I##J*(paq->d1beta1 + paq->d2beta2 + paq->d3beta3) \
  )

#define BSSN_DT_AIJ(I, J) ( \
    exp(-4.0*paq->phi)*( paq->alpha*(paq->ricciTF##I##J - 8.0*PI*paq->STF##I##J) - paq->D##I##D##J##aTF ) \
    + paq->alpha*((paq->K + 2.0*paq->theta)*paq->A##I##J - 2.0*( \
        paq->gammai11*paq->A1##I*paq->A1##J + paq->gammai12*paq->A1##I*paq->A2##J + paq->gammai13*paq->A1##I*paq->A3##J \
        + paq->gammai21*paq->A2##I*paq->A1##J + paq->gammai22*paq->A2##I*paq->A2##J + paq->gammai23*paq->A2##I*paq->A3##J \
        + paq->gammai31*paq->A3##I*paq->A1##J + paq->gammai32*paq->A3##I*paq->A2##J + paq->gammai33*paq->A3##I*paq->A3##J \
      )) \
    + paq->beta1*derivative(paq->i, paq->j, paq->k, 1, A##I##J##_a) \
    + paq->beta2*derivative(paq->i, paq->j, paq->k, 2, A##I##J##_a) \
    + paq->beta3*derivative(paq->i, paq->j, paq->k, 3, A##I##J##_a) \
    + paq->A##I##1*paq->d##J##beta1 + paq->A##I##2*paq->d##J##beta2 + paq->A##I##3*paq->d##J##beta3 \
    + paq->A##J##1*paq->d##I##beta1 + paq->A##J##2*paq->d##I##beta2 + paq->A##J##3*paq->d##I##beta3 \
    - (2.0/3.0)*paq->A##I##J*(paq->d1beta1 + paq->d2beta2 + paq->d3beta3) \
  )

#define BSSN_DT_GAMMAI(I) (BSSN_DT_GAMMAI_NOSHIFT(I) + BSSN_DT_GAMMAI_SHIFT(I))

#define BSSN_DT_GAMMAI_NOSHIFT(I) ( \
    - 2.0*(paq->Acont##I##1*paq->d1a + paq->Acont##I##2*paq->d2a + paq->Acont##I##3*paq->d3a) \
    + 2.0*paq->alpha*( \
        paq->G##I##11*paq->Acont11 + paq->G##I##22*paq->Acont22 + paq->G##I##33*paq->Acont33 \
          + 2.0*(paq->G##I##12*paq->Acont12 + paq->G##I##13*paq->Acont13 + paq->G##I##23*paq->Acont23) \
        - (1.0/3.0) * ( \
            2.0*(paq->gammai##I##1*paq->d1K + paq->gammai##I##2*paq->d2K + paq->gammai##I##3*paq->d3K) \
            + paq->gammai##I##1*paq->d1theta + paq->gammai##I##2*paq->d2theta + paq->gammai##I##3*paq->d3theta \
          ) \
        - 2.0*paq->alpha*Z4c_K1_DAMPING_AMPLITUDE*( \
            paq->Gamma##I - paq->Gammad##I \
          ) \
        - 8.0*PI*(paq->gammai##I##1*paq->S1 + paq->gammai##I##2*paq->S2 + paq->gammai##I##3*paq->S3) \
        + 6.0 * (paq->Acont##I##1*paq->d1phi + paq->Acont##I##2*paq->d2phi + paq->Acont##I##3*paq->d3phi) \
      ) \
    )

#if USE_BSSN_SHIFT
#define BSSN_DT_GAMMAI_SHIFT(I) ( \
    + paq->beta1*derivative(paq->i, paq->j, paq->k, 1, Gamma##I##_a) \
    + paq->beta2*derivative(paq->i, paq->j, paq->k, 2, Gamma##I##_a) \
    + paq->beta3*derivative(paq->i, paq->j, paq->k, 3, Gamma##I##_a) \
    - paq->Gamma1*paq->d1beta##I + paq->Gamma2*paq->d2beta##I + paq->Gamma3*paq->d3beta##I \
    + (2.0/3.0) * paq->Gamma##I * (paq->d1beta1 + paq->d2beta2 + paq->d3beta3) \
    + (1.0/3.0) * ( \
        paq->gammai##I##1*double_derivative(paq->i, paq->j, paq->k, 1, 1, beta1_a) + paq->gammai##I##1*double_derivative(paq->i, paq->j, paq->k, 2, 1, beta2_a) + paq->gammai##I##1*double_derivative(paq->i, paq->j, paq->k, 3, 1, beta3_a) +  \
        paq->gammai##I##2*double_derivative(paq->i, paq->j, paq->k, 1, 2, beta1_a) + paq->gammai##I##2*double_derivative(paq->i, paq->j, paq->k, 2, 2, beta2_a) + paq->gammai##I##2*double_derivative(paq->i, paq->j, paq->k, 3, 2, beta3_a) +  \
        paq->gammai##I##3*double_derivative(paq->i, paq->j, paq->k, 1, 3, beta1_a) + paq->gammai##I##3*double_derivative(paq->i, paq->j, paq->k, 2, 3, beta2_a) + paq->gammai##I##3*double_derivative(paq->i, paq->j, paq->k, 3, 3, beta3_a) \
      ) \
    + ( \
        paq->gammai11*double_derivative(paq->i, paq->j, paq->k, 1, 1, beta##I##_a) + paq->gammai22*double_derivative(paq->i, paq->j, paq->k, 2, 2, beta##I##_a) + paq->gammai33*double_derivative(paq->i, paq->j, paq->k, 3, 3, beta##I##_a) \
        + 2.0*(paq->gammai12*double_derivative(paq->i, paq->j, paq->k, 1, 2, beta##I##_a) + paq->gammai13*double_derivative(paq->i, paq->j, paq->k, 1, 3, beta##I##_a) + paq->gammai23*double_derivative(paq->i, paq->j, paq->k, 2, 3, beta##I##_a)) \
      ) \
  )
#else
#define BSSN_DT_GAMMAI_SHIFT(I) 0.0
#endif

#define BSSN_MI(I) exp(6.0*paq->phi)*( \
    - 2.0/3.0*derivative(paq->i, paq->j, paq->k, I, DIFFK_a) \
    - 8*PI*(paq->gammai##I##1*paq->S1 + paq->gammai##I##2*paq->S2+ paq->gammai##I##3*paq->S3) \
    - 2.0/3.0*2.0*paq->d##I##theta \
    + 6.0*( \
      paq->gammai11*paq->A1##I*paq->d1phi + paq->gammai21*paq->A2##I*paq->d1phi + paq->gammai31*paq->A3##I*paq->d1phi \
      + paq->gammai12*paq->A1##I*paq->d2phi + paq->gammai22*paq->A2##I*paq->d2phi + paq->gammai32*paq->A3##I*paq->d2phi \
      + paq->gammai13*paq->A1##I*paq->d3phi + paq->gammai23*paq->A2##I*paq->d3phi + paq->gammai33*paq->A3##I*paq->d3phi \
    ) + ( \
      /* (gamma^jk D_j A_ki) */ \
      paq->gammai11*derivative(paq->i, paq->j, paq->k, 1, A1##I##_a) + paq->gammai12*derivative(paq->i, paq->j, paq->k, 2, A1##I##_a) + paq->gammai13*derivative(paq->i, paq->j, paq->k, 3, A1##I##_a) \
      + paq->gammai21*derivative(paq->i, paq->j, paq->k, 1, A2##I##_a) + paq->gammai22*derivative(paq->i, paq->j, paq->k, 2, A2##I##_a) + paq->gammai23*derivative(paq->i, paq->j, paq->k, 3, A2##I##_a) \
      + paq->gammai31*derivative(paq->i, paq->j, paq->k, 1, A3##I##_a) + paq->gammai32*derivative(paq->i, paq->j, paq->k, 2, A3##I##_a) + paq->gammai33*derivative(paq->i, paq->j, paq->k, 3, A3##I##_a) \
      - paq->Gamma1*paq->A1##I - paq->Gamma2*paq->A2##I - paq->Gamma3*paq->A3##I \
      - paq->GL11##I*paq->Acont11 - paq->GL21##I*paq->Acont21 - paq->GL31##I*paq->Acont31 \
      - paq->GL12##I*paq->Acont12 - paq->GL22##I*paq->Acont22 - paq->GL32##I*paq->Acont32 \
      - paq->GL13##I*paq->Acont13 - paq->GL23##I*paq->Acont23 - paq->GL33##I*paq->Acont33 \
    ) \
  )

#define BSSN_MI_SCALE(I) exp(6.0*paq->phi)*( \
    fabs(2.0/3.0*derivative(paq->i, paq->j, paq->k, I, DIFFK_a)) \
    + fabs(8*PI*(paq->gammai##I##1*paq->S1 + paq->gammai##I##2*paq->S2+ paq->gammai##I##3*paq->S3)) \
    + 6.0*fabs( \
      paq->gammai11*paq->A1##I*paq->d1phi + paq->gammai21*paq->A2##I*paq->d1phi + paq->gammai31*paq->A3##I*paq->d1phi \
      + paq->gammai12*paq->A1##I*paq->d2phi + paq->gammai22*paq->A2##I*paq->d2phi + paq->gammai32*paq->A3##I*paq->d2phi \
      + paq->gammai13*paq->A1##I*paq->d3phi + paq->gammai23*paq->A2##I*paq->d3phi + paq->gammai33*paq->A3##I*paq->d3phi \
    ) + fabs( \
      /* (gamma^jk D_j A_ki) */ \
      paq->gammai11*derivative(paq->i, paq->j, paq->k, 1, A1##I##_a) + paq->gammai12*derivative(paq->i, paq->j, paq->k, 2, A1##I##_a) + paq->gammai13*derivative(paq->i, paq->j, paq->k, 3, A1##I##_a) \
      + paq->gammai21*derivative(paq->i, paq->j, paq->k, 1, A2##I##_a) + paq->gammai22*derivative(paq->i, paq->j, paq->k, 2, A2##I##_a) + paq->gammai23*derivative(paq->i, paq->j, paq->k, 3, A2##I##_a) \
      + paq->gammai31*derivative(paq->i, paq->j, paq->k, 1, A3##I##_a) + paq->gammai32*derivative(paq->i, paq->j, paq->k, 2, A3##I##_a) + paq->gammai33*derivative(paq->i, paq->j, paq->k, 3, A3##I##_a) \
      - paq->Gamma1*paq->A1##I - paq->Gamma2*paq->A2##I - paq->Gamma3*paq->A3##I \
      - paq->GL11##I*paq->Acont11 - paq->GL21##I*paq->Acont21 - paq->GL31##I*paq->Acont31 \
      - paq->GL12##I*paq->Acont12 - paq->GL22##I*paq->Acont22 - paq->GL32##I*paq->Acont32 \
      - paq->GL13##I*paq->Acont13 - paq->GL23##I*paq->Acont23 - paq->GL33##I*paq->Acont33 \
    ) \
  )

/*
 * Full metric calcs
 */

// metric derivs

#define SET_DKM00(k) paq->d##k##m00 = DKM00(k);
#define DKM00(k) \
      -2.0*paq->alpha*paq->d##k##a \
      + paq->d##k##m11*paq->beta1*paq->beta1 + paq->d##k##m22*paq->beta2*paq->beta2 + paq->d##k##m33*paq->beta3*paq->beta3 \
      + 2.0*(paq->d##k##m12*paq->beta1*paq->beta2 + paq->d##k##m13*paq->beta1*paq->beta3 + paq->d##k##m23*paq->beta2*paq->beta3) \
      + 2.0*exp(4.0*paq->phi)*( \
          paq->gamma11*paq->d##k##beta1*paq->beta1 + paq->gamma12*paq->d##k##beta2*paq->beta1 + paq->gamma13*paq->d##k##beta3*paq->beta1 \
          + paq->gamma12*paq->d##k##beta1*paq->beta2 + paq->gamma22*paq->d##k##beta2*paq->beta2 + paq->gamma23*paq->d##k##beta3*paq->beta2 \
          + paq->gamma13*paq->d##k##beta1*paq->beta3 + paq->gamma23*paq->d##k##beta2*paq->beta3 + paq->gamma33*paq->d##k##beta3*paq->beta3 \
        );


#define SET_DKM0I(k, i) paq->d##k##m0##i = DKM0I(k, i);
#define DKM0I(k, i) \
      exp(4.0*paq->phi)*( \
        paq->gamma1##i * paq->d##k##beta##1 + paq->gamma2##i * paq->d##k##beta##2 + paq->gamma3##i * paq->d##k##beta##3 \
        + (paq->d##k##g1##i + 4.0*paq->d##k##phi*paq->gamma1##i) * paq->beta##1   \
        + (paq->d##k##g2##i + 4.0*paq->d##k##phi*paq->gamma2##i) * paq->beta##2 \
        + (paq->d##k##g3##i + 4.0*paq->d##k##phi*paq->gamma3##i) * paq->beta##3 \
      );

#define SET_DKMIJ(k, i, j) paq->d##k##m##i##j = DKMIJ(k, i, j);
#define DKMIJ(k, i, j) exp(4.0*paq->phi)*(paq->d##k##g##i##j + 4.0*paq->d##k##phi*paq->gamma##i##j);


// metric

#define SET_M00() paq->m00 = M00();
#define M00() \
      -paq->alpha*paq->alpha + exp(4.0*paq->phi)*( \
        paq->gamma11*paq->beta1*paq->beta1 + paq->gamma22*paq->beta2*paq->beta2 + paq->gamma33*paq->beta3*paq->beta3 \
        + 2.0*(paq->gamma12*paq->beta1*paq->beta2 + paq->gamma13*paq->beta1*paq->beta3 + paq->gamma23*paq->beta2*paq->beta3) \
      );

#define SET_M0I(i) paq->m0##i = M0I(i);
#define M0I(i) exp(4.0*paq->phi)*(paq->gamma1##i*paq->beta1 + paq->gamma2##i*paq->beta2 + paq->gamma3##i*paq->beta3);

#define SET_MIJ(i, j) paq->m##i##j = MIJ(i, j);
#define MIJ(i, j) exp(4.0*paq->phi)*(paq->gamma##i##j);

// inverse metric

#define SET_Mi00() paq->mi00 = Mi00();
#define Mi00() -1.0/paq->alpha/paq->alpha;

#define SET_Mi0I(i) paq->mi0##i = Mi0I(i);
#define Mi0I(i) 1.0/paq->alpha/paq->alpha*paq->beta##i;

#define SET_MiIJ(i, j) paq->mi##i##j = MiIJ(i, j);
#define MiIJ(i, j) exp(-4.0*paq->phi)*(paq->gamma##i##j) - 1.0/paq->alpha/paq->alpha*paq->beta##i*paq->beta##j;


/*
 * Conversion to ADM quantities for raytracing
 */
#define BSSN_RP_DG(I,J,L) \
  P*(4.0*paq->gamma##I##J*paq->d##L##phi+ paq->d##L##g##I##J)

#define BSSN_RP_DDG(I,J,L,M) \
  P*( \
    16.0*paq->gamma##I##J*paq->d##M##phi*paq->d##L##phi \
    + 4.0*(paq->d##L##g##I##J*paq->d##M##phi + paq->d##M##g##I##J*paq->d##L##phi + paq->gamma##I##J*paq->d##L##d##M##phi) \
    + paq->d##L##d##M##g##I##J \
  )

#define BSSN_RP_K(I,J) \
  P*(paq->A##I##J + 1.0/3.0*paq->gamma##I##J*paq->K)

#define BSSN_RP_DK(I,J,L) \
  P*( \
    4.0*(paq->A##I##J + 1.0/3.0*paq->gamma##I##J*paq->K)*paq->d##L##phi + derivative(paq->i, paq->j, paq->k, L, A##I##J##_a) \
    + 1.0/3.0*paq->K*paq->d##L##g##I##J + 1.0/3.0*paq->gamma##I##J*paq->d##L##K \
  )

#define BSSN_RP_GAMMA(I,J,K) \
  paq->G##I##J##K + 2.0*( \
      (I==J ? 1.0 : 0.0)*paq->d##K##phi \
      + (I==K ? 1.0 : 0.0)*paq->d##J##phi \
      - paq->gamma##J##K*(paq->gammai##I##1*paq->d1phi + paq->gammai##I##2*paq->d2phi + paq->gammai##I##3*paq->d3phi) \
    )

#define BSSN_RP_GAMMAL(I,J,K) P*( \
  paq->GL##I##J##K + 2.0*( \
      paq->gamma##I##J*paq->d##K##phi \
      + paq->gamma##I##K*paq->d##J##phi \
      - paq->gamma##J##K*paq->d##I##phi \
  ) )

/*
 * Enforce standard ordering of indexes for tensor components
 */

// actual fields:
#define gamma21 gamma12
#define gamma31 gamma13
#define gamma32 gamma23
#define gammai21 gammai12
#define gammai31 gammai13
#define gammai32 gammai23
#define A21 A12
#define A31 A13
#define A32 A23

#define DIFFgamma21 DIFFgamma12
#define DIFFgamma31 DIFFgamma13
#define DIFFgamma32 DIFFgamma23
#define DIFFgammai21 DIFFgammai12
#define DIFFgammai31 DIFFgammai13
#define DIFFgammai32 DIFFgammai23

#define DIFFgamma21_a DIFFgamma12_a
#define DIFFgamma31_a DIFFgamma13_a
#define DIFFgamma32_a DIFFgamma23_a
#define A21_a A12_a
#define A31_a A13_a
#define A32_a A23_a

#define A21_adj A12_adj
#define A31_adj A13_adj
#define A32_adj A23_adj

#define A21_adj_ext A12_adj_ext
#define A31_adj_ext A13_adj_ext
#define A32_adj_ext A23_adj_ext

// local variables:
// ricci tensor
#define ricciTF21 ricciTF12
#define ricciTF31 ricciTF13
#define ricciTF32 ricciTF23
#define ricci21 ricci12
#define ricci31 ricci13
#define ricci32 ricci23

// covariant double-derivatives of phi
#define D2D1phi D1D2phi
#define D3D1phi D1D3phi
#define D3D2phi D2D3phi

// covariant double-derivatives of alpha
#define D2D1aTF D1D2aTF
#define D3D1aTF D1D3aTF
#define D3D2aTF D2D3aTF

// Inverse ext. curvature
#define Acont21 Acont12
#define Acont31 Acont13
#define Acont32 Acont23

// christoffel symbols
#define G121 G112
#define G131 G113
#define G132 G123
#define G221 G212
#define G231 G213
#define G232 G223
#define G321 G312
#define G331 G313
#define G332 G323

// lower christoffel symbols
#define GL121 GL112
#define GL131 GL113
#define GL132 GL123
#define GL221 GL212
#define GL231 GL213
#define GL232 GL223
#define GL321 GL312
#define GL331 GL313
#define GL332 GL323

// Metric derivatives
#define d1g21 d1g12
#define d1g31 d1g13
#define d1g32 d1g23
#define d2g21 d2g12
#define d2g31 d2g13
#define d2g32 d2g23
#define d3g21 d3g12
#define d3g31 d3g13
#define d3g32 d3g23

// second derivatives of the metric
// bad metric indices
#define d1d1g21 d1d1g12
#define d1d1g31 d1d1g13
#define d1d1g32 d1d1g23
#define d1d2g21 d1d2g12
#define d1d2g31 d1d2g13
#define d1d2g32 d1d2g23
#define d1d3g21 d1d3g12
#define d1d3g31 d1d3g13
#define d1d3g32 d1d3g23
#define d2d2g21 d2d2g12
#define d2d2g31 d2d2g13
#define d2d2g32 d2d2g23
#define d2d3g21 d2d3g12
#define d2d3g31 d2d3g13
#define d2d3g32 d2d3g23
#define d3d3g21 d3d3g12
#define d3d3g31 d3d3g13
#define d3d3g32 d3d3g23
// bad derivative indices
#define d2d1g11 d1d2g11
#define d2d1g12 d1d2g12
#define d2d1g13 d1d2g13
#define d2d1g22 d1d2g22
#define d2d1g23 d1d2g23
#define d2d1g33 d1d2g33
#define d3d1g11 d1d3g11
#define d3d1g12 d1d3g12
#define d3d1g13 d1d3g13
#define d3d1g22 d1d3g22
#define d3d1g23 d1d3g23
#define d3d1g33 d1d3g33
#define d3d2g11 d2d3g11
#define d3d2g12 d2d3g12
#define d3d2g13 d2d3g13
#define d3d2g22 d2d3g22
#define d3d2g23 d2d3g23
#define d3d2g33 d2d3g33
// bad in both indices
#define d2d1g21 d1d2g12
#define d2d1g31 d1d2g13
#define d2d1g32 d1d2g23
#define d3d1g21 d1d3g12
#define d3d1g31 d1d3g13
#define d3d1g32 d1d3g23
#define d3d2g21 d2d3g12
#define d3d2g31 d2d3g13
#define d3d2g32 d2d3g23

// Full Metric
#define m10 m01
#define m20 m02
#define m30 m03
#define m21 m12
#define m31 m13
#define m32 m23

// Full Metric derivatives
#define d1m10 d1m01
#define d1m20 d1m02
#define d1m30 d1m03
#define d1m21 d1m12
#define d1m31 d1m13
#define d1m32 d1m23
#define d2m10 d2m01
#define d2m20 d2m02
#define d2m30 d2m03
#define d2m21 d2m12
#define d2m31 d2m13
#define d2m32 d2m23
#define d3m10 d3m01
#define d3m20 d3m02
#define d3m30 d3m03
#define d3m21 d3m12
#define d3m31 d3m13
#define d3m32 d3m23

// source terms
#define S21 S12
#define S31 S13
#define S32 S23

#define STF21 STF12
#define STF31 STF13
#define STF32 STF23

#endif
