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


#define BSSN_ZERO_FIELD(field, reg, idx) \
  field->_array##reg[idx] = 0.0;

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


// Initialize all fields
#define BSSN_RK_INITIALIZE_FIELD(field) \
  field->stepInit();

#define BSSN_RK_INITIALIZE \
  BSSN_APPLY_TO_FIELDS(BSSN_RK_INITIALIZE_FIELD)


// Evolve all fields
#define BSSN_RK_EVOLVE_PT_FIELD(field) \
  field->_array_c[bd->idx] = ev_##field(bd);

#define BSSN_RK_EVOLVE_PT \
  BSSN_APPLY_TO_FIELDS(BSSN_RK_EVOLVE_PT_FIELD)


// Finalize an RK step for all BSSN fields
#define BSSN_FINALIZE_K_FIELD(field, n) \
  field->K##n##Finalize();

#define BSSN_FINALIZE_K(n) \
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_FINALIZE_K_FIELD, n)


// Initialize all fields
#define BSSN_SET_DT_FIELD(field, dt) \
  field->setDt(dt)

#define BSSN_SET_DT(dt) \
  BSSN_APPLY_TO_FIELDS_ARGS(BSSN_SET_DT_FIELD, dt)


/*
 * Aux. variable calculations
 */

#define BSSN_CALCULATE_CHRISTOFFEL(I, J, K) bd->G##I##J##K = 0.5*( \
    bd->gammai##I##1 * (bd->d##J##g##K##1 + bd->d##K##g##J##1 - bd->d1g##J##K) + \
    bd->gammai##I##2 * (bd->d##J##g##K##2 + bd->d##K##g##J##2 - bd->d2g##J##K) + \
    bd->gammai##I##3 * (bd->d##J##g##K##3 + bd->d##K##g##J##3 - bd->d3g##J##K) \
  )

#define BSSN_CALCULATE_CHRISTOFFEL_LOWER(I, J, K) bd->GL##I##J##K = 0.5*( \
    bd->d##J##g##K##I + bd->d##K##g##J##I - bd->d##I##g##J##K \
  )

#define BSSN_CALCULATE_DGAMMA(I, J, K) bd->d##I##g##J##K = derivative(bd->i, bd->j, bd->k, I, DIFFgamma##J##K->_array_a);

#define BSSN_CALCULATE_ACONT(I, J) bd->Acont##I##J = ( \
    bd->gammai##I##1*bd->gammai##J##1*bd->A11 + bd->gammai##I##2*bd->gammai##J##1*bd->A21 + bd->gammai##I##3*bd->gammai##J##1*bd->A31 \
    + bd->gammai##I##1*bd->gammai##J##2*bd->A12 + bd->gammai##I##2*bd->gammai##J##2*bd->A22 + bd->gammai##I##3*bd->gammai##J##2*bd->A32 \
    + bd->gammai##I##1*bd->gammai##J##3*bd->A13 + bd->gammai##I##2*bd->gammai##J##3*bd->A23 + bd->gammai##I##3*bd->gammai##J##3*bd->A33 \
  );

// needs the gamma*ldlphi vars defined:
// not actually trace free yet!
#define BSSN_CALCULATE_DIDJALPHA(I, J) bd->D##I##D##J##aTF = double_derivative(bd->i, bd->j, bd->k, I, J, DIFFalpha->_array_a) - ( \
    (bd->G1##I##J + 2.0*( (1==I)*bd->d##J##phi + (1==J)*bd->d##I##phi - bd->gamma##I##J*gammai1ldlphi))*bd->d1a + \
    (bd->G2##I##J + 2.0*( (2==I)*bd->d##J##phi + (2==J)*bd->d##I##phi - bd->gamma##I##J*gammai2ldlphi))*bd->d2a + \
    (bd->G3##I##J + 2.0*( (3==I)*bd->d##J##phi + (3==J)*bd->d##I##phi - bd->gamma##I##J*gammai3ldlphi))*bd->d3a \
  );


#define BSSN_CALCULATE_RICCI_UNITARY_TERM1(K, L, I, J) \
  bd->gammai##K##L*bd->d##K##d##L##g##I##J

#define BSSN_CALCULATE_RICCI_UNITARY_TERM2(K, I, J) \
  bd->gamma##K##I*derivative(bd->i, bd->j, bd->k, J, Gamma##K->_array_a)

#define BSSN_CALCULATE_RICCI_UNITARY_TERM3(K, I, J) \
  bd->Gammad##K*bd->GL##I##J##K

#if EXCLUDE_SECOND_ORDER_SMALL
#define BSSN_CALCULATE_RICCI_UNITARY_TERM4(K, L, M, I, J) 0
#else
#define BSSN_CALCULATE_RICCI_UNITARY_TERM4(K, L, M, I, J) \
  bd->gammai##L##M*( \
    bd->G##K##L##I*bd->GL##J##K##M + bd->G##K##L##J*bd->GL##I##K##M \
    /* symmetrize below because symmetry */ \
    + 0.5*(bd->G##K##I##M*bd->GL##K##L##J + bd->G##K##J##M*bd->GL##K##L##I) \
  )
#endif

#define BSSN_CALCULATE_RICCI_UNITARY(I, J) bd->Uricci##I##J = ( \
    - 0.5*( \
      COSMO_SUMMATION_2_ARGS(BSSN_CALCULATE_RICCI_UNITARY_TERM1, I, J) \
    ) \
    + 0.5*( \
      COSMO_SUMMATION_1_ARGS(BSSN_CALCULATE_RICCI_UNITARY_TERM2, I, J) \
      + COSMO_SUMMATION_1_ARGS(BSSN_CALCULATE_RICCI_UNITARY_TERM2, J, I) \
    ) \
    + 0.5*( \
      COSMO_SUMMATION_1_ARGS(BSSN_CALCULATE_RICCI_UNITARY_TERM3, I, J) \
      + COSMO_SUMMATION_1_ARGS(BSSN_CALCULATE_RICCI_UNITARY_TERM3, J, I) \
    ) \
    + COSMO_SUMMATION_3_ARGS(BSSN_CALCULATE_RICCI_UNITARY_TERM4, I, J) \
  );

#define BSSN_CALCULATE_DIDJGAMMA_PERMS(I, J)           \
  bd->d##I##d##J##g11 = double_derivative(bd->i, bd->j, bd->k, I, J, DIFFgamma11->_array_a); \
  bd->d##I##d##J##g12 = double_derivative(bd->i, bd->j, bd->k, I, J, DIFFgamma12->_array_a); \
  bd->d##I##d##J##g13 = double_derivative(bd->i, bd->j, bd->k, I, J, DIFFgamma13->_array_a); \
  bd->d##I##d##J##g22 = double_derivative(bd->i, bd->j, bd->k, I, J, DIFFgamma22->_array_a); \
  bd->d##I##d##J##g23 = double_derivative(bd->i, bd->j, bd->k, I, J, DIFFgamma23->_array_a); \
  bd->d##I##d##J##g33 = double_derivative(bd->i, bd->j, bd->k, I, J, DIFFgamma33->_array_a)


/*
 * Evolution equations for indexed components
 */

#if USE_BSSN_SHIFT
#define BSSN_DT_DIFFGAMMAIJ(I, J) ( \
    - 2.0*bd->alpha*bd->A##I##J \
    + bd->gamma##I##1*bd->d##J##beta1 + bd->gamma##I##2*bd->d##J##beta2 + bd->gamma##I##3*bd->d##J##beta3 \
    + bd->gamma##J##1*bd->d##I##beta1 + bd->gamma##J##2*bd->d##I##beta2 + bd->gamma##J##3*bd->d##I##beta3 \
    - (2.0/3.0)*bd->gamma##I##J*(bd->d1beta1 + bd->d2beta2 + bd->d3beta3) \
  )
#else
#define BSSN_DT_DIFFGAMMAIJ(I, J) ( \
    - 2.0*bd->alpha*bd->A##I##J \
  )
#endif

#if EXCLUDE_SECOND_ORDER_SMALL
# define BSSN_DT_AIJ_SECOND_ORDER_AA(I, J) 0
#else
# define BSSN_DT_AIJ_SECOND_ORDER_AA(I, J) ( \
      bd->gammai11*bd->A1##I*bd->A1##J + bd->gammai12*bd->A1##I*bd->A2##J + bd->gammai13*bd->A1##I*bd->A3##J \
      + bd->gammai21*bd->A2##I*bd->A1##J + bd->gammai22*bd->A2##I*bd->A2##J + bd->gammai23*bd->A2##I*bd->A3##J \
      + bd->gammai31*bd->A3##I*bd->A1##J + bd->gammai32*bd->A3##I*bd->A2##J + bd->gammai33*bd->A3##I*bd->A3##J \
    )
#endif

#if EXCLUDE_SECOND_ORDER_SMALL
# define BSSN_DT_AIJ_SECOND_ORDER_KA(I, J) ( \
    (bd->K_FRW + 2.0*bd->theta)*bd->A##I##J \
  )
#else
# define BSSN_DT_AIJ_SECOND_ORDER_KA(I, J) ( \
    (bd->K + 2.0*bd->theta)*bd->A##I##J \
  )
#endif

#if USE_BSSN_SHIFT
#define BSSN_DT_AIJ(I, J) ( \
    exp(-4.0*bd->phi)*( bd->alpha*(bd->ricciTF##I##J - 8.0*PI*bd->STF##I##J) - bd->D##I##D##J##aTF ) \
    + bd->alpha*(BSSN_DT_AIJ_SECOND_ORDER_KA(I,J) - 2.0*BSSN_DT_AIJ_SECOND_ORDER_AA(I,J)) \
    + bd->beta1*derivative(bd->i, bd->j, bd->k, 1, A##I##J->_array_a) \
    + bd->beta2*derivative(bd->i, bd->j, bd->k, 2, A##I##J->_array_a) \
    + bd->beta3*derivative(bd->i, bd->j, bd->k, 3, A##I##J->_array_a) \
    + bd->A##I##1*bd->d##J##beta1 + bd->A##I##2*bd->d##J##beta2 + bd->A##I##3*bd->d##J##beta3 \
    + bd->A##J##1*bd->d##I##beta1 + bd->A##J##2*bd->d##I##beta2 + bd->A##J##3*bd->d##I##beta3 \
    - (2.0/3.0)*bd->A##I##J*(bd->d1beta1 + bd->d2beta2 + bd->d3beta3) \
  )
#else
#define BSSN_DT_AIJ(I, J) ( \
    exp(-4.0*bd->phi)*( bd->alpha*(bd->ricciTF##I##J - 8.0*PI*bd->STF##I##J) - bd->D##I##D##J##aTF ) \
    + bd->alpha*(BSSN_DT_AIJ_SECOND_ORDER_KA(I,J) - 2.0*BSSN_DT_AIJ_SECOND_ORDER_AA(I,J)) \
  )
#endif

#define BSSN_DT_GAMMAI(I) (BSSN_DT_GAMMAI_NOSHIFT(I) + BSSN_DT_GAMMAI_SHIFT(I))

#if EXCLUDE_SECOND_ORDER_SMALL
#define BSSN_DT_GAMMAI_SECOND_ORDER(I) 0
#else
#define BSSN_DT_GAMMAI_SECOND_ORDER(I) ( \
  bd->G##I##11*bd->Acont11 + bd->G##I##22*bd->Acont22 + bd->G##I##33*bd->Acont33 \
  + 2.0*(bd->G##I##12*bd->Acont12 + bd->G##I##13*bd->Acont13 + bd->G##I##23*bd->Acont23) \
  + 6.0*(bd->Acont##I##1*bd->d1phi + bd->Acont##I##2*bd->d2phi + bd->Acont##I##3*bd->d3phi) \
)
#endif

#define BSSN_DT_GAMMAI_NOSHIFT(I) ( \
    - 2.0*(bd->Acont##I##1*bd->d1a + bd->Acont##I##2*bd->d2a + bd->Acont##I##3*bd->d3a) \
    + 2.0*bd->alpha*( \
        + BSSN_DT_GAMMAI_SECOND_ORDER(I) \
        - (1.0/3.0) * ( \
            2.0*(bd->gammai##I##1*bd->d1K + bd->gammai##I##2*bd->d2K + bd->gammai##I##3*bd->d3K) \
            + bd->gammai##I##1*bd->d1theta + bd->gammai##I##2*bd->d2theta + bd->gammai##I##3*bd->d3theta \
          ) \
        - 2.0*bd->alpha*Z4c_K1_DAMPING_AMPLITUDE*( \
            bd->Gamma##I - bd->Gammad##I \
          ) \
        - 8.0*PI*(bd->gammai##I##1*bd->S1 + bd->gammai##I##2*bd->S2 + bd->gammai##I##3*bd->S3) \
      ) \
    )

#if USE_BSSN_SHIFT
#define BSSN_DT_GAMMAI_SHIFT(I) ( \
    + bd->beta1*derivative(bd->i, bd->j, bd->k, 1, Gamma##I->_array_a) \
    + bd->beta2*derivative(bd->i, bd->j, bd->k, 2, Gamma##I->_array_a) \
    + bd->beta3*derivative(bd->i, bd->j, bd->k, 3, Gamma##I->_array_a) \
    - bd->Gamma1*bd->d1beta##I + bd->Gamma2*bd->d2beta##I + bd->Gamma3*bd->d3beta##I \
    + (2.0/3.0) * bd->Gamma##I * (bd->d1beta1 + bd->d2beta2 + bd->d3beta3) \
    + (1.0/3.0) * ( \
        bd->gammai##I##1*double_derivative(bd->i, bd->j, bd->k, 1, 1, beta1->_array_a) + bd->gammai##I##1*double_derivative(bd->i, bd->j, bd->k, 2, 1, beta2->_array_a) + bd->gammai##I##1*double_derivative(bd->i, bd->j, bd->k, 3, 1, beta3->_array_a) +  \
        bd->gammai##I##2*double_derivative(bd->i, bd->j, bd->k, 1, 2, beta1->_array_a) + bd->gammai##I##2*double_derivative(bd->i, bd->j, bd->k, 2, 2, beta2->_array_a) + bd->gammai##I##2*double_derivative(bd->i, bd->j, bd->k, 3, 2, beta3->_array_a) +  \
        bd->gammai##I##3*double_derivative(bd->i, bd->j, bd->k, 1, 3, beta1->_array_a) + bd->gammai##I##3*double_derivative(bd->i, bd->j, bd->k, 2, 3, beta2->_array_a) + bd->gammai##I##3*double_derivative(bd->i, bd->j, bd->k, 3, 3, beta3->_array_a) \
      ) \
    + ( \
        bd->gammai11*double_derivative(bd->i, bd->j, bd->k, 1, 1, beta##I->_array_a) + bd->gammai22*double_derivative(bd->i, bd->j, bd->k, 2, 2, beta##I->_array_a) + bd->gammai33*double_derivative(bd->i, bd->j, bd->k, 3, 3, beta##I->_array_a) \
        + 2.0*(bd->gammai12*double_derivative(bd->i, bd->j, bd->k, 1, 2, beta##I->_array_a) + bd->gammai13*double_derivative(bd->i, bd->j, bd->k, 1, 3, beta##I->_array_a) + bd->gammai23*double_derivative(bd->i, bd->j, bd->k, 2, 3, beta##I->_array_a)) \
      ) \
  )
#else
#define BSSN_DT_GAMMAI_SHIFT(I) 0.0
#endif

// Momentum constraint with index lowered using the conformal metric:
// M_I == \bar{gamma}_{IJ} M^J
#define BSSN_MI(I) exp(6.0*bd->phi)*( \
    - 2.0/3.0*bd->d##I##K \
    /* Note: S_I was lowered with the full metric, not conformal. */ \
    - 8*PI*(bd->S##I) \
    - 2.0/3.0*2.0*bd->d##I##theta \
    + 6.0*( \
      bd->gammai11*bd->A1##I*bd->d1phi + bd->gammai21*bd->A2##I*bd->d1phi + bd->gammai31*bd->A3##I*bd->d1phi \
      + bd->gammai12*bd->A1##I*bd->d2phi + bd->gammai22*bd->A2##I*bd->d2phi + bd->gammai32*bd->A3##I*bd->d2phi \
      + bd->gammai13*bd->A1##I*bd->d3phi + bd->gammai23*bd->A2##I*bd->d3phi + bd->gammai33*bd->A3##I*bd->d3phi \
    ) + ( \
      /* (gamma^jk D_j A_ki) */ \
      bd->gammai11*derivative(bd->i, bd->j, bd->k, 1, A1##I->_array_a) + bd->gammai12*derivative(bd->i, bd->j, bd->k, 2, A1##I->_array_a) + bd->gammai13*derivative(bd->i, bd->j, bd->k, 3, A1##I->_array_a) \
      + bd->gammai21*derivative(bd->i, bd->j, bd->k, 1, A2##I->_array_a) + bd->gammai22*derivative(bd->i, bd->j, bd->k, 2, A2##I->_array_a) + bd->gammai23*derivative(bd->i, bd->j, bd->k, 3, A2##I->_array_a) \
      + bd->gammai31*derivative(bd->i, bd->j, bd->k, 1, A3##I->_array_a) + bd->gammai32*derivative(bd->i, bd->j, bd->k, 2, A3##I->_array_a) + bd->gammai33*derivative(bd->i, bd->j, bd->k, 3, A3##I->_array_a) \
      - bd->Gamma1*bd->A1##I - bd->Gamma2*bd->A2##I - bd->Gamma3*bd->A3##I \
      - bd->GL11##I*bd->Acont11 - bd->GL21##I*bd->Acont21 - bd->GL31##I*bd->Acont31 \
      - bd->GL12##I*bd->Acont12 - bd->GL22##I*bd->Acont22 - bd->GL32##I*bd->Acont32 \
      - bd->GL13##I*bd->Acont13 - bd->GL23##I*bd->Acont23 - bd->GL33##I*bd->Acont33 \
    ) \
  )

#define BSSN_MI_SCALE(I) exp(6.0*bd->phi)*( \
    fabs(2.0/3.0*derivative(bd->i, bd->j, bd->k, I, DIFFK->_array_a)) \
    + fabs(8*PI*(bd->S##I)) \
    + 6.0*fabs( \
      bd->gammai11*bd->A1##I*bd->d1phi + bd->gammai21*bd->A2##I*bd->d1phi + bd->gammai31*bd->A3##I*bd->d1phi \
      + bd->gammai12*bd->A1##I*bd->d2phi + bd->gammai22*bd->A2##I*bd->d2phi + bd->gammai32*bd->A3##I*bd->d2phi \
      + bd->gammai13*bd->A1##I*bd->d3phi + bd->gammai23*bd->A2##I*bd->d3phi + bd->gammai33*bd->A3##I*bd->d3phi \
    ) + fabs( \
      /* (gamma^jk D_j A_ki) */ \
      bd->gammai11*derivative(bd->i, bd->j, bd->k, 1, A1##I->_array_a) + bd->gammai12*derivative(bd->i, bd->j, bd->k, 2, A1##I->_array_a) + bd->gammai13*derivative(bd->i, bd->j, bd->k, 3, A1##I->_array_a) \
      + bd->gammai21*derivative(bd->i, bd->j, bd->k, 1, A2##I->_array_a) + bd->gammai22*derivative(bd->i, bd->j, bd->k, 2, A2##I->_array_a) + bd->gammai23*derivative(bd->i, bd->j, bd->k, 3, A2##I->_array_a) \
      + bd->gammai31*derivative(bd->i, bd->j, bd->k, 1, A3##I->_array_a) + bd->gammai32*derivative(bd->i, bd->j, bd->k, 2, A3##I->_array_a) + bd->gammai33*derivative(bd->i, bd->j, bd->k, 3, A3##I->_array_a) \
      - bd->Gamma1*bd->A1##I - bd->Gamma2*bd->A2##I - bd->Gamma3*bd->A3##I \
      - bd->GL11##I*bd->Acont11 - bd->GL21##I*bd->Acont21 - bd->GL31##I*bd->Acont31 \
      - bd->GL12##I*bd->Acont12 - bd->GL22##I*bd->Acont22 - bd->GL32##I*bd->Acont32 \
      - bd->GL13##I*bd->Acont13 - bd->GL23##I*bd->Acont23 - bd->GL33##I*bd->Acont33 \
    ) \
  )

/*
 * Full metric calcs
 */

// metric derivs

#define SET_DKM00(k) bd->d##k##m00 = DKM00(k);
#define DKM00(k) \
      -2.0*bd->alpha*bd->d##k##a \
      + bd->d##k##m11*bd->beta1*bd->beta1 + bd->d##k##m22*bd->beta2*bd->beta2 + bd->d##k##m33*bd->beta3*bd->beta3 \
      + 2.0*(bd->d##k##m12*bd->beta1*bd->beta2 + bd->d##k##m13*bd->beta1*bd->beta3 + bd->d##k##m23*bd->beta2*bd->beta3) \
      + 2.0*exp(4.0*bd->phi)*( \
          bd->gamma11*bd->d##k##beta1*bd->beta1 + bd->gamma12*bd->d##k##beta2*bd->beta1 + bd->gamma13*bd->d##k##beta3*bd->beta1 \
          + bd->gamma12*bd->d##k##beta1*bd->beta2 + bd->gamma22*bd->d##k##beta2*bd->beta2 + bd->gamma23*bd->d##k##beta3*bd->beta2 \
          + bd->gamma13*bd->d##k##beta1*bd->beta3 + bd->gamma23*bd->d##k##beta2*bd->beta3 + bd->gamma33*bd->d##k##beta3*bd->beta3 \
        );


#define SET_DKM0I(k, i) bd->d##k##m0##i = DKM0I(k, i);
#define DKM0I(k, i) \
      exp(4.0*bd->phi)*( \
        bd->gamma1##i * bd->d##k##beta##1 + bd->gamma2##i * bd->d##k##beta##2 + bd->gamma3##i * bd->d##k##beta##3 \
        + (bd->d##k##g1##i + 4.0*bd->d##k##phi*bd->gamma1##i) * bd->beta##1   \
        + (bd->d##k##g2##i + 4.0*bd->d##k##phi*bd->gamma2##i) * bd->beta##2 \
        + (bd->d##k##g3##i + 4.0*bd->d##k##phi*bd->gamma3##i) * bd->beta##3 \
      );

#define SET_DKMIJ(k, i, j) bd->d##k##m##i##j = DKMIJ(k, i, j);
#define DKMIJ(k, i, j) exp(4.0*bd->phi)*(bd->d##k##g##i##j + 4.0*bd->d##k##phi*bd->gamma##i##j);


// metric

#define SET_M00() bd->m00 = M00();
#define M00() \
      -bd->alpha*bd->alpha + exp(4.0*bd->phi)*( \
        bd->gamma11*bd->beta1*bd->beta1 + bd->gamma22*bd->beta2*bd->beta2 + bd->gamma33*bd->beta3*bd->beta3 \
        + 2.0*(bd->gamma12*bd->beta1*bd->beta2 + bd->gamma13*bd->beta1*bd->beta3 + bd->gamma23*bd->beta2*bd->beta3) \
      );

#define SET_M0I(i) bd->m0##i = M0I(i);
#define M0I(i) exp(4.0*bd->phi)*(bd->gamma1##i*bd->beta1 + bd->gamma2##i*bd->beta2 + bd->gamma3##i*bd->beta3);

#define SET_MIJ(i, j) bd->m##i##j = MIJ(i, j);
#define MIJ(i, j) exp(4.0*bd->phi)*(bd->gamma##i##j);

// inverse metric

#define SET_Mi00() bd->mi00 = Mi00();
#define Mi00() -1.0/bd->alpha/bd->alpha;

#define SET_Mi0I(i) bd->mi0##i = Mi0I(i);
#define Mi0I(i) 1.0/bd->alpha/bd->alpha*bd->beta##i;

#define SET_MiIJ(i, j) bd->mi##i##j = MiIJ(i, j);
#define MiIJ(i, j) exp(-4.0*bd->phi)*(bd->gamma##i##j) - 1.0/bd->alpha/bd->alpha*bd->beta##i*bd->beta##j;


/*
 * Conversion to ADM quantities for raytracing
 */
#define BSSN_RP_DG(I,J,L) \
  P*(4.0*bd->gamma##I##J*bd->d##L##phi+ bd->d##L##g##I##J)

#define BSSN_RP_DDG(I,J,L,M) \
  P*( \
    16.0*bd->gamma##I##J*bd->d##M##phi*bd->d##L##phi \
    + 4.0*(bd->d##L##g##I##J*bd->d##M##phi + bd->d##M##g##I##J*bd->d##L##phi + bd->gamma##I##J*bd->d##L##d##M##phi) \
    + bd->d##L##d##M##g##I##J \
  )

#define BSSN_RP_K(I,J) \
  P*(bd->A##I##J + 1.0/3.0*bd->gamma##I##J*bd->K)

#define BSSN_RP_DK(I,J,L) \
  P*( \
    4.0*(bd->A##I##J + 1.0/3.0*bd->gamma##I##J*bd->K)*bd->d##L##phi + derivative(bd->i, bd->j, bd->k, L, A##I##J->_array_a) \
    + 1.0/3.0*bd->K*bd->d##L##g##I##J + 1.0/3.0*bd->gamma##I##J*bd->d##L##K \
  )

#define BSSN_RP_GAMMA(I,J,K) \
  bd->G##I##J##K + 2.0*( \
      (I==J ? 1.0 : 0.0)*bd->d##K##phi \
      + (I==K ? 1.0 : 0.0)*bd->d##J##phi \
      - bd->gamma##J##K*(bd->gammai##I##1*bd->d1phi + bd->gammai##I##2*bd->d2phi + bd->gammai##I##3*bd->d3phi) \
    )

#define BSSN_RP_GAMMAL(I,J,K) P*( \
  bd->GL##I##J##K + 2.0*( \
      bd->gamma##I##J*bd->d##K##phi \
      + bd->gamma##I##K*bd->d##J##phi \
      - bd->gamma##J##K*bd->d##I##phi \
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

#define DIFFgamma21 DIFFgamma12
#define DIFFgamma31 DIFFgamma13
#define DIFFgamma32 DIFFgamma23
#define A21 A12
#define A31 A13
#define A32 A23

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
