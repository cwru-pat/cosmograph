#ifndef COSMO_DEFINES
#define COSMO_DEFINES

#define N 256
#define dt 0.1

#define RESTRICT __restrict__

#define INDEX(i,j,k) (((i+N)%N)*N*N + ((j+N)%N)*N + (k+N)%N)

#define LOOP3(i,j,k) \
  for(idx_t i=0; i<N; ++i) \
    for(idx_t j=0; j<N; ++j) \
      for(idx_t k=0; k<N; ++k)

#define INTERNAL_LOOP3(i,j,k) \
  for(idx_t i=1; i<N-1; ++i) \
    for(idx_t j=1; j<N-1; ++j) \
      for(idx_t k=1; k<N-1; ++k)


#define RK4_ARRAY_CREATE(name) \
        real_t * name, * name##_p, * name_##i;

#define RK4_ARRAY_ALLOC(name) \
        name       = new real_t[N*N*N]; \
        name##_p   = new real_t[N*N*N]; \
        name##_i   = new real_t[N*N*N];

#define RK4_ARRAY_DELETE(name) \
        delete [] name;        \
        delete [] name##_p;    \
        delete [] name##_i;

// RK4 has a diagonal tableau, so we only need to compute
// the coefficients one at a time ("_i" arrays), given the
// values in the previous step ("_p" arrays).
#define RK4_ARRAY_ADDMAP(name)          \
        fields[#name]   = name;         \
        fields[#name "_p"] = name##_p;  \
        fields[#name "_i"] = name##_i;

#define RK4_ARRAY_CYCLE(name) \
        std::swap(name, name_p);

#define BSSN_APPLY_TO_FIELDS(function)  \
        function(gamma11);              \
        function(gamma12);              \
        function(gamma13);              \
        function(gamma22);              \
        function(gamma23);              \
        function(gamma33);              \
        function(gammai11);             \
        function(gammai12);             \
        function(gammai13);             \
        function(gammai22);             \
        function(gammai23);             \
        function(gammai33);             \
        function(phi);                  \
        function(A11);                  \
        function(A12);                  \
        function(A13);                  \
        function(A22);                  \
        function(A23);                  \
        function(A33);                  \
        function(K);                    \
        function(Gamma1);               \
        function(Gamma2);               \
        function(Gamma3);               \
        function(beta1);                \
        function(beta2);                \
        function(beta3);                \
        function(alpha);

#define BSSN_APPLY_TO_IJ_PERMS(function, more_args) \
        function(1, 1, more_args);                  \
        function(1, 2, more_args);                  \
        function(1, 3, more_args);                  \
        function(2, 2, more_args);                  \
        function(2, 3, more_args);                  \
        function(3, 3, more_args);

// standard ordering of indexes for tensor components
#define gamma21 gamma12
#define gamma31 gamma13
#define gamma32 gamma23
#define gammai21 gammai12
#define gammai31 gammai13
#define gammai32 gammai23
#define A21 A12
#define A31 A13
#define A32 A23
        
#define ricciTF21 ricciTF12
#define ricciTF31 ricciTF13
#define ricciTF32 ricciTF23

#define D2D1a D1D2a
#define D3D1a D1D3a
#define D3D2a D2D3a

#define G121 G112
#define G131 G113
#define G132 G123
#define G221 G212
#define G231 G213
#define G232 G223
#define G321 G312
#define G331 G313
#define G332 G323

#define d1g21 d1g12
#define d1g31 d1g13
#define d1g32 d1g23
#define d2g21 d2g12
#define d2g31 d2g13
#define d2g32 d2g23
#define d3g21 d3g12
#define d3g31 d3g13
#define d3g32 d3g23



#define BSSN_GAMMA_EQUATION(I, J, idx) ( \
    /* -2 \alpha A_{ij} */ \
  - 2*alpha[idx]*A ## I ## J ## [idx] \
    /* \beta^k * d_k \gamma_{ij} */ \
  + beta1[idx]*der(gamma ## I ## J, 1, idx) + beta2[idx]*der(gamma ## I ## J, 2, idx) + beta3[idx]*der(gamma ## I ## J, 3, idx) \
    /* \gamma_{ik} d_j beta^k */ \
  + gamma ## I ## 1[idx] *der(beta1, J, idx) + gamma ## I ## 2[idx] *der(beta2, J, idx) + gamma ## I ## 3[idx] *der(beta3, J, idx) \
    /* \gamma_{jk} d_i beta */ \
  + gamma ## J ## 1[idx] *der(beta1, I, idx) + gamma ## J ## 2[idx] *der(beta2, I, idx) + gamma ## J ## 3[idx] *der(beta3, I, idx) \
    /* -2/3 \gamma_{ij} d_k \beta^k */ \
  - 2.0/3.0 * gamma ## I ## J ## [idx] * ( der(beta1, 1, idx) + der(beta2, 2, idx) + der(beta3, 3, idx) ) \
)

#define BSSN_A_EQUATION(I, J, idx) ( \
  exp(-4.0*phi[idx])*( \
    - D##I##D##J##alpha \
    + alpha[idx]*(ricciTF ## I ## J ## [idx] + /* tracefree source */) \
  ) \
  + alpha[idx] * ( \
      /* K A_{ij} */ \
    K * A ## I ## J ## [idx] \
      /* A_{il} A^l_j = gamma^{lm} A_{il} A_{mj} */ \
    - 2.0 * ( \
      gamma11[idx]*A##I##1[idx]*A##J##1[idx] + gamma12[idx]*A##I##1[idx]*A##J##2[idx] + gamma13[idx]*A##I##1[idx]*A##J##3[idx] \
      + gamma21[idx]*A##I##2[idx]*A##J##1[idx] + gamma22[idx]*A##I##2[idx]*A##J##2[idx] + gamma23[idx]*A##I##2[idx]*A##J##3[idx] \
      + gamma31[idx]*A##I##3[idx]*A##J##1[idx] + gamma32[idx]*A##I##3[idx]*A##J##2[idx] + gamma33[idx]*A##I##3[idx]*A##J##3[idx] \
    ) \
  ) \
    /* \beta^k * d_k \A_{ij} */ \
  + beta1[idx]*der(A ## I ## J, 1, idx) + beta2[idx]*der(A ## I ## J, 2, idx) + beta3[idx]*der(A ## I ## J, 3, idx) \
    /* \A_{ik} d_j beta^k */ \
  + A ## I ## 1[idx] *der(beta1, J, idx) + A ## I ## 2[idx] *der(beta2, J, idx) + A ## I ## 3[idx] *der(beta3, J, idx) \
    /* \A_{jk} d_i beta */ \
  + A ## J ## 1[idx] *der(beta1, I, idx) + A ## J ## 2[idx] *der(beta2, I, idx) + A ## J ## 3[idx] *der(beta3, I, idx) \
    /* -2/3 \A_{ij} d_k \beta^k */ \
  - 2.0/3.0 * A ## I ## J ## [idx] * ( der(beta1, 1, idx) + der(beta2, 2, idx) + der(beta3, 3, idx) ) \
)


#endif
