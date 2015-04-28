#ifndef BSSN_MACROS
#define BSSN_MACROS

/*
 * applying functions to lots of vars
 */

#define BSSN_APPLY_TO_FIELDS(function)  \
  function(gamma11);                    \
  function(gamma12);                    \
  function(gamma13);                    \
  function(gamma22);                    \
  function(gamma23);                    \
  function(gamma33);                    \
  function(phi);                        \
  function(A11);                        \
  function(A12);                        \
  function(A13);                        \
  function(A22);                        \
  function(A23);                        \
  function(A33);                        \
  function(K);                          \
  function(Gamma1);                     \
  function(Gamma2);                     \
  function(Gamma3);

#define BSSN_LIST_FIELDS(function)  \
  function(gamma11),                \
  function(gamma12),                \
  function(gamma13),                \
  function(gamma22),                \
  function(gamma23),                \
  function(gamma33),                \
  function(phi),                    \
  function(A11),                    \
  function(A12),                    \
  function(A13),                    \
  function(A22),                    \
  function(A23),                    \
  function(A33),                    \
  function(K),                      \
  function(Gamma1),                 \
  function(Gamma2),                 \
  function(Gamma3)

#define BSSN_APPLY_TO_SOURCES(function) \
  function(r);                          \
  function(S);                          \
  function(S1);                         \
  function(S2);                         \
  function(S3);                         \
  function(STF11);                      \
  function(STF12);                      \
  function(STF13);                      \
  function(STF22);                      \
  function(STF23);                      \
  function(STF33);

#define BSSN_LIST_SOURCES(function) \
  function(r),                      \
  function(S),                      \
  function(S1),                     \
  function(S2),                     \
  function(S3),                     \
  function(STF11),                  \
  function(STF12),                  \
  function(STF13),                  \
  function(STF22),                  \
  function(STF23),                  \
  function(STF33)

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

#define BSSN_SWAP_ARRAYS(reg_prefix_1, reg_prefix_2)                               \
  std::swap(gamma11##reg_prefix_1,          gamma11##reg_prefix_2);                \
  std::swap(fields["gamma11"#reg_prefix_1], fields["gamma11"#reg_prefix_2]);       \
  std::swap(gamma12##reg_prefix_1,          gamma12##reg_prefix_2);                \
  std::swap(fields["gamma12"#reg_prefix_1], fields["gamma12"#reg_prefix_2]);       \
  std::swap(gamma13##reg_prefix_1,          gamma13##reg_prefix_2);                \
  std::swap(fields["gamma13"#reg_prefix_1], fields["gamma13"#reg_prefix_2]);       \
  std::swap(gamma22##reg_prefix_1,          gamma22##reg_prefix_2);                \
  std::swap(fields["gamma22"#reg_prefix_1], fields["gamma22"#reg_prefix_2]);       \
  std::swap(gamma23##reg_prefix_1,          gamma23##reg_prefix_2);                \
  std::swap(fields["gamma23"#reg_prefix_1], fields["gamma23"#reg_prefix_2]);       \
  std::swap(gamma33##reg_prefix_1,          gamma33##reg_prefix_2);                \
  std::swap(fields["gamma33"#reg_prefix_1], fields["gamma33"#reg_prefix_2]);       \
  std::swap(phi##reg_prefix_1,              phi##reg_prefix_2);                    \
  std::swap(fields["phi"#reg_prefix_1],     fields["phi"#reg_prefix_2]);           \
  std::swap(A11##reg_prefix_1,              A11##reg_prefix_2);                    \
  std::swap(fields["A11"#reg_prefix_1],     fields["A11"#reg_prefix_2]);           \
  std::swap(A12##reg_prefix_1,              A12##reg_prefix_2);                    \
  std::swap(fields["A12"#reg_prefix_1],     fields["A12"#reg_prefix_2]);           \
  std::swap(A13##reg_prefix_1,              A13##reg_prefix_2);                    \
  std::swap(fields["A13"#reg_prefix_1],     fields["A13"#reg_prefix_2]);           \
  std::swap(A22##reg_prefix_1,              A22##reg_prefix_2);                    \
  std::swap(fields["A22"#reg_prefix_1],     fields["A22"#reg_prefix_2]);           \
  std::swap(A23##reg_prefix_1,              A23##reg_prefix_2);                    \
  std::swap(fields["A23"#reg_prefix_1],     fields["A23"#reg_prefix_2]);           \
  std::swap(A33##reg_prefix_1,              A33##reg_prefix_2);                    \
  std::swap(fields["A33"#reg_prefix_1],     fields["A33"#reg_prefix_2]);           \
  std::swap(K##reg_prefix_1,                K##reg_prefix_2);                      \
  std::swap(fields["K"#reg_prefix_1],       fields["K"#reg_prefix_2]);             \
  std::swap(Gamma1##reg_prefix_1,           Gamma1##reg_prefix_2);                 \
  std::swap(fields["Gamma1"#reg_prefix_1],  fields["Gamma1"#reg_prefix_2]);        \
  std::swap(Gamma2##reg_prefix_1,           Gamma2##reg_prefix_2);                 \
  std::swap(fields["Gamma2"#reg_prefix_1],  fields["Gamma2"#reg_prefix_2]);        \
  std::swap(Gamma3##reg_prefix_1,           Gamma3##reg_prefix_2);                 \
  std::swap(fields["Gamma3"#reg_prefix_1],  fields["Gamma3"#reg_prefix_2]);

#define BSSN_COPY_ARRAYS(reg_prefix_from, reg_prefix_to)        \
  std::copy(gamma11##reg_prefix_from,  gamma11##reg_prefix_from + POINTS,  gamma11##reg_prefix_to  ); \
  std::copy(gamma12##reg_prefix_from,  gamma12##reg_prefix_from + POINTS,  gamma12##reg_prefix_to  ); \
  std::copy(gamma13##reg_prefix_from,  gamma13##reg_prefix_from + POINTS,  gamma13##reg_prefix_to  ); \
  std::copy(gamma22##reg_prefix_from,  gamma22##reg_prefix_from + POINTS,  gamma22##reg_prefix_to  ); \
  std::copy(gamma23##reg_prefix_from,  gamma23##reg_prefix_from + POINTS,  gamma23##reg_prefix_to  ); \
  std::copy(gamma33##reg_prefix_from,  gamma33##reg_prefix_from + POINTS,  gamma33##reg_prefix_to  ); \
  std::copy(phi##reg_prefix_from,      phi##reg_prefix_from + POINTS,      phi##reg_prefix_to      ); \
  std::copy(A11##reg_prefix_from,      A11##reg_prefix_from + POINTS,      A11##reg_prefix_to      ); \
  std::copy(A12##reg_prefix_from,      A12##reg_prefix_from + POINTS,      A12##reg_prefix_to      ); \
  std::copy(A13##reg_prefix_from,      A13##reg_prefix_from + POINTS,      A13##reg_prefix_to      ); \
  std::copy(A22##reg_prefix_from,      A22##reg_prefix_from + POINTS,      A22##reg_prefix_to      ); \
  std::copy(A23##reg_prefix_from,      A23##reg_prefix_from + POINTS,      A23##reg_prefix_to      ); \
  std::copy(A33##reg_prefix_from,      A33##reg_prefix_from + POINTS,      A33##reg_prefix_to      ); \
  std::copy(K##reg_prefix_from,        K##reg_prefix_from + POINTS,        K##reg_prefix_to        ); \
  std::copy(Gamma1##reg_prefix_from,   Gamma1##reg_prefix_from + POINTS,   Gamma1##reg_prefix_to   ); \
  std::copy(Gamma2##reg_prefix_from,   Gamma2##reg_prefix_from + POINTS,   Gamma2##reg_prefix_to   ); \
  std::copy(Gamma3##reg_prefix_from,   Gamma3##reg_prefix_from + POINTS,   Gamma3##reg_prefix_to   );

#define BSSN_ZERO_ARRAY(reg, index) \
  gamma11##reg[idx] = 0;            \
  gamma12##reg[idx] = 0;            \
  gamma13##reg[idx] = 0;            \
  gamma22##reg[idx] = 0;            \
  gamma23##reg[idx] = 0;            \
  gamma33##reg[idx] = 0;            \
  phi##reg[idx]     = 0;            \
  A11##reg[idx]     = 0;            \
  A12##reg[idx]     = 0;            \
  A13##reg[idx]     = 0;            \
  A22##reg[idx]     = 0;            \
  A23##reg[idx]     = 0;            \
  A33##reg[idx]     = 0;            \
  K##reg[idx]       = 0;            \
  Gamma1##reg[idx]  = 0;            \
  Gamma2##reg[idx]  = 0;            \
  Gamma3##reg[idx]  = 0;

// arr_c[idx] = arr_p[idx] + dt*mult*evfn(arr_a);
#define BSSN_COMPUTE_RK_STEP(mult) \
  gamma11_c[paq->idx] = gamma11_p[paq->idx] + dt*mult*ev_gamma11(paq); \
  gamma12_c[paq->idx] = gamma12_p[paq->idx] + dt*mult*ev_gamma12(paq); \
  gamma13_c[paq->idx] = gamma13_p[paq->idx] + dt*mult*ev_gamma13(paq); \
  gamma22_c[paq->idx] = gamma22_p[paq->idx] + dt*mult*ev_gamma22(paq); \
  gamma23_c[paq->idx] = gamma23_p[paq->idx] + dt*mult*ev_gamma23(paq); \
  gamma33_c[paq->idx] = gamma33_p[paq->idx] + dt*mult*ev_gamma33(paq); \
  phi_c[paq->idx]     = phi_p[paq->idx]     + dt*mult*ev_phi(paq);     \
  A11_c[paq->idx]     = A11_p[paq->idx]     + dt*mult*ev_A11(paq);     \
  A12_c[paq->idx]     = A12_p[paq->idx]     + dt*mult*ev_A12(paq);     \
  A13_c[paq->idx]     = A13_p[paq->idx]     + dt*mult*ev_A13(paq);     \
  A22_c[paq->idx]     = A22_p[paq->idx]     + dt*mult*ev_A22(paq);     \
  A23_c[paq->idx]     = A23_p[paq->idx]     + dt*mult*ev_A23(paq);     \
  A33_c[paq->idx]     = A33_p[paq->idx]     + dt*mult*ev_A33(paq);     \
  K_c[paq->idx]       = K_p[paq->idx]       + dt*mult*ev_K(paq);       \
  Gamma1_c[paq->idx]  = Gamma1_p[paq->idx]  + dt*mult*ev_Gamma1(paq);  \
  Gamma2_c[paq->idx]  = Gamma2_p[paq->idx]  + dt*mult*ev_Gamma2(paq);  \
  Gamma3_c[paq->idx]  = Gamma3_p[paq->idx]  + dt*mult*ev_Gamma3(paq);

#define BSSN_ADD_C_TO_F(mult) \
  gamma11_f[paq->idx] += mult*gamma11_c[paq->idx];  \
  gamma12_f[paq->idx] += mult*gamma12_c[paq->idx];  \
  gamma13_f[paq->idx] += mult*gamma13_c[paq->idx];  \
  gamma22_f[paq->idx] += mult*gamma22_c[paq->idx];  \
  gamma23_f[paq->idx] += mult*gamma23_c[paq->idx];  \
  gamma33_f[paq->idx] += mult*gamma33_c[paq->idx];  \
  phi_f[paq->idx]     += mult*phi_c[paq->idx];      \
  A11_f[paq->idx]     += mult*A11_c[paq->idx];      \
  A12_f[paq->idx]     += mult*A12_c[paq->idx];      \
  A13_f[paq->idx]     += mult*A13_c[paq->idx];      \
  A22_f[paq->idx]     += mult*A22_c[paq->idx];      \
  A23_f[paq->idx]     += mult*A23_c[paq->idx];      \
  A33_f[paq->idx]     += mult*A33_c[paq->idx];      \
  K_f[paq->idx]       += mult*K_c[paq->idx];        \
  Gamma1_f[paq->idx]  += mult*Gamma1_c[paq->idx];   \
  Gamma2_f[paq->idx]  += mult*Gamma2_c[paq->idx];   \
  Gamma3_f[paq->idx]  += mult*Gamma3_c[paq->idx];

// arr_f = (1.0/3.0)*(arr_f - arr_p) + (1.0/6.0)*evfn(arr_a)
#define BSSN_FINAL_RK4_STEP() \
  gamma11_f[paq->idx] = (1.0/3.0)*(gamma11_f[paq->idx] - gamma11_p[paq->idx]) + (1.0/6.0)*dt*ev_gamma11(paq); \
  gamma12_f[paq->idx] = (1.0/3.0)*(gamma12_f[paq->idx] - gamma12_p[paq->idx]) + (1.0/6.0)*dt*ev_gamma12(paq); \
  gamma13_f[paq->idx] = (1.0/3.0)*(gamma13_f[paq->idx] - gamma13_p[paq->idx]) + (1.0/6.0)*dt*ev_gamma13(paq); \
  gamma22_f[paq->idx] = (1.0/3.0)*(gamma22_f[paq->idx] - gamma22_p[paq->idx]) + (1.0/6.0)*dt*ev_gamma22(paq); \
  gamma23_f[paq->idx] = (1.0/3.0)*(gamma23_f[paq->idx] - gamma23_p[paq->idx]) + (1.0/6.0)*dt*ev_gamma23(paq); \
  gamma33_f[paq->idx] = (1.0/3.0)*(gamma33_f[paq->idx] - gamma33_p[paq->idx]) + (1.0/6.0)*dt*ev_gamma33(paq); \
  phi_f[paq->idx]     = (1.0/3.0)*(phi_f[paq->idx]     - phi_p[paq->idx])     + (1.0/6.0)*dt*ev_phi(paq);     \
  A11_f[paq->idx]     = (1.0/3.0)*(A11_f[paq->idx]     - A11_p[paq->idx])     + (1.0/6.0)*dt*ev_A11(paq);     \
  A12_f[paq->idx]     = (1.0/3.0)*(A12_f[paq->idx]     - A12_p[paq->idx])     + (1.0/6.0)*dt*ev_A12(paq);     \
  A13_f[paq->idx]     = (1.0/3.0)*(A13_f[paq->idx]     - A13_p[paq->idx])     + (1.0/6.0)*dt*ev_A13(paq);     \
  A22_f[paq->idx]     = (1.0/3.0)*(A22_f[paq->idx]     - A22_p[paq->idx])     + (1.0/6.0)*dt*ev_A22(paq);     \
  A23_f[paq->idx]     = (1.0/3.0)*(A23_f[paq->idx]     - A23_p[paq->idx])     + (1.0/6.0)*dt*ev_A23(paq);     \
  A33_f[paq->idx]     = (1.0/3.0)*(A33_f[paq->idx]     - A33_p[paq->idx])     + (1.0/6.0)*dt*ev_A33(paq);     \
  K_f[paq->idx]       = (1.0/3.0)*(K_f[paq->idx]       - K_p[paq->idx])       + (1.0/6.0)*dt*ev_K(paq);       \
  Gamma1_f[paq->idx]  = (1.0/3.0)*(Gamma1_f[paq->idx]  - Gamma1_p[paq->idx])  + (1.0/6.0)*dt*ev_Gamma1(paq);  \
  Gamma2_f[paq->idx]  = (1.0/3.0)*(Gamma2_f[paq->idx]  - Gamma2_p[paq->idx])  + (1.0/6.0)*dt*ev_Gamma2(paq);  \
  Gamma3_f[paq->idx]  = (1.0/3.0)*(Gamma3_f[paq->idx]  - Gamma3_p[paq->idx])  + (1.0/6.0)*dt*ev_Gamma3(paq);

/*
 * Aux. variable calculations
 */

#define BSSN_COMPUTE_LOCAL_GAMMAI_PF(IJ, f1, f2, f3, f4) \
  paq->gammai##IJ = paq->gamma##f1*paq->gamma##f2 - paq->gamma##f3*paq->gamma##f4; \
  paq->gammai##IJ##_adj[0][1][1] = paq->gamma##f1##_adj[0][1][1]*paq->gamma##f2##_adj[0][1][1] - paq->gamma##f3##_adj[0][1][1]*paq->gamma##f4##_adj[0][1][1]; \
  paq->gammai##IJ##_adj[1][0][1] = paq->gamma##f1##_adj[1][0][1]*paq->gamma##f2##_adj[1][0][1] - paq->gamma##f3##_adj[1][0][1]*paq->gamma##f4##_adj[1][0][1]; \
  paq->gammai##IJ##_adj[1][1][0] = paq->gamma##f1##_adj[1][1][0]*paq->gamma##f2##_adj[1][1][0] - paq->gamma##f3##_adj[1][1][0]*paq->gamma##f4##_adj[1][1][0]; \
  paq->gammai##IJ##_adj[1][1][1] = paq->gamma##f1##_adj[1][1][1]*paq->gamma##f2##_adj[1][1][1] - paq->gamma##f3##_adj[1][1][1]*paq->gamma##f4##_adj[1][1][1]; \
  paq->gammai##IJ##_adj[1][1][2] = paq->gamma##f1##_adj[1][1][2]*paq->gamma##f2##_adj[1][1][2] - paq->gamma##f3##_adj[1][1][2]*paq->gamma##f4##_adj[1][1][2]; \
  paq->gammai##IJ##_adj[1][2][1] = paq->gamma##f1##_adj[1][2][1]*paq->gamma##f2##_adj[1][2][1] - paq->gamma##f3##_adj[1][2][1]*paq->gamma##f4##_adj[1][2][1]; \
  paq->gammai##IJ##_adj[2][1][1] = paq->gamma##f1##_adj[2][1][1]*paq->gamma##f2##_adj[2][1][1] - paq->gamma##f3##_adj[2][1][1]*paq->gamma##f4##_adj[2][1][1];

// extended points around faces
#define BSSN_COMPUTE_LOCAL_GAMMAI_F2(IJ, f1, f2, f3, f4) \
  paq->gammai##IJ##_adj_ext[0][0] = paq->gamma##f1##_adj_ext[0][0]*paq->gamma##f2##_adj_ext[0][0] - paq->gamma##f3##_adj_ext[0][0]*paq->gamma##f4##_adj_ext[0][0]; \
  paq->gammai##IJ##_adj_ext[0][1] = paq->gamma##f1##_adj_ext[0][1]*paq->gamma##f2##_adj_ext[0][1] - paq->gamma##f3##_adj_ext[0][1]*paq->gamma##f4##_adj_ext[0][1]; \
  paq->gammai##IJ##_adj_ext[1][0] = paq->gamma##f1##_adj_ext[1][0]*paq->gamma##f2##_adj_ext[1][0] - paq->gamma##f3##_adj_ext[1][0]*paq->gamma##f4##_adj_ext[1][0]; \
  paq->gammai##IJ##_adj_ext[1][1] = paq->gamma##f1##_adj_ext[1][1]*paq->gamma##f2##_adj_ext[1][1] - paq->gamma##f3##_adj_ext[1][1]*paq->gamma##f4##_adj_ext[1][1]; \
  paq->gammai##IJ##_adj_ext[2][0] = paq->gamma##f1##_adj_ext[2][0]*paq->gamma##f2##_adj_ext[2][0] - paq->gamma##f3##_adj_ext[2][0]*paq->gamma##f4##_adj_ext[2][0]; \
  paq->gammai##IJ##_adj_ext[2][1] = paq->gamma##f1##_adj_ext[2][1]*paq->gamma##f2##_adj_ext[2][1] - paq->gamma##f3##_adj_ext[2][1]*paq->gamma##f4##_adj_ext[2][1];

#define BSSN_CALCULATE_CHRISTOFFEL(I, J, K) paq->G##I##J##K = 0.5*( \
    paq->gammai##I##1 * (paq->d##J##g##K##1 + paq->d##K##g##J##1 - paq->d1g##J##K) + \
    paq->gammai##I##2 * (paq->d##J##g##K##2 + paq->d##K##g##J##2 - paq->d2g##J##K) + \
    paq->gammai##I##3 * (paq->d##J##g##K##3 + paq->d##K##g##J##3 - paq->d3g##J##K) \
  )

#define BSSN_CALCULATE_CHRISTOFFEL_LOWER(I, J, K) paq->GL##I##J##K = 0.5*( \
    paq->d##J##g##K##I + paq->d##K##g##J##I - paq->d##I##g##J##K \
  )

#define BSSN_CALCULATE_DGAMMAI(I, J, K) paq->d##I##gi##J##K = der_ext(paq->gammai##J##K##_adj, paq->gammai##J##K##_adj_ext, I);

#define BSSN_CALCULATE_DGAMMA(I, J, K) paq->d##I##g##J##K = der_ext(paq->gamma##J##K##_adj, paq->gamma##J##K##_adj_ext, I);

#define BSSN_CALCULATE_ACONT(I, J) paq->Acont##I##J = ( \
    paq->gammai##I##1*paq->gammai##J##1*paq->A11 + paq->gammai##I##2*paq->gammai##J##1*paq->A21 + paq->gammai##I##3*paq->gammai##J##1*paq->A31 \
    + paq->gammai##I##1*paq->gammai##J##2*paq->A12 + paq->gammai##I##2*paq->gammai##J##2*paq->A22 + paq->gammai##I##3*paq->gammai##J##2*paq->A32 \
    + paq->gammai##I##1*paq->gammai##J##3*paq->A13 + paq->gammai##I##2*paq->gammai##J##3*paq->A23 + paq->gammai##I##3*paq->gammai##J##3*paq->A33 \
  );

// unitary piece only:
// calculated using http://relativity.livingreviews.org/Articles/lrr-2011-6/fulltext.html
#define BSSN_CALCULATE_RICCITF_UNITARY(I, J) paq->ricciTF##I##J = ( \
    - 0.5*( \
      paq->gammai11*paq->d1d1g##I##J + paq->gammai22*paq->d2d2g##I##J + paq->gammai33*paq->d3d3g##I##J \
      + 2.0*(paq->gammai12*paq->d1d2g##I##J + paq->gammai13*paq->d1d3g##I##J + paq->gammai23*paq->d2d3g##I##J) \
    ) \
    + 0.5*( \
      paq->gamma1##I*der_ext(paq->Gamma1_adj, paq->Gamma1_adj_ext, J) + paq->gamma2##I*der_ext(paq->Gamma2_adj, paq->Gamma2_adj_ext, J) + paq->gamma3##I*der_ext(paq->Gamma3_adj, paq->Gamma3_adj_ext, J) + \
      paq->gamma1##J*der_ext(paq->Gamma1_adj, paq->Gamma1_adj_ext, I) + paq->gamma2##J*der_ext(paq->Gamma2_adj, paq->Gamma2_adj_ext, I) + paq->gamma3##J*der_ext(paq->Gamma3_adj, paq->Gamma3_adj_ext, I) \
    ) \
    - 0.5*( \
      paq->d1g##I##1*paq->d1gi##J##1 + paq->d1g##I##2*paq->d2gi##J##1 + paq->d1g##I##3*paq->d3gi##J##1 \
        + paq->d2g##I##1*paq->d1gi##J##2 + paq->d2g##I##2*paq->d2gi##J##2 + paq->d2g##I##3*paq->d3gi##J##2 \
        + paq->d3g##I##1*paq->d1gi##J##3 + paq->d3g##I##2*paq->d2gi##J##3 + paq->d3g##I##3*paq->d3gi##J##3 \
      + paq->d1g##J##1*paq->d1gi##I##1 + paq->d1g##J##2*paq->d2gi##I##1 + paq->d1g##J##3*paq->d3gi##I##1 \
        + paq->d2g##J##1*paq->d1gi##I##2 + paq->d2g##J##2*paq->d2gi##I##2 + paq->d2g##J##3*paq->d3gi##I##2 \
        + paq->d3g##J##1*paq->d1gi##I##3 + paq->d3g##J##2*paq->d2gi##I##3 + paq->d3g##J##3*paq->d3gi##I##3 \
      - paq->Gamma1*paq->d1g##I##J - paq->Gamma2*paq->d2g##I##J - paq->Gamma3*paq->d3g##I##J \
    ) \
    - ( \
        paq->G1##I##1*paq->G1##J##1 + paq->G1##I##2*paq->G2##J##1 + paq->G1##I##3*paq->G3##J##1 \
      + paq->G2##I##1*paq->G1##J##2 + paq->G2##I##2*paq->G2##J##2 + paq->G2##I##3*paq->G3##J##2 \
      + paq->G3##I##1*paq->G1##J##3 + paq->G3##I##2*paq->G2##J##3 + paq->G3##I##3*paq->G3##J##3 \
    ) \
  );

#define BSSN_CALCULATE_RICCITF_UNITARY_ALT(I, J) paq->ricciTF##I##J = ( \
    - 0.5*( \
      paq->gammai11*paq->d1d1g##I##J + paq->gammai22*paq->d2d2g##I##J + paq->gammai33*paq->d3d3g##I##J \
      + 2.0*(paq->gammai12*paq->d1d2g##I##J + paq->gammai13*paq->d1d3g##I##J + paq->gammai23*paq->d2d3g##I##J) \
    ) \
    + 0.5*( \
      paq->gamma1##I*der_ext(paq->Gamma1_adj, paq->Gamma1_adj_ext, J) + paq->gamma2##I*der_ext(paq->Gamma2_adj, paq->Gamma2_adj_ext, J) + paq->gamma3##I*der_ext(paq->Gamma3_adj, paq->Gamma3_adj_ext, J) + \
      paq->gamma1##J*der_ext(paq->Gamma1_adj, paq->Gamma1_adj_ext, I) + paq->gamma2##J*der_ext(paq->Gamma2_adj, paq->Gamma2_adj_ext, I) + paq->gamma3##J*der_ext(paq->Gamma3_adj, paq->Gamma3_adj_ext, I) \
    ) \
    + 0.5*( \
      paq->Gamma1*paq->GL##I##J##1 + paq->Gamma2*paq->GL##I##J##2 + paq->Gamma3*paq->GL##I##J##3 \
      + paq->Gamma1*paq->GL##J##I##1 + paq->Gamma2*paq->GL##J##I##2 + paq->Gamma3*paq->GL##J##I##3 \
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
  paq->d##I##d##J##g11 = dder_ext(paq->gamma11_adj, paq->gamma11_adj_ext, I, J); \
  paq->d##I##d##J##g12 = dder_ext(paq->gamma12_adj, paq->gamma12_adj_ext, I, J); \
  paq->d##I##d##J##g13 = dder_ext(paq->gamma13_adj, paq->gamma13_adj_ext, I, J); \
  paq->d##I##d##J##g22 = dder_ext(paq->gamma22_adj, paq->gamma22_adj_ext, I, J); \
  paq->d##I##d##J##g23 = dder_ext(paq->gamma23_adj, paq->gamma23_adj_ext, I, J); \
  paq->d##I##d##J##g33 = dder_ext(paq->gamma33_adj, paq->gamma33_adj_ext, I, J)


/*
 * Evolution equations for indexed components
 */

#define BSSN_DT_GAMMAIJ(I, J) ( \
    - 2.0*paq->A##I##J \
  )

#define BSSN_DT_AIJ(I, J) ( \
    exp(-4.0*paq->phi)*(paq->ricciTF##I##J - 8.0*PI*paq->STF##I##J) \
    + (paq->K*paq->A##I##J - 2.0*( \
        paq->gammai11*paq->A1##I*paq->A1##J + paq->gammai12*paq->A1##I*paq->A2##J + paq->gammai13*paq->A1##I*paq->A3##J \
        + paq->gammai21*paq->A2##I*paq->A1##J + paq->gammai22*paq->A2##I*paq->A2##J + paq->gammai23*paq->A2##I*paq->A3##J \
        + paq->gammai31*paq->A3##I*paq->A1##J + paq->gammai32*paq->A3##I*paq->A2##J + paq->gammai33*paq->A3##I*paq->A3##J \
      )) \
  )

#define BSSN_DT_GAMMAI(I) ( \
    + 2.0*( \
        paq->G##I##11*paq->Acont11 + paq->G##I##22*paq->Acont22 + paq->G##I##33*paq->Acont33 \
          + 2.0*(paq->G##I##12*paq->Acont12 + paq->G##I##13*paq->Acont13 + paq->G##I##23*paq->Acont23) \
        - (2.0/3.0) * (paq->gammai##I##1*der_ext(paq->K_adj, paq->K_adj_ext, 1) + paq->gammai##I##2*der_ext(paq->K_adj, paq->K_adj_ext, 2) + paq->gammai##I##3*der_ext(paq->K_adj, paq->K_adj_ext, 3)) \
        - 8.0*PI*(paq->gammai##I##1*paq->S1 + paq->gammai##I##2*paq->S2 + paq->gammai##I##3*paq->S3) \
        + 6.0 * (paq->Acont##I##1*paq->d1phi + paq->Acont##I##2*paq->d2phi + paq->Acont##I##3*paq->d3phi) \
      ) \
  )

#define BSSN_MI(I) ( \
    /* http://arxiv.org/pdf/1109.5782v2.pdf */ \
    -2.0/3.0*der_ext(paq->K_adj, paq->K_adj_ext, I) - 8*PI*paq->S##I \
    + 6.0*( \
      paq->gamma11*paq->A1##I*paq->d1phi + paq->gamma21*paq->A2##I*paq->d1phi + paq->gamma31*paq->A3##I*paq->d1phi \
      + paq->gamma12*paq->A1##I*paq->d2phi + paq->gamma22*paq->A2##I*paq->d2phi + paq->gamma32*paq->A3##I*paq->d2phi \
      + paq->gamma13*paq->A1##I*paq->d3phi + paq->gamma23*paq->A2##I*paq->d3phi + paq->gamma33*paq->A3##I*paq->d3phi \
    ) + ( \
      /* (gamma^jk D_j A_ki) */ \
      paq->gammai11*der(paq->A1##I##_adj, 1) + paq->gammai12*der(paq->A1##I##_adj, 2) + paq->gammai13*der(paq->A1##I##_adj, 3) \
      + paq->gammai21*der(paq->A2##I##_adj, 1) + paq->gammai22*der(paq->A2##I##_adj, 2) + paq->gammai23*der(paq->A2##I##_adj, 3) \
      + paq->gammai31*der(paq->A3##I##_adj, 1) + paq->gammai32*der(paq->A3##I##_adj, 2) + paq->gammai33*der(paq->A3##I##_adj, 3) \
      - paq->Gamma1*paq->A1##I - paq->Gamma2*paq->A2##I - paq->Gamma3*paq->A3##I \
      - paq->GL11##I*paq->Acont11 - paq->GL22##I*paq->Acont22 - paq->GL33##I*paq->Acont33 \
      - 2.0*(paq->GL12##I*paq->Acont12 + paq->GL13##I*paq->Acont13 + paq->GL23##I*paq->Acont23) \
    ) - 2.0*paq->d##I##phi*( \
      paq->gammai11*paq->A11 + paq->gammai22*paq->A22 + paq->gammai33*paq->A33 \
      + 2.0*(paq->gammai12*paq->A12 + paq->gammai13*paq->A13 + paq->gammai23*paq->A23) \
    ) \
  )

#define BSSN_MI_MAG(I) ( \
    fabs(2.0/3.0*der_ext(paq->K_adj, paq->K_adj_ext, I) - 8*PI*paq->S##I) \
    + 6.0*fabs( \
      paq->gamma11*paq->A1##I*paq->d1phi + paq->gamma21*paq->A2##I*paq->d1phi + paq->gamma31*paq->A3##I*paq->d1phi \
      + paq->gamma12*paq->A1##I*paq->d2phi + paq->gamma22*paq->A2##I*paq->d2phi + paq->gamma32*paq->A3##I*paq->d2phi \
      + paq->gamma13*paq->A1##I*paq->d3phi + paq->gamma23*paq->A2##I*paq->d3phi + paq->gamma33*paq->A3##I*paq->d3phi \
    ) + fabs( \
      /* (gamma^jk D_j A_ki) */ \
      paq->gammai11*der(paq->A1##I##_adj, 1) + paq->gammai12*der(paq->A1##I##_adj, 2) + paq->gammai13*der(paq->A1##I##_adj, 3) \
      + paq->gammai21*der(paq->A2##I##_adj, 1) + paq->gammai22*der(paq->A2##I##_adj, 2) + paq->gammai23*der(paq->A2##I##_adj, 3) \
      + paq->gammai31*der(paq->A3##I##_adj, 1) + paq->gammai32*der(paq->A3##I##_adj, 2) + paq->gammai33*der(paq->A3##I##_adj, 3) \
      - paq->Gamma1*paq->A1##I - paq->Gamma2*paq->A2##I - paq->Gamma3*paq->A3##I \
      - paq->GL11##I*paq->Acont11 - paq->GL22##I*paq->Acont22 - paq->GL33##I*paq->Acont33 \
      - 2.0*(paq->GL12##I*paq->Acont12 + paq->GL13##I*paq->Acont13 + paq->GL23##I*paq->Acont23) \
    ) + fabs(2.0*paq->d##I##phi*( \
      paq->gammai11*paq->A11 + paq->gammai22*paq->A22 + paq->gammai33*paq->A33 \
      + 2.0*(paq->gammai12*paq->A12 + paq->gammai13*paq->A13 + paq->gammai23*paq->A23) \
    )) \
  )

/*
 * Full metric calcs
 */

// metric derivs

#define SET_DKM00(k) paq->d##k##m00 = DKM00(k);
#define DKM00(k) 0.0;

#define SET_DKM0I(k, i) paq->d##k##m0##i = DKM0I(k, i);
#define DKM0I(k, i) \
      exp(4.0*paq->phi)*( \
        (paq->d##k##g1##i + 4.0*paq->d##k##phi*paq->gamma1##i) \
        + (paq->d##k##g2##i + 4.0*paq->d##k##phi*paq->gamma2##i) \
        + (paq->d##k##g3##i + 4.0*paq->d##k##phi*paq->gamma3##i) \
      );

#define SET_DKMIJ(k, i, j) paq->d##k##m##i##j = DKMIJ(k, i, j);
#define DKMIJ(k, i, j) exp(4.0*paq->phi)*(paq->d##k##g##i##j + 4.0*paq->d##k##phi*paq->gamma##i##j);


// metric

#define SET_M00() paq->m00 = M00();
#define M00() (\
      -1.0 + exp(4.0*paq->phi) \
    );

#define SET_M0I(i) paq->m0##i = M0I(i);
#define M0I(i) 0.0;

#define SET_MIJ(i, j) paq->m##i##j = MIJ(i, j);
#define MIJ(i, j) exp(4.0*paq->phi)*(paq->gamma##i##j);

// inverse metric

#define SET_Mi00() paq->mi00 = Mi00();
#define Mi00() -1.0;

#define SET_Mi0I(i) paq->mi0##i = Mi0I(i);
#define Mi0I(i) 0.0;

#define SET_MiIJ(i, j) paq->mi##i##j = MiIJ(i, j);
#define MiIJ(i, j) exp(-4.0*paq->phi)*(paq->gamma##i##j);



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

#define A21_adj A12_adj
#define A31_adj A13_adj
#define A32_adj A23_adj

// local variables:
// ricci tensor
#define ricciTF21 ricciTF12
#define ricciTF31 ricciTF13
#define ricciTF32 ricciTF23

// covariant double-derivatives of phi
#define D2D1phi D1D2phi
#define D3D1phi D1D3phi
#define D3D2phi D2D3phi

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

// inverse metric derivatives
#define d1gi21 d1gi12
#define d1gi31 d1gi13
#define d1gi32 d1gi23
#define d2gi21 d2gi12
#define d2gi31 d2gi13
#define d2gi32 d2gi23
#define d3gi21 d3gi12
#define d3gi31 d3gi13
#define d3gi32 d3gi23

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
