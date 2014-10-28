#ifndef HYDRO_MACROS
#define HYDRO_MACROS

// No \tau term needed; P = w \rho EOS.

#define HYDRO_APPLY_TO_FIELDS(function) \
    /* conserved quantities */          \
    function(UD);                       \
    function(US1);                      \
    function(US2);                      \
    function(US3);                      \
    /*function(Ut);*/

#define HYDRO_APPLY_TO_FLUXES(function) \
    /* fluxes */                        \
    function(FD);                       \
    function(FS1);                      \
    function(FS2);                      \
    function(FS3);                      \
    /*function(Ft);*/

#define HYDRO_APPLY_TO_SOURCES(function) \
    /* Sources */                        \
    /* function(SD); = 0 */              \
    function(SS1);                       \
    function(SS2);                       \
    function(SS3);                       \
    /*function(St);*/

#define HYDRO_APPLY_TO_PRIMITIVES(function) \
    /* primitives */                        \
    function(r);                            \
    function(v1);                           \
    function(v2);                           \
    function(v3);

#define HYDRO_APPLY_TO_FLUXES_INT(function)     \
    /* fluxes @ interfaces */                   \
    function(FD_int);                           \
    function(FS1_int);                          \
    function(FS2_int);                          \
    function(FS3_int);


// T^ab * g_ab,j source term in flud EOM

#define TABGAB_J(j)  (                                                       \
      paq->d##j##m00 * ( uu_fac + w_EOS*paq->m00 )                           \
      + paq->d##j##m11 * ( uu_fac*v1_a[idx]*v1_a[idx] + w_EOS*paq->m11 )     \
      + paq->d##j##m22 * ( uu_fac*v2_a[idx]*v2_a[idx] + w_EOS*paq->m22 )     \
      + paq->d##j##m33 * ( uu_fac*v3_a[idx]*v3_a[idx] + w_EOS*paq->m33 )     \
      + 2.0*(                                                                \
          + paq->d##j##m01 * ( uu_fac*v1_a[idx] + w_EOS*paq->m01 )           \
          + paq->d##j##m02 * ( uu_fac*v2_a[idx] + w_EOS*paq->m02 )           \
          + paq->d##j##m03 * ( uu_fac*v3_a[idx] + w_EOS*paq->m03 )           \
          + paq->d##j##m12 * ( uu_fac*v1_a[idx]*v2_a[idx] + w_EOS*paq->m12 ) \
          + paq->d##j##m13 * ( uu_fac*v1_a[idx]*v3_a[idx] + w_EOS*paq->m13 ) \
          + paq->d##j##m23 * ( uu_fac*v2_a[idx]*v3_a[idx] + w_EOS*paq->m23 ) \
        )                                                                    \
    )


#endif
