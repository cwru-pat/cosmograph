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

#endif
