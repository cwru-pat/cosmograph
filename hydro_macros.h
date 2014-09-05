#ifndef HYDRO_MACROS
#define HYDRO_MACROS

#define HYDRO_APPLY_TO_FIELDS(function) \
  /* conserved variables */             \
    function(U0);                       \
    function(U1);                       \
    function(U2);                       \
    function(U3);                       \
    function(U4);                       \
  /* flux variables */                  \
    function(F0_1);                     \
    function(F1_1);                     \
    function(F2_1);                     \
    function(F3_1);                     \
    function(F4_1);                     \
    function(F0_2);                     \
    function(F1_2);                     \
    function(F2_2);                     \
    function(F3_2);                     \
    function(F4_2);                     \
    function(F0_3);                     \
    function(F1_3);                     \
    function(F2_3);                     \
    function(F3_3);                     \
    function(F4_3);                     \
  /* source terms */                    \
    function(S1);                       \
    function(S2);                       \
    function(S3);                       \
    function(S4);

#endif
