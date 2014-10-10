#ifndef HYDRO_MACROS
#define HYDRO_MACROS

#define HYDRO_APPLY_TO_FIELDS(function) \
    /* conserved quantities */          \
    function(UD);                       \
    function(US1);                      \
    function(US2);                      \
    function(US3);                      \
    function(Ut);                       \
    /* fluxes */                        \
    function(FS11);                     \
    function(FS12);                     \
    function(FS13);                     \
    function(FS21);                     \
    function(FS22);                     \
    function(FS23);                     \
    function(FS31);                     \
    function(FS32);                     \
    function(FS33);                     \
    /* Sources */                       \
    function(SD);                       \
    function(SS1);                      \
    function(SS2);                      \
    function(SS3);                      \
    function(St);

#endif

