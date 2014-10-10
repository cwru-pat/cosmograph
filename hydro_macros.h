#ifndef HYDRO_MACROS
#define HYDRO_MACROS

#define HYDRO_APPLY_TO_FIELDS(function) \
    /* conserved quantities */          \
    function(UD);                       \
    function(US1);                      \
    function(US2);                      \
    function(US3);                      \
    function(Ut);                       \
    /* primitives */                    \
    function(Ut);                       \

#define HYDRO_APPLY_TO_FLUXES(function) \
    /* fluxes */                        \
    function(FD);                       \
    function(FS1);                      \
    function(FS2);                      \
    function(FS3);                      \
    function(Ft);

#define HYDRO_APPLY_TO_SOURCES(function) \
    /* Sources */                        \
    function(SD);                        \
    function(SS1);                       \
    function(SS2);                       \
    function(SS3);                       \
    function(St);

#endif

