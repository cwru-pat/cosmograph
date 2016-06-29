#ifndef COSMO_MACROS
#define COSMO_MACROS

/* Cosmological run parameters: */
// Code runs much slower with global N (slow indexing?)
#define N 32
#define NX N
#define NY N
#define NZ N
#define POINTS ((NX)*(NY)*(NZ))
// physical box size (eg., in hubble units)
#define H_LEN_FRAC 0.5
// #define dx (H_LEN_FRAC/(1.0*N))
// #define dt (0.1*dx)
/****/


/* Stability test parameters: 
#define R 1
#define N 128
#define NX (50*R)
// NY, NZ should be at least the # points in the stencil
#define NY 10
#define NZ 10
#define POINTS ((NX)*(NY)*(NZ))
// box size in hubble units
#define H_LEN_FRAC 0.5
#define dx (1.0/NX)
#define dt (0.1*dx)
****/

#define USE_REFERENCE_FRW false
#define NORMALIZE_GAMMAIJ_AIJ true

// Numerical Error Damping strength parameters
#define KO_ETA 1.0
#define BS_H_DAMPING_AMPLITUDE 1.0
#define JM_K_DAMPING_AMPLITUDE 0.0

// not thoroughly tested / debugged:
#define USE_Z4c_DAMPING false
#if USE_Z4c_DAMPING
  #define Z4c_K1_DAMPING_AMPLITUDE 0.05
  #define Z4c_K2_DAMPING_AMPLITUDE 0.0
#else
  #define Z4c_K1_DAMPING_AMPLITUDE 0.0
  #define Z4c_K2_DAMPING_AMPLITUDE 0.0
#endif

// unfinished: #define USE_CONFORMAL_SYNC_ALPHA false
#define USE_HARMONIC_ALPHA false // dust sims require sync. gauge

#define USE_BSSN_SHIFT true

// Stencil order
#define STENCIL_ORDER 8

#define STENCIL_CONCATENATOR(function, order) function ## order
#define STENCIL_EVALUATOR(function, order) STENCIL_CONCATENATOR(function, order)
#define STENCIL_ORDER_FUNCTION(function) STENCIL_EVALUATOR(function, STENCIL_ORDER)

// WENO "epsilon" parameter
#define EPS 0.0001

#define PI (4.0*atan(1.0))
#define SIGN(x) (((x) < 0.0) ? -1 : ((x) > 0.0))
#define pw2(x) ((x)*(x))
#define C_RE(c) ((c)[0])
#define C_IM(c) ((c)[1])
#define ROUND_2_IDXT(f) ((idx_t)(f >= 0.0 ? (f + 0.5) : (f - 0.5)))

#define RESTRICT __restrict__

// standard index, implementing periodic boundary conditions
#define INDEX(i,j,k) ( ((i+NX)%(NX))*(NY)*(NZ) + ((j+NY)%(NY))*(NZ) + (k+NZ)%(NZ) )
// indexing without periodicity
#define NP_INDEX(i,j,k) ((NZ)*(NY)*(i) + (NZ)*(j) + (k))
// indexing of flux arrays; indexes cell boundaries with 'd' = 1,2,3, using periodic BCs
#define F_INDEX(i,j,k,d) ( ((i+NX)%(NX))*(NY)*(NZ)*3 + ((j+NY)%(NY))*(NZ)*3 + (k+NZ)%(NZ)*3 + (d+2)%3 )
// non-periodic indexing of flux arrays; indexes cell boundaries with 'd' = 1,2,3
#define F_NP_INDEX(i,j,k,d) ( (NZ)*(NY)*(i)*3 + (NZ)*(j)*3 + (k)*3 + (d+2)%3 )

// map spatial (i,j) to array index
#define aIDX(i,j) ( i <= j  ? (7-i)*i/2 - 4 + j : (7-j)*j/2 - 4 + i )
// map spatial (i) to array index
#define iIDX(i) ( i - 1 )


#define LOOP3(i,j,k) \
  for(i=0; i<NX; ++i) \
    for(j=0; j<NY; ++j) \
      for(k=0; k<NZ; ++k)

#define AREA_LOOP(j,k) \
  for(j=0; j<NY; ++j) \
    for(k=0; k<NZ; ++k)

#define INTERNAL_LOOP3(i,j,k) \
  for(i=1; i<NX-1; ++i) \
    for(j=1; j<NY-1; ++j) \
      for(k=1; k<NZ-1; ++k)

#define DECLARE_REAL_T(name) real_t name

// structure to store 27 immediately adjacent points
#define DECLARE_ADJACENT_REAL_T(name) real_t name##_adj[3][3][3]
// structure to store "faces" 2 points away (6 of them; _ext[dimension (x,y,z)][direction (-,+)])
#define DECLARE_ADJ_ADJACENT_REAL_T(name) real_t name##_adj_ext[3][2]


// RK4 method, using 4 "registers".  One for the "_p"revious step data, one
// for the data being "_a"ctively used for calculation, one for the
// Runge-Kutta "_c"oefficient being calculated, and lastly the "_f"inal
// result of the calculation.
#define RK4_ARRAY_ADDMAP(name) \
        fields[#name "_a"] = & name->_array_a; \
        fields[#name "_c"] = & name->_array_c; \
        fields[#name "_p"] = & name->_array_p; \
        fields[#name "_f"] = & name->_array_f

#define RK4_ARRAY_CREATE(name) \
        register_t * name

#define RK4_ARRAY_ALLOC(name) \
        name = new register_t(); \
        name->init(NX, NY, NZ, dt)

#define RK4_ARRAY_DELETE(name) \
        name->~RK4Register()

#define RK4_SET_LOCAL_VALUES(name) \
    bd->name = name->_array_a[bd->idx];

// A GEN1 method; just declares one register.
// Sets up an "_a" (active) register.
#define GEN1_ARRAY_ADDMAP(name)         \
        fields[#name "_a"] = &name##_a

#define GEN1_ARRAY_CREATE(name) \
        arr_t name##_a

#define GEN1_ARRAY_ALLOC(name) \
        name##_a.init(NX, NY, NZ)

#define GEN1_ARRAY_DELETE(name) \
        name##_a.~CosmoArray()

#define GEN1_SET_LOCAL_VALUES(name) \
    bd->name = name##_a[bd->idx];

#endif
