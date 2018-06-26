#ifndef COSMO_MACROS
#define COSMO_MACROS


/************************************************/
/* Variables that can be changed at compilation */
/* time using, eg, the -D option for gcc.       */
/************************************************/

// simulation size
// Code runs much slower with global N (slow indexing?)
#ifndef COSMO_N
  #define COSMO_N 16
#endif
#ifndef NX
  #define NX COSMO_N
#endif
#ifndef NY
  #define NY COSMO_N
#endif
#ifndef NZ
  #define NZ COSMO_N
#endif
#define POINTS ((NX)*(NY)*(NZ))

// physical box size (in units of the initial Hubble^-1 scale)
// eg; L = H_LEN_FRAC = N*dx
#ifndef H_LEN_FRAC
  #define H_LEN_FRAC 0.5
#endif

// compile using reference integrator?
#ifndef USE_REFERENCE_FRW
  #define USE_REFERENCE_FRW true
#endif

// Stencil order
#ifndef STENCIL_ORDER
  #define STENCIL_ORDER 8
#endif

// evolve shift as well? (if not, assumed to be zero)
#ifndef USE_BSSN_SHIFT
  #define USE_BSSN_SHIFT false
#endif

// Gamma-driver gauge settings (must turn on bssn_shift as well)
#ifndef USE_GAMMA_DRIVER
  #define USE_GAMMA_DRIVER false
#endif

// Gamma-driver gauge settings (must turn on bssn_shift as well)
#ifndef USE_GENERALIZED_NEWTON
  #define USE_GENERALIZED_NEWTON false
#endif


// Optionally exclude some second-order terms
#ifndef EXCLUDE_SECOND_ORDER_SMALL
  #define EXCLUDE_SECOND_ORDER_SMALL false
#endif
#ifndef EXCLUDE_SECOND_ORDER_FRW
  #define EXCLUDE_SECOND_ORDER_FRW false
#endif

// Optionally compile without raytracing, multigrid classes
#ifndef USE_MULTIGRID
  #define USE_MULTIGRID true
#endif
#ifndef USE_COSMOTRACE
  #define USE_COSMOTRACE true
#endif

//Potential types
#ifndef USE_COSMO_CONST_POTENTIAL
  #define USE_COSMO_CONST_POTENTIAL true
  #ifndef COSMO_CONST
    #define COSMO_CONST 0.0003
  #endif
#endif

/*****************************************/
/* Additional variable/macro definitions */
/*****************************************/

// not really tested:
#ifndef USE_Z4c_DAMPING
#  define USE_Z4c_DAMPING false
#endif
#if USE_Z4c_DAMPING
#  define Z4c_K1_DAMPING_AMPLITUDE 0.5
#  define Z4c_K2_DAMPING_AMPLITUDE 0.1
#else
#  define Z4c_K1_DAMPING_AMPLITUDE 0.0
#  define Z4c_K2_DAMPING_AMPLITUDE 0.0
#endif

#define STENCIL_CONCATENATOR(function, order) function ## order
#define STENCIL_EVALUATOR(function, order) STENCIL_CONCATENATOR(function, order)
#define STENCIL_ORDER_FUNCTION(function) STENCIL_EVALUATOR(function, STENCIL_ORDER)

#define PI (4.0*atan(1.0))
#define SIGN(x) (((x) < 0.0) ? -1 : ((x) > 0.0))
#define pw2(x) ((x)*(x))
#define C_RE(c) ((c)[0])
#define C_IM(c) ((c)[1])
#define ROUND_2_IDXT(f) ((idx_t)(f >= 0.0 ? (f + 0.5) : (f - 0.5)))

#define RESTRICT __restrict__

// standard index, implementing periodic boundary conditions
#define INDEX(i,j,k) ( ((i+4*NX)%(NX))*(NY)*(NZ) + ((j+4*NY)%(NY))*(NZ) + (k+4*NZ)%(NZ) )
// indexing without periodicity
#define NP_INDEX(i,j,k) ((NZ)*(NY)*(i) + (NZ)*(j) + (k))
// indexing of flux arrays; indexes cell boundaries with 'd' = 1,2,3, using periodic BCs
#define F_INDEX(i,j,k,d) ( ((i+NX)%(NX))*(NY)*(NZ)*3 + ((j+NY)%(NY))*(NZ)*3 + (k+NZ)%(NZ)*3 + (d+2)%3 )
// non-periodic indexing of flux arrays; indexes cell boundaries with 'd' = 1,2,3
#define F_NP_INDEX(i,j,k,d) ( (NZ)*(NY)*(i)*3 + (NZ)*(j)*3 + (k)*3 + (d+2)%3 )

// index with variable grid resolution
#define H_INDEX(i,j,k,nx,ny,nz) (((i+nx)%(nx))*(ny)*(nz) + ((j+ny)%(ny))*(nz) + (k+(nz))%(nz))

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
        delete name

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

#define GEN1_ARRAY_DELETE(name)

#define GEN1_SET_LOCAL_VALUES(name) \
    bd->name = name##_a[bd->idx];


// macros for summing
#define COSMO_SUMMATION_1(MACRO) \
  ( MACRO(1) + MACRO(2) + MACRO(3) )

#define COSMO_SUMMATION_1_ARGS(MACRO, ...) \
  ( MACRO(1, __VA_ARGS__) + MACRO(2, __VA_ARGS__) + MACRO(3, __VA_ARGS__) )

#define COSMO_SUMMATION_2(MACRO) ( \
  MACRO(1, 1) + MACRO(1, 2) + MACRO(1, 3) \
  + MACRO(2, 1) + MACRO(2, 2) + MACRO(2, 3) \
  + MACRO(3, 1) + MACRO(3, 2) + MACRO(3, 3) \
  )

#define COSMO_SUMMATION_2_ARGS(MACRO, ...) ( \
  MACRO(1, 1, __VA_ARGS__) + MACRO(1, 2, __VA_ARGS__) + MACRO(1, 3, __VA_ARGS__) \
  + MACRO(2, 1, __VA_ARGS__) + MACRO(2, 2, __VA_ARGS__) + MACRO(2, 3, __VA_ARGS__) \
  + MACRO(3, 1, __VA_ARGS__) + MACRO(3, 2, __VA_ARGS__) + MACRO(3, 3, __VA_ARGS__) \
  )

#define COSMO_SUMMATION_3(MACRO) ( \
  MACRO(1, 1, 1) + MACRO(1, 1, 2) + MACRO(1, 1, 3) \
  + MACRO(1, 2, 1) + MACRO(1, 2, 2) + MACRO(1, 2, 3) \
  + MACRO(1, 3, 1) + MACRO(1, 3, 2) + MACRO(1, 3, 3) \
  + MACRO(2, 1, 1) + MACRO(2, 1, 2) + MACRO(2, 1, 3) \
  + MACRO(2, 2, 1) + MACRO(2, 2, 2) + MACRO(2, 2, 3) \
  + MACRO(2, 3, 1) + MACRO(2, 3, 2) + MACRO(2, 3, 3) \
  + MACRO(3, 1, 1) + MACRO(3, 1, 2) + MACRO(3, 1, 3) \
  + MACRO(3, 2, 1) + MACRO(3, 2, 2) + MACRO(3, 2, 3) \
  + MACRO(3, 3, 1) + MACRO(3, 3, 2) + MACRO(3, 3, 3) \
  )

#define COSMO_SUMMATION_3_ARGS(MACRO, ...) ( \
  MACRO(1, 1, 1, __VA_ARGS__) + MACRO(1, 1, 2, __VA_ARGS__) + MACRO(1, 1, 3, __VA_ARGS__) \
  + MACRO(1, 2, 1, __VA_ARGS__) + MACRO(1, 2, 2, __VA_ARGS__) + MACRO(1, 2, 3, __VA_ARGS__) \
  + MACRO(1, 3, 1, __VA_ARGS__) + MACRO(1, 3, 2, __VA_ARGS__) + MACRO(1, 3, 3, __VA_ARGS__) \
  + MACRO(2, 1, 1, __VA_ARGS__) + MACRO(2, 1, 2, __VA_ARGS__) + MACRO(2, 1, 3, __VA_ARGS__) \
  + MACRO(2, 2, 1, __VA_ARGS__) + MACRO(2, 2, 2, __VA_ARGS__) + MACRO(2, 2, 3, __VA_ARGS__) \
  + MACRO(2, 3, 1, __VA_ARGS__) + MACRO(2, 3, 2, __VA_ARGS__) + MACRO(2, 3, 3, __VA_ARGS__) \
  + MACRO(3, 1, 1, __VA_ARGS__) + MACRO(3, 1, 2, __VA_ARGS__) + MACRO(3, 1, 3, __VA_ARGS__) \
  + MACRO(3, 2, 1, __VA_ARGS__) + MACRO(3, 2, 2, __VA_ARGS__) + MACRO(3, 2, 3, __VA_ARGS__) \
  + MACRO(3, 3, 1, __VA_ARGS__) + MACRO(3, 3, 2, __VA_ARGS__) + MACRO(3, 3, 3, __VA_ARGS__) \
  )

#endif
