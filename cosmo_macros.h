#ifndef COSMO_DEFINES
#define COSMO_DEFINES

/* Cosmological run parameters: */
#define N 64
#define NX N
#define NY N
#define NZ N
#define POINTS ((NX)*(NY)*(NZ))
// box size in hubble units
#define H_LEN_FRAC 0.5
#define dx (H_LEN_FRAC/(1.0*N))
#define dt (0.1*dx)
/****/

/* Stability test parameters: 
#define R 1
#define N 50
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

#define USE_REFERENCE_FRW true
#define NORMALIZE_GAMMAIJ_AIJ true

// Numerical Error Damping strength parameters
#define KO_ETA 1.0
#define BS_H_DAMPING_AMPLITUDE 1.0
#define JM_K_DAMPING_AMPLITUDE 0.0

#define USE_Z4c_DAMPING false
#if USE_Z4c_DAMPING
  #define Z4c_K1_DAMPING_AMPLITUDE 0.05
  #define Z4c_K2_DAMPING_AMPLITUDE 0.0
#else
  #define Z4c_K1_DAMPING_AMPLITUDE 0.0
  #define Z4c_K2_DAMPING_AMPLITUDE 0.0
#endif

#define USE_CONFORMAL_SYNC_ALPHA false
#define USE_HARMONIC_ALPHA false

#define USE_BSSN_SHIFT false

// Stencil order
#define STENCIL_ORDER 8

#define STENCIL_CONCATENATOR(function, order) function ## order
#define STENCIL_EVALUATOR(function, order) STENCIL_CONCATENATOR(function, order)
#define STENCIL_ORDER_FUNCTION(function) STENCIL_EVALUATOR(function, STENCIL_ORDER)

#define USE_WEDGE_INTEGRATOR false

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
// FFT indexing
#define FFT_INDEX(i,j,k) ((NZ/2+1)*NY*((i+NX)%NX) + (NZ/2+1)*((j+NY)%NY) + ((k+NZ)%NZ))
// FFT indexing without periodicity
#define FFT_NP_INDEX(i,j,k) ((NZ/2+1)*NY*(i) + (NZ/2+1)*(j) + (k))

#define LOOP3(i,j,k) \
  for(i=0; i<NX; ++i) \
    for(j=0; j<NY; ++j) \
      for(k=0; k<NZ; ++k)

#define PARALLEL_LOOP3(i,j,k) \
  _Pragma("omp parallel for default(shared) private(i, j, k)") \
  LOOP3(i,j,k)

#define AREA_LOOP(j,k) \
  for(j=0; j<NY; ++j) \
    for(k=0; k<NZ; ++k)

#define PARALLEL_AREA_LOOP(j,k) \
  _Pragma("omp parallel for default(shared) private(j, k)") \
  AREA_LOOP(j,k)

#define INTERNAL_LOOP3(i,j,k) \
  for(i=1; i<NX-1; ++i) \
    for(j=1; j<NY-1; ++j) \
      for(k=1; k<NZ-1; ++k)

#define DECLARE_REAL_T(name) real_t name

// structure to store 27 immediately adjacent points
#define DECLARE_ADJACENT_REAL_T(name) real_t name##_adj[3][3][3]
// structure to store "faces" 2 points away (6 of them; _ext[dimension (x,y,z)][direction (-,+)])
#define DECLARE_ADJ_ADJACENT_REAL_T(name) real_t name##_adj_ext[3][2]

#define SET_LOCAL_VALUES_P(name) \
    paq->name = name##_a[paq->idx];

#if USE_WEDGE_INTEGRATOR

    #define AREA (NY*NZ)
    #define STENCIL_SUPPORT_POINTS (STENCIL_ORDER + 1)

    #define WEDGE_SLICE_1_LEN ( 3*STENCIL_SUPPORT_POINTS - 2 )
    #define WEDGE_SLICE_2_LEN ( 2*STENCIL_SUPPORT_POINTS - 1 )
    #define WEDGE_SLICE_3_LEN ( STENCIL_SUPPORT_POINTS )
    
    #define WEDGE_TAIL_LEN  (3*STENCIL_SUPPORT_POINTS + 1)/2;
    #define WEDGE_AFTER_LEN 3*(STENCIL_SUPPORT_POINTS - 1)/2;

    #if STENCIL_SUPPORT_POINTS > NX
        #error "Grid in NX direction must be larger than stencil base!"
    #endif

    // RK4 method, using a single register by means of a "wedge" integrator.
     #define RK4_ARRAY_ADDMAP(name)          \
            fields[#name "_p"]  = name##_p;  \
            fields[#name "_K1"] = name##_K1; \
            fields[#name "_K2"] = name##_K2; \
            fields[#name "_K3"] = name##_K3; \
            fields[#name "_t"]  = name##_t;  \
            fields[#name "_r"]  = name##_r;  \
            fields[#name "_a"]  = name##_a;

    #define RK4_ARRAY_CREATE(name) \
            real_t * name##_p, * name##_K1, * name##_K2, * name##_K3, * name##_t, * name##_r , * name##_a

    #define RK4_ARRAY_ALLOC(name)                             \
            name##_p    = new real_t[POINTS];                 \
            name##_K1   = new real_t[AREA*WEDGE_SLICE_1_LEN]; \
            name##_K2   = new real_t[AREA*WEDGE_SLICE_2_LEN]; \
            name##_K3   = new real_t[AREA*WEDGE_SLICE_3_LEN]; \
            name##_t    = new real_t[AREA*WEDGE_TAIL_LEN];    \
            name##_r    = new real_t[AREA*WEDGE_AFTER_LEN];   \
            name##_a = name##_p;

    #define RK4_ARRAY_DELETE(name) \
            delete [] name##_p;    \
            delete [] name##_K1;   \
            delete [] name##_K2;   \
            delete [] name##_K3;   \
            delete [] name##_t;    \
            delete [] name##_r;

#else
    // RK4 method, using 4 "registers".  One for the "_p"revious step data, one
    // for the data being "_a"ctively used for calculation, one for the
    // Runge-Kutta "_c"oefficient being calculated, and lastly the "_f"inal
    // result of the calculation.
     #define RK4_ARRAY_ADDMAP(name)         \
            fields[#name "_a"] = name##_a;  \
            fields[#name "_c"] = name##_c;  \
            fields[#name "_p"] = name##_p;  \
            fields[#name "_f"] = name##_f

    #define RK4_ARRAY_CREATE(name) \
            periodicArray<idx_t, real_t> * name##_a, * name##_c, * name##_p, * name##_f \
                    , * name##_t, * name##_r, * name##_K1, * name##_K2, * name##_K3

    #define RK4_ARRAY_ALLOC(name) \
            name##_a   = new periodicArray<idx_t, real_t> (NX, NY, NZ); \
            name##_c   = new periodicArray<idx_t, real_t> (NX, NY, NZ); \
            name##_p   = new periodicArray<idx_t, real_t> (NX, NY, NZ); \
            name##_f   = new periodicArray<idx_t, real_t> (NX, NY, NZ)

    #define RK4_ARRAY_DELETE(name) \
            delete name##_a;    \
            delete name##_c;    \
            delete name##_p;    \
            delete name##_f

#endif

// A GEN1 method; just declares one register.
// Sets up an "_a" (active) register.
#define GEN1_ARRAY_ADDMAP(name)         \
        fields[#name "_a"] = name##_a

#define GEN1_ARRAY_CREATE(name) \
        periodicArray<idx_t, real_t> * name##_a

#define GEN1_ARRAY_ALLOC(name) \
        name##_a   = new periodicArray<idx_t, real_t> (NX, NY, NZ);

#define GEN1_ARRAY_DELETE(name) \
        delete name##_a

#endif
