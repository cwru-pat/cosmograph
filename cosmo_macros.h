#ifndef COSMO_DEFINES
#define COSMO_DEFINES

#define N 128
#define POINTS (N*N*N)
// box size in hubble units
#define H_LEN_FRAC 0.5
#define dx (H_LEN_FRAC/(1.0*N))
#define dt (0.1*dx)

// WENO "epsilon" parameter
#define EPS 0.0001

// Numerical Error Damping strength parameters
#define KO_eta 0.0
#define BS_H_DAMPING_AMPLITUDE 40.0
#define JM_K_DAMPING_AMPLITUDE 0.0

#define Z4c_DAMPING 0
#define Z4c_K1_DAMPING_AMPLITUDE 0.0
#define Z4c_K2_DAMPING_AMPLITUDE 0.0

#define PI (4.0*atan(1.0))
#define SIGN(x) (((x) < 0.0) ? -1 : ((x) > 0.0))
#define pw2(x) ((x)*(x))
#define C_RE(c) ((c)[0])
#define C_IM(c) ((c)[1])
#define ROUND_2_IDXT(f) ((idx_t)(f >= 0.0 ? (f + 0.5) : (f - 0.5)))

#define RESTRICT __restrict__

// standard index, implementing periodic boundary conditions
#define INDEX(i,j,k) ( ((i+N)%N)*N*N + ((j+N)%N)*N + (k+N)%N )
// indexing without periodicity
#define NP_INDEX(i,j,k) (N*N*(i) + N*(j) + (k))
// indexing of flux arrays; indexes cell boundaries with 'd' = 1,2,3, using periodic BCs
#define F_INDEX(i,j,k,d) (((i+N)%N)*N*N*3 + ((j+N)%N)*N*3 + ((k+N)%N)*3 + (d+2)%3 )
// non-periodic indexing of flux arrays; indexes cell boundaries with 'd' = 1,2,3
#define F_NP_INDEX(i,j,k,d) ((i)*N*N*3 + (j)*N*3 + (k)*3 + (d+2)%3 )
// FFT indexing
#define FFT_INDEX(i,j,k) ((N/2+1)*N*((i+N)%N) + (N/2+1)*((j+N)%N) + ((k+N)%N))
// FFT indexing without periodicity
#define FFT_NP_INDEX(i,j,k) ((N/2+1)*N*(i) + (N/2+1)*(j) + (k))

#define LOOP3(i,j,k) \
  for(i=0; i<N; ++i) \
    for(j=0; j<N; ++j) \
      for(k=0; k<N; ++k)

#define INTERNAL_LOOP3(i,j,k) \
  for(i=1; i<N-1; ++i) \
    for(j=1; j<N-1; ++j) \
      for(k=1; k<N-1; ++k)

#define DECLARE_REAL_T(name) real_t name

// structure to store 27 immediately adjacent points
#define DECLARE_ADJACENT_REAL_T(name) real_t name##_adj[3][3][3]
// structure to store "faces" 2 points away (6 of them; _ext[dimension (x,y,z)][direction (-,+)])
#define DECLARE_ADJ_ADJACENT_REAL_T(name) real_t name##_adj_ext[3][2]

// indexes for accessing adjacent values
#define SET_LOCAL_INDEXES \
    idx_t idx001 = INDEX(paq->i-1, paq->j-1, paq->k   );   \
    idx_t idx010 = INDEX(paq->i-1, paq->j,   paq->k-1 );   \
    idx_t idx011 = INDEX(paq->i-1, paq->j,   paq->k   );   \
    idx_t idx012 = INDEX(paq->i-1, paq->j,   paq->k+1 );   \
    idx_t idx021 = INDEX(paq->i-1, paq->j+1, paq->k   );   \
    idx_t idx100 = INDEX(paq->i,   paq->j-1, paq->k-1 );   \
    idx_t idx101 = INDEX(paq->i,   paq->j-1, paq->k   );   \
    idx_t idx102 = INDEX(paq->i,   paq->j-1, paq->k+1 );   \
    idx_t idx110 = INDEX(paq->i,   paq->j,   paq->k-1 );   \
    idx_t idx111 = INDEX(paq->i,   paq->j,   paq->k   );   \
    idx_t idx112 = INDEX(paq->i,   paq->j,   paq->k+1 );   \
    idx_t idx120 = INDEX(paq->i,   paq->j+1, paq->k-1 );   \
    idx_t idx121 = INDEX(paq->i,   paq->j+1, paq->k   );   \
    idx_t idx122 = INDEX(paq->i,   paq->j+1, paq->k+1 );   \
    idx_t idx201 = INDEX(paq->i+1, paq->j-1, paq->k   );   \
    idx_t idx210 = INDEX(paq->i+1, paq->j,   paq->k-1 );   \
    idx_t idx211 = INDEX(paq->i+1, paq->j,   paq->k   );   \
    idx_t idx212 = INDEX(paq->i+1, paq->j,   paq->k+1 );   \
    idx_t idx221 = INDEX(paq->i+1, paq->j+1, paq->k   );

// Points, Faces, and Edges only
#define SET_LOCAL_VALUES_PFE(name) \
    paq->name = name##_a[paq->idx];                \
    paq->name##_adj[0][0][1] = name##_a[idx001];   \
    paq->name##_adj[0][1][0] = name##_a[idx010];   \
    paq->name##_adj[0][1][1] = name##_a[idx011];   \
    paq->name##_adj[0][1][2] = name##_a[idx012];   \
    paq->name##_adj[0][2][1] = name##_a[idx021];   \
    paq->name##_adj[1][0][0] = name##_a[idx100];   \
    paq->name##_adj[1][0][1] = name##_a[idx101];   \
    paq->name##_adj[1][0][2] = name##_a[idx102];   \
    paq->name##_adj[1][1][0] = name##_a[idx110];   \
    paq->name##_adj[1][1][1] = name##_a[idx111];   \
    paq->name##_adj[1][1][2] = name##_a[idx112];   \
    paq->name##_adj[1][2][0] = name##_a[idx120];   \
    paq->name##_adj[1][2][1] = name##_a[idx121];   \
    paq->name##_adj[1][2][2] = name##_a[idx122];   \
    paq->name##_adj[2][0][1] = name##_a[idx201];   \
    paq->name##_adj[2][1][0] = name##_a[idx210];   \
    paq->name##_adj[2][1][1] = name##_a[idx211];   \
    paq->name##_adj[2][1][2] = name##_a[idx212];   \
    paq->name##_adj[2][2][1] = name##_a[idx221];

// Points and Faces only
#define SET_LOCAL_VALUES_PF(name) \
    paq->name = name##_a[paq->idx];                \
    paq->name##_adj[0][1][1] = name##_a[idx011];   \
    paq->name##_adj[1][0][1] = name##_a[idx101];   \
    paq->name##_adj[1][1][0] = name##_a[idx110];   \
    paq->name##_adj[1][1][1] = name##_a[idx111];   \
    paq->name##_adj[1][1][2] = name##_a[idx112];   \
    paq->name##_adj[1][2][1] = name##_a[idx121];   \
    paq->name##_adj[2][1][1] = name##_a[idx211];   \

// Point only
#define SET_LOCAL_VALUES_P(name) \
    paq->name = name##_a[paq->idx];

// Extended faces
// _ext[dimension (x,y,z)][direction (-,+)])
#define SET_LOCAL_VALUES_F2(name) \
    paq->name##_adj_ext[0][0] = name##_a[INDEX(paq->i-2, paq->j  , paq->k  )];   \
    paq->name##_adj_ext[1][0] = name##_a[INDEX(paq->i  , paq->j-2, paq->k  )];   \
    paq->name##_adj_ext[2][0] = name##_a[INDEX(paq->i  , paq->j  , paq->k-2)];   \
    paq->name##_adj_ext[0][1] = name##_a[INDEX(paq->i+2, paq->j  , paq->k  )];   \
    paq->name##_adj_ext[1][1] = name##_a[INDEX(paq->i  , paq->j+2, paq->k  )];   \
    paq->name##_adj_ext[2][1] = name##_a[INDEX(paq->i  , paq->j  , paq->k+2)];

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
        real_t * name##_a, * name##_c, * name##_p, * name##_f

#define RK4_ARRAY_ALLOC(name) \
        name##_a   = new real_t[N*N*N]; \
        name##_c   = new real_t[N*N*N]; \
        name##_p   = new real_t[N*N*N]; \
        name##_f   = new real_t[N*N*N]

#define RK4_ARRAY_DELETE(name) \
        delete [] name##_a;    \
        delete [] name##_c;    \
        delete [] name##_p;    \
        delete [] name##_f


// A GEN2 method; any method needing 2 registers.
// Sets up a "_f" (final) and "_a" (active) register.
#define GEN2_ARRAY_ADDMAP(name)         \
        fields[#name "_f"] = name##_f;  \
        fields[#name "_a"] = name##_a

#define GEN2_ARRAY_CREATE(name) \
        real_t * name##_f, * name##_a

#define GEN2_ARRAY_ALLOC(name) \
        name##_f   = new real_t[N*N*N]; \
        name##_a   = new real_t[N*N*N]

#define GEN2_ARRAY_DELETE(name) \
        delete [] name##_f;    \
        delete [] name##_a


// A GEN1 method; just declares one register.
// Sets up an "_a" (active) register.
#define GEN1_ARRAY_ADDMAP(name)         \
        fields[#name "_a"] = name##_a

#define GEN1_ARRAY_CREATE(name) \
        real_t * name##_a

#define GEN1_ARRAY_ALLOC(name) \
        name##_a   = new real_t[N*N*N]

#define GEN1_ARRAY_DELETE(name) \
        delete [] name##_a


// A FLUX method; just declares one register.
// Sets up a "_a" (active) register.
#define FLUX_ARRAY_ADDMAP(name)         \
        fields[#name "_a"] = name##_a

#define FLUX_ARRAY_CREATE(name) \
        real_t * name##_a

#define FLUX_ARRAY_ALLOC(name) \
        name##_a   = new real_t[N*N*N*3]

#define FLUX_ARRAY_DELETE(name) \
        delete [] name##_a

#define ACC_DEF_SIM_FIELDS() \
    real_t * const phi_a = bssnSim.fields["phi_a"];

#define ACC_SIM_FIELDS bssnSim.fields["phi_a"][0:POINTS]


#endif
