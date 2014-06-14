#ifndef COSMO_DEFINES
#define COSMO_DEFINES

#define N 256
#define POINTS (N*N*N)
#define dt 0.1

#define RESTRICT __restrict__

#define INDEX(i,j,k) (((i+N)%N)*N*N + ((j+N)%N)*N + (k+N)%N)

#define LOOP3(i,j,k) \
  for(idx_t i=0; i<N; ++i) \
    for(idx_t j=0; j<N; ++j) \
      for(idx_t k=0; k<N; ++k)

#define INTERNAL_LOOP3(i,j,k) \
  for(idx_t i=1; i<N-1; ++i) \
    for(idx_t j=1; j<N-1; ++j) \
      for(idx_t k=1; k<N-1; ++k)


#define RK4_ARRAY_CREATE(name) \
        real_t * name##_a, * name##_k, * name##_p, * name##_f

#define RK4_ARRAY_ALLOC(name) \
        name##_a   = new real_t[N*N*N]; \
        name##_k   = new real_t[N*N*N]; \
        name##_p   = new real_t[N*N*N]; \
        name##_f   = new real_t[N*N*N]

#define RK4_ARRAY_DELETE(name) \
        delete [] name##_a;    \
        delete [] name##_k;    \
        delete [] name##_p;    \
        delete [] name##_f

#define DECLARE_REAL_T(name) real_t name

#define SET_LOCAL_VALUES(name) paq->name = name##_a[paq->idx]

// RK4 method, using 4 "registers".  One for the "_p"revious step data, one
// for the data being "_a"ctively used for calculation, one for the
// Runge-"_k"utta coefficient being calculated, and lastly the "_f"inal
// result of the calculation.
 #define RK4_ARRAY_ADDMAP(name)          \
        fields[#name "_a"] = name##_a;  \
        fields[#name "_k"] = name##_k;  \
        fields[#name "_p"] = name##_p;  \
        fields[#name "_f"] = name##_f

#endif
