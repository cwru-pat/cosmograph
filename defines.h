#ifndef COSMO_DEFINES
#define COSMO_DEFINES

#define N 256
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
        real_t name, name_p, name_i;

#define RK4_ARRAY_ALLOC(name) \
        name       = new real_t[N*N*N]; \
        name_p     = new real_t[N*N*N]; \
        name_i     = new real_t[N*N*N];

#define RK4_ARRAY_DELETE(name) \
        delete [] name;        \
        delete [] name_p;      \
        delete [] name_i;

// RK4 has a diagonal tableau, so we only need to compute
// the coefficients one at a time ("_i" arrays), given the
// values in the previous step ("_p" arrays).
#define RK4_ARRAY_ADDMAP(name)        \
        fields["name"]   = name;      \
        fields["name_p"] = name_p;    \
        fields["name_i"] = name_i;

#define RK4_ARRAY_CYCLE(name) \
        std::swap(name, name_p);

#define BSSN_APPLY_TO_FIELDS(function)  \
        function(gammaxx);              \
        function(gammaxy);              \
        function(gammaxy);              \
        function(gammayy);              \
        function(gammayz);              \
        function(gammazz);              \
        function(gammaxz);              \
        function(W);                    \
        function(Axx);                  \
        function(Axy);                  \
        function(Ayy);                  \
        function(Ayz);                  \
        function(Azz);                  \
        function(Axz);                  \
        function(K);                    \
        function(Gammax);               \
        function(Gammay);               \
        function(Gammaz);


#endif
