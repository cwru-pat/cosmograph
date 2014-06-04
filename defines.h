#ifndef COSMO_DEFINES
#define COSMO_DEFINES

#define N 256
#define dt 0.1

#define INDEX(i,j,k) (((i+N)%N)*N*N + ((j+N)%N)*N + (k+N)%N)

#define LOOP3(i,j,k) for(i=0; i<N; i++) for(j=0; j<N; j++) for(k=0; k<N; k++)

#define INTERNAL_LOOP3(i,j,k) for(i=1; i<N-1; i++) for(j=1; j<N-1; j++) for(k=1; k<N-1; k++)

#endif
