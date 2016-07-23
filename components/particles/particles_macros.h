#ifndef COSMO_PARTICLES_MACROS
#define COSMO_PARTICLES_MACROS

#define PARTICLES_INTERPOLATION(var) \
  linearInterpolation( \
    corner_pp_in[0][0][0].var, corner_pp_in[0][0][1].var, corner_pp_in[0][1][0].var, corner_pp_in[0][1][1].var, \
    corner_pp_in[1][0][0].var, corner_pp_in[1][0][1].var, corner_pp_in[1][1][0].var, corner_pp_in[1][1][1].var, \
    x_d \
  )

#define PARTICLES_ROUND(val) ((idx_t)( (val) + 0.5))

#define DER(field) (derivative(x_idx, y_idx, z_idx, a+1, field))

#define PARTICLES_PARALLEL_LOOP(pr) \
  typename particle_vec::iterator pr; \
  _Pragma("omp parallel for default(shared) private(pr)") \
  for(pr = particles->begin(); pr < particles->end(); ++pr) \

#endif
