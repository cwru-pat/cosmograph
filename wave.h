#ifndef COSMO_WAVE
#define COSMO_WAVE

/** BSSN class **/
class WAVE
{
  /* arrays for storing fields */

  // wave equation fields
  real_t *phi;
  real_t *phi_p;
  real_t *www;
  real_t *www_p;

public:
  void init()
  {
    phi = (real_t *) malloc(N*N*N*sizeof(real_t));
    phi_p = (real_t *) malloc(N*N*N*sizeof(real_t));

    www = (real_t *) malloc(N*N*N*sizeof(real_t));
    www_p = (real_t *) malloc(N*N*N*sizeof(real_t));

    idx_t i,j,k;
    LOOP3(i, j, k)
    {
      phi[INDEX(i,j,k)] = 0;
      phi_p[INDEX(i,j,k)] = 0;
      www[INDEX(i,j,k)] = 0;
      www_p[INDEX(i,j,k)] = 0;
    }
  }

  void step()
  {
    real_t lap_phi_p, lap_www_p;
    idx_t index, i, j, k;

    // evolve fields
    INTERNAL_LOOP3(i, j, k)
    {
      index = INDEX(i,j,k);

      lap_phi_p = lap_stencil(i, j, k, phi);
      lap_www_p = lap_stencil(i, j, k, www);

      phi[index] = phi_p[index] + dt*www_p[index] + dt*dt/2.0*lap_phi_p;
      www[index] = www_p[index] + dt*lap_phi_p + dt*dt/2.0*lap_www_p;
    }


    /* switch field with field_p */
    real_t * swap = phi_p;
    phi_p = phi;
    phi = swap;

    swap = www_p;
    www_p = www;
    www = swap;
  }

  real_t derivative_stencil(idx_t i, idx_t j, idx_t k, int d, real_t *field)
  {
    switch (d) {
      case 1:
        return field[INDEX(i+1,j,k)] - field[INDEX(i-1,j,k)];
        break;
      case 2:
        return field[INDEX(i,j+1,k)] - field[INDEX(i,j-1,k)];
        break;
      case 3:
        return field[INDEX(i,j,k+1)] - field[INDEX(i,j,k-1)];
        break;
    }
  }

  real_t lap_stencil(idx_t i, idx_t j, idx_t k, real_t *field)
  {
    return (
        - 6.0*field[INDEX(i,j,k)]
        + field[INDEX(i+1,j,k)] + field[INDEX(i-1,j,k)]
        + field[INDEX(i,j+1,k)] + field[INDEX(i,j-1,k)]
        + field[INDEX(i,j,k+1)] + field[INDEX(i,j,k-1)]
      );
  }

};

#endif
