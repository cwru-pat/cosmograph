
#include "Wave.h"

namespace cosmo
{

inline real_t derivative_stencil(const idx_t &i, const idx_t &j,
    const idx_t &k, const int &d, const real_t *field)
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

inline real_t lap_stencil(const idx_t &i, const idx_t &j, const idx_t &k, real_t *field)
{
  return (
      - 6.0*field[INDEX(i,j,k)]
      + field[INDEX(i+1,j,k)] + field[INDEX(i-1,j,k)]
      + field[INDEX(i,j+1,k)] + field[INDEX(i,j-1,k)]
      + field[INDEX(i,j,k+1)] + field[INDEX(i,j,k-1)]
    );
}


void Wave::step()
{
  idx_t i, j, k;

  // evolve fields
  INTERNAL_LOOP3(i, j, k)
  {
    const idx_t index = INDEX(i,j,k);

    const real_t lap_phi_p = lap_stencil(i, j, k, phi);
    const real_t lap_www_p = lap_stencil(i, j, k, www);

    phi[index] = phi_p[index] + dt*www_p[index] + dt*dt/2.0*lap_phi_p;
    www[index] = www_p[index] + dt*lap_phi_p + dt*dt/2.0*lap_www_p;
  }

  //step_boundary()

  /* swap buffers */
  std::swap(phi, phi_p);
  std::swap(www, www_p);
}

#if 0
void Wave::step_boundary()
{
  real_t lap_phi_p, lap_www_p;
  idx_t index, i, j, k;
  /* evolve fields */
  INTERNAL_LOOP3(i, j, k)
  {

  }
}
#endif

/* Wave constructor */
Wave::Wave()
{
  phi   = new real_t[N*N*N];
  phi_p = new real_t[N*N*N];
  www   = new real_t[N*N*N];
  www_p = new real_t[N*N*N];

  idx_t i,j,k;
  LOOP3(i, j, k)
  {
    phi[i] = 0;
    phi_p[i] = 0;
    www[i] = 0;
    www_p[i] = 0;
  }
}

/* Wave destructor */
Wave::~Wave()
{
  delete [] phi;
  delete [] phi_p;
  delete [] www;
  delete [] www_p;
}

}
