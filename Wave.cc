
#include "Wave.h"
#include "globals.h"

namespace cosmo
{


void Wave::step()
{
  idx_t index;
  real_t lap_www_p, lap_phi_p;

  // evolve fields
  LOOP3(i, j, k)
  {
    index = INDEX(i,j,k);

    lap_phi_p = lap_stencil(i, j, k, phi_p);
    //lap_www_p = lap_stencil(i, j, k, www_p);

    phi[index] = phi_p[index] + dt*www_p[index]/* + dt*dt/2.0*lap_phi_p*/;
    www[index] = www_p[index] + dt*lap_phi_p/* + dt*dt/2.0*lap_www_p*/;
  }

  /* swap buffers */
  std::swap(phi, phi_p);
  std::swap(www, www_p);
}

/* Initial conditions */
void Wave::init()
{
  LOOP3(i, j, k)
  {
    phi[INDEX(i,j,k)] = 0;
    phi_p[INDEX(i,j,k)] = exp( -pow((real_t) i - N/2,2.0)/100 );;
    www[INDEX(i,j,k)] = 0;
    www_p[INDEX(i,j,k)] = 0;
  }
}

void Wave::dump_strip(std::string field, int axis, idx_t n1, idx_t n2)
{
  io_dump_strip(fields[field], axis, n1, n2);
}

/* Wave constructor */
Wave::Wave()
{
  phi   = new real_t[N*N*N];
  phi_p = new real_t[N*N*N];
  www   = new real_t[N*N*N];
  www_p = new real_t[N*N*N];

  fields["www"]   = www;
  fields["phi"]   = phi;
  fields["phi_p"] = phi_p;
  fields["www_p"] = www_p;
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
