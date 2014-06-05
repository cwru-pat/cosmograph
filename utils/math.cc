
#include "math.h"

namespace cosmo
{

inline real_t lap_stencil(idx_t i, idx_t j, idx_t k, real_t *field)
{
return (
    - 6.0*field[INDEX(i,j,k)]
    + field[INDEX(i+1,j,k)] + field[INDEX(i-1,j,k)]
    + field[INDEX(i,j+1,k)] + field[INDEX(i,j-1,k)]
    + field[INDEX(i,j,k+1)] + field[INDEX(i,j,k-1)]
  );
}

inline real_t derivative_stencil(idx_t i, idx_t j, idx_t k, int d, real_t *field)
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

}
