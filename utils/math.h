#ifndef COSMO_UTILS_MATH_H
#define COSMO_UTILS_MATH_H

#include "../cosmo.h"

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

inline real_t derivative_stencil(idx_t i, idx_t j, idx_t k, int d,
    real_t *field)
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

  /* XXX */
  return 0;
}

inline real_t double_derivative_stencil(idx_t i, idx_t j, idx_t k, int d,
    real_t *field)
{
  switch (d) {
    case 1:
      return field[INDEX(i+1,j,k)] + field[INDEX(i-1,j,k)] - 2.0*field[INDEX(i,j,k)];
      break;
    case 2:
      return field[INDEX(i,j+1,k)] + field[INDEX(i,j-1,k)] - 2.0*field[INDEX(i,j,k)];
      break;
    case 3:
      return field[INDEX(i,j,k+1)] + field[INDEX(i,j,k-1)] - 2.0*field[INDEX(i,j,k)];
      break;
  }

  /* XXX */
  return 0;
}

inline real_t double_derivative(idx_t i, idx_t j, idx_t k, int d1, int d2,
    real_t *field)
{
  if(d1 == d2) {
    return double_derivative_stencil(i, j, k, d1, field);
  } else {
    switch (d1) {
      case 1:
        return derivative_stencil(i+1, j, k, d2, field) - derivative_stencil(i-1, j, k, d2, field);
        break;
      case 2:
        return derivative_stencil(i, j+1, k, d2, field) - derivative_stencil(i, j-1, k, d2, field);
        break;
      case 3:
        return derivative_stencil(i, j, k+1, d2, field) - derivative_stencil(i, j, k-1, d2, field);
        break;
    }
  }

  /* XXX */
  return 0;
}

inline real_t average(real_t *field)
{
  // note this may have poor precision for large datasets
  real_t sum = 0.0; 
  LOOP3(i, j, k)
  {
    sum += field[NP_INDEX(i,j,k)];
  }
  return sum/POINTS;
}

inline real_t average_lap(real_t *field)
{
  // note this may have poor precision for large datasets
  idx_t i, j, k;
  real_t sum = 0.0;
  LOOP3(i, j, k)
  {
    sum += lap_stencil(i, j, k, field)/dx/dx;
  }
  return sum/POINTS;
}

inline real_t standard_deviation(real_t *field)
{
  // note this may have poor precision for large datasets
  idx_t i, j, k;
  real_t sum = 0.0;
  real_t avg = average(field);
  LOOP3(i, j, k)
  {
    sum += pw2(avg - field[NP_INDEX(i,j,k)]);
  }
  return sqrt(sum/(POINTS-1));
}


}

#endif
