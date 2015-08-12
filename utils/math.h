#ifndef COSMO_UTILS_MATH_H
#define COSMO_UTILS_MATH_H

#include "../cosmo.h"

namespace cosmo
{

inline real_t derivative(idx_t i, idx_t j, idx_t k, int d,
    real_t *field)
{
  switch (d) {
    case 1:
      return (
        1.0/280.0*field[INDEX(i-4,j,k)]
        - 4.0/105.0*field[INDEX(i-3,j,k)]
        + 1.0/5.0*field[INDEX(i-2,j,k)]
        - 4.0/5.0*field[INDEX(i-1,j,k)]
        + 4.0/5.0*field[INDEX(i+1,j,k)]
        - 1.0/5.0*field[INDEX(i+2,j,k)]
        + 4.0/105.0*field[INDEX(i+3,j,k)]
        - 1.0/280.0*field[INDEX(i+4,j,k)]
      )/dx;
      break;
    case 2:
      return (
        1.0/280.0*field[INDEX(i,j-4,k)]
        - 4.0/105.0*field[INDEX(i,j-3,k)]
        + 1.0/5.0*field[INDEX(i,j-2,k)]
        - 4.0/5.0*field[INDEX(i,j-1,k)]
        + 4.0/5.0*field[INDEX(i,j+1,k)]
        - 1.0/5.0*field[INDEX(i,j+2,k)]
        + 4.0/105.0*field[INDEX(i,j+3,k)]
        - 1.0/280.0*field[INDEX(i,j+4,k)]
      )/dx;
      break;
    case 3:
      return (
        1.0/280.0*field[INDEX(i,j,k-4)]
        - 4.0/105.0*field[INDEX(i,j,k-3)]
        + 1.0/5.0*field[INDEX(i,j,k-2)]
        - 4.0/5.0*field[INDEX(i,j,k-1)]
        + 4.0/5.0*field[INDEX(i,j,k+1)]
        - 1.0/5.0*field[INDEX(i,j,k+2)]
        + 4.0/105.0*field[INDEX(i,j,k+3)]
        - 1.0/280.0*field[INDEX(i,j,k+4)]
      )/dx;
      break;
  }

  /* XXX */
  return 0;
}

inline real_t mixed_derivative_stencil(idx_t i, idx_t j, idx_t k, int d1, int d2, real_t *field)
{
  if( (d1 == 1 && d2 == 2) || (d1 == 2 && d2 == 1) ) {
    return (
      135.0*(
        - field[INDEX(i+1,j-1,k)] + field[INDEX(i+1,j+1,k)]
        + field[INDEX(i-1,j-1,k)] - field[INDEX(i-1,j+1,k)]
      ) - 27.0/2.0*(
        - field[INDEX(i+2,j-2,k)] + field[INDEX(i+2,j+2,k)]
        + field[INDEX(i-2,j-2,k)] - field[INDEX(i-2,j+2,k)]
      ) + (
        - field[INDEX(i+3,j-3,k)] + field[INDEX(i+3,j+3,k)]
        + field[INDEX(i-3,j-3,k)] - field[INDEX(i-3,j+3,k)]
      )
    )/360.0/dx/dx;
  }

  if( (d1 == 1 && d2 == 3) || (d1 == 3 && d2 == 1) ) {
    return (
      135.0*(
        - field[INDEX(i+1,j,k-1)] + field[INDEX(i+1,j,k+1)]
        + field[INDEX(i-1,j,k-1)] - field[INDEX(i-1,j,k+1)]
      ) - 27.0/2.0*(
        - field[INDEX(i+2,j,k-2)] + field[INDEX(i+2,j,k+2)]
        + field[INDEX(i-2,j,k-2)] - field[INDEX(i-2,j,k+2)]
      ) + (
        - field[INDEX(i+3,j,k-3)] + field[INDEX(i+3,j,k+3)]
        + field[INDEX(i-3,j,k-3)] - field[INDEX(i-3,j,k+3)]
      )
    )/360.0/dx/dx;
  }

  if( (d1 == 3 && d2 == 2) || (d1 == 2 && d2 == 3) ) {
    return (
      135.0*(
        - field[INDEX(i,j+1,k-1)] + field[INDEX(i,j+1,k+1)]
        + field[INDEX(i,j-1,k-1)] - field[INDEX(i,j-1,k+1)]
      ) - 27.0/2.0*(
        - field[INDEX(i,j+2,k-2)] + field[INDEX(i,j+2,k+2)]
        + field[INDEX(i,j-2,k-2)] - field[INDEX(i,j-2,k+2)]
      ) + (
        - field[INDEX(i,j+3,k-3)] + field[INDEX(i,j+3,k+3)]
        + field[INDEX(i,j-3,k-3)] - field[INDEX(i,j-3,k+3)]
      )
    )/360.0/dx/dx;
  }

  /* XXX */
  return 0;
}

inline real_t double_derivative_stencil(idx_t i, idx_t j, idx_t k, int d,
    real_t *field)
{
  switch (d) {
    case 1:
      return (
          - 1.0/560.0*field[INDEX(i-4,j,k)]
          + 8.0/315.0*field[INDEX(i-3,j,k)]
          - 1.0/5.0*field[INDEX(i-2,j,k)]
          + 8.0/5.0*field[INDEX(i-1,j,k)]
          - 205.0/72.0*field[INDEX(i-0,j,k)]
          + 8.0/5.0*field[INDEX(i+1,j,k)]
          - 1.0/5.0*field[INDEX(i+2,j,k)]
          + 8.0/315.0*field[INDEX(i+3,j,k)]
          - 1.0/560.0*field[INDEX(i+4,j,k)]
        )/dx/dx;
      break;
    case 2:
      return (
          - 1.0/560.0*field[INDEX(i,j-4,k)]
          + 8.0/315.0*field[INDEX(i,j-3,k)]
          - 1.0/5.0*field[INDEX(i,j-2,k)]
          + 8.0/5.0*field[INDEX(i,j-1,k)]
          - 205.0/72.0*field[INDEX(i,j-0,k)]
          + 8.0/5.0*field[INDEX(i,j+1,k)]
          - 1.0/5.0*field[INDEX(i,j+2,k)]
          + 8.0/315.0*field[INDEX(i,j+3,k)]
          - 1.0/560.0*field[INDEX(i,j+4,k)]
        )/dx/dx;
      break;
    case 3:
      return (
          - 1.0/560.0*field[INDEX(i,j,k-4)]
          + 8.0/315.0*field[INDEX(i,j,k-3)]
          - 1.0/5.0*field[INDEX(i,j,k-2)]
          + 8.0/5.0*field[INDEX(i,j,k-1)]
          - 205.0/72.0*field[INDEX(i,j,k-0)]
          + 8.0/5.0*field[INDEX(i,j,k+1)]
          - 1.0/5.0*field[INDEX(i,j,k+2)]
          + 8.0/315.0*field[INDEX(i,j,k+3)]
          - 1.0/560.0*field[INDEX(i,j,k+4)]
        )/dx/dx;
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
    return mixed_derivative_stencil(i, j, k, d1, d2, field);
  }

  /* XXX */
  return 0;
}

inline real_t average(real_t *field)
{
  // note this may have poor precision for large datasets
  real_t sum = 0.0; 
  idx_t i, j, k;
  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:sum)
  LOOP3(i, j, k)
  {
    sum += field[NP_INDEX(i,j,k)];
  }
  return sum/POINTS;
}

inline real_t volume_average(real_t *phi)
{
  real_t sum = 0.0; 
  idx_t i, j, k;
  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:sum)
  LOOP3(i, j, k)
  {
    sum += exp(6.0*phi[NP_INDEX(i,j,k)]);
  }
  return sum/POINTS;
}

inline real_t conformal_average(real_t *field, real_t *phi)
{
  real_t sum = 0.0; 
  idx_t i, j, k;
  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:sum)
  LOOP3(i, j, k)
  {
    sum += exp(6.0*phi[NP_INDEX(i,j,k)])*field[NP_INDEX(i,j,k)];
  }
  real_t vol = volume_average(phi);
  return sum/POINTS/vol;
}

inline real_t standard_deviation(real_t *field, real_t avg)
{
  // note this may have poor precision for large datasets
  idx_t i, j, k;
  real_t sum = 0.0;
  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:sum)
  LOOP3(i, j, k)
  {
    sum += pw2(avg - field[NP_INDEX(i,j,k)]);
  }
  return sqrt(sum/(POINTS-1));
}

inline real_t standard_deviation(real_t *field)
{
  real_t avg = average(field);
  return standard_deviation(field, avg);
}

inline real_t conformal_standard_deviation(real_t *field, real_t *phi, real_t avg)
{
  real_t sum = 0.0; 
  idx_t i, j, k;

  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:sum)
  LOOP3(i, j, k)
  {
    sum += exp(6.0*phi[NP_INDEX(i,j,k)])*pw2(avg - field[NP_INDEX(i,j,k)]);
  }
  real_t vol = volume_average(phi);
  return sqrt(sum/(POINTS-1)/vol);
}

inline real_t conformal_standard_deviation(real_t *field, real_t *phi)
{
  real_t avg = conformal_average(field, phi);
  return conformal_standard_deviation(field, phi, avg);
}

inline idx_t numNaNs(real_t *field)
{
  idx_t i, j, k;
  idx_t NaNs = 0;
  
  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:NaNs)
  LOOP3(i,j,k)
  {
    real_t val = field[NP_INDEX(i,j,k)];
    union { float val; uint32_t x; } u = { (float) val };
    if((u.x << 1) > 0xff000000u)
    {
      NaNs += 1;
    }
  }

  return NaNs;
}

}

#endif
