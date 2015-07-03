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

inline real_t derivative(idx_t i, idx_t j, idx_t k, int d,
    real_t *field)
{
  switch (d) {
    case 1:
      return (-1.0*field[INDEX(i+2,j,k)] + 8.0*field[INDEX(i+1,j,k)] - 8.0*field[INDEX(i-1,j,k)] + field[INDEX(i-2,j,k)])/12.0/dx;
      break;
    case 2:
      return (-1.0*field[INDEX(i,j+2,k)] + 8.0*field[INDEX(i,j+1,k)] - 8.0*field[INDEX(i,j-1,k)] + field[INDEX(i,j-2,k)])/12.0/dx;
      break;
    case 3:
      return (-1.0*field[INDEX(i,j,k+2)] + 8.0*field[INDEX(i,j,k+1)] - 8.0*field[INDEX(i,j,k-1)] + field[INDEX(i,j,k-2)])/12.0/dx;
      break;
  }

  /* XXX */
  return 0;
}

inline real_t mixed_derivative_stencil(idx_t i, idx_t j, idx_t k, int d1, int d2, real_t *field)
{
  if( (d1 == 1 && d2 == 2) || (d1 == 2 && d2 == 1) ) {
    return -(
       1.00*field[INDEX(i-2,j+2,k)] + -8.00*field[INDEX(i-1,j+2,k)] +  8.00*field[INDEX(i+1,j+2,k)] + -1.00*field[INDEX(i+2,j+2,k)] +
      -8.00*field[INDEX(i-2,j+1,k)] +  64.0*field[INDEX(i-1,j+1,k)] + -64.0*field[INDEX(i+1,j+1,k)] +  8.00*field[INDEX(i+2,j+1,k)] +
       8.00*field[INDEX(i-2,j-1,k)] + -64.0*field[INDEX(i-1,j-1,k)] +  64.0*field[INDEX(i+1,j-1,k)] + -8.00*field[INDEX(i+2,j-1,k)] +
      -1.00*field[INDEX(i-2,j-2,k)] +  8.00*field[INDEX(i-1,j-2,k)] + -8.00*field[INDEX(i+1,j-2,k)] +  1.00*field[INDEX(i+2,j-2,k)]
    )/12.0/12.0/dx/dx;
  }

  if( (d1 == 1 && d2 == 3) || (d1 == 3 && d2 == 1) ) {
    return -(
       1.00*field[INDEX(i-2,j,k+2)] + -8.00*field[INDEX(i-1,j,k+2)] +  8.00*field[INDEX(i+1,j,k+2)] + -1.00*field[INDEX(i+2,j,k+2)] +
      -8.00*field[INDEX(i-2,j,k+1)] +  64.0*field[INDEX(i-1,j,k+1)] + -64.0*field[INDEX(i+1,j,k+1)] +  8.00*field[INDEX(i+2,j,k+1)] +
       8.00*field[INDEX(i-2,j,k-1)] + -64.0*field[INDEX(i-1,j,k-1)] +  64.0*field[INDEX(i+1,j,k-1)] + -8.00*field[INDEX(i+2,j,k-1)] +
      -1.00*field[INDEX(i-2,j,k-2)] +  8.00*field[INDEX(i-1,j,k-2)] + -8.00*field[INDEX(i+1,j,k-2)] +  1.00*field[INDEX(i+2,j,k-2)]
    )/12.0/12.0/dx/dx;
  }

  if( (d1 == 3 && d2 == 2) || (d1 == 2 && d2 == 3) ) {
    return -(
       1.00*field[INDEX(i,j+2,k-2)] + -8.00*field[INDEX(i,j+2,k-1)] +  8.00*field[INDEX(i,j+2,k+1)] + -1.00*field[INDEX(i,j+2,k+2)] +
      -8.00*field[INDEX(i,j+1,k-2)] +  64.0*field[INDEX(i,j+1,k-1)] + -64.0*field[INDEX(i,j+1,k+1)] +  8.00*field[INDEX(i,j+1,k+2)] +
       8.00*field[INDEX(i,j-1,k-2)] + -64.0*field[INDEX(i,j-1,k-1)] +  64.0*field[INDEX(i,j-1,k+1)] + -8.00*field[INDEX(i,j-1,k+2)] +
      -1.00*field[INDEX(i,j-2,k-2)] +  8.00*field[INDEX(i,j-2,k-1)] + -8.00*field[INDEX(i,j-2,k+1)] +  1.00*field[INDEX(i,j-2,k+2)]
    )/144.0/dx/dx;
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
          1.0*field[INDEX(i+4,j,k)] +
          -16.0*field[INDEX(i+3,j,k)] +
          64.0*field[INDEX(i+2,j,k)] +
          16.0*field[INDEX(i+1,j,k)] +
          -130.0*field[INDEX(i+0,j,k)] +
          16.0*field[INDEX(i-1,j,k)] +
          64.0*field[INDEX(i-2,j,k)] +
          -16.0*field[INDEX(i-3,j,k)] +
          1.0*field[INDEX(i-4,j,k)]
        )/12.0/12.0/dx/dx;
      break;
    case 2:
      return (
          1.0*field[INDEX(i,j+4,k)] +
          -16.0*field[INDEX(i,j+3,k)] +
          64.0*field[INDEX(i,j+2,k)] +
          16.0*field[INDEX(i,j+1,k)] +
          -130.0*field[INDEX(i,j+0,k)] +
          16.0*field[INDEX(i,j-1,k)] +
          64.0*field[INDEX(i,j-2,k)] +
          -16.0*field[INDEX(i,j-3,k)] +
          1.0*field[INDEX(i,j-4,k)]
        )/12.0/12.0/dx/dx;
      break;
    case 3:
      return (
          1.0*field[INDEX(i,j,k+4)] +
          -16.0*field[INDEX(i,j,k+3)] +
          64.0*field[INDEX(i,j,k+2)] +
          16.0*field[INDEX(i,j,k+1)] +
          -130.0*field[INDEX(i,j,k+0)] +
          16.0*field[INDEX(i,j,k-1)] +
          64.0*field[INDEX(i,j,k-2)] +
          -16.0*field[INDEX(i,j,k-3)] +
          1.0*field[INDEX(i,j,k-4)]
        )/12.0/12.0/dx/dx;
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
