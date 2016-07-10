#ifndef COSMO_UTILS_MATH_H
#define COSMO_UTILS_MATH_H

#include "../cosmo_macros.h"
#include "../cosmo_types.h"
#include "../cosmo_includes.h"
#include "../cosmo_globals.h"

namespace cosmo
{

/**
 * @brief Kreiss-Oliger dissipation for 2nd-order stencils
 * @details
 * Works only for 2nd-order accurate (Odx2) derivatives
 * accuracy a = 2r - 2
 * r = (a + 2)/2
 * for a = 2, r = 2
 * Q = (-1)^r * (dx)^(2r-1) D_+^r D_- ^r / 2^(2r)
 *   = dx^3 / 2^6 D_+^2 D_-^2
 *   = dx^3 / 64 * [stencil: 1, -4, 6, -4, 1]
 * TODO: higher order?
 * 
 * @param i gridpoint in x-dir
 * @param j gridpoint in y-dir
 * @param k gridpoint in z-dir
 * @param field releavnt field
 * @return dissipation factor
 */
inline real_t KO_dissipation_Q(idx_t i, idx_t j, idx_t k, arr_t & field)
{
  #if STENCIL_ORDER == 2
    real_t stencil = (
        1.0*field[INDEX(i-2,j,k)] + 1.0*field[INDEX(i,j-2,k)] + 1.0*field[INDEX(i,j,k-2)]
      - 4.0*field[INDEX(i-1,j,k)] - 4.0*field[INDEX(i,j-1,k)] - 4.0*field[INDEX(i,j,k-1)]
      + 6.0*field[INDEX(i  ,j,k)] + 6.0*field[INDEX(i,j  ,k)] + 6.0*field[INDEX(i,j,k  )]
      - 4.0*field[INDEX(i+1,j,k)] - 4.0*field[INDEX(i,j+1,k)] - 4.0*field[INDEX(i,j,k+1)]
      + 1.0*field[INDEX(i+2,j,k)] + 1.0*field[INDEX(i,j+2,k)] + 1.0*field[INDEX(i,j,k+2)]
    )/dx/dx;
    real_t dissipation = KO_ETA*pow(dx, 3.0)/64.0*stencil;
    return dissipation;
  #else
    return 0.0;
  #endif

  // for a = 6, r = 4
  // Q = (dx)^(7) D_+^4 D_- ^4 / 2^(8)
  //   = (dx)^(7)/256 * [stencil: ]
  // ...
}

inline real_t derivative_Odx2(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field)
{
  switch (d) {
    case 1:
      return ((
        - 1.0/2.0*field[INDEX(i-1,j,k)]
        + 1.0/2.0*field[INDEX(i+1,j,k)]
      )/dx);
      break;
    case 2:
      return ((
        - 1.0/2.0*field[INDEX(i,j-1,k)]
        + 1.0/2.0*field[INDEX(i,j+1,k)]
      )/dx);
      break;
    case 3:
      return ((
        - 1.0/2.0*field[INDEX(i,j,k-1)]
        + 1.0/2.0*field[INDEX(i,j,k+1)]
      )/dx);
      break;
  }

  /* XXX */
  return 0;
}

inline real_t derivative_Odx4(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field)
{
  switch (d) {
    case 1:
      return (
        + 1.0/12.0*field[INDEX(i-2,j,k)]
        - 2.0/3.0*field[INDEX(i-1,j,k)]
        + 2.0/3.0*field[INDEX(i+1,j,k)]
        - 1.0/12.0*field[INDEX(i+2,j,k)]
      )/dx;
      break;
    case 2:
      return (
        + 1.0/12.0*field[INDEX(i,j-2,k)]
        - 2.0/3.0*field[INDEX(i,j-1,k)]
        + 2.0/3.0*field[INDEX(i,j+1,k)]
        - 1.0/12.0*field[INDEX(i,j+2,k)]
      )/dx;
      break;
    case 3:
      return (
        + 1.0/12.0*field[INDEX(i,j,k-2)]
        - 2.0/3.0*field[INDEX(i,j,k-1)]
        + 2.0/3.0*field[INDEX(i,j,k+1)]
        - 1.0/12.0*field[INDEX(i,j,k+2)]
      )/dx;
      break;
  }

  /* XXX */
  return 0;
}

inline real_t derivative_Odx6(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field)
{
  switch (d) {
    case 1:
      return (
        - 1.0/60.0*field[INDEX(i-3,j,k)]
        + 3.0/20.0*field[INDEX(i-2,j,k)]
        - 3.0/4.0*field[INDEX(i-1,j,k)]
        + 3.0/4.0*field[INDEX(i+1,j,k)]
        - 3.0/20.0*field[INDEX(i+2,j,k)]
        + 1.0/60.0*field[INDEX(i+3,j,k)]
      )/dx;
      break;
    case 2:
      return (
        - 1.0/60.0*field[INDEX(i,j-3,k)]
        + 3.0/20.0*field[INDEX(i,j-2,k)]
        - 3.0/4.0*field[INDEX(i,j-1,k)]
        + 3.0/4.0*field[INDEX(i,j+1,k)]
        - 3.0/20.0*field[INDEX(i,j+2,k)]
        + 1.0/60.0*field[INDEX(i,j+3,k)]
      )/dx;
      break;
    case 3:
      return (
        - 1.0/60.0*field[INDEX(i,j,k-3)]
        + 3.0/20.0*field[INDEX(i,j,k-2)]
        - 3.0/4.0*field[INDEX(i,j,k-1)]
        + 3.0/4.0*field[INDEX(i,j,k+1)]
        - 3.0/20.0*field[INDEX(i,j,k+2)]
        + 1.0/60.0*field[INDEX(i,j,k+3)]
      )/dx;
      break;
  }

  /* XXX */
  return 0;
}

inline real_t derivative_Odx8(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field)
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

inline real_t mixed_derivative_stencil_Odx2(idx_t i, idx_t j, idx_t k, int d1, int d2, arr_t & field)
{
  if( (d1 == 1 && d2 == 2) || (d1 == 2 && d2 == 1) ) {
    return (
      - field[INDEX(i+1,j-1,k)] + field[INDEX(i+1,j+1,k)]
      + field[INDEX(i-1,j-1,k)] - field[INDEX(i-1,j+1,k)]
    )/4.0/dx/dx;
  }

  if( (d1 == 1 && d2 == 3) || (d1 == 3 && d2 == 1) ) {
    return (
      - field[INDEX(i+1,j,k-1)] + field[INDEX(i+1,j,k+1)]
      + field[INDEX(i-1,j,k-1)] - field[INDEX(i-1,j,k+1)]
    )/4.0/dx/dx;
  }

  if( (d1 == 3 && d2 == 2) || (d1 == 2 && d2 == 3) ) {
    return (
      - field[INDEX(i,j+1,k-1)] + field[INDEX(i,j+1,k+1)]
      + field[INDEX(i,j-1,k-1)] - field[INDEX(i,j-1,k+1)]
    )/4.0/dx/dx;
  }

  /* XXX */
  return 0;
}

inline real_t mixed_derivative_stencil_Odx4(idx_t i, idx_t j, idx_t k, int d1, int d2, arr_t & field)
{
  if( (d1 == 1 && d2 == 2) || (d1 == 2 && d2 == 1) ) {
    return (
      (
        - field[INDEX(i+1,j-1,k)] + field[INDEX(i+1,j+1,k)]
        + field[INDEX(i-1,j-1,k)] - field[INDEX(i-1,j+1,k)]
      ) - 1.0/16.0*(
        - field[INDEX(i+2,j-2,k)] + field[INDEX(i+2,j+2,k)]
        + field[INDEX(i-2,j-2,k)] - field[INDEX(i-2,j+2,k)]
      )
    )/3.0/dx/dx;
  }

  if( (d1 == 1 && d2 == 3) || (d1 == 3 && d2 == 1) ) {
    return (
      (
        - field[INDEX(i+1,j,k-1)] + field[INDEX(i+1,j,k+1)]
        + field[INDEX(i-1,j,k-1)] - field[INDEX(i-1,j,k+1)]
      ) - 1.0/16.0*(
        - field[INDEX(i+2,j,k-2)] + field[INDEX(i+2,j,k+2)]
        + field[INDEX(i-2,j,k-2)] - field[INDEX(i-2,j,k+2)]
      )
    )/3.0/dx/dx;
  }

  if( (d1 == 3 && d2 == 2) || (d1 == 2 && d2 == 3) ) {
    return (
      (
        - field[INDEX(i,j+1,k-1)] + field[INDEX(i,j+1,k+1)]
        + field[INDEX(i,j-1,k-1)] - field[INDEX(i,j-1,k+1)]
      ) - 1.0/16.0*(
        - field[INDEX(i,j+2,k-2)] + field[INDEX(i,j+2,k+2)]
        + field[INDEX(i,j-2,k-2)] - field[INDEX(i,j-2,k+2)]
      )
    )/3.0/dx/dx;
  }

  /* XXX */
  return 0;
}

inline real_t mixed_derivative_stencil_Odx6(idx_t i, idx_t j, idx_t k, int d1, int d2, arr_t & field)
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

inline real_t mixed_derivative_stencil_Odx8(idx_t i, idx_t j, idx_t k, int d1,
 int d2, arr_t & field)
{
  if( (d1 == 1 && d2 == 2) || (d1 == 2 && d2 == 1) ) {
    return (
      2.0/5.0*(
        - field[INDEX(i+1,j-1,k)] + field[INDEX(i+1,j+1,k)]
        + field[INDEX(i-1,j-1,k)] - field[INDEX(i-1,j+1,k)]
      ) - 1.0/20.0*(
        - field[INDEX(i+2,j-2,k)] + field[INDEX(i+2,j+2,k)]
        + field[INDEX(i-2,j-2,k)] - field[INDEX(i-2,j+2,k)]
      ) + 2.0/315.0*(
        - field[INDEX(i+3,j-3,k)] + field[INDEX(i+3,j+3,k)]
        + field[INDEX(i-3,j-3,k)] - field[INDEX(i-3,j+3,k)]
      ) - 1.0/2240.0*(
        - field[INDEX(i+4,j-4,k)] + field[INDEX(i+4,j+4,k)]
        + field[INDEX(i-4,j-4,k)] - field[INDEX(i-4,j+4,k)]
      )
    )/dx/dx;
  }

  if( (d1 == 1 && d2 == 3) || (d1 == 3 && d2 == 1) ) {
    return (
      2.0/5.0*(
        - field[INDEX(i+1,j,k-1)] + field[INDEX(i+1,j,k+1)]
        + field[INDEX(i-1,j,k-1)] - field[INDEX(i-1,j,k+1)]
      ) - 1.0/20.0*(
        - field[INDEX(i+2,j,k-2)] + field[INDEX(i+2,j,k+2)]
        + field[INDEX(i-2,j,k-2)] - field[INDEX(i-2,j,k+2)]
      ) + 2.0/315.0*(
        - field[INDEX(i+3,j,k-3)] + field[INDEX(i+3,j,k+3)]
        + field[INDEX(i-3,j,k-3)] - field[INDEX(i-3,j,k+3)]
      ) - 1.0/2240.0*(
        - field[INDEX(i+4,j,k-4)] + field[INDEX(i+4,j,k+4)]
        + field[INDEX(i-4,j,k-4)] - field[INDEX(i-4,j,k+4)]
      )
    )/dx/dx;
  }

  if( (d1 == 3 && d2 == 2) || (d1 == 2 && d2 == 3) ) {
    return (
      2.0/5.0*(
        - field[INDEX(i,j+1,k-1)] + field[INDEX(i,j+1,k+1)]
        + field[INDEX(i,j-1,k-1)] - field[INDEX(i,j-1,k+1)]
      ) - 1.0/20.0*(
        - field[INDEX(i,j+2,k-2)] + field[INDEX(i,j+2,k+2)]
        + field[INDEX(i,j-2,k-2)] - field[INDEX(i,j-2,k+2)]
      ) + 2.0/315.0*(
        - field[INDEX(i,j+3,k-3)] + field[INDEX(i,j+3,k+3)]
        + field[INDEX(i,j-3,k-3)] - field[INDEX(i,j-3,k+3)]
      ) - 1.0/2240.0*(
        - field[INDEX(i,j+4,k-4)] + field[INDEX(i,j+4,k+4)]
        + field[INDEX(i,j-4,k-4)] - field[INDEX(i,j-4,k+4)]
      )
    )/dx/dx;
  }

  /* XXX */
  return 0;
}

inline real_t double_derivative_stencil_Odx2(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field)
{
  switch (d) {
    case 1:
      return (
          field[INDEX(i-1,j,k)]
          - 2.0*field[INDEX(i-0,j,k)]
          + field[INDEX(i+1,j,k)]
        )/dx/dx;
      break;
    case 2:
      return (
          field[INDEX(i,j-1,k)]
          - 2.0*field[INDEX(i,j-0,k)]
          + field[INDEX(i,j+1,k)]
        )/dx/dx;
      break;
    case 3:
      return (
          field[INDEX(i,j,k-1)]
          - 2.0*field[INDEX(i,j,k-0)]
          + field[INDEX(i,j,k+1)]
        )/dx/dx;
      break;
  }

  /* XXX */
  return 0;
}

inline real_t double_derivative_stencil_Odx4(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field)
{
  switch (d) {
    case 1:
      return (
          - 1.0/12.0*field[INDEX(i-2,j,k)]
          + 4.0/3.0*field[INDEX(i-1,j,k)]
          - 5.0/2.0*field[INDEX(i-0,j,k)]
          + 4.0/3.0*field[INDEX(i+1,j,k)]
          - 1.0/12.0*field[INDEX(i+2,j,k)]
        )/dx/dx;
      break;
    case 2:
      return (
          - 1.0/12.0*field[INDEX(i,j-2,k)]
          + 4.0/3.0*field[INDEX(i,j-1,k)]
          - 5.0/2.0*field[INDEX(i,j-0,k)]
          + 4.0/3.0*field[INDEX(i,j+1,k)]
          - 1.0/12.0*field[INDEX(i,j+2,k)]
        )/dx/dx;
      break;
    case 3:
      return (
          - 1.0/12.0*field[INDEX(i,j,k-2)]
          + 4.0/3.0*field[INDEX(i,j,k-1)]
          - 5.0/2.0*field[INDEX(i,j,k-0)]
          + 4.0/3.0*field[INDEX(i,j,k+1)]
          - 1.0/12.0*field[INDEX(i,j,k+2)]
        )/dx/dx;
      break;
  }

  /* XXX */
  return 0;
}

inline real_t double_derivative_stencil_Odx6(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field)
{
  switch (d) {
    case 1:
      return (
          1.0/90.0*field[INDEX(i-3,j,k)]
          - 3.0/20.0*field[INDEX(i-2,j,k)]
          + 3.0/2.0*field[INDEX(i-1,j,k)]
          - 49.0/18.0*field[INDEX(i-0,j,k)]
          + 3.0/2.0*field[INDEX(i+1,j,k)]
          - 3.0/20.0*field[INDEX(i+2,j,k)]
          + 1.0/90.0*field[INDEX(i+3,j,k)]
        )/dx/dx;
      break;
    case 2:
      return (
          1.0/90.0*field[INDEX(i,j-3,k)]
          - 3.0/20.0*field[INDEX(i,j-2,k)]
          + 3.0/2.0*field[INDEX(i,j-1,k)]
          - 49.0/18.0*field[INDEX(i,j-0,k)]
          + 3.0/2.0*field[INDEX(i,j+1,k)]
          - 3.0/20.0*field[INDEX(i,j+2,k)]
          + 1.0/90.0*field[INDEX(i,j+3,k)]
        )/dx/dx;
      break;
    case 3:
      return (
          1.0/90.0*field[INDEX(i,j,k-3)]
          - 3.0/20.0*field[INDEX(i,j,k-2)]
          + 3.0/2.0*field[INDEX(i,j,k-1)]
          - 49.0/18.0*field[INDEX(i,j,k-0)]
          + 3.0/2.0*field[INDEX(i,j,k+1)]
          - 3.0/20.0*field[INDEX(i,j,k+2)]
          + 1.0/90.0*field[INDEX(i,j,k+3)]
        )/dx/dx;
      break;
  }

  /* XXX */
  return 0;
}

inline real_t double_derivative_stencil_Odx8(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field)
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

/**
 * @brief Compute a derivative using a stencil order defined by a
 * preprocessor directive
 * 
 * @param i x-index
 * @param j x-index
 * @param k x-index
 * @param d direction of derivative
 * @param field field to differentiate
 * @return derivative
 */
inline real_t derivative(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field)
{
  return STENCIL_ORDER_FUNCTION(derivative_Odx)(i, j, k, d, field);
}

/**
 * @brief Compute a mixed derivative using a stencil order defined by a
 * preprocessor directive
 * 
 * @param i x-index
 * @param j x-index
 * @param k x-index
 * @param d1 direction of derivative in one direction
 * @param d2 direction of derivative in another direction
 * @param field field to differentiate
 * @return derivative
 */
inline real_t mixed_derivative_stencil(idx_t i, idx_t j, idx_t k, int d1, int d2, arr_t & field)
{
  return STENCIL_ORDER_FUNCTION(mixed_derivative_stencil_Odx)(i, j, k, d1, d2, field);
}

/**
 * @brief Compute a second-order derivative using a stencil order defined by a
 * preprocessor directive
 * 
 * @param i x-index
 * @param j x-index
 * @param k x-index
 * @param d direction of 2nd order derivative to compute
 * @param field field to differentiate
 * @return derivative
 */
inline real_t double_derivative_stencil(idx_t i, idx_t j, idx_t k, int d,
    arr_t & field)
{
  return STENCIL_ORDER_FUNCTION(double_derivative_stencil_Odx)(i, j, k, d, field);
}

/**
 * @brief A more generic function for 2nd derivs; calls
 * either @double_derivative_stencil or @mixed_derivative_stencil
 * 
 * @param i x-index
 * @param j x-index
 * @param k x-index
 * @param d1 direction of first derivative
 * @param d2 direction of second derivative
 * @param field field to differentiate
 * @return derivative
 */
inline real_t double_derivative(idx_t i, idx_t j, idx_t k, int d1, int d2,
    arr_t & field)
{
  if(d1 == d2) {
    return double_derivative_stencil(i, j, k, d1, field);
  } else {
    return mixed_derivative_stencil(i, j, k, d1, d2, field);
  }

  /* XXX */
  return 0;
}

/**
 * @brief Computes the laplacian of a field
 * @details sums up double_derivatives
 * 
 * @brief A more generic function for 2nd derivs; calls
 * either @double_derivative_stencil or @mixed_derivative_stencil
 * 
 * @param i x-index
 * @param j x-index
 * @param k x-index
 * @return laplacian
 */
inline real_t laplacian(idx_t i, idx_t j, idx_t k, arr_t & field)
{
  return (
    double_derivative(i, j, k, 1, 1, field)
    + double_derivative(i, j, k, 2, 2, field)
    + double_derivative(i, j, k, 3, 3, field)
  );
}

/**
 * @brief Compute the average of a field
 * 
 * @param field field to average
 * @return average
 */
inline real_t average(arr_t & field)
{
  // note this may have poor precision for large datasets
  real_t sum = 0.0; 
  idx_t i=0, j=0, k=0;
  
  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:sum)
  LOOP3(i, j, k)
  {
    sum += field[NP_INDEX(i,j,k)];
  }
  return sum/POINTS;
}

/**
 * @brief Compute the average volume of a simulation
 * 
 * @param DIFFphi difference variable conformal factor field
 * @param phi_FRW reference FRW conformal factor
 * 
 * @return [description]
 */
inline real_t volume_average(arr_t & DIFFphi, real_t phi_FRW)
{
  real_t sum = 0.0; 
  idx_t i=0, j=0, k=0;

  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:sum)
  LOOP3(i, j, k)
  {
    sum += exp(6.0*(DIFFphi[NP_INDEX(i,j,k)] + phi_FRW));
  }
  return sum/POINTS;
}

/**
 * @brief Compute the volume-weighted average of a field
 * 
 * @param field field to average
 * @param DIFFphi difference variable conformal factor field
 * @param phi_FRW reference FRW conformal factor
 * @return volume-weighted average
 */
inline real_t conformal_average(arr_t & field, arr_t & DIFFphi, real_t phi_FRW)
{
  real_t sum = 0.0; 
  idx_t i=0, j=0, k=0;

  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:sum)
  LOOP3(i, j, k)
  {
    sum += exp(6.0*(DIFFphi[NP_INDEX(i,j,k)] + phi_FRW))*field[NP_INDEX(i,j,k)];
  }
  real_t vol = volume_average(DIFFphi, phi_FRW);
  return sum/POINTS/vol;
}

/**
 * @brief Compute the standard deviation of a field
 * 
 * @param field field to analyze
 * @param avg pre-computed average of a field
 * @return standard deviation
 */
inline real_t standard_deviation(arr_t & field, real_t avg)
{
  // note this may have poor precision for large datasets
  idx_t i=0, j=0, k=0;
  real_t sum = 0.0;
  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:sum)
  LOOP3(i, j, k)
  {
    sum += pw2(avg - field[NP_INDEX(i,j,k)]);
  }
  return sqrt(sum/(POINTS-1));
}

/**
 * @brief Compute the standard deviation of a field
 * 
 * @param field field to analyze
 * @return standard deviation
 */
inline real_t standard_deviation(arr_t & field)
{
  real_t avg = average(field);
  return standard_deviation(field, avg);
}

/**
 * @brief Compute the volume-weighted standard deviation of a field
 * 
 * @param field field to analyze
 * @param DIFFphi difference variable conformal factor field
 * @param phi_FRW reference FRW conformal factor
 * @param avg pre-computed volume-weighted average
 * @return volume-weighted standard deviation
 */
inline real_t conformal_standard_deviation(arr_t & field, arr_t & DIFFphi, real_t phi_FRW, real_t avg)
{
  real_t sum = 0.0; 
  idx_t i=0, j=0, k=0;

  #pragma omp parallel for default(shared) private(i, j, k) reduction(+:sum)
  LOOP3(i, j, k)
  {
    sum += exp(6.0*(DIFFphi[NP_INDEX(i,j,k)] + phi_FRW))*pw2(avg - field[NP_INDEX(i,j,k)]);
  }
  real_t vol = volume_average(DIFFphi, phi_FRW);
  return sqrt(sum/(POINTS-1)/vol);
}

/**
 * @brief Compute the volume-weighted standard deviation of a field
 * 
 * @param field field to analyze
 * @param DIFFphi difference variable conformal factor field
 * @param phi_FRW reference FRW conformal factor
 * @return volume-weighted standard deviation
 */
inline real_t conformal_standard_deviation(arr_t & field, arr_t & DIFFphi, real_t phi_FRW)
{
  real_t avg = conformal_average(field, DIFFphi, phi_FRW);
  return conformal_standard_deviation(field, DIFFphi, phi_FRW, avg);
}

/**
 * @brief Compute the maximum value of a field
 */
inline real_t max(arr_t & field)
{
  idx_t i=0, j=0, k=0;
  real_t max_val = field[0];
  LOOP3(i, j, k)
  {
    if(field[INDEX(i,j,k)] > max_val) {
      max_val = field[INDEX(i,j,k)];
    }
  }
  return max_val;
}

/**
 * @brief Compute the minimum value of a field
 */
inline real_t min(arr_t & field)
{
  idx_t i=0, j=0, k=0;
  real_t min_val = field[0];
  LOOP3(i, j, k)
  {
    if(field[INDEX(i,j,k)] < min_val) {
      min_val = field[INDEX(i,j,k)];
    }
  }
  return min_val;
}

/**
 * @brief compute the number of NAN values in a field
 * @details may have portability issues
 */
inline idx_t numNaNs(arr_t & field)
{
  idx_t i=0, j=0, k=0;
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
