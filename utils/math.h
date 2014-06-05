#ifndef COSMO_UTILS_MATH_H
#define COSMO_UTILS_MATH_H

#include "../cosmo.h"

namespace cosmo
{

inline real_t lap_stencil(idx_t i, idx_t j, idx_t k, real_t *field);

inline real_t derivative_stencil(idx_t i, idx_t j, idx_t k, int d, real_t *field);

}

#endif
