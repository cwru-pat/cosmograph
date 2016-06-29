#ifndef COSMO_TYPES
#define COSMO_TYPES

#include "utils/Array.h"
#include "utils/RK4Register.h"
#include <string>
#include <map>

namespace cosmo
{

// changing this affects FFTs:
typedef double real_t;
// see http://www.fftw.org/doc/Precision.html

typedef long int idx_t;

typedef CosmoArray<idx_t, real_t> arr_t;

typedef RK4Register<idx_t, real_t> register_t;

typedef std::map <std::string, arr_t *> map_t;

} /* namespace cosmo */

#endif
