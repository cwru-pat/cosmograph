#ifndef COSMO_TYPES
#define COSMO_TYPES

#include "utils/Array.h"
#include "utils/RK4Register.h"
#include <string>
#include <map>

namespace cosmo
{

#if USE_LONG_DOUBLES
// changing this affects FFTs:
typedef long double real_t; /**< real type; changing this may require changes to HDF5 and FFTW functionality */
// see http://www.fftw.org/doc/Precision.html
#else
typedef double real_t;
#endif

typedef long int idx_t; /**< indexing type, must be long enough to support large arrays */

typedef CosmoArray<idx_t, real_t> arr_t; /**< base array type */

typedef RK4Register<idx_t, real_t> register_t; /**< RK4 group of registers (_a, _p, _c, _f types) */

typedef std::map <std::string, arr_t *> map_t; /**< Map type; maps strings to array references */

} /* namespace cosmo */

#endif
