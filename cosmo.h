#ifndef COSMO_H
#define COSMO_H

/* C includes */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <zlib.h>

/* C++ includes */
#include <iostream>
#include <iomanip>
#include <map>
#include <fstream>
#include <algorithm>
#include <string>
#include <sstream>
#include <random>

#include <omp.h>

#include "cosmo_macros.h"

/* independent components */
#include "utils/periodicArray.h"
#include "utils/ConfigParser.h"
#include "utils/Timer.h"
#include "utils/reference_frw.h"

namespace cosmo
{

// changing real_t affects FFTs (see http://www.fftw.org/doc/Precision.html)
typedef double real_t;
typedef long int idx_t;
typedef periodicArray<idx_t, real_t> arr_t;

} /* namespace cosmo */

#include "utils/Fourier.h"
#include "utils/math.h"

#include "io.h"
#include "ICs.h"
#include "bssn.h"
#include "bssn_data.h"
#include "static.h"

#endif
