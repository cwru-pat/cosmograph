#ifndef COSMO_H
#define COSMO_H

/* C includes */
#include <zlib.h>
#include <omp.h>

/* standard C++ includes */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <map>
#include <fstream>
#include <algorithm>
#include <string>
#include <sstream>
#include <random>


#include "cosmo_macros.h"

namespace cosmo
{

// changing this affects FFTs:
typedef double real_t;
// see http://www.fftw.org/doc/Precision.html

typedef long int idx_t;

} /* namespace cosmo */

#include "utils/ConfigParser.h"
#include "utils/Timer.h"
#include "utils/Fourier.h"
#include "utils/math.h"
#include "utils/reference_frw.h"

#include "cosmotrace/raytrace.h"

#include "io.h"
#include "ICs.h"
#include "bssn.h"
#include "bssn_data.h"
#include "static.h"
#include "particles.h"

#endif
