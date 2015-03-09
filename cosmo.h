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
#include "utils/io.h"
#include "utils/math.h"

#include "bssn.h"
#include "bssn_data.h"
#include "hydro.h"
#include "hydro_data.h"
#include "lambda.h"
#include "ICs.h"
#include "frw.h"

#endif
