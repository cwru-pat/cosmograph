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

#include "cosmo_macros.h"

namespace cosmo
{

typedef float real_t;
typedef long int idx_t;

} /* namespace cosmo */

#include "utils/math.h"
#include "utils/io.h"
#include "utils/ConfigParser.h"
#include "utils/Timer.h"
#include "Wave.h"

#include "bssn.h"
#include "bssn_data.h"
#include "hydro.h"

#endif
