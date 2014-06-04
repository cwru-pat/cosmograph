#ifndef COSMO_H
#define COSMO_H

/* C includes */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>

/* C++ includes */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <map>

namespace cosmo
{

#include "defines.h"

typedef float real_t;
typedef unsigned int idx_t;

#include "utils/math.h"

#include "utils/Timer.h"
TimerManager _timer;

#include "wave.h"


} /* namespace cosmo */

#endif
