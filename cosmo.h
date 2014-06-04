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
typedef float real_t;
typedef unsigned int idx_t;

#include "utils/Timer.h"

TimerManager _timer;

} /* namespace cosmo */

#endif
