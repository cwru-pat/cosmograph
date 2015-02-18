#ifndef COSMO_ICS
#define COSMO_ICS

#include "cosmo.h"
#include "globals.h"

#include "ICs_data.h"

namespace cosmo
{

void set_BH_ICs(std::map <std::string, real_t *> & bssn_fields);

} // namespace cosmo

#endif
