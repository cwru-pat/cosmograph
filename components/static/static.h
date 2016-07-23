#ifndef COSMO_STATIC
#define COSMO_STATIC

#include "../../cosmo_types.h"
#include "../../cosmo_macros.h"
#include "../../utils/Array.h"
#include "../../utils/FRW.h"

namespace cosmo
{

/** Static matter class **/
class Static
{
  /* Fluid field */
  // just a density variable
  GEN1_ARRAY_CREATE(DIFFD);

public:
  map_t fields;

  Static();
  ~Static();

  void addBSSNSrc(map_t & bssn_fields, FRW<real_t> *frw);

  void init();
};

}

#endif
