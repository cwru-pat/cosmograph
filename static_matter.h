#ifndef COSMO_STATIC
#define COSMO_STATIC

#include "cosmo.h"
#include "globals.h"

namespace cosmo
{

/** Static matter class **/
class Static
{
  /* Fluid fields */
  GEN1_ARRAY_CREATE(r);

public:
  std::map <std::string, real_t *> fields;

  Static();
  ~Static();

  void addBSSNSrc(std::map <std::string, real_t *> & bssn_fields);

  void init();
};

}

#endif
