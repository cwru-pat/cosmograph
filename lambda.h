#ifndef COSMO_LAMBDA
#define COSMO_LAMBDA

#include "cosmo.h"
#include "globals.h"

namespace cosmo
{

/** Lambda class **/
class Lambda
{
  real_t lambda; /* CC energy density*/ 

public:

  Lambda(real_t L);
  ~Lambda();

  void addBSSNSrc(std::map <std::string, real_t *> & bssn_fields);

};

}

#endif
