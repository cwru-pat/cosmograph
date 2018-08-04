#ifndef COSMO_LAMBDA
#define COSMO_LAMBDA

#include "../../cosmo_types.h"
#include "../bssn/bssn.h"

namespace cosmo
{

/** Lambda class **/
class Lambda
{
  real_t lambda; /* CC energy density*/ 

public:

  Lambda();
  ~Lambda();

  void setLambda(real_t lambda_in);
  void addBSSNSource(BSSN *bssn);

};

}

#endif
