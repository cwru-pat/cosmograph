#ifndef COSMO_DUST_LAMBDA_SIM_H
#define COSMO_DUST_LAMBDA_SIM_H

#include "dust.h"

namespace cosmo
{

/**
 * derived class based on DustSim class (dust.h)
 */
class DustLambdaSim : public DustSim
{
public:
  DustLambdaSim();
  ~DustLambdaSim(){}

  void init();
  void setICs();
  // use DustSim implementation of runStep;
  // only ICs need to be different.
};

} /* namespace cosmo */

#endif
