#ifndef COSMO_DUST_SIM_H
#define COSMO_DUST_SIM_H

#include "sim.h"
#include "../static/static.h"


namespace cosmo
{

/**
 * derived class based on CosmoSim class (sim.h)
 */
class DustSim : public CosmoSim
{
protected:
  Static * staticSim;

public:
  DustSim();
  ~DustSim(){}

  void init();
  void setICs();
  void initDustStep();
  void outputDustStep();
  void runDustStep();
  void runStep();
};

} /* namespace cosmo */

#endif
