#ifndef COSMO_VACUUM_SIM_H
#define COSMO_VACUUM_SIM_H

#include "sim.h"

namespace cosmo
{

/**
 * derived class based on CosmoSim class (sim.h)
 */
class VacuumSim : public CosmoSim
{
public:
  VacuumSim(){}
  ~VacuumSim(){}

  void init();
  void setICs();
  void initVacuumStep();
  void outputVacuumStep();
  void runVacuumStep();
  void runStep();
};

} /* namespace cosmo */

#endif
