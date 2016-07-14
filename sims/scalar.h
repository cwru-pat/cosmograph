#ifndef COSMO_SCALAR_SIM_H
#define COSMO_SCALAR_SIM_H

#include "sim.h"
#include "../scalar/scalar.h"
#include "../elliptic_solver/multigrid.h"

namespace cosmo
{

/**
 * derived class based on CosmoSim class (sim.h)
 */
class ScalarSim : public CosmoSim
{
protected:
  Scalar * scalarSim;

public:
  ScalarSim(){}
  ~ScalarSim(){}

  void init();
  void setScalarWaveICs();
  void setScalarMultigridICs();
  void initScalarStep();
  void outputScalarStep();
  void runScalarStep();
  void runStep();
  void setAnalyticScalarTestICs();
};

} /* namespace cosmo */

#endif
