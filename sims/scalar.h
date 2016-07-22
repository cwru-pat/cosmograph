/** @file scalar.h
 * @brief Functionality for scalar simulations implemented via ScalarSim class.
 */

#ifndef COSMO_SCALAR_SIM_H
#define COSMO_SCALAR_SIM_H

#include "sim.h"
#include "../scalar/scalar.h"

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
  void setICs();
  void initScalarStep();
  void outputScalarStep();
  void runScalarStep();
  void runStep();
};

} /* namespace cosmo */

#endif
