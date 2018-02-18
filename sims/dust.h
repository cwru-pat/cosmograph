#ifndef COSMO_DUST_SIM_H
#define COSMO_DUST_SIM_H

#include "sim.h"
#include "../components/static/static.h"


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
  ~DustSim()
  {
    std::cout << "Cleaning up...";
    delete iodata;
    delete bssnSim;
    delete fourier;
    if(use_bardeen)
    {
      delete bardeen;
    }
    std::cout << "done.\n";
    std::cout << std::flush;
  }

  void init();
  void setICs();
  void initDustStep();
  void outputDustStep();
  void runDustStep();
  void runStep();
};

} /* namespace cosmo */

#endif
