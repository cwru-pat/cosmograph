#ifndef COSMO_PARTICLE_SIM_H
#define COSMO_PARTICLE_SIM_H

#include "sim.h"
#include "../components/particles/particles.h"


namespace cosmo
{

/**
 * derived class based on CosmoSim class (sim.h)
 */
class ParticleSim : public CosmoSim
{
protected:
  Particles * particles;

public:
  ParticleSim() {}
  ~ParticleSim()
  {
    std::cout << "Cleaning up...";

    delete particles;
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
  void initParticleStep();
  void outputParticleStep();
  void runParticleStep();
  void runStep();
};

} /* namespace cosmo */

#endif
