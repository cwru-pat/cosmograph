#ifndef COSMO_PARTICLE_SIM_H
#define COSMO_PARTICLE_SIM_H

#include "sim.h"
#include "../particles/particles.h"


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
  ParticleSim();
  ~ParticleSim(){}

  void init();
  void setICs();
  void initParticleStep();
  void outputParticleStep();
  void runParticleStep();
  void runStep();
};

} /* namespace cosmo */

#endif
