#include "particles.h"
#include "../components/particles/particle_ic.h"

namespace cosmo
{

void ParticleSim::init()
{
  _timer["init"].start();

  // initialize base class
  simInit();

  iodata->log("Running 'particles' type simulation.");
  particles = new Particles();

  _timer["init"].stop();
}

void ParticleSim::setICs()
{
  if(_config("ic_type", "") == "vectorpert")
  {
    particle_ic_set_vectorpert(bssnSim, particles, iodata);
  }
  else if(_config("ic_type", "") == "sinusoid")
  {
    particle_ic_set_sinusoid(bssnSim, particles, iodata);
  }
  else if(_config("ic_type", "") == "sinusoid_to_compare")
  {
    particle_ic_set_sinusoid_to_compare(bssnSim, particles, iodata);
  }
  else
  {
    particle_ic_set_random(bssnSim, particles, fourier, iodata);
  }
}

void ParticleSim::initParticleStep()
{
  _timer["RK_steps"].start();
    bssnSim->stepInit();
    particles->stepInit(bssnSim->fields);
    bssnSim->clearSrc();
    particles->addParticlesToBSSNSrc(bssnSim);
  _timer["RK_steps"].stop();
}

void ParticleSim::outputParticleStep()
{
  _timer["output"].start();
    prepBSSNOutput();
    io_bssn_fields_snapshot(iodata, step, bssnSim->fields);
    io_bssn_fields_powerdump(iodata, step, bssnSim->fields, fourier);
    io_bssn_dump_statistics(iodata, step, bssnSim->fields, bssnSim->frw);
    io_bssn_constraint_violation(iodata, step, bssnSim);
    io_print_particles(iodata, step, particles);
    if(step == 0)
    {
      outputStateInformation();
    }
  _timer["output"].stop();
}

void ParticleSim::runParticleStep()
{
  _timer["RK_steps"].start();
    // First RK step
    bssnSim->RKEvolve();
    particles->RK1Step(bssnSim->fields);
    bssnSim->K1Finalize();
    particles->regSwap_c_a();

    // Second RK step source
    bssnSim->clearSrc();
    particles->addParticlesToBSSNSrc(bssnSim);
    // Second RK step
    bssnSim->RKEvolve();
    particles->RK2Step(bssnSim->fields);
    bssnSim->K2Finalize();
    particles->regSwap_c_a();

    // Third RK step source
    bssnSim->clearSrc();
    particles->addParticlesToBSSNSrc(bssnSim);
    // Third RK step
    bssnSim->RKEvolve();
    particles->RK3Step(bssnSim->fields);
    bssnSim->K3Finalize();
    particles->regSwap_c_a();

    // Fourth RK step source
    bssnSim->clearSrc();
    particles->addParticlesToBSSNSrc(bssnSim);
    // Fourth RK step
    bssnSim->RKEvolve();
    particles->RK4Step(bssnSim->fields);
    bssnSim->K4Finalize();
    particles->stepTerm();
    
    // "current" data should be in the _p array.
  _timer["RK_steps"].stop();
}

void ParticleSim::runStep()
{
  runCommonStepTasks();

  initParticleStep();
  outputParticleStep();
  runParticleStep();
}

} /* namespace cosmo */
