#include "dust.h"
#include "../static/static_ic.h"

namespace cosmo
{

DustSim::DustSim()
{
  // just check to make sure we can use this class.
  if(USE_HARMONIC_ALPHA) {
    iodata->log("Error - not using synchronous gauge! You must use it for dust sims.");
    iodata->log("Please change this setting in cosmo_macros.h and recompile.");
    throw -1;
  }
}

void DustSim::init()
{
  _timer["init"].start();

  // initialize base class
  simInit();

  iodata->log("Initializing 'dust' type simulation.");
  staticSim = new Static();
  staticSim->init();

  _timer["init"].stop();
}

/**
 * @brief      Set conformally flat initial conditions for a w=0
 *  fluid in synchronous gauge.
 */
void DustSim::setICs()
{
  _timer["ICs"].start();

  iodata->log("Setting dust initial conditions (ICs).");
  dust_ic_set_random(bssnSim, staticSim, fourier, iodata);
  iodata->log("Finished setting ICs.");
  
  _timer["ICs"].stop();
}

void DustSim::initDustStep()
{
  _timer["RK_steps"].start();
    bssnSim->stepInit();
    bssnSim->clearSrc();
    staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
  _timer["RK_steps"].stop();
}

void DustSim::outputDustStep()
{
  _timer["output"].start();
    prepBSSNOutput();
    io_bssn_fields_snapshot(iodata, step, bssnSim->fields);
    io_bssn_fields_powerdump(iodata, step, bssnSim->fields, fourier);
    io_bssn_dump_statistics(iodata, step, bssnSim->fields, bssnSim->frw);
    io_bssn_constraint_violation(iodata, step, bssnSim);
  _timer["output"].stop();
}

void DustSim::runDustStep()
{
  _timer["RK_steps"].start();
    // First RK step
    // source already set in initDustStep() (used for output)
    bssnSim->RKEvolve();
    bssnSim->K1Finalize();

    // Second RK step
    bssnSim->clearSrc();
    staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
    bssnSim->RKEvolve();
    bssnSim->K2Finalize();

    // Third RK step
    bssnSim->clearSrc();
    staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
    bssnSim->RKEvolve();
    bssnSim->K3Finalize();

    // Fourth RK step
    bssnSim->clearSrc();
    staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
    bssnSim->RKEvolve();
    bssnSim->K4Finalize();

    // "current" data should be in the _p array.
  _timer["RK_steps"].stop();
}

void DustSim::runStep()
{
  runCommonStepTasks();

  initDustStep();
  outputDustStep();
  runDustStep();
}

} /* namespace cosmo */
