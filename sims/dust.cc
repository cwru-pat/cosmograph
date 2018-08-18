#include "dust.h"
#include "../components/static/static_ic.h"

namespace cosmo
{

DustSim::DustSim()
{
  // just check to make sure we can use this class.
  if(_config("lapse", "") != "" && _config("lapse", "") != "static")
  {
    iodata->log("Error - not using synchronous gauge! You must use it for dust sims.");
    iodata->log("Please change this setting in the config file and re-run.");
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
  lambda = new Lambda();
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

  if(_config("ic_type", "") == "shell")
  {
    dust_ic_set_sphere(bssnSim, staticSim, iodata);
  }
  else if(_config("ic_type", "") == "sinusoid")
  {
    dust_ic_set_sinusoid(bssnSim, staticSim, lambda, fourier, iodata);
  }
  else
  {
    iodata->log("Creating gaussian random field.");
    dust_ic_set_random(bssnSim, staticSim, fourier, iodata);
  }
  iodata->log("Finished setting ICs.");
  
  _timer["ICs"].stop();
}

void DustSim::initDustStep()
{
  _timer["RK_steps"].start();
    bssnSim->stepInit();
    bssnSim->clearSrc();
    staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
    lambda->addBSSNSource(bssnSim);
  _timer["RK_steps"].stop();
}

void DustSim::outputDustStep()
{
  _timer["output"].start();
    prepBSSNOutput();
    if(use_bardeen)
      io_svt_violation(iodata, step, bardeen, t);

    io_bssn_fields_snapshot(iodata, step, bssnSim->fields);
    io_bssn_fields_powerdump(iodata, step, bssnSim->fields, fourier);
    io_bssn_dump_statistics(iodata, step, bssnSim->fields, bssnSim->frw);
    io_bssn_constraint_violation(iodata, step, bssnSim);
    if(step == 0)
    {
      outputStateInformation();
    }
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
    lambda->addBSSNSource(bssnSim);
    bssnSim->RKEvolve();
    bssnSim->K2Finalize();

    // Third RK step
    bssnSim->clearSrc();
    staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
    lambda->addBSSNSource(bssnSim);
    bssnSim->RKEvolve();
    bssnSim->K3Finalize();

    // Fourth RK step
    bssnSim->clearSrc();
    staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
    lambda->addBSSNSource(bssnSim);
    bssnSim->RKEvolve();
    bssnSim->K4Finalize();

    // "current" data should be in the _p array.
  _timer["RK_steps"].stop();
}

void DustSim::runStep()
{
  initDustStep();
  runCommonStepTasks();
  outputDustStep();
  runDustStep();
}

} /* namespace cosmo */
