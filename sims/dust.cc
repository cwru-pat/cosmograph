#include "dust.h"
#include "../components/static/static_ic.h"
#include "../components/phase_space_sheet/sheets_ic.h"

namespace cosmo
{

DustSim::DustSim()
{
  // just check to make sure we can use this class.
  if(_config("lapse", "") != "" && _config("lapse", "") != "Static" && _config("lapse", "") != "ConformalFLRW")
  {
    iodata->log("Error - not using synchronous gauge! You must use it for dust sims.");
    iodata->log("Please change this setting in the config file and re-run.");
    throw -1;
  }

  take_ray_step = false;
  raysheet_flip_step = std::stoi(_config("raysheet_flip_step", "-1"));
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
  raySheet = new Sheet();
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
  else if(_config("ic_type", "") == "semianalytic")
  {
    dust_ic_set_semianalytic(bssnSim, staticSim, lambda, fourier, iodata);
  }
  else if(_config("ic_type", "") == "sinusoid_3d")
  {
    dust_ic_set_sinusoid_3d(bssnSim, staticSim, lambda, fourier, iodata);
  }
  else
  {
    iodata->log("Setting cosmological ICs.");
    dust_ic_set_random(bssnSim, staticSim, lambda, fourier, iodata);
  }

  // Set raytracing initial conditions
  sheets_ic_rays(bssnSim, raySheet, iodata);

  iodata->log("Finished setting ICs.");
  
  _timer["ICs"].stop();
}

void DustSim::initDustStep()
{
  _timer["RK_steps"].start();
    bssnSim->stepInit();
    bssnSim->clearSrc();
    if(take_ray_step) raySheet->stepInit();
    staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
    lambda->addBSSNSource(bssnSim);
  _timer["RK_steps"].stop();

  arr_t & DIFFr_a = *bssnSim->fields["DIFFr_a"];
  real_t rho_tot_avg = average(DIFFr_a);
  real_t rho_L = lambda->getLambda();
  real_t Omega_L_flip = std::stod(_config("raysheet_flip_omega_L", "0.0"));
  if(Omega_L_flip > 0.0 && rho_L/rho_tot_avg > Omega_L_flip)
  {
    take_ray_step = 1;
    num_steps = 2*step;
    iodata->log("\nFlipping sign of dt @ step = " + std::to_string(step) );
    iodata->log("--Omega_L was " + std::to_string(rho_L/rho_tot_avg) );
    iodata->log("--Setting final number of simulation steps to " + std::to_string(num_steps) );
    dt = -std::abs(dt);
    bssnSim->setDt(dt);
    raySheet->setDt(dt);
  }
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
    if(take_ray_step)
      io_raysheet_dump(iodata, step, raySheet, bssnSim, lambda);
  _timer["output"].stop();
}

void DustSim::runDustStep()
{
  _timer["RK_steps"].start();
    // First RK step
    // source already set in initDustStep() (used for output)
    bssnSim->RKEvolve();
    if(take_ray_step) raySheet->RKStep(bssnSim);
    bssnSim->K1Finalize();
    if(take_ray_step) raySheet->K1Finalize();

    // Second RK step
    bssnSim->clearSrc();
    staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
    lambda->addBSSNSource(bssnSim);
    bssnSim->RKEvolve();
    if(take_ray_step) raySheet->RKStep(bssnSim);
    bssnSim->K2Finalize();
    if(take_ray_step) raySheet->K2Finalize();

    // Third RK step
    bssnSim->clearSrc();
    staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
    lambda->addBSSNSource(bssnSim);
    bssnSim->RKEvolve();
    if(take_ray_step) raySheet->RKStep(bssnSim);
    bssnSim->K3Finalize();
    if(take_ray_step) raySheet->K3Finalize();

    // Fourth RK step
    bssnSim->clearSrc();
    staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
    lambda->addBSSNSource(bssnSim);
    bssnSim->RKEvolve();
    if(take_ray_step) raySheet->RKStep(bssnSim);
    bssnSim->K4Finalize();
    if(take_ray_step) raySheet->K4Finalize();

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
