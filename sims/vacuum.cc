#include "vacuum.h"
#include "../components/bssn/bssn_ic.h"

namespace cosmo
{

void VacuumSim::init()
{
  _timer["init"].start();

  // initialize base class
  simInit();

  iodata->log("Running 'vacuum' type simulation.");
  _timer["init"].stop();
}

/**
 * @brief      Set vacuum initial conditions
 *
 * @param[in]  map to BSSN fields
 * @param      initialized IOData
 */
void VacuumSim::setICs()
{
  _timer["ICs"].start();
  iodata->log("Setting initial conditions (ICs).");

  iodata->log("Setting initial conditions (ICs).");

  if(_config["ic_type"] == "stability")
  {
    bssn_ic_awa_stability(bssnSim);
  }
  else if(_config["ic_type"] == "linear_wave")
  {
    bssn_ic_awa_linear_wave(bssnSim);
  }
  else if(_config["ic_type"] == "linear_wave_desitter")
  {
    bssn_ic_awa_linear_wave_desitter(bssnSim);
  }
  else
  {
    iodata->log("IC type not recognized!");
    throw -1;
  }
  iodata->log("Finished setting ICs.");
  _timer["ICs"].stop();
}

void VacuumSim::initVacuumStep()
{
  _timer["RK_steps"].start();
    bssnSim->stepInit();
  _timer["RK_steps"].stop();
}

void VacuumSim::outputVacuumStep()
{
  _timer["output"].start();
    prepBSSNOutput();
    io_bssn_fields_snapshot(iodata, step, bssnSim->fields);
    io_bssn_fields_powerdump(iodata, step, bssnSim->fields, fourier);
    io_bssn_dump_statistics(iodata, step, bssnSim->fields, bssnSim->frw);
    io_bssn_constraint_violation(iodata, step, bssnSim);
  _timer["output"].stop();
}

void VacuumSim::runVacuumStep()
{
  _timer["RK_steps"].start();
    // Full RK step minus init()
    bssnSim->step();
  _timer["RK_steps"].stop();
}

void VacuumSim::runStep()
{
  runCommonStepTasks();

  initVacuumStep();
  outputVacuumStep();
  runVacuumStep();
}

} /* namespace cosmo */
