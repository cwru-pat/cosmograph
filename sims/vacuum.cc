#include "vacuum.h"

namespace cosmo
{

void VacuumSim::init()
{
  _timer["init"].start();

  // initialize base class
  simInit();

  iodata->log("Running 'vacuum' type simulation.");
  // TODO: Set vacuum ICs (eg, AwA test)
  setVacuumICs();

  _timer["init"].stop();
}

/**
 * @brief      Set vacuum initial conditions
 *
 * @param[in]  map to BSSN fields
 * @param      initialized IOData
 */
void VacuumSim::setVacuumICs()
{
  // meh
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
