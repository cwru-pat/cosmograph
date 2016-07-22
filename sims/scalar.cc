#include "scalar.h"
#include "../scalar/scalar_ic.h"

namespace cosmo
{

void ScalarSim::init()
{
  _timer["init"].start();

  // initialize base class
  simInit();

  iodata->log("Running 'scalar' type simulation.");
  if(!USE_HARMONIC_ALPHA) {
    iodata->log("Warning - not using harmonic gauge! You may want to use it.");
  }

  scalarSim = new Scalar();

  _timer["init"].stop();
}

void ScalarSim::setICs()
{
  _timer["ICs"].start();
  iodata->log("Setting initial conditions (ICs).");

  if(_config["scalar_ic_type"] == "wave")
  {
    scalar_ic_set_wave(bssnSim, scalarSim);
  }
  else if(_config["scalar_ic_type"] == "Lambda")
  {
    scalar_ic_set_Lambda(bssnSim, scalarSim);
  }
  else if(_config["scalar_ic_type"] == "semianalytic_test")
  {
    scalar_ic_set_semianalytic_test(bssnSim, scalarSim);
  }
  else if(_config["scalar_ic_type"] == "multigrid")
  {
    scalar_ic_set_multigrid(bssnSim, scalarSim);
  }
  else
  {
    iodata->log("IC type not recognized!");
    throw -1;
  }

  iodata->log("Finished setting ICs.");
  _timer["ICs"].stop();
}

void ScalarSim::initScalarStep()
{
  _timer["RK_steps"].start();
    bssnSim->stepInit();
    scalarSim->stepInit();
    bssnSim->clearSrc();
    scalarSim->addBSSNSource(bssnSim);
  _timer["RK_steps"].stop();
}

void ScalarSim::outputScalarStep()
{
  _timer["output"].start();
    prepBSSNOutput();
    io_bssn_fields_snapshot(iodata, step, bssnSim->fields);
    io_bssn_fields_powerdump(iodata, step, bssnSim->fields, fourier);
    io_bssn_dump_statistics(iodata, step, bssnSim->fields, bssnSim->frw);
    io_bssn_constraint_violation(iodata, step, bssnSim);
    io_scalar_snapshot(iodata, step, scalarSim);
  _timer["output"].stop();
}

void ScalarSim::runScalarStep()
{
  idx_t i=0, j=0, k=0;
  BSSNData b_data;
  _timer["RK_steps"].start();

    // First RK step
    #pragma omp parallel for default(shared) private(i, j, k, b_data)
    LOOP3(i,j,k)
    {
      bssnSim->RKEvolvePt(i, j, k, &b_data);
      scalarSim->RKEvolvePt(&b_data);
    }
    bssnSim->K1Finalize();
    scalarSim->K1Finalize();

    // Second RK step
    bssnSim->clearSrc();
    scalarSim->addBSSNSource(bssnSim);

    #pragma omp parallel for default(shared) private(i, j, k, b_data)
    LOOP3(i,j,k)
    {
      bssnSim->RKEvolvePt(i, j, k, &b_data);
      scalarSim->RKEvolvePt(&b_data);
    }
    bssnSim->K2Finalize();
    scalarSim->K2Finalize();

    // Third RK step
    bssnSim->clearSrc();
    scalarSim->addBSSNSource(bssnSim);

    #pragma omp parallel for default(shared) private(i, j, k, b_data)
    LOOP3(i,j,k)
    {
      bssnSim->RKEvolvePt(i, j, k, &b_data);
      scalarSim->RKEvolvePt(&b_data);
    }
    bssnSim->K3Finalize();
    scalarSim->K3Finalize();

    // Fourth RK step
    bssnSim->clearSrc();
    scalarSim->addBSSNSource(bssnSim);
    #pragma omp parallel for default(shared) private(i, j, k, b_data)
    LOOP3(i,j,k)
    {
      bssnSim->RKEvolvePt(i, j, k, &b_data);
      scalarSim->RKEvolvePt(&b_data);
    }
    bssnSim->K4Finalize();
    scalarSim->K4Finalize();

    // "current" data should be in the _p array.
  _timer["RK_steps"].stop();
}

void ScalarSim::runStep()
{
  runCommonStepTasks();
  
  initScalarStep();
  outputScalarStep();
  runScalarStep();
}

} /* namespace cosmo */
