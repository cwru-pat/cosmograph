#include "dust_lambda.h"

namespace cosmo
{

DustLambdaSim::DustLambdaSim()
{
  // Reference metric is needed for Lambda
  if(!USE_REFERENCE_FRW)
  {
    iodata->log("Error - not using reference metric! You must use it for lambda sims.");
    throw -1;
  }
}

void DustLambdaSim::init()
{
  _timer["init"].start();

  // initialize base class
  simInit();

  iodata->log("Running lambda + 'dust' type simulation.");
  staticSim = new Static();
  staticSim->init();

  iodata->log("Creating initial conditions.");
  setDustICs();
  addLambdaICs();

  _timer["init"].stop();
}

/**
 * @brief Add in a lambda component to a dust simulation.
 */
void DustLambdaSim::addLambdaICs()
{
  ICsData icd = cosmo_get_ICsData();
  auto & frw = bssnSim->frw;
  
  // Add in lambda
  real_t rho_matter = icd.rho_K_matter;
  real_t rho_lambda = icd.rho_K_lambda;
  real_t rho_FRW = rho_lambda + rho_matter;
  real_t K_frw = -sqrt(24.0*PI*rho_FRW);

  frw->set_K(K_frw);
  // dust already added; just need lambda
  frw->addFluid(rho_lambda, -1.0 /* w=-1 */);
}

} /* namespace cosmo */
