#include "cosmo_includes.h"
#include "cosmo_types.h"
#include "globals.h"

#include "sims/sim.h"
#include "sims/dust.h"
#include "sims/dust_lambda.h"
#include "sims/particles.h"
#include "sims/scalar.h"
#include "sims/vacuum.h"

using namespace std;
using namespace cosmo;

/* global definitions */
TimerManager _timer;
ConfigParser _config;
#ifndef dt
  real_t dt;
#endif
#ifndef dx
  real_t dx;
#endif

int main(int argc, char **argv)
{
  // If not compiled in, set dt, dx
  // TODO: set other things (and check for performance hits)
  //       set these in config file?
  //       big performance hit setting N dynamically, should deal with BCs separately.
  #ifndef dx
    dx = H_LEN_FRAC/(1.0*N);
  #endif
  #ifndef dt
    dt = 0.1*dx;
  #endif

  // read in config file
  if(argc != 2)
  {
    std::cout << "Error: please supply exactly one config filename as an argument.\n";
    return EXIT_FAILURE;
  }
  // TODO: eliminate global _config; feed directly into constructor.
  _config.parse(argv[1]);

  // number of threads
  omp_set_num_threads(stoi(_config["omp_num_threads"]));

  // Create simulation according to simulation_type
  CosmoSim * cosmoSim;
  std::string simulation_type = _config["simulation_type"];
  if( simulation_type == "dust" )
  {
    cosmoSim = new DustSim();
  }
  else if( simulation_type == "dust_lambda" )
  {
    cosmoSim = new DustLambdaSim();
  }
  else if( simulation_type == "particles" )
  {
    cosmoSim = new ParticleSim();
  }
  else if( simulation_type == "scalar" )
  {
    cosmoSim = new ScalarSim();
  }
  else if( simulation_type == "vacuum" )
  {
    cosmoSim = new VacuumSim();
  }
  else
  {
    std::cerr << "Invalid simulation type specified.";
    throw 2;
  }

  // Initialize simulation 
  cosmoSim->init();

  // Run simulation
  cosmoSim->run();

  // Clean up
  cosmoSim->~CosmoSim();

  return EXIT_SUCCESS;
}
