#include "cosmo_includes.h"
#include "cosmo_types.h"
#include "globals.h"

#include "cosmosim.h"

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
  _config.parse(argv[1]);

  // number of threads
  omp_set_num_threads(stoi(_config["omp_num_threads"]));

  // Create and run simulation
  CosmoSim * cosmoSim;
// TODO: eliminate global _config; feed directly into constructor.
  cosmoSim = new CosmoSim();
  cosmoSim->init();
  cosmoSim->run();

  cosmoSim.~cosmoSim();

  return EXIT_SUCCESS;
}
