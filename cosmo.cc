#include "cosmo_includes.h"
#include "cosmo_types.h"
#include "cosmo_globals.h"

#include "sims/sim.h"
#include "sims/dust.h"
#include "sims/dust_lambda.h"
#include "sims/particles.h"
#include "sims/scalar.h"
#include "sims/vacuum.h"
#include "sims/sheets.h"

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

#ifndef Lx
  real_t Lx;
#endif
#ifndef Ly
  real_t Ly;
#endif
#ifndef Lz
  real_t Lz;
#endif

#ifndef Nx
  int Nx;
#endif
#ifndef Ny
  int Ny;
#endif
#ifndef Nz
  int Nz;
#endif


int main(int argc, char **argv)
{
  // read in config file
  if(argc != 2)
  {
    std::cout << "Error: please supply exactly one config filename as an argument.\n";
    return EXIT_FAILURE;
  }
  // TODO: eliminate global _config; feed directly into constructor.
  _config.parse(argv[1]);

  // If not compiled in, set dt, dx
  // TODO: set other things (and check for performance hits)
  //       set these in config file?
  //       big performance hit setting N dynamically, should deal with BCs separately.
  
  
  Lx = stod(_config( "Lx", "1" ));
  Ly = stod(_config( "Ly", "1" ));
  Lz = stod(_config( "Lz", "1" ));

  Nx = stoi(_config( "Nx", "32" ));
  Ny = stoi(_config( "Ny", "32" ));
  Nz = stoi(_config( "Nz", "32" ));

  if( fabs(H_LEN_FRAC/(1.0*COSMO_N) - Lx / (double)Nx) > 1E-9
      || fabs(H_LEN_FRAC/(1.0*COSMO_N) - Ly / (double)Ny) > 1E-9
      || fabs(H_LEN_FRAC/(1.0*COSMO_N) - Lz / (double)Nz) > 1E-9)
  {
    std::cout << "Macro setting is not equal to parameter setting!\n";
    return EXIT_FAILURE;
  }
    
  
  dx = stold(_config( "dx", stringify(H_LEN_FRAC/(1.0*COSMO_N)) ));
  dt = stold(_config( "dt_frac", "0.1" ))*dx;

  
  
  // Set number of threads if specified
  int num_threads = stoi(_config("omp_num_threads", "1"));
  if(num_threads >= 1)
    omp_set_num_threads(num_threads);

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
  else if( simulation_type == "sheets")
  {
   
    cosmoSim = new SheetSim();
  }
  else
  {
    std::cerr << "Invalid simulation type specified.";
    throw 2;
  }

  // Initialize simulation 
  cosmoSim->init();

  // Generate initial conditions
  cosmoSim->setICs();

  // Run simulation
  cosmoSim->run();

  delete cosmoSim;
  return EXIT_SUCCESS;
}
