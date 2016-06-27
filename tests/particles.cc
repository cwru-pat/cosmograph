// g++ -g --std=c++11 -lrt particles.cc ../bssn/bssn.cc ../particles/particles.cc ../utils/Timer.cc ../utils/ConfigParser.cc

#include "../cosmo_includes.h"
#include "../cosmo_types.h"
#include "../cosmo_globals.h"

#include "../particles/particles.h"
#include "../bssn/bssn.h"

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

int main()
{
  std::cout << "Creating BSSN class..." << std::endl;
  // initialize GR sim; particles needs to reference this
  // but - don't evolve, remaining Minkowski/flat
  BSSN bssnSim;
  bssnSim.init();

  // make a particle
  std::cout << "Creating Particles class..." << std::endl;
  Particles particles;
  particles.init(5);

  std::cout << "Integrating Particle EOMs..." << std::endl;

  // initialize particle sim
  particles.stepInit(bssnSim.fields);

  // First RK step
  particles.RK1Step(bssnSim.fields);
  particles.regSwap_c_a();

  // Second RK step
  particles.RK2Step(bssnSim.fields);
  particles.regSwap_c_a();

  // Third RK step
  particles.RK3Step(bssnSim.fields);
  particles.regSwap_c_a();

  // Fourth RK step
  particles.RK4Step(bssnSim.fields);

  // Wrap up
  particles.stepTerm();

  std::cout << "Testing adding particle SRC to BSSN sim..." << std::endl;
  particles.addParticlesToBSSNSrc(bssnSim.fields);

  std::cout << "Done!" << std::endl;
  exit(EXIT_SUCCESS);
}
