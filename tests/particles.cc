// g++ -g --std=c++11 -lrt particles.cc ../components/bssn/bssn.cc ../components/bssn/BSSNGaugeHandler.cc ../components/particles/particles.cc ../utils/Timer.cc ../utils/ConfigParser.cc

#include "../cosmo_includes.h"
#include "../cosmo_types.h"
#include "../cosmo_globals.h"

#include "../components/particles/particles.h"
#include "../components/bssn/bssn.h"

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
  dt = 0.5;
  dx = 1.0;

  std::cout << "Creating BSSN class..." << std::endl;
  // initialize GR sim; particles needs to reference this
  // but - don't evolve, remaining Minkowski/flat
  BSSN bssnSim(&_config);

  // make a particle
  std::cout << "Creating Particles class..." << std::endl;
  Particles particles;
  particles.init(1);

  std::cout << "Integrating Particle EOMs..." << std::endl;
  for(int s=0; s<10000; s++)
  {
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
  }

  // Don't source for now
  // std::cout << "Testing adding particle SRC to BSSN sim..." << std::endl;
  // particles.addParticlesToBSSNSrc(bssnSim.fields);

  std::cout << "Done!" << std::endl;
  exit(EXIT_SUCCESS);
}
