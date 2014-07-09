
#include "cosmo.h"
#include "globals.h"

using namespace std;
using namespace cosmo;

/* global definitions */
TimerManager _timer;
ConfigParser _config;

int main(int argc, char **argv)
{
  _timer["MAIN"].start();

  // read in config file
  if(argc != 2)
  {
    cout << "Error: please supply exactly one config filename as an argument.\n";
    return EXIT_FAILURE;
  }
  else
  {
    _config.parse(argv[1]);
  }

  // Create simulation
  BSSN bssnSim;

  // initial conditions
  _timer["init"].start();
  bssnSim.init();
  _timer["init"].stop();
  _timer["output"].start();
  //io_dump_strip(bssnSim.fields["alpha_p"], 1, 1, 1);
  _timer["output"].stop();

  // evolve simulation
  _timer["loop"].start();
  for(idx_t i=0; i < 10; ++i) {
    cout << "  phi_p = "     << bssnSim.fields["phi_p"][0]
         << "; K_p = "       << bssnSim.fields["K_p"][0]
         << "; A11_p = "     << bssnSim.fields["A11_p"][0]
         << "; gamma11_p = " << bssnSim.fields["gamma11_p"][0]
         << "; Gamma1_p = "  << bssnSim.fields["Gamma1_p"][0]
         << " \n";

    bssnSim.step();

  }
  _timer["loop"].stop();

  _timer["MAIN"].stop();

  cout << endl << _timer << endl;

  return EXIT_SUCCESS;
}
