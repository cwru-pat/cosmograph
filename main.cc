
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
  Hydro hydroSim(&bssnSim); // one fluid

  // initial conditions
  _timer["init"].start();
  bssnSim.init();
  hydroSim.init();
  _timer["init"].stop();


  // evolve simulation
  _timer["loop"].start();
  for(idx_t i=0; i < 50; ++i) {

    cout << "  gamma11_p = " << bssnSim.fields["gamma11_p"][0]
         << "; phi_p = " << bssnSim.fields["phi_p"][0]
         << "; K_p = " << bssnSim.fields["K_p"][0]
         << " \n";

    // Run RK steps explicitly here (ties together BSSN + Hydro stuff).
    // See bssn class or hydro class for more comments.
    bssnSim.stepInit();

    // First RK step
    LOOP3(i, j, k)
    {
      bssnSim.K1CalcPt(i, j, k);
    }
    bssnSim.regSwap_c_a();

    // Second RK step
    LOOP3(i, j, k)
    {
      bssnSim.K2CalcPt(i, j, k);
    }
    bssnSim.regSwap_c_a();

    // Third RK step
    LOOP3(i, j, k)
    {
      bssnSim.K3CalcPt(i, j, k);
    }
    bssnSim.regSwap_c_a();

    // Fourth RK step
    LOOP3(i, j, k)
    {
      bssnSim.K4CalcPt(i, j, k);
    }

    // Wrap up
    bssnSim.stepTerm();

  }
  _timer["loop"].stop();

  _timer["MAIN"].stop();

  cout << endl << _timer << endl;

  return EXIT_SUCCESS;
}
