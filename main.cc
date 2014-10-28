
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
  BSSNData b_paq; // data structure associated with bssn sim

  Hydro hydroSim; // one fluid
  HydroData h_paq;

  // initial conditions
  _timer["init"].start();
  bssnSim.init();
  hydroSim.init();
  _timer["init"].stop();


  // evolve simulation
  _timer["loop"].start();
  for(idx_t i=0; i < 5; ++i) {

    // "final" calculation always ends up in "_p" register.
    cout << "step # " << i
         << "; phi_p = " << bssnSim.fields["phi_p"][0]
         << "; K_p = " << bssnSim.fields["K_p"][0]
         << " \n";

    // Run RK steps explicitly here (ties together BSSN + Hydro stuff).
    // See bssn class or hydro class for more comments.
    bssnSim.stepInit();
    hydroSim.init();

    // First RK step & Set Hydro Vars
    LOOP3(i, j, k)
    {
      bssnSim.K1CalcPt(i, j, k, &b_paq);

      // hydro acts upon _a arrays, takes data from existing data in b_paq
      // need to set full metric components first; this calculation is only
      // done when explicitly called.
      bssnSim.set_full_metric_der(&b_paq);
      hydroSim.setQuantitiesCell(&b_paq, &h_paq);
    }
    bssnSim.regSwap_c_a();

    // Subsequent BSSN steps
      // Second RK step
      LOOP3(i, j, k)
      {
        bssnSim.K2CalcPt(i, j, k, &b_paq);
      }
      bssnSim.regSwap_c_a();
      // Third RK step
      LOOP3(i, j, k)
      {
        bssnSim.K3CalcPt(i, j, k, &b_paq);
      }
      bssnSim.regSwap_c_a();
      // Fourth RK step
      LOOP3(i, j, k)
      {
        bssnSim.K4CalcPt(i, j, k, &b_paq);
      }

    // Subsequent hydro step
      LOOP3(i, j, k)
      {
        hydroSim.setAllFluxInt(i, j, k);
      }
      LOOP3(i, j, k)
      {
        hydroSim.evolveFluid(i, j, k);
      }

    // Wrap up
    // clear out source term in bssn calc for replacement
    bssnSim.clearSrc();
    // add in new evolved source
    hydroSim.addBSSNSrc(bssnSim.fields);
    // bssn _f <-> _p
    bssnSim.stepTerm();

  }
  _timer["loop"].stop();

  _timer["MAIN"].stop();

  cout << endl << _timer << endl;

  return EXIT_SUCCESS;
}
