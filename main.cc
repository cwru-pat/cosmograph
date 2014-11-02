
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
  BSSNData b_paq = {0}; // data structure associated with bssn sim

  Hydro hydroSim (1.0/3.0); // fluid with some w_EOS
  HydroData h_paq = {0};

  // initial conditions
  _timer["init"].start();
  bssnSim.init();
  hydroSim.init();
  _timer["init"].stop();

  ofstream outFile;
  outFile.open("vals.dat");

  // evolve simulation
  _timer["loop"].start();
  for(idx_t s=0; s < 1000; ++s) {

    real_t tmp = bssnSim.fields["phi_p"][0];
    outFile << tmp << "\n";
    outFile.flush();

    // Run RK steps explicitly here (ties together BSSN + Hydro stuff).
    // See bssn class or hydro class for more comments.

    // Init arrays and calculate source term for next step
      // _p is copied to _a here, which hydro uses
      bssnSim.stepInit();
      // clear existing data
      bssnSim.clearSrc();
      // add hydro source to bssn sim
      hydroSim.addBSSNSrc(bssnSim.fields);

    // First RK step & Set Hydro Vars
    LOOP3(i, j, k)
    {
      bssnSim.K1CalcPt(i, j, k, &b_paq);

      // hydro takes data from existing data in b_paq (data in _a register)
      // need to set full metric components first; this calculation is only
      // done when explicitly called.
      bssnSim.set_full_metric_der(&b_paq);
      bssnSim.set_full_metric(&b_paq);
      hydroSim.setQuantitiesCell(&b_paq, &h_paq);
    }

DUMPSTUFF("postinit")

    // reset source using new metric
    bssnSim.clearSrc();
    // add hydro source to bssn sim
    hydroSim.addBSSNSrc(bssnSim.fields);

    bssnSim.regSwap_c_a();

    // Subsequent BSSN steps
      // Second RK step
      LOOP3(i, j, k)
      {
        bssnSim.K2CalcPt(i, j, k, &b_paq);
      }

      // reset source using new metric
      bssnSim.clearSrc();
      // add hydro source to bssn sim
      hydroSim.addBSSNSrc(bssnSim.fields);

      bssnSim.regSwap_c_a();
      // Third RK step
      LOOP3(i, j, k)
      {
        bssnSim.K3CalcPt(i, j, k, &b_paq);
      }

      // reset source using new metric
      bssnSim.clearSrc();
      // add hydro source to bssn sim
      hydroSim.addBSSNSrc(bssnSim.fields);

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
      // bssn _f <-> _p
      bssnSim.stepTerm();
      // hydro _a <-> _f
      hydroSim.stepTerm();

  }
  _timer["loop"].stop();

 outFile.close(); 

  _timer["MAIN"].stop();

  cout << endl << _timer << endl;

  return EXIT_SUCCESS;
}
