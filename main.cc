
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

  Lambda lambdaSim (0.0);

  // initial conditions
  _timer["init"].start();
  bssnSim.init();
  hydroSim.init();
  hydroSim.addBSSNSrc(bssnSim.fields);
  _timer["init"].stop();

  ofstream outFile;
  outFile.open("vals.dat");

  // evolve simulation
  _timer["loop"].start();
  for(idx_t i=0; i < 10000; ++i) {

    real_t tmp = bssnSim.fields["phi_p"][0];
    outFile << tmp << "\n";
    outFile.flush();

    // "final" calculation always ends up in "_p" register.
    cout << "step # " << i
         << ": phi_p = " << bssnSim.fields["phi_p"][10]
         << "\t K_p = " << bssnSim.fields["K_p"][10]
         << "\t A12_a = " << bssnSim.fields["A12_a"][10]
         << "\t UD_a = " << hydroSim.fields["UD_a"][10]
         << "\t r_a = " << bssnSim.fields["r_a"][10]
         << "\t S_a = " << bssnSim.fields["S_a"][10]
         << " \n";

    // Run RK steps explicitly here (ties together BSSN + Hydro stuff).
    // See bssn class or hydro class for more comments.
    bssnSim.stepInit();

    // First RK step & Set Hydro Vars
    LOOP3(i, j, k)
    {
      bssnSim.K1CalcPt(i, j, k, &b_paq);

      // hydro acts upon _a arrays, takes data from existing data in b_paq
      // need to set full metric components first; this calculation is only
      // done when explicitly called.
      bssnSim.set_full_metric_der(&b_paq);
      bssnSim.set_full_metric(&b_paq);
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

    // Calculate new source term for next step of bssn simulation
      // clear existing data
      bssnSim.clearSrc();
      // add in new evolved source
      hydroSim.addBSSNSrc(bssnSim.fields);

    // Add in 'Lambda' term
      lambdaSim.addBSSNSrc(bssnSim.fields);

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
