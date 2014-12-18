
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
  int steps;
  real_t rho_K, peak_amplitude;

  // read in config file
  if(argc != 2)
  {
    cout << "Error: please supply exactly one config filename as an argument.\n";
    return EXIT_FAILURE;
  }
  else
  {
    _config.parse(argv[1]);
    peak_amplitude = (real_t) stold(_config["peak_amplitude"]);
    rho_K = (real_t) stold(_config["rho_K"]);
    steps = stoi(_config["steps"]);
  }

  // Create simulation
  _timer["init"].start();
    // Fluid fields
    Hydro hydroSim (0.0/3.0); // fluid with some w_EOS
    HydroData h_paq = {0};
    hydroSim.init();
    // DE
    Lambda lambdaSim (0.0);

    // GR Fields
    BSSN bssnSim;
    BSSNData b_paq = {0}; // data structure associated with bssn sim
    bssnSim.init();

    // generic reusable fourier class for N^3 arrays
    Fourier fourier;
    fourier.Initialize(N, hydroSim.fields["UD_a"] /* just any N^3 array for planning */);

    // Determine initial conditions.
    ICsData i_paq = {0};
    i_paq.peak_k = 2.0/((real_t) N);
    i_paq.peak_amplitude = peak_amplitude; // figure out units here
    // Note that this is going to be the conformal density, not physical density.
    // This specifies the spectrum of fluctuations around the mean, but not the mean.
    set_gaussian_random_field(hydroSim.fields["UD_a"], &fourier, &i_paq);
    // Set physical density fluctuations and metric using UD_a
    set_physical_from_conformal(bssnSim.fields, hydroSim.fields, &fourier);
    // Set a background (roughly, an average) density, and extrinsic curvature
    set_density_and_K(bssnSim.fields, hydroSim.fields, rho_K);
  _timer["init"].stop();


  // evolve simulation
  _timer["loop"].start();
  for(idx_t s=0; s < steps; ++s) {

    // output simulation information
    io_dump_quantities(bssnSim.fields, hydroSim.fields, _config["outfile"]);

    // Run RK steps explicitly here (ties together BSSN + Hydro stuff).
    // See bssn class or hydro class for more comments.

    // Init arrays and calculate source term for next step
      // _p is copied to _a here, which hydro uses
      bssnSim.stepInit();
      // clear existing data
      bssnSim.clearSrc();
      // add hydro source to bssn sim
      hydroSim.addBSSNSrc(bssnSim.fields);
      lambdaSim.addBSSNSrc(bssnSim.fields);

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

    // reset source using new metric
    bssnSim.clearSrc();
    // add hydro source to bssn sim
    hydroSim.addBSSNSrc(bssnSim.fields);
    lambdaSim.addBSSNSrc(bssnSim.fields);

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
      lambdaSim.addBSSNSrc(bssnSim.fields);

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
      lambdaSim.addBSSNSrc(bssnSim.fields);

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

  _timer["MAIN"].stop();

  cout << endl << _timer << endl;

  return EXIT_SUCCESS;
}
