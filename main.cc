
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
  idx_t steps, slice_output_interval, grid_output_interval;
  idx_t i, j, k, s;

  real_t rho_K_matter, rho_K_lambda, peak_amplitude;

  // read in config file
  if(argc != 2)
  {
    cout << "Error: please supply exactly one config filename as an argument.\n";
    return EXIT_FAILURE;
  }
  else
  {
    _config.parse(argv[1]);
    peak_amplitude = (real_t) stold(_config["peak_amplitude"]); // fluctuation amplitude
    rho_K_matter = (real_t) stold(_config["rho_K_matter"]); // background density
    rho_K_lambda = (real_t) stold(_config["rho_K_lambda"]); // DE density
    steps = stoi(_config["steps"]);
    slice_output_interval = stoi(_config["slice_output_interval"]);
    grid_output_interval = stoi(_config["grid_output_interval"]);
    omp_set_num_threads(stoi(_config["omp_num_threads"]));
  }

  // IO init
  IOData iodata;
  iodata.output_dir = _config["output_dir"];
  io_init(&iodata);

  // Create simulation
  std::cout << "Creating initial conditions...\n";
  _timer["init"].start();
    // Fluid fields
    Hydro hydroSim (0.0/3.0); // fluid with some w_EOS
    HydroData h_paq = {0};
    hydroSim.init();
    // DE
    Lambda lambdaSim (rho_K_lambda);

    // GR Fields
    BSSN bssnSim;
    BSSNData b_paq = {0}; // data structure associated with bssn sim
    bssnSim.init();

    // generic reusable fourier class for N^3 arrays
    Fourier fourier;
    fourier.Initialize(N, hydroSim.fields["UD_a"] /* just any N^3 array for planning */);

    // Determine initial conditions.
    ICsData i_paq = {0};
    i_paq.peak_k = 1.0/((real_t) N);
    i_paq.peak_amplitude = peak_amplitude; // figure out units here
    // Note that this is going to be the conformal density, not physical density.
    // This specifies the spectrum of fluctuations around the mean, but not the mean.
    set_gaussian_random_field(hydroSim.fields["UD_a"], &fourier, &i_paq);
    // Set physical density fluctuations and metric using UD_a
    set_physical_from_conformal(bssnSim.fields, hydroSim.fields, &fourier);
    // Set a background (roughly, an average) density, and extrinsic curvature
    set_matter_density_and_K(bssnSim.fields, hydroSim.fields, rho_K_matter);
    set_lambda_K(bssnSim.fields, rho_K_lambda);
  _timer["init"].stop();

  // evolve simulation
  std::cout << "Running Simulation...\n";
  _timer["loop"].start();
  for(s=0; s < steps; ++s) {

    // output simulation information
    _timer["output"].start();
    io_dump_quantities(bssnSim.fields, hydroSim.fields, _config["outfile"], &iodata);
    if(s%slice_output_interval == 0)
    {
      io_dump_2dslice(bssnSim.fields["phi_p"], "phi_slice." + to_string(s), &iodata);
      io_dump_2dslice(hydroSim.fields["UD_f"], "UD_slice."  + to_string(s), &iodata);
    }
    if(s%grid_output_interval == 0)
    {
      io_dump_3dslice(bssnSim.fields["gamma11_p"], "gamma11." + to_string(s), &iodata);
      io_dump_3dslice(bssnSim.fields["gamma12_p"], "gamma12." + to_string(s), &iodata);
      io_dump_3dslice(bssnSim.fields["gamma13_p"], "gamma13." + to_string(s), &iodata);
      io_dump_3dslice(bssnSim.fields["gamma22_p"], "gamma22." + to_string(s), &iodata);
      io_dump_3dslice(bssnSim.fields["gamma23_p"], "gamma23." + to_string(s), &iodata);
      io_dump_3dslice(bssnSim.fields["gamma33_p"], "gamma33." + to_string(s), &iodata);
      io_dump_3dslice(bssnSim.fields["phi_p"],     "phi."     + to_string(s), &iodata);
      io_dump_3dslice(bssnSim.fields["A11_p"],     "A11."     + to_string(s), &iodata);
      io_dump_3dslice(bssnSim.fields["A12_p"],     "A12."     + to_string(s), &iodata);
      io_dump_3dslice(bssnSim.fields["A13_p"],     "A13."     + to_string(s), &iodata);
      io_dump_3dslice(bssnSim.fields["A22_p"],     "A22."     + to_string(s), &iodata);
      io_dump_3dslice(bssnSim.fields["A23_p"],     "A23."     + to_string(s), &iodata);
      io_dump_3dslice(bssnSim.fields["A33_p"],     "A33."     + to_string(s), &iodata);
      io_dump_3dslice(bssnSim.fields["K_p"],       "K."       + to_string(s), &iodata);
      io_dump_3dslice(bssnSim.fields["ricci_a"],   "ricci."   + to_string(s), &iodata);
      io_dump_3dslice(hydroSim.fields["UD_f"],     "UD."      + to_string(s), &iodata);
      io_dump_3dslice(hydroSim.fields["US1_f"],    "US1."     + to_string(s), &iodata);
      io_dump_3dslice(hydroSim.fields["US2_f"],    "US2."     + to_string(s), &iodata);
      io_dump_3dslice(hydroSim.fields["US3_f"],    "US3."     + to_string(s), &iodata);
    }
    _timer["output"].stop();

    // Run RK steps explicitly here (ties together BSSN + Hydro stuff).
    // See bssn class or hydro class for more comments.
    _timer["RK_steps"].start();

    // Init arrays and calculate source term for next step
      // _p is copied to _a here, which hydro uses
      bssnSim.stepInit();
      // clear existing data
      bssnSim.clearSrc();
      // add hydro source to bssn sim
      hydroSim.addBSSNSrc(bssnSim.fields);
      lambdaSim.addBSSNSrc(bssnSim.fields);

    // First RK step & Set Hydro Vars
    #pragma omp parallel for default(shared) private(i, j, k, b_paq, h_paq)
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
      #pragma omp parallel for default(shared) private(i, j, k, b_paq)
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
      #pragma omp parallel for default(shared) private(i, j, k, b_paq)
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
      #pragma omp parallel for default(shared) private(i, j, k, b_paq)
      LOOP3(i, j, k)
      {
        bssnSim.K4CalcPt(i, j, k, &b_paq);
      }


    // Subsequent hydro step
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k)
      {
        hydroSim.setAllFluxInt(i, j, k);
      }
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k)
      {
        hydroSim.evolveFluid(i, j, k);
      }

    // Wrap up
      // bssn _f <-> _p
      bssnSim.stepTerm();
      // hydro _a <-> _f
      hydroSim.stepTerm();
    _timer["RK_steps"].stop();
  }
  _timer["loop"].stop();

  _timer["MAIN"].stop();

  cout << endl << _timer << endl;

  return EXIT_SUCCESS;
}
