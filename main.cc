
#include "cosmo.h"
#include "globals.h"
#include <cmath>
#include <cfloat>

using namespace std;
using namespace cosmo;

/* global definitions */
TimerManager _timer;
ConfigParser _config;

int main(int argc, char **argv)
{
  _timer["MAIN"].start();
  idx_t i, j, k, s, steps;

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

  // IO init - will use this for logging.
  IOData iodata;
  io_init(&iodata, _config["output_dir"]);
  // save a copy of config.txt
  io_config_backup(&iodata, argv[1]);

  steps = stoi(_config["steps"]);
  omp_set_num_threads(stoi(_config["omp_num_threads"]));

  // Create simulation
  _timer["init"].start();
    LOG(iodata.log, "Creating initial conditions...\n");

    // Trial FRW class
    FRW<real_t> frw (0.0, 0.0);
    frw.addFluid(0.5, 0.0);

    // Fluid fields
    Hydro hydroSim (0.0/3.0); // fluid with some w_EOS
    HydroData h_paq = {0};
    hydroSim.init();
    // DE
    Lambda lambdaSim;

    // GR Fields
    BSSN bssnSim;
    BSSNData b_paq = {0}; // data structure associated with bssn sim
    bssnSim.init();

    // generic reusable fourier class for N^3 arrays
    Fourier fourier;
    fourier.Initialize(N, hydroSim.fields["UD_a"] /* just any N^3 array for planning */);

    // Hydro ICs
    LOG(iodata.log, "Using shock initial conditions...\n");
    set_shock_ICs(bssnSim.fields, hydroSim.fields, &fourier, &iodata);

  _timer["init"].stop();

  // evolve simulation
  LOG(iodata.log, "Running simulation...\n");
  _timer["loop"].start();
  for(s=0; s < steps; ++s) {

    // output simulation information
    _timer["output"].start();
      io_show_progress(s, steps);
      io_data_dump(bssnSim.fields, hydroSim.fields, &iodata, s, &fourier);
      io_dump_strip(hydroSim.fields["US1_a"], 1, 15, 15, &iodata);
      io_dump_strip(hydroSim.fields["UD_a"], 1, 31, 31, &iodata);
    _timer["output"].stop();

    // Run RK steps explicitly here (ties together BSSN + Hydro stuff).
    // See bssn class or hydro class for more comments.
    _timer["Hydro_steps"].start();

    // No BSSN evolution here- just flat metric. Set Hydro Vars.
    #pragma omp parallel for default(shared) private(i, j, k, b_paq, h_paq)
    LOOP3(i, j, k)
    {
      bssnSim.set_paq_values(i, j, k, &b_paq);
      bssnSim.set_full_metric_der(&b_paq);
      bssnSim.set_full_metric(&b_paq);
      hydroSim.setQuantitiesCell(&b_paq, &h_paq);
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
      // hydro _a <-> _f
      hydroSim.stepTerm();
    _timer["Hydro_steps"].stop();
  }
  _timer["loop"].stop();

  _timer["MAIN"].stop();

  LOG(iodata.log, endl << _timer << endl);

  return EXIT_SUCCESS;
}
