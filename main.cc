
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
  idx_t i, j, k, s, steps;
  real_t total_hamiltonian_constraint = 0.0;

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
  std::cout << "Creating initial conditions...\n";
  _timer["init"].start();
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

    if(_config["ic_type"] == "flat")
    {
      // "flat dynamic" initial conditions:
      std::cout << "Using flat initial conditions...\n";
      set_flat_dynamic_ICs(bssnSim.fields, hydroSim.fields, &fourier);
    }
    else
    {
      // default to "conformal" initial conditions:
      std::cout << "Using conformal initial conditions...\n";
      set_conformal_ICs(bssnSim.fields, hydroSim.fields, &fourier);
    }

  _timer["init"].stop();

  // evolve simulation
  std::cout << "Running Simulation...\n";
  _timer["loop"].start();
  for(s=0; s < steps; ++s) {

    // output simulation information
    _timer["output"].start();
      io_show_progress(s, steps);
      io_data_dump(bssnSim.fields, hydroSim.fields, &iodata, s, &fourier);
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

    // First RK step, Set Hydro Vars, & calc. constraint
    total_hamiltonian_constraint = 0.0;
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

    _timer["output"].start();
      if(s%iodata.meta_output_interval == 0)
      {
        LOOP3(i,j,k)
          total_hamiltonian_constraint += fabs(
              bssnSim.hamiltonianConstraintCalc(NP_INDEX(i,j,k))/bssnSim.hamiltonianConstraintMag(NP_INDEX(i,j,k))
            );
        io_dump_data(total_hamiltonian_constraint/POINTS, &iodata, "avg_H_violation");
      }
    _timer["output"].stop();
  }
  _timer["loop"].stop();


  _timer["MAIN"].stop();

  cout << endl << _timer << endl;

  return EXIT_SUCCESS;
}
