
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
    // GR Fields
    BSSN bssnSim;
    bssnSim.clearSrc();
    BSSNData b_paq = {0}; // data structure associated with bssn sim
    bssnSim.init();

    // Determine initial conditions.
    ICsData i_paq = {0};
    set_BH_ICs(bssnSim.fields);
  _timer["init"].stop();

  // evolve simulation
  std::cout << "Running Simulation...\n";
  _timer["loop"].start();
  for(s=0; s < steps; ++s) {

    // output simulation information
    _timer["output"].start();
    std::cout << "Running step " << s << ".  Phi is:" << bssnSim.fields["phi_p"][INDEX(N/2, N/2, N/2)] << "\n";
    std::cout << "                Alpha is:" << bssnSim.fields["alpha_p"][INDEX(1, N/2, N/2)]
              << ", " << bssnSim.fields["alpha_p"][INDEX(14, N/2, N/2)]
              << "\n";
    io_dump_strip(bssnSim.fields["phi_p"], 1, N/2, N/2, &iodata);
    if(s%slice_output_interval == 0)
    {
      io_dump_2dslice(bssnSim.fields["phi_p"], "phi_slice." + to_string(s), &iodata);
    }
    if(s%grid_output_interval == 0)
    {
      // io_dump_3dslice(bssnSim.fields["phi_p"],     "phi."     + to_string(s), &iodata);
    }
    _timer["output"].stop();

    // Run RK steps explicitly here (ties together BSSN + Hydro stuff).
    // See bssn class or hydro class for more comments.
    _timer["RK_steps"].start();

    // Init arrays and calculate source term for next step
      // _p is copied to _a here, which hydro uses
      bssnSim.stepInit();

    // First RK step
    #pragma omp parallel for default(shared) private(i, j, k, b_paq)
    LOOP3(i, j, k)
    {
      bssnSim.K1CalcPt(i, j, k, &b_paq);
    }
    bssnSim.apply_boundaries();
    bssnSim.regSwap_c_a();

    // Subsequent BSSN steps
      // Second RK step
      #pragma omp parallel for default(shared) private(i, j, k, b_paq)
      LOOP3(i, j, k)
      {
        bssnSim.K2CalcPt(i, j, k, &b_paq);
      }
      bssnSim.apply_boundaries();
      bssnSim.regSwap_c_a();

      // Third RK step
      #pragma omp parallel for default(shared) private(i, j, k, b_paq)
      LOOP3(i, j, k)
      {
        bssnSim.K3CalcPt(i, j, k, &b_paq);
      }
      bssnSim.apply_boundaries();
      bssnSim.regSwap_c_a();

      // Fourth RK step
      #pragma omp parallel for default(shared) private(i, j, k, b_paq)
      LOOP3(i, j, k)
      {
        bssnSim.K4CalcPt(i, j, k, &b_paq);
      }
      bssnSim.apply_boundaries();

    // Wrap up
      // bssn _f <-> _p
      bssnSim.stepTerm();
    _timer["RK_steps"].stop();
  }
  _timer["loop"].stop();

  _timer["MAIN"].stop();

  cout << endl << _timer << endl;

  return EXIT_SUCCESS;
}
