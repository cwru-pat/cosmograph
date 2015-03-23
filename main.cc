
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

    if(_config["ic_type"] == "flat")
    {
      // "flat dynamic" initial conditions:
      LOG(iodata.log, "Using flat initial conditions...\n");
      set_flat_dynamic_ICs(bssnSim.fields, hydroSim.fields, &fourier, &iodata);
    }
    else
    {
      // default to "conformal" initial conditions:
      LOG(iodata.log, "Using conformal initial conditions...\n");
      set_conformal_ICs(bssnSim.fields, hydroSim.fields, &fourier, &iodata);
    }

  _timer["init"].stop();

  // evolve simulation
  LOG(iodata.log, "Running simulation...\n");
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
      real_t total_hamiltonian_constraint = 0.0,
             mean_hamiltonian_constraint = 0.0,
             stdev_hamiltonian_constraint = 0.0,
             total_momentum_constraint_1 = 0.0,
             total_momentum_constraint_2 = 0.0,
             total_momentum_constraint_3 = 0.0;
      if(s%iodata.meta_output_interval == 0)
      {
        idx_t isNaN = 0;
        #pragma omp parallel for default(shared) private(i, j, k, b_paq) \
         reduction(+:total_hamiltonian_constraint, mean_hamiltonian_constraint, total_momentum_constraint_1, total_momentum_constraint_2, total_momentum_constraint_3, isNaN)
        LOOP3(i,j,k)
        {
          real_t violation_fraction = bssnSim.hamiltonianConstraintCalc(i,j,k)/bssnSim.hamiltonianConstraintMag(i,j,k);
          total_hamiltonian_constraint += fabs(violation_fraction);
          mean_hamiltonian_constraint += violation_fraction/POINTS;

          union { float violation_fraction; uint32_t x; } u = { (float) violation_fraction };
          if((u.x << 1) > 0xff000000u && isNaN == 0)
          {
            LOG(iodata.log, "NAN at (" << i << "," << j << "," << k << ")!\n");
            isNaN += 1;
          }

          bssnSim.set_paq_values(i, j, k, &b_paq);
          total_momentum_constraint_1 += fabs(bssnSim.momentumConstraintCalc(&b_paq, 1)/bssnSim.momentumConstraintMag(&b_paq, 1));
          total_momentum_constraint_2 += fabs(bssnSim.momentumConstraintCalc(&b_paq, 2)/bssnSim.momentumConstraintMag(&b_paq, 2));
          total_momentum_constraint_3 += fabs(bssnSim.momentumConstraintCalc(&b_paq, 3)/bssnSim.momentumConstraintMag(&b_paq, 3));
        }
        if(isNaN > 0)
        {
          LOG(iodata.log, "NAN detected!\n");
          _timer["output"].stop();
          break;
        }
        #pragma omp parallel for default(shared) private(i, j, k, b_paq) reduction(+:stdev_hamiltonian_constraint)
        LOOP3(i,j,k)
        {
          stdev_hamiltonian_constraint += pw2(bssnSim.hamiltonianConstraintCalc(i,j,k)/bssnSim.hamiltonianConstraintMag(i,j,k) - mean_hamiltonian_constraint);
        }
        io_dump_data(mean_hamiltonian_constraint, &iodata, "avg_H_violation");
        io_dump_data(sqrt(stdev_hamiltonian_constraint/(POINTS-1.0)), &iodata, "std_H_violation");
        std::cout << "\n" << "H-viol:            " << mean_hamiltonian_constraint
                  << "\n" << "H-viol std:        " << sqrt(stdev_hamiltonian_constraint/(POINTS-1.0))
                  << "\n" << "M-viol L1 average: " << 1.0/POINTS*sqrt(pw2(total_momentum_constraint_1)+pw2(total_momentum_constraint_2)+pw2(total_momentum_constraint_3))
                  << "\n";
      }
    _timer["output"].stop();
  }
  _timer["loop"].stop();

  _timer["MAIN"].stop();

  LOG(iodata.log, endl << _timer << endl);

  return EXIT_SUCCESS;
}
