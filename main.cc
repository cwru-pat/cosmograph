
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

    // Fluid fields
    // Static matter (w=0)
    Static staticSim;
    staticSim.init();
    // DE
    Lambda lambdaSim;

    // GR Fields
    BSSN bssnSim;
    bssnSim.init();

    // generic reusable fourier class for N^3 arrays
    Fourier fourier;
    fourier.Initialize(N, staticSim.fields["D_a"] /* just any N^3 array for planning */);

    // "conformal" initial conditions:
    LOG(iodata.log, "Using conformal initial conditions...\n");
    set_conformal_ICs(bssnSim.fields, staticSim.fields, &fourier, &iodata);

    // Trial FRW class
    staticSim.addBSSNSrc(bssnSim.fields);
    real_t frw_rho = average(bssnSim.fields["r_a"]);
    FRW<real_t> frw (0.0, -sqrt(24.0*PI*frw_rho) /* K */);
    frw.addFluid(frw_rho /* rho */, 0.0 /* 'w' */);

  _timer["init"].stop();

  // evolve simulation
  LOG(iodata.log, "Running simulation...\n");
  _timer["loop"].start();
  for(s=0; s < steps; ++s) {

    _timer["Reference FRW"].start();
    // LOG(iodata.log,    "\n"
    //                 << frw.get_phi()
    //                 << "\n"
    //                 << bssnSim.fields["phi_a"][0]
    //                 << "\n"
    //               );
    int subiters = 100;
    for(int p=0; p<subiters; ++p)
    {
      frw.step(dt/((real_t) subiters));
    }
    _timer["Reference FRW"].stop();

    // output simulation information
    _timer["output"].start();
      io_show_progress(s, steps);
      io_data_dump(bssnSim.fields, staticSim.fields, &iodata, s, &fourier);
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
      staticSim.addBSSNSrc(bssnSim.fields);
      lambdaSim.addBSSNSrc(bssnSim.fields);

    // First RK step, Set Hydro Vars, & calc. constraint
    #pragma omp parallel for default(shared) private(i, j, k)
    LOOP3(i, j, k)
    {
      BSSNData b_paq = {0}; // data structure associated with bssn sim
      bssnSim.K1CalcPt(i, j, k, &b_paq);
    }

    // reset source using new metric
    bssnSim.clearSrc();
    // add hydro source to bssn sim
    staticSim.addBSSNSrc(bssnSim.fields);
    lambdaSim.addBSSNSrc(bssnSim.fields);

    bssnSim.regSwap_c_a();

    // Subsequent BSSN steps
      // Second RK step
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k)
      {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.K2CalcPt(i, j, k, &b_paq);
      }

      // reset source using new metric
      bssnSim.clearSrc();
      // add hydro source to bssn sim
      staticSim.addBSSNSrc(bssnSim.fields);
      lambdaSim.addBSSNSrc(bssnSim.fields);

      bssnSim.regSwap_c_a();
      // Third RK step
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k)
      {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.K3CalcPt(i, j, k, &b_paq);
      }

      // reset source using new metric
      bssnSim.clearSrc();
      // add hydro source to bssn sim
      staticSim.addBSSNSrc(bssnSim.fields);
      lambdaSim.addBSSNSrc(bssnSim.fields);

      bssnSim.regSwap_c_a();

      // Fourth RK step
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k)
      {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.K4CalcPt(i, j, k, &b_paq);
      }

    // Wrap up
      // bssn _f <-> _p
      bssnSim.stepTerm();

    _timer["RK_steps"].stop();
    _timer["output"].start();
      real_t total_hamiltonian_constraint = 0.0,
             mean_hamiltonian_constraint = 0.0,
             stdev_hamiltonian_constraint = 0.0;
      if(s%iodata.meta_output_interval == 0)
      {
        idx_t isNaN = 0;
        bssnSim.stepInit();
        bssnSim.clearSrc();
        staticSim.addBSSNSrc(bssnSim.fields);
        #pragma omp parallel for default(shared) private(i, j, k) \
         reduction(+:total_hamiltonian_constraint, mean_hamiltonian_constraint, isNaN)
        LOOP3(i,j,k)
        {
          BSSNData b_paq = {0};
          bssnSim.set_paq_values(i, j, k, &b_paq);
          real_t ham = (b_paq.ricci + 2.0/3.0*pw2(b_paq.K) - 16.0*PI*b_paq.rho);
          real_t ham2 = sqrt(pw2(b_paq.ricci) + pw2(2.0/3.0*pw2(b_paq.K)) + pw2(16.0*PI*b_paq.rho));
          real_t violation_fraction = ham/ham2;
          total_hamiltonian_constraint += fabs(violation_fraction);
          mean_hamiltonian_constraint += violation_fraction/POINTS;

          union { float violation_fraction; uint32_t x; } u = { (float) violation_fraction };
          if((u.x << 1) > 0xff000000u && isNaN == 0)
          {
            LOG(iodata.log, "NAN at (" << i << "," << j << "," << k << ")!\n");
            isNaN += 1;
          }
        }
        if(isNaN > 0)
        {
          LOG(iodata.log, "\nNAN detected!\n");
          _timer["output"].stop();
          break;
        }
        #pragma omp parallel for default(shared) private(i, j, k) reduction(+:stdev_hamiltonian_constraint)
        LOOP3(i,j,k)
        {
          BSSNData b_paq = {0};
          bssnSim.set_paq_values(i, j, k, &b_paq);
          real_t ham = (b_paq.ricci + 2.0/3.0*pw2(b_paq.K) - 16.0*PI*b_paq.rho);
          real_t ham2 = sqrt(pw2(b_paq.ricci) + pw2(2.0/3.0*pw2(b_paq.K)) + pw2(16.0*PI*b_paq.rho));
          real_t violation_fraction = ham/ham2;
          stdev_hamiltonian_constraint += pw2(violation_fraction - mean_hamiltonian_constraint);
        }
        stdev_hamiltonian_constraint = sqrt(stdev_hamiltonian_constraint/(POINTS-1.0));
        io_dump_data(mean_hamiltonian_constraint, &iodata, "avg_H_violation");
        io_dump_data(stdev_hamiltonian_constraint, &iodata, "std_H_violation");
      }
    _timer["output"].stop();
  }
  _timer["loop"].stop();

  _timer["output"].start();
  LOG(iodata.log, "\nAverage conformal factor reached " << average(bssnSim.fields["phi_p"]) << "\n");
  LOG(iodata.log, "Ending simulation.\n");
  _timer["output"].stop();

  _timer["MAIN"].stop();

  LOG(iodata.log, endl << _timer << endl);

  return EXIT_SUCCESS;
}
