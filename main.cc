
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
  idx_t i=0, j=0, k=0, s=0, steps=0;

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

    // GR Fields
    BSSN bssnSim;
    bssnSim.init();

    // generic reusable fourier class for NX*NY*NZ arrays
    Fourier fourier;
    fourier.Initialize(NX, NY, NZ, staticSim.fields["DIFFD_a"] /* just any array for planning */);

    // Standard FRW spacetime integrator
    FRW<real_t> frw (0.0, 0.0);

    // "conformal" initial conditions:
    LOG(iodata.log, "Using conformal initial conditions...\n");
    set_conformal_ICs(bssnSim.fields, staticSim.fields, &fourier, &iodata, &frw);

  _timer["init"].stop();

  // evolve simulation
  LOG(iodata.log, "Running simulation...\n");

  _timer["loop"].start();
  for(s=0; s < steps; ++s) {

    // initialize data for RK step in first loop
    // Init arrays and calculate source term for next step
      // _p is copied to _a here (which matter sectors use)
      bssnSim.stepInit();
      // clear existing source data
      bssnSim.clearSrc();
      staticSim.addBSSNSrc(bssnSim.fields, &frw);

    // output simulation information
    // these generally output any data in the _a registers (which should 
    // be identical to _p at this point).
    _timer["output"].start();
      // set_paq_values calculates ricci_a and AijAij_a data, needed for output
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i,j,k)
      {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.set_paq_values(i,j,k,&b_paq, &frw);
        // Additionally set KD (killing vector "Delta" quantities)
        bssnSim.set_KillingDelta(i, j, k, &b_paq);
      }
      io_show_progress(s, steps);
      io_data_dump(bssnSim.fields, staticSim.fields, &iodata, s, &fourier, &frw);
    _timer["output"].stop();

    // Run RK steps explicitly here (ties together BSSN + Hydro stuff).
    // See bssn class or hydro class for more comments.
    _timer["RK_steps"].start();


    // Pre-set detgamma and set gammai in _a register
    #pragma omp parallel for default(shared) private(i, j, k)
    LOOP3(i, j, k) {
      bssnSim.set_detgamma(i, j, k);
      bssnSim.set_DIFFgammai_values(i, j, k);
    }
    // FRW simulation should be in the correct state here
    // First RK step, Set Hydro Vars, & calc. constraint
    #pragma omp parallel for default(shared) private(i, j, k)
    LOOP3(i, j, k)
    {
      BSSNData b_paq = {0}; // data structure associated with bssn sim
      bssnSim.K1CalcPt(i, j, k, &b_paq, &frw);
    }

    // Intermediate RK step is now in _c register, move to _a for use in next step.
    bssnSim.regSwap_c_a();
    // reset source using new metric
    bssnSim.clearSrc();
    staticSim.addBSSNSrc(bssnSim.fields, &frw);

    // Subsequent BSSN steps
      // FRW simulation should progress to t -> t + dt/2
      for(int d=0; d<FRW_SUBSTEPS/2; ++d)
      {
        frw.step(dt/( (real_t) FRW_SUBSTEPS) );
      }

      // Second RK step
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k) {
        bssnSim.set_detgamma(i, j, k);
        bssnSim.set_DIFFgammai_values(i, j, k);
      }
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k)
      {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.K2CalcPt(i, j, k, &b_paq, &frw);
      }

      // Intermediate RK step is now in _c register, move to _a for use in next step.
      bssnSim.regSwap_c_a();
      // reset source using new metric
      bssnSim.clearSrc();
      staticSim.addBSSNSrc(bssnSim.fields, &frw);

      // Third RK step
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k) {
        bssnSim.set_detgamma(i, j, k);
        bssnSim.set_DIFFgammai_values(i, j, k);
      }
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k)
      {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.K3CalcPt(i, j, k, &b_paq, &frw);
      }

      // Intermediate RK step is now in _c register, move to _a for use in next step.
      bssnSim.regSwap_c_a();
      // reset source using new metric
      bssnSim.clearSrc();
      staticSim.addBSSNSrc(bssnSim.fields, &frw);

      // FRW simulation should progress to t -> t + dt
      // (progress dt/2 further)
      for(int d=0; d<FRW_SUBSTEPS/2; ++d)
      {
        frw.step(dt/( (real_t) FRW_SUBSTEPS) );
      }

      // Fourth RK step
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k) {
        bssnSim.set_detgamma(i, j, k);
        bssnSim.set_DIFFgammai_values(i, j, k);
      }
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k)
      {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.K4CalcPt(i, j, k, &b_paq, &frw);
      }

    // Wrap up
      // bssn _f <-> _p
      bssnSim.stepTerm();
      // "current" data is in the _p array.

    _timer["RK_steps"].stop();


    _timer["meta_output_interval"].start();

      if(s%iodata.meta_output_interval == 0)
      {
        idx_t isNaN = 0;
        real_t H_calcs[7], M_calcs[7];

        // make sure clean data is in _a array
        // Set values at points
        bssnSim.stepInit(); // copy _p to _a
        bssnSim.clearSrc();
        staticSim.addBSSNSrc(bssnSim.fields, &frw);
        #pragma omp parallel for default(shared) private(i, j, k)
        LOOP3(i, j, k) {
          bssnSim.set_detgamma(i, j, k);
          bssnSim.set_DIFFgammai_values(i, j, k);
        }
        #pragma omp parallel for default(shared) private(i, j, k)
        LOOP3(i, j, k) {
          BSSNData b_paq = {0}; // data structure associated with bssn sim
          bssnSim.set_paq_values(i,j,k,&b_paq, &frw);
        }

        // Constraint Violation Calculations
        bssnSim.setHamiltonianConstraintCalcs(H_calcs, &frw, false);
        io_dump_data(H_calcs[4], &iodata, "H_violations"); // mean(H/[H])
        io_dump_data(H_calcs[5], &iodata, "H_violations"); // stdev(H/[H])
        io_dump_data(H_calcs[6], &iodata, "H_violations"); // mean(H/[H])
        io_dump_data(H_calcs[2], &iodata, "H_violations"); // max(H)

        bssnSim.setMomentumConstraintCalcs(M_calcs, &frw);
        io_dump_data(M_calcs[4], &iodata, "M_violations"); // mean(M/[M])
        io_dump_data(M_calcs[5], &iodata, "M_violations"); // stdev(M/[M])
        io_dump_data(M_calcs[6], &iodata, "M_violations"); // mean(M/[M])
        io_dump_data(M_calcs[2], &iodata, "M_violations"); // max(M)

        if(numNaNs(bssnSim.fields["DIFFphi_a"]) > 0)
        {
          LOG(iodata.log, "\nNAN detected!\n");
          _timer["meta_output_interval"].stop();
          break;
        }
      }

    _timer["meta_output_interval"].stop();

  }
  _timer["loop"].stop();

  _timer["output"].start();
  io_dump_3dslice(bssnSim.fields["dk0_slice_phi_a"], "dk0_slice_phi", &iodata);
  LOG(iodata.log, "\nAverage conformal factor reached " << average(bssnSim.fields["DIFFphi_p"]) << "\n");
  LOG(iodata.log, "Ending simulation.\n");
  _timer["output"].stop();

  _timer["MAIN"].stop();

  LOG(iodata.log, endl << _timer << endl);

  return EXIT_SUCCESS;
}
