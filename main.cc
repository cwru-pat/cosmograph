
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
    // DE
    Lambda lambdaSim;

    // GR Fields
    BSSN bssnSim;
    bssnSim.init();

    // generic reusable fourier class for NX*NY*NZ arrays
    Fourier fourier;
    fourier.Initialize(NX, NY, NZ, staticSim.fields["D_a"] /* just any array for planning */);

    // "conformal" initial conditions:
    LOG(iodata.log, "Using conformal initial conditions...\n");
    set_conformal_ICs(bssnSim.fields, staticSim.fields, &fourier, &iodata);

    // Trial FRW class
    staticSim.addBSSNSrc(bssnSim.fields);
    real_t frw_rho = average(bssnSim.fields["r_a"]);
    real_t frw_K0  = -sqrt(24.0*PI*frw_rho);
    FRW<real_t> frw (0.0, frw_K0 /* K */);
    frw.addFluid(frw_rho /* rho */, 0.0 /* 'w' */);

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
      staticSim.addBSSNSrc(bssnSim.fields);
      lambdaSim.addBSSNSrc(bssnSim.fields);

    // output simulation information
    // these generally output any data in the _a registers (which should 
    // be identical to _p at this point).
    _timer["output"].start();
      // set_paq_values calculates ricci_a and AijAij_a data, needed for output
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i,j,k)
      {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.set_paq_values(i,j,k,&b_paq);
        // Additionally set KD (killing vector "Delta" quantities)
        bssnSim.set_KillingDelta(i, j, k, &b_paq);
      }
      io_show_progress(s, steps);
      io_data_dump(bssnSim.fields, staticSim.fields, &iodata, s, &fourier);
    _timer["output"].stop();

    // Run RK steps explicitly here (ties together BSSN + Hydro stuff).
    // See bssn class or hydro class for more comments.
    _timer["RK_steps"].start();


    // Pre-set detgamma and set gammai in _a register
    #pragma omp parallel for default(shared) private(i, j, k)
    LOOP3(i, j, k) {
      bssnSim.set_detgamma(i,j,k);
      bssnSim.set_gammai_values(i, j, k);
    }
    // First RK step, Set Hydro Vars, & calc. constraint
    #pragma omp parallel for default(shared) private(i, j, k)
    LOOP3(i, j, k)
    {
      BSSNData b_paq = {0}; // data structure associated with bssn sim
      bssnSim.K1CalcPt(i, j, k, &b_paq);
    }

    // Intermediate RK step is now in _c register, move to _a for use in next step.
    bssnSim.regSwap_c_a();
    // reset source using new metric
    bssnSim.clearSrc();
    staticSim.addBSSNSrc(bssnSim.fields);
    lambdaSim.addBSSNSrc(bssnSim.fields);


    // Subsequent BSSN steps
      // Second RK step
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k) {
        bssnSim.set_detgamma(i,j,k);
        bssnSim.set_gammai_values(i, j, k);
      }
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k)
      {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.K2CalcPt(i, j, k, &b_paq);
      }

      // Intermediate RK step is now in _c register, move to _a for use in next step.
      bssnSim.regSwap_c_a();
      // reset source using new metric
      bssnSim.clearSrc();
      staticSim.addBSSNSrc(bssnSim.fields);
      lambdaSim.addBSSNSrc(bssnSim.fields);

      // Third RK step
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k) {
        bssnSim.set_detgamma(i,j,k);
        bssnSim.set_gammai_values(i, j, k);
      }
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k)
      {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.K3CalcPt(i, j, k, &b_paq);
      }

      // Intermediate RK step is now in _c register, move to _a for use in next step.
      bssnSim.regSwap_c_a();
      // reset source using new metric
      bssnSim.clearSrc();
      staticSim.addBSSNSrc(bssnSim.fields);
      lambdaSim.addBSSNSrc(bssnSim.fields);

      // Fourth RK step
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k) {
        bssnSim.set_detgamma(i,j,k);
        bssnSim.set_gammai_values(i, j, k);
      }
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k)
      {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.K4CalcPt(i, j, k, &b_paq);
      }

    // Wrap up
      // bssn _f <-> _p
      bssnSim.stepTerm();
      // "current" data is in the _p array.

    _timer["RK_steps"].stop();


    _timer["meta_output_interval"].start();

      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k) {
        idx_t idx = INDEX(i,j,k);
        real_t K_ref = frw_K0/2.0;
        // K starts negative, gets closer to 0 monatonically
        // check to see when simulation has crossed a reference K
        real_t K2 = bssnSim.fields["K_p"][idx]; // current step is in _p array
        real_t K1 = bssnSim.fields["K_f"][idx]; // previous step is in _f array after swap

        if( K_ref >= K1 && K_ref <= K2 )
        {
          // linear interpolation to determine "time" of K_ref
          real_t ref_time = (K_ref - K1) / (K2 - K1);
          // get phi value at this "time"
          real_t phi2 = bssnSim.fields["phi_p"][idx];
          real_t phi1 = bssnSim.fields["phi_f"][idx];
          real_t phi_at_ref = (phi2 - phi1)*ref_time + phi1;
          // store phi value on the slice
          bssnSim.fields["dk0_slice_phi_a"][idx] = phi_at_ref;
        }
      }

      if(s%iodata.meta_output_interval == 0)
      {
        idx_t isNaN = 0;

        // make sure clean data is in _a array
        // Set values at points
        bssnSim.stepInit(); // copy _p to _a
        bssnSim.clearSrc();
        staticSim.addBSSNSrc(bssnSim.fields);
        lambdaSim.addBSSNSrc(bssnSim.fields);
        #pragma omp parallel for default(shared) private(i, j, k)
        LOOP3(i, j, k) {
          bssnSim.set_detgamma(i,j,k);
          bssnSim.set_gammai_values(i, j, k);
        }
        #pragma omp parallel for default(shared) private(i, j, k)
        LOOP3(i, j, k) {
          BSSNData b_paq = {0}; // data structure associated with bssn sim
          bssnSim.set_paq_values(i,j,k,&b_paq);
        }

        // Constraint Violation Calculations
        real_t mean_hamiltonian_constraint = bssnSim.hamiltonianConstraintMagMean();
        real_t stdev_hamiltonian_constraint = bssnSim.hamiltonianConstraintMagStDev(mean_hamiltonian_constraint);
        real_t max_hamiltonian_constraint = bssnSim.hamiltonianConstraintMagMax();
        real_t mean_momentum_constraint = bssnSim.momentumConstraintMagMean();
        real_t stdev_momentum_constraint = bssnSim.momentumConstraintMagStDev(mean_momentum_constraint);
        real_t max_momentum_constraint = bssnSim.momentumConstraintMagMax();
        io_dump_data(mean_hamiltonian_constraint, &iodata, "H_violations");
        io_dump_data(stdev_hamiltonian_constraint, &iodata, "H_violations");
        io_dump_data(max_hamiltonian_constraint, &iodata, "H_violations");
        io_dump_data(mean_momentum_constraint, &iodata, "M_violations");
        io_dump_data(stdev_momentum_constraint, &iodata, "M_violations");
        io_dump_data(max_momentum_constraint, &iodata, "M_violations");

        // LOG(iodata.log,
        //   "\nFractional Mean / St. Dev Momentum constraint violation is: " <<
        //   mean_momentum_constraint << " / " << stdev_momentum_constraint <<
        //   "\nFractional Mean / St. Dev Hamiltonian constraint violation is: " <<
        //   mean_hamiltonian_constraint << " / " << stdev_hamiltonian_constraint <<
        //   "\n"
        // );

        if(numNaNs(bssnSim.fields["phi_a"]) > 0)
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
  LOG(iodata.log, "\nAverage conformal factor reached " << average(bssnSim.fields["phi_p"]) << "\n");
  LOG(iodata.log, "Ending simulation.\n");
  _timer["output"].stop();

  _timer["MAIN"].stop();

  LOG(iodata.log, endl << _timer << endl);

  return EXIT_SUCCESS;
}
