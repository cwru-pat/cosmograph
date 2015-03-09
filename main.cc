
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
  idx_t steps, slice_output_interval, grid_output_interval,
        spec_output_interval, meta_output_interval;
  idx_t i, j, k, s;

  real_t rho_K_matter, rho_K_lambda, rho_K_lambda_frac, peak_amplitude,
         peak_amplitude_frac, length_scale, ic_spec_cut_frac;
  real_t total_hamiltonian_constraint = 0;

  // read in config file
  if(argc != 2)
  {
    cout << "Error: please supply exactly one config filename as an argument.\n";
    return EXIT_FAILURE;
  }
  else
  {
    _config.parse(argv[1]);
    
    length_scale = (real_t) stold(_config["length_scale"]); // volume in hubble units
    rho_K_matter = 3.0/PI/8.0*pw2(length_scale/(N*dx)); // matter density satisfies FRW equation

    // power spectrum amplitude as a fraction of the density
    peak_amplitude_frac = (real_t) stold(_config["peak_amplitude_frac"]); // fluctuation amplitude
    peak_amplitude = rho_K_matter*peak_amplitude_frac;

    rho_K_lambda_frac = (real_t) stold(_config["rho_K_lambda_frac"]); // DE density
    rho_K_lambda = rho_K_lambda_frac*rho_K_matter;

    ic_spec_cut_frac = (real_t) stold(_config["ic_spec_cut_frac"]); // power spectrum cutoff parameter

    steps = stoi(_config["steps"]);
    slice_output_interval = stoi(_config["slice_output_interval"]);
    grid_output_interval = stoi(_config["grid_output_interval"]);
    meta_output_interval = stoi(_config["meta_output_interval"]);
    spec_output_interval = stoi(_config["spec_output_interval"]);
    omp_set_num_threads(stoi(_config["omp_num_threads"]));
  }

  // IO init
  IOData iodata;
  iodata.output_dir = _config["output_dir"]; // gets a '/' appended if needed
  io_init(&iodata);
  // save a copy of config.txt
  ifstream source(argv[1], ios::binary);
  ofstream dest(iodata.output_dir + "config.txt", ios::binary);
  dest << source.rdbuf();
  source.close();
  dest.close();

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
    /* (peak scale in hubble units) * (to pixel scale) */
    i_paq.peak_k = (1.0/0.07)*(length_scale/((real_t) N));
    i_paq.peak_amplitude = peak_amplitude; // figure out units here
    i_paq.ic_spec_cut = N*ic_spec_cut_frac; // cut spectrum off around p ~ ic_spec_cut
                                            // (max is p ~ sqrt(2.5)*N )

    // 1) Either "conformal" initial conditions:
    set_conformal_ICs(bssnSim.fields, hydroSim.fields,
        &fourier, &i_paq, rho_K_matter, rho_K_lambda);
    // 2) or "flat" initial conditions:
    // set_flat_dynamic_ICs(bssnSim.fields, hydroSim.fields,
    //     &fourier, &i_paq, rho_K_matter, rho_K_lambda);

  _timer["init"].stop();

  // evolve simulation
  std::cout << "Running Simulation...\n";
  _timer["loop"].start();
  for(s=0; s < steps; ++s) {

    // output simulation information
    _timer["output"].start();
    io_show_progress(s, steps);
    if(s%slice_output_interval == 0)
    {
      io_dump_2dslice(bssnSim.fields["K_p"], "K_slice." + to_string(s), &iodata);
      io_dump_2dslice(bssnSim.fields["phi_p"], "phi_slice." + to_string(s), &iodata);
      io_dump_2dslice(hydroSim.fields["UD_a"], "UD_slice."  + to_string(s), &iodata);
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
      io_dump_3dslice(bssnSim.fields["AijAij_a"],  "AijAij."   + to_string(s), &iodata);
      io_dump_3dslice(hydroSim.fields["UD_a"],     "UD."      + to_string(s), &iodata);
      io_dump_3dslice(hydroSim.fields["US1_a"],    "US1."     + to_string(s), &iodata);
      io_dump_3dslice(hydroSim.fields["US2_a"],    "US2."     + to_string(s), &iodata);
      io_dump_3dslice(hydroSim.fields["US3_a"],    "US3."     + to_string(s), &iodata);
    }
    if(s%spec_output_interval == 0)
    {
      fourier.powerDump(bssnSim.fields["phi_p"], &iodata);
    }
    if(s%meta_output_interval == 0)
    {
      // some average values
      io_dump_quantities(bssnSim.fields, hydroSim.fields, _config["dump_file"], &iodata);
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
      if(s%meta_output_interval == 0)
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
