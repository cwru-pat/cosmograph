
#include "main.h"

using namespace std;
using namespace cosmo;

/* global definitions */
TimerManager _timer;
ConfigParser _config;
#ifndef dt
  real_t dt;
#endif
#ifndef dx
  real_t dx;
#endif

int main(int argc, char **argv)
{
  _timer["MAIN"].start();
  idx_t i=0, j=0, k=0, s=0, steps=0;

  // If not compiled in, set dt, dx (TODO: set in config file?)
  #ifndef dt
    dt = 0.1*dx;
  #endif
  #ifndef dx
    dx = H_LEN_FRAC/(1.0*N);
  #endif

  // read in config file
  if(argc != 2)
  {
    cout << "Error: please supply exactly one config filename as an argument.\n";
    return EXIT_FAILURE;
  }
  _config.parse(argv[1]);
  steps = stoi(_config["steps"]);
  omp_set_num_threads(stoi(_config["omp_num_threads"]));

  // IO init - will use this for logging.
  IOData iodata;
  io_init(&iodata, _config["output_dir"]);
  // save a copy of config.txt
  io_config_backup(&iodata, argv[1]);

  // Create simulation
  _timer["init"].start();
    LOG(iodata.log, "Creating initial conditions...\n");

    // Static matter (w=0)
    Static staticSim;
    staticSim.init();

    // GR Fields
    BSSN bssnSim;
    bssnSim.init();

    // Generic reusable fourier class for NX*NY*NZ arrays
    Fourier fourier;
    fourier.Initialize(NX, NY, NZ,
      staticSim.fields["DIFFD_a"] /* just any array for planning */);

    // vector of rays (rayracing code)
    std::vector<RayTrace<real_t, idx_t> *> rays;
    init_ray_vector(&rays, 1000);

    // set initial conditions
    set_ICs(bssnSim.fields, staticSim.fields, &fourier, &iodata, bssnSim.frw);

    // Particles
    Particles particles;
    // init_particle_vector(&particles, staticSim.fields, bssnSim.fields);

  // evolve simulation
  LOG(iodata.log, "Running simulation...\n");
  _timer["loop"].start();
  for(s=0; s < steps; ++s) {
    // check for NAN in a field every step
    if(numNaNs(bssnSim.fields["DIFFphi_a"]) > 0)
    {
      LOG(iodata.log, "\nNAN detected!\n");
      break;
    }

    // initialize data for RK step in first loop
    // Init arrays and calculate source term for next step
      // _p is copied to _a here (which matter sectors reference)
      bssnSim.stepInit();
      // reset source data
      bssnSim.clearSrc();
      staticSim.addBSSNSrc(bssnSim.fields, bssnSim.frw);
      particles.stepInit(bssnSim.fields);

    // output simulation information
    // these generally output any data in the _a registers (which should 
    // be identical to _p at this point).
    _timer["output"].start();
      call_io_routines(&bssnSim, &staticSim, &iodata, s, steps,
                       &fourier, bssnSim.frw, &rays);
    _timer["output"].stop();

    // Run RK steps explicitly here (ties together BSSN + Hydro stuff).
    // See bssn class or hydro class for more comments.
    _timer["RK_steps"].start();
      // FRW simulation & particles should be in a correct state here

      /**
       * Schematic writeup of Particle (p_) and bssn (b_) RK4 computations.
       * (r_ variable is bssn src)
       * 
       * step Init:
       * b_p; p_f = 0;
       * b_a = b_p;
       * p_p; p_a = p_c = p_p; p_f = 0
       * r_a = r(p_a)
       * 
       * RK1 step:
       * b_c = b_p + dt/2 * f(b_a, r_a);
       * b_f += b_c
       * b_c <-> b_a;
       * r_a = 0;
       * p_c = p_p + dt/2 * f( p_a, b_c )
       * p_f += p_c
       * r_a = r(p_c)
       * p_a <-> p_c
       * 
       * RK2 step:
       * b_c = b_p + dt/2 * f(b_a, r_a);
       * b_f += 2 b_c
       * b_c <-> b_a
       * r_a = 0
       * p_c = p_p + dt/2 * f( p_a, b_c )
       * p_f += 2 p_c
       * r_a = r(p_c)
       * p_a <-> p_c
       * 
       * RK3 step:
       * b_c = b_p + dt * f( b_a, r_a )
       * b_f += b_c
       * b_c <-> b_a
       * r_a = 0
       * p_c = p_p + dt * f( p_a, b_c )
       * p_f += p_c
       * r_a = r(p_c)
       * p_a <-> p_c
       * 
       * RK4 step:
       * b_f = 1/3 * (b_f - b_p) + dt/6 * f(b_a)
       * TODO: need to set b_c for vv
       * p_c = p_p + dt/2 * f( p_a, b_c )
       * p_f += p_c
       * p_a <-> p_c
       * 
       * (p_f = 5p_p + 1/2 K1 + K2 + K3 + 1/2 K4)
       * p_f = p_f / 3 - 2/3 p_p      // finalize
       * p_f <-> p_p
       **/

      // First RK step, Set Hydro Vars, & calc. constraint
      bssnSim.K1Calc();
      // reset source using new metric
      bssnSim.clearSrc();
      staticSim.addBSSNSrc(bssnSim.fields, bssnSim.frw);
      particles.RK1Step(bssnSim.fields);

      // Second RK step
      bssnSim.K2Calc();
      // reset source using new metric
      bssnSim.clearSrc();
      staticSim.addBSSNSrc(bssnSim.fields, bssnSim.frw);
      particles.RK2Step(bssnSim.fields);

      // Third RK step
      bssnSim.K3Calc();
      // reset source using new metric
      bssnSim.clearSrc();
      staticSim.addBSSNSrc(bssnSim.fields, bssnSim.frw);
      particles.RK3Step(bssnSim.fields);

      // Fourth RK step
      bssnSim.K4Calc();
      particles.RK4Step(bssnSim.fields);

      // Wrap up
        // bssn _f <-> _p
        bssnSim.stepTerm();
        particles.stepTerm();
        // "current" data is in the _p array.
    _timer["RK_steps"].stop();


    // evolve any light rays
    _timer["Raytrace_step"].start();
      for(RayTrace<real_t, idx_t> * ray : rays)
      {
        // set primitives from BSSN sim
        bssnSim.setRaytracePrimitives(ray);
        // evolve ray
        ray->setDerivedQuantities();
        ray->evolveRay();
      }
    _timer["Raytrace_step"].stop();      

  }
  _timer["loop"].stop();

  _timer["output"].start();
    LOG(iodata.log, "\nAverage conformal factor reached " << average(bssnSim.fields["DIFFphi_p"]) << "\n");
    LOG(iodata.log, "Ending simulation.\n");
  _timer["output"].stop();

  _timer["MAIN"].stop();

  LOG(iodata.log, endl << _timer << endl);

  return EXIT_SUCCESS;
}


/**
 * @brief Output simulation information
 * @details Output various quantities
 */
void cosmo::call_io_routines(BSSN * bssnSim, Static * staticSim,
                      IOData *iodata, idx_t step, idx_t steps,
                      Fourier *fourier, FRW<real_t> *frw,
                      std::vector<RayTrace<real_t, idx_t> *> * rays)
{
  idx_t i, j, k;

  // set_paq_values calculates ricci_a and AijAij_a data, needed for output
  // and subsequent constraint calculations
  #pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    BSSNData b_paq = {0}; // data structure associated with bssn sim
    bssnSim->set_paq_values(i, j, k, &b_paq);
    // Additionally set KD (killing vector "Delta" quantities)
    bssnSim->set_KillingDelta(i, j, k, &b_paq);
  }

  // dump field snapshots
  io_fields_snapshot(bssnSim->fields, staticSim->fields, iodata, step, fourier, bssnSim->frw);

  // power spectra
  if(step % iodata->spec_output_interval == 0)
  {
    fourier->powerDump(bssnSim->fields["DIFFphi_a"], iodata);
    fourier->powerDump(bssnSim->fields["DIFFr_a"], iodata);
  }

  // compute & output additional information  
  _timer["meta_output_interval"].start();
    if(step % iodata->meta_output_interval == 0)
    {
      idx_t isNaN = 0;
      real_t H_calcs[7], M_calcs[7];

      // Constraint Violation Calculations
      bssnSim->setHamiltonianConstraintCalcs(H_calcs, false);
      io_dump_value(H_calcs[4], iodata, "H_violations"); // mean(H/[H])
      io_dump_value(H_calcs[5], iodata, "H_violations"); // stdev(H/[H])
      io_dump_value(H_calcs[6], iodata, "H_violations"); // max(H/[H])
      io_dump_value(H_calcs[2], iodata, "H_violations"); // max(H)

      bssnSim->setMomentumConstraintCalcs(M_calcs);
      io_dump_value(M_calcs[4], iodata, "M_violations"); // mean(M/[M])
      io_dump_value(M_calcs[5], iodata, "M_violations"); // stdev(M/[M])
      io_dump_value(M_calcs[6], iodata, "M_violations"); // max(M/[M])
      io_dump_value(M_calcs[2], iodata, "M_violations"); // max(M)

      // g_11 along a 1-d strip
      io_dump_strip(bssnSim->fields["DIFFgamma11_a"], 1, 0, 0, iodata);

      // statistical information
      io_dump_statistics(bssnSim->fields, staticSim->fields, iodata, frw);
    }
  _timer["meta_output_interval"].stop();

  // rayrace information
  RaytraceData<real_t> tmp_rd = {0};
  for(RayTrace<real_t, idx_t> * ray : *rays)
  {
    tmp_rd = ray->getRaytraceData();
    io_dump_value(tmp_rd.E, iodata, "ray_functions");
    io_dump_value(tmp_rd.x[0], iodata, "ray_functions");
    io_dump_value(tmp_rd.x[1], iodata, "ray_functions");
    io_dump_value(tmp_rd.x[2], iodata, "ray_functions");
    io_dump_value(ray->RicciLensingScalarSum(), iodata, "ray_functions");
    io_dump_value(ray->WeylLensingScalarSum_Re(), iodata, "ray_functions");
    io_dump_value(ray->WeylLensingScalarSum_Im(), iodata, "ray_functions");
  }

  io_show_progress(step, steps);
}
