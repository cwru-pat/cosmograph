#ifndef COSMO_SIM_H
#define COSMO_SIM_H

#include "cosmo_includes.h"
#include "cosmo_types.h"
#include "globals.h"

#include "utils/Fourier.h"
#include "utils/reference_frw.h"

#include "ICs.h"
#include "io.h"
#include "bssn.h"
#include "static.h"

namespace cosmo
{

class CosmoSim
{
  idx_t step;
  idx_t num_steps;
  bool dt_flip;
  idx_t dt_flip_step;

  std::string simulation_type;
  IOData iodata;
  Fourier * fourier;
  
  BSSN * bssnSim;
  Static * staticSim;
  Particles * particles;

  int verbosity;

  bool ray_integrate;
  idx_t ray_flip_step;
  std::vector<RayTrace<real_t, idx_t> *> rays;

public:
  CosmoSim()
  {
    // Initialize iodata struct first
    io_init(&iodata, _config["output_dir"]);
    // save a copy of config.txt
    io_config_backup(&iodata, _config.getFileName());

    // fix number of simulation steps
    step = 0;
    num_steps = stoi(_config["steps"]);

    // integrating any light rays?
    if( stoi(_config["ray_integrate"]) == 1 )
    {
      ray_integrate = true;
      ray_flip_step = stoi(_config["ray_flip_step"]);
    }
    else
    {
      ray_integrate = false;
    }

    // Store simulation type
    simulation_type = _config["simulation_type"];
  }

  ~CosmoSim()
  {
    LOG(iodata.log, std::endl << _timer << std::endl);
  }

  /**
   * @brief      Initialize individual simulation class instances
   */
  void init()
  {
    _timer["init"].start();

    // Initialize simulation

    // GR fields
    bssnSim = new BSSN();
    bssnSim->init();

    // FFT helper
    fourier = new Fourier();
    fourier->Initialize(NX, NY, NZ,
      bssnSim->fields["DIFFphi_a"]->_array /* arbitrary array for planning */);

    // set ICs according to simulation_type
    if( simulation_type == "dust" )
    {
      LOG(iodata.log, "Running 'dust' type simulation.\n");
      staticSim = new Static();
      staticSim->init();
      LOG(iodata.log, "Creating initial conditions.\n");
      ICs_set_dust(bssnSim->fields, staticSim->fields, fourier, &iodata, bssnSim->frw);
    }
    else if( simulation_type == "particles" )
    {
      LOG(iodata.log, "Running 'particles' type simulation.\n");
      particles = new Particles();
      ICs_set_particle(particles, bssnSim->fields, fourier, &iodata);
    }
    else if( simulation_type == "vacuum" )
    {
      LOG(iodata.log, "Running 'vacuum' type simulation.\n");
// TODO: Set vacuum ICs (eg, AwA test)
      ICs_set_vacuum(bssnSim->fields, &iodata);
    }
    else
    {
      LOG(iodata.log, "Invalid simulation type specified.\n");
      throw 2;
    }

    if(ray_integrate)
    {
      int rays_num = std::stoi(_config["rays_num"]);
      init_ray_vector(&rays, rays_num);
    }

    _timer["init"].stop();
  }

  /**
   * @brief      Run the simulation.
   */
  void run()
  {
    LOG(iodata.log, "Running simulation...\n");

    _timer["loop"].start();
    while(step <= num_steps)
    {
      runStep();
      step++;
    }
    _timer["loop"].stop();

    LOG(iodata.log, "\nEnding simulation.\n");
    LOG(iodata.log, "Average conformal factor reached "
      << average(*bssnSim->fields["DIFFphi_p"]) << "\n");
  }

  /**
   * @brief Run a step of a simulation.
   * @details Run an RK4 step of a simulation, perform any I/O or reductions
   *  required.
   *  
   *  Schematic writeup of Particle (p_) and bssn (b_) RK4 computations:
   *  (r_ variable is bssn src)
   *  
   *  step Init:
   *  b_p; p_f = 0;
   *  b_a = b_p;
   *  p_p; p_a = p_c = p_p; p_f = 0
   *  r_a = r(p_a)
   *  
   *  (Output content in _a registers)
   *  
   *  RK1 step:
   *  b_c = b_p + dt/2 * f(b_a, r_a);
   *  b_f += b_c
   *  b_c <-> b_a;
   *  r_a = 0;
   *  p_c = p_p + dt/2 * f( p_a, b_c )
   *  p_f += p_c
   *  p_a <-> p_c
   *  r_a = r(p_a)
   *  
   *  RK2 step:
   *  b_c = b_p + dt/2 * f(b_a, r_a);
   *  b_f += 2 b_c
   *  b_c <-> b_a
   *  r_a = 0
   *  p_c = p_p + dt/2 * f( p_a, b_c )
   *  p_f += 2 p_c
   *  p_a <-> p_c
   *  r_a = r(p_a)
   *  
   *  RK3 step:
   *  b_c = b_p + dt * f( b_a, r_a )
   *  b_f += b_c
   *  b_c <-> b_a
   *  r_a = 0
   *  p_c = p_p + dt * f( p_a, b_c )
   *  p_f += p_c
   *  p_a <-> p_c
   *  r_a = r(p_a)
   *  
   *  RK4 step:
   *  b_f = 1/3 * (b_f - b_p) + dt/6 * f(b_a)
   *  p_c = p_p + dt/2 * f( p_a, b_c )
   *  p_f += p_c
   *  p_a <-> p_c
   *  
   *  Finalize:
   *  b_f <-> b_p
   *  (p_f = 5p_p + 1/2 K1 + K2 + K3 + 1/2 K4)
   *  p_f = p_f / 3 - 2/3 p_p
   *  p_f <-> p_p
   *  
   *  _p registers now contain "correct"
   */
  void runStep()
  {
    // check for NAN every step
    if(simNumNaNs() > 0)
    {
      LOG(iodata.log, "\nNAN detected!\n");
      throw 10;
    }

    // progress bar in terminal
    io_show_progress(step, num_steps);

    // Evolve light rays when integrating backwards
    if(ray_integrate)
    {
      if(step == ray_flip_step) {
        LOG(iodata.log, "\nFlipping sign of dt @ step = " << step << "\n");
        dt = -dt;
      }
      if(step >= ray_flip_step) {
        outputRayTraceStep();
        runRayTraceStep();
      }
    }

    // evolve BSSN + matter sources
    // call simulation-specific output routines
    if( simulation_type == "dust" )
    {
      initDustStep();
      outputDustStep();
      runDustStep();
    }
    else if( simulation_type == "particles" )
    {
      initParticleStep();
      outputParticleStep();
      runParticleStep();
    }
    else if( simulation_type == "scalar" )
    {

    }
    else if( simulation_type == "vacuum" )
    {
      initVacuumStep();
      outputVacuumStep();
      runVacuumStep();
    }

  }

  void initParticleStep()
  {
    _timer["RK_steps"].start();
      bssnSim->stepInit();
      bssnSim->clearSrc();
      particles->stepInit(bssnSim->fields);
    _timer["RK_steps"].stop();
  }

  void outputParticleStep()
  {
    _timer["output"].start();
      prepBSSNOutput();
      io_bssn_fields_snapshot(&iodata, step, bssnSim->fields);
      io_bssn_fields_powerdump(&iodata, step, bssnSim->fields, fourier);
      io_bssn_dump_statistics(&iodata, step, bssnSim->fields, bssnSim->frw);
      io_bssn_constraint_violation(&iodata, step, bssnSim);
    _timer["output"].stop();
  }

  void runParticleStep()
  {
    _timer["RK_steps"].start();
      // First RK step
      bssnSim->K1Calc();
      bssnSim->clearSrc();
      particles->RK1Step(bssnSim->fields);
      bssnSim->regSwap_c_a();

      // Second RK step
      bssnSim->K2Calc();
      bssnSim->clearSrc();
      particles->RK2Step(bssnSim->fields);
      bssnSim->regSwap_c_a();

      // Third RK step
      bssnSim->K3Calc();
      bssnSim->clearSrc();
      particles->RK3Step(bssnSim->fields);
      bssnSim->regSwap_c_a();

      // Fourth RK step
      bssnSim->K4Calc();
      particles->RK4Step(bssnSim->fields);

      // Wrap up
      bssnSim->stepTerm();
      particles->stepTerm();
      // "current" data should be in the _p array.
    _timer["RK_steps"].stop();
  }

  void initDustStep()
  {
    COSMOSIM_COUT << "Initializing dust step... " << std::flush;
    _timer["RK_steps"].start();
      bssnSim->stepInit();
      bssnSim->clearSrc();
      staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
    _timer["RK_steps"].stop();
    COSMOSIM_COUT << "done.\n";
  }

  void outputDustStep()
  {
    COSMOSIM_COUT << "Output dust step... " << std::flush;
    _timer["output"].start();
      prepBSSNOutput();
      io_bssn_fields_snapshot(&iodata, step, bssnSim->fields);
      io_bssn_fields_powerdump(&iodata, step, bssnSim->fields, fourier);
      io_bssn_dump_statistics(&iodata, step, bssnSim->fields, bssnSim->frw);
      io_bssn_constraint_violation(&iodata, step, bssnSim);
    _timer["output"].stop();
    COSMOSIM_COUT << "done.\n";
  }

  void runDustStep()
  {
    COSMOSIM_COUT << "Running dust step... " << std::flush;
    _timer["RK_steps"].start();
      // First RK step
      COSMOSIM_COUT << "RK1; " << std::flush;
      bssnSim->K1Calc();
      bssnSim->clearSrc();
      staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
      bssnSim->regSwap_c_a();

      // Second RK step
      COSMOSIM_COUT << "RK2; " << std::flush;
      bssnSim->K2Calc();
      bssnSim->clearSrc();
      staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
      bssnSim->regSwap_c_a();

      // Third RK step
      COSMOSIM_COUT << "RK3; " << std::flush;
      bssnSim->K3Calc();
      bssnSim->clearSrc();
      staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
      bssnSim->regSwap_c_a();

      // Fourth RK step
      COSMOSIM_COUT << "RK4; " << std::flush;
      bssnSim->K4Calc();

      // Wrap up
      bssnSim->stepTerm();
      // "current" data should be in the _p array.
    _timer["RK_steps"].stop();
    COSMOSIM_COUT << "done.\n";
  }

  void initVacuumStep()
  {
    _timer["RK_steps"].start();
      bssnSim->stepInit();
    _timer["RK_steps"].stop();
  }

  void outputVacuumStep()
  {
    _timer["output"].start();
      prepBSSNOutput();
      io_bssn_fields_snapshot(&iodata, step, bssnSim->fields);
      io_bssn_fields_powerdump(&iodata, step, bssnSim->fields, fourier);
      io_bssn_dump_statistics(&iodata, step, bssnSim->fields, bssnSim->frw);
      io_bssn_constraint_violation(&iodata, step, bssnSim);
    _timer["output"].stop();
  }

  void runVacuumStep()
  {
    _timer["RK_steps"].start();
      // Full RK step minus init()
      bssnSim->step();
    _timer["RK_steps"].stop();
  }

  void runRayTraceStep()
  {
    // evolve any light rays
    _timer["Raytrace_step"].start();
    auto ray = rays.begin();

    #pragma omp parallel for default(shared) private(ray)
    for(ray = rays.begin(); ray < rays.end(); ++ray)
    {
      // set primitives from BSSN sim
      bssnSim->setRaytracePrimitives(*ray);
      // evolve ray
      (*ray)->setDerivedQuantities();
      (*ray)->evolveRay();
    }    
    _timer["Raytrace_step"].stop();
  }

  void outputRayTraceStep()
  {
    _timer["output"].start();
    io_raytrace_dump(&iodata, step, &rays);
    _timer["output"].stop();
  }

  void prepBSSNOutput()
  {
    idx_t i, j, k;

    #pragma omp parallel for default(shared) private(i, j, k)
    LOOP3(i,j,k)
    {
      // set_paq_values calculates ricci_a and AijAij_a data, needed for output
      // and potentially subsequent Killing calculations
      BSSNData b_paq = {0}; // data structure associated with bssn sim
      bssnSim->set_paq_values(i, j, k, &b_paq);
      
      // Additionally set KD (killing vector "Delta" quantities)
      bssnSim->set_KillingDelta(i, j, k, &b_paq);
    }
  }

  idx_t simNumNaNs()
  {
    // check for NAN in a field
    return numNaNs(*bssnSim->fields["DIFFphi_a"]);
  }

  void setVerbosity(int verbosity_in)
  {
    verbosity = verbosity_in;
  }

};

} /* namespace cosmo */

#endif
