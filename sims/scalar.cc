#include "scalar.h"

namespace cosmo
{

void ScalarSim::init()
{
  _timer["init"].start();

  // initialize base class
  simInit();

  iodata->log("Running 'scalar' type simulation.");
  if(!USE_HARMONIC_ALPHA) {
    iodata->log("Warning - not using harmonic gauge! You may want to use it.");
  }

  scalarSim = new Scalar();
  iodata->log("Creating initial conditions.");
  setScalarMultigridICs();
  iodata->log("Finished setting ICs.");

  _timer["init"].stop();
}

/**
 * @brief small-amplitude wave on a flat metric
 *  use a stationary gaussian wave packet for now
 */
void ScalarSim::setScalarWaveICs()
{
  // BSSN is already initialized to flat, just initialize scalar fields
  arr_t & phi = scalarSim->phi._array_p; // gaussian 
  arr_t & Pi = scalarSim->Pi._array_p; // Pi = 0
  arr_t & psi1 = scalarSim->psi1._array_p; // derivative of phi in x-dir
  arr_t & psi2 = scalarSim->psi3._array_p; // derivative of phi in y-dir
  arr_t & psi3 = scalarSim->psi2._array_p; // derivative of phi in z-dir

  arr_t & K_p = *bssnSim->fields["DIFFK_p"]; // extrinsic curvature
  arr_t & K_a = *bssnSim->fields["DIFFK_a"]; // extrinsic curvature

  // iterators
  idx_t i, j, k;
  // gaussian parameters
  real_t amplitude = 1.0e-10;
  real_t sigx = (real_t) NX / 10.0;
  real_t sigy = (real_t) NY / 10.0;
  real_t sigz = (real_t) NZ / 10.0;

  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    phi[INDEX(i,j,k)] = amplitude*exp(
        -pw2(i - ((real_t) NX-1)/2.0)/pw2(sigx)
        // plane-wave only, so commented:
        // -pw2(j - ((real_t) NY-1)/2.0)/pw2(sigy)
        // -pw2(k - ((real_t) NZ-1)/2.0)/pw2(sigz)
      );
  }

  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    psi1[INDEX(i,j,k)] = derivative(i, j, k, 1, phi);
    psi2[INDEX(i,j,k)] = derivative(i, j, k, 2, phi);
    psi3[INDEX(i,j,k)] = derivative(i, j, k, 3, phi);
  }

  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    K_a[INDEX(i,j,k)] = K_p[INDEX(i,j,k)] = -std::sqrt(24*PI*scalarSim->V(phi[INDEX(i,j,k)]));
  }

  return;
}

/**
 * @brief Use the multigrid solver to solve for metric factors given
 * a particular scalar field implementation.
 */
void ScalarSim::setScalarMultigridICs()
{
  idx_t i, j, k;

  // Choose a configuration for the scalar fields first:
  arr_t & phi = scalarSim->phi._array_p; // field
  arr_t & psi1 = scalarSim->psi1._array_p; // derivative of phi in x-dir
  arr_t & psi2 = scalarSim->psi3._array_p; // derivative of phi in y-dir
  arr_t & psi3 = scalarSim->psi2._array_p; // derivative of phi in z-dir

  arr_t & Pi = scalarSim->Pi._array_p; // time-derivative of field phi

  std::random_device rd;
  std::mt19937 gen(7.0 /*rd()*/);
  std::uniform_real_distribution<real_t> dist(0, 2.0*PI);

  // cutoff @ "ic_spec_cut"; maybe initialize this field
  // according to some power spectrum?
  real_t n_max = std::stoi(_config["n_max"]);
  real_t phi_0 = std::stod(_config["phi_0"]);
  real_t delta = std::stod(_config["delta_phi"]);

  // background value
  LOOP3(i,j,k)
    phi[INDEX(i,j,k)] = phi_0;

  // sum over different modes
  for(int n = -n_max; n <= n_max; ++n)
  {
    if(n != 0)
    {
      // random phases
      real_t x_phase = dist(gen),
             y_phase = dist(gen),
             z_phase = dist(gen);

      #pragma omp parallel for default(shared) private(i,j,k)
      LOOP3(i,j,k)
      {
        // some simusoid modes
        phi[INDEX(i,j,k)] += delta*(
                cos(2.0*PI*((real_t) n/NX)*i + x_phase )
                 + cos(2.0*PI*((real_t) n/NY)*j + y_phase )
                 + cos(2.0*PI*((real_t) n/NZ)*k + z_phase )
            );
      }
    }
  }

  // initialize psi according to values in phi
  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    psi1[INDEX(i,j,k)] = derivative(i, j, k, 1, phi);
    psi2[INDEX(i,j,k)] = derivative(i, j, k, 2, phi);
    psi3[INDEX(i,j,k)] = derivative(i, j, k, 3, phi);
  }

  // PI is zero for now


  // compute background/average K
  real_t K_src = 0;
  LOOP3(i, j, k)
  {
    K_src += 3.0*scalarSim->V(phi[INDEX(i,j,k)])
      +3.0/2.0*(
        pw2(psi1[INDEX(i,j,k)]) + pw2(psi2[INDEX(i,j,k)]) + pw2(psi3[INDEX(i,j,k)])
      );
  }
  K_src = -std::sqrt(K_src/NX/NY/NZ);

  arr_t & K_p = *bssnSim->fields["DIFFK_p"]; // extrinsic curvature
  arr_t & K_a = *bssnSim->fields["DIFFK_a"]; // extrinsic curvature

  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    K_a[INDEX(i,j,k)] = K_p[INDEX(i,j,k)] = K_src;
  }


  // solve for BSSN fields using multigrid class:
  FASMultigrid multigrid (N, N*dx, 4);

  idx_t u_exp[2] = { 1, 5 };
  multigrid.build_rho(2, u_exp);
  LOOP3(i, j, k)
  {
    real_t value = PI*(
        pw2(psi1[INDEX(i,j,k)]) + pw2(psi2[INDEX(i,j,k)]) + pw2(psi3[INDEX(i,j,k)])
      );
    multigrid.setPolySrcAtPt(i, j, k, 0, value);

    value = PI*scalarSim->V(phi[INDEX(i,j,k)]) - K_a[INDEX(i,j,k)]*K_a[INDEX(i,j,k)]/12.0;
    multigrid.setPolySrcAtPt(i, j, k, 1, value);
  }
  multigrid.initializeRhoHeirarchy();

  if(std::stoi(_config["debug_multigrid"]))
  {
    multigrid.printSourceStrip(0, 5);
    multigrid.printSourceStrip(1, 5);
  }

  // set initial guess and solve using 3 V-cycles
  multigrid.setTrialSolution(0);
  multigrid.VCycles(std::stoi(_config["num_v_cycles"]));

  if(std::stoi(_config["debug_multigrid"]))
  {
    multigrid.printSolutionStrip(1);
    multigrid.printSolutionStrip(5);
  }

  // set bssn phi using multigrid solution
  arr_t & phi_p = *bssnSim->fields["DIFFphi_p"];
  arr_t & phi_a = *bssnSim->fields["DIFFphi_a"];

  real_t * u = multigrid.getSolution();
  
  #pragma omp parallel for
  LOOP3(i, j, k)
  {
    phi_p[INDEX(i,j,k)] = std::log(std::abs(u[INDEX(i,j,k)]));
    phi_a[INDEX(i,j,k)] = std::log(std::abs(u[INDEX(i,j,k)]));
  }

  // TODO: Why doesn't this work?
  // multigrid.~FASMultigrid();

  std::cout << "Finished setting ICs!\n" << std::flush;

  return;
}

void ScalarSim::initScalarStep()
{
  _timer["RK_steps"].start();
    bssnSim->stepInit();
    scalarSim->stepInit();
    bssnSim->clearSrc();
    scalarSim->addBSSNSource(bssnSim);
  _timer["RK_steps"].stop();
}

void ScalarSim::outputScalarStep()
{
  _timer["output"].start();
    prepBSSNOutput();
    io_bssn_fields_snapshot(iodata, step, bssnSim->fields);
    io_bssn_fields_powerdump(iodata, step, bssnSim->fields, fourier);
    io_bssn_dump_statistics(iodata, step, bssnSim->fields, bssnSim->frw);
    io_bssn_constraint_violation(iodata, step, bssnSim);
    io_scalar_snapshot(iodata, step, scalarSim);
  _timer["output"].stop();
}

void ScalarSim::runScalarStep()
{
  idx_t i=0, j=0, k=0;
  BSSNData b_data;
  _timer["RK_steps"].start();

    // First RK step
    #pragma omp parallel for default(shared) private(i, j, k, b_data)
    LOOP3(i,j,k)
    {
      bssnSim->RKEvolvePt(i, j, k, &b_data);
      scalarSim->RKEvolvePt(&b_data);
    }
    bssnSim->K1Finalize();
    scalarSim->K1Finalize();

    // Second RK step
    bssnSim->clearSrc();
    scalarSim->addBSSNSource(bssnSim);
    #pragma omp parallel for default(shared) private(i, j, k, b_data)
    LOOP3(i,j,k)
    {
      bssnSim->RKEvolvePt(i, j, k, &b_data);
      scalarSim->RKEvolvePt(&b_data);
    }
    bssnSim->K2Finalize();
    scalarSim->K2Finalize();

    // Third RK step
    bssnSim->clearSrc();
    scalarSim->addBSSNSource(bssnSim);
    #pragma omp parallel for default(shared) private(i, j, k, b_data)
    LOOP3(i,j,k)
    {
      bssnSim->RKEvolvePt(i, j, k, &b_data);
      scalarSim->RKEvolvePt(&b_data);
    }
    bssnSim->K3Finalize();
    scalarSim->K3Finalize();

    // Fourth RK step
    bssnSim->clearSrc();
    scalarSim->addBSSNSource(bssnSim);
    #pragma omp parallel for default(shared) private(i, j, k, b_data)
    LOOP3(i,j,k)
    {
      bssnSim->RKEvolvePt(i, j, k, &b_data);
      scalarSim->RKEvolvePt(&b_data);
    }
    bssnSim->K4Finalize();
    scalarSim->K4Finalize();

    // "current" data should be in the _p array.
  _timer["RK_steps"].stop();
}

void ScalarSim::runStep()
{
  runCommonStepTasks();
  
  initScalarStep();
  outputScalarStep();
  runScalarStep();
}

} /* namespace cosmo */
