#include "particles.h"

namespace cosmo
{

ParticleSim::ParticleSim()
{
  if(USE_REFERENCE_FRW)
  {
    iodata->log("Error: USE_REFERENCE_FRW must be false for particle sims!");
  }
}

void ParticleSim::init()
{
  _timer["init"].start();

  // initialize base class
  simInit();

  iodata->log("Running 'particles' type simulation.");
  particles = new Particles();
  setParticleICs();

  _timer["init"].stop();
}

/**
 * @brief Initialize particles from gaussian random field data
 * @details Initialize particles from gaussian random field data, so particle
 *  masses are that needed to recreate the corresponding density field.
 * 
 * @param particles reference to Particles class containing particles
 * @param bssn_fields map to fields from bssn class
 * @param      fourier       { parameter_description }
 * @param      iod           { parameter_description }
 */
void ParticleSim::setParticleICs()
{
  idx_t i, j, k;
  ICsData icd = cosmo_get_ICsData();
  real_t rho_FRW = icd.rho_K_matter;
  arr_t & DIFFr = *bssnSim->fields["DIFFr_a"];
  arr_t & DIFFphi_p = *bssnSim->fields["DIFFphi_p"];
  arr_t & DIFFphi_a = *bssnSim->fields["DIFFphi_a"];
  arr_t & DIFFphi_c = *bssnSim->fields["DIFFphi_c"];
  arr_t & DIFFphi_f = *bssnSim->fields["DIFFphi_f"];
  iodata->log( "Generating ICs with peak at k = " + stringify(icd.peak_k) );
  iodata->log( "Generating ICs with peak amp. = " + stringify(icd.peak_amplitude) );

  // the conformal factor in front of metric is the solution to
  // d^2 exp(\phi) = -2*pi exp(5\phi) * \rho
  // generate gaussian random field 1 + xi = exp(phi) (use phi_p as a proxy):
  set_gaussian_random_field(DIFFphi_p._array, fourier, &icd);

  // rho = -lap(phi)/xi^5/2pi
  LOOP3(i, j, k) {
    Particle<real_t> particle = {0};

    // Particle position in grid
    particle.X[0] = i*dx;
    particle.X[1] = j*dx;
    particle.X[2] = k*dx;

    // Particle mass
    DIFFr[NP_INDEX(i,j,k)] = rho_FRW - 0.5/PI/(
      pow(1.0 + DIFFphi_p[NP_INDEX(i,j,k)], 5.0)
    )*(
      double_derivative(i, j, k, 1, 1, DIFFphi_p)
      + double_derivative(i, j, k, 2, 2, DIFFphi_p)
      + double_derivative(i, j, k, 3, 3, DIFFphi_p)
    );
    real_t rho = rho_FRW + DIFFr[INDEX(i,j,k)];
    real_t rootdetg = std::exp(6.0*DIFFphi_p[INDEX(i,j,k)]);
    particle.M = rho*dx*dx*dx*rootdetg;

    particles->addParticle( particle );
  }

  // phi = ln(xi)
  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k) {
    idx_t idx = NP_INDEX(i,j,k);
    DIFFphi_a[idx] = log1p(DIFFphi_p[idx]);
    DIFFphi_f[idx] = log1p(DIFFphi_p[idx]);
    DIFFphi_p[idx] = log1p(DIFFphi_p[idx]);
  }

  // Make sure min density value > 0
  // Set conserved density variable field
  real_t min = icd.rho_K_matter;
  real_t max = min;
  LOOP3(i,j,k)
  {
    real_t rho = DIFFr[NP_INDEX(i,j,k)];

    if(rho < min)
    {
      min = rho;
    }
    if(rho > max)
    {
      max = rho;
    }
    if(rho != rho)
    {
      iodata->log("Error: NaN energy density.");
      throw -1;
    }
  }

  iodata->log( "Minimum fluid density: " + stringify(min) );
  iodata->log( "Maximum fluid density: " + stringify(max) );
  iodata->log( "Average fluctuation density: " + stringify(average(DIFFr)) );
  iodata->log( "Std.dev fluctuation density: " + stringify(standard_deviation(DIFFr)) );
  if(min < 0.0)
  {
    iodata->log( "Error: negative density in some regions.");
    throw -1;
  }

  arr_t & DIFFK_p = *bssnSim->fields["DIFFK_p"];
  arr_t & DIFFK_a = *bssnSim->fields["DIFFK_a"];
  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    DIFFK_a[idx] = -sqrt(24.0*PI*rho_FRW);
    DIFFK_p[idx] = -sqrt(24.0*PI*rho_FRW);
  }
}

void ParticleSim::initParticleStep()
{
  _timer["RK_steps"].start();
    bssnSim->stepInit();
    particles->stepInit(bssnSim->fields);
    bssnSim->clearSrc();
    particles->addParticlesToBSSNSrc(bssnSim->fields);
  _timer["RK_steps"].stop();
}

void ParticleSim::outputParticleStep()
{
  _timer["output"].start();
    prepBSSNOutput();
    io_bssn_fields_snapshot(iodata, step, bssnSim->fields);
    io_bssn_fields_powerdump(iodata, step, bssnSim->fields, fourier);
    io_bssn_dump_statistics(iodata, step, bssnSim->fields, bssnSim->frw);
    io_bssn_constraint_violation(iodata, step, bssnSim);
  _timer["output"].stop();
}

void ParticleSim::runParticleStep()
{
  _timer["RK_steps"].start();
    // First RK step
    bssnSim->K1Calc();
    particles->RK1Step(bssnSim->fields);
    bssnSim->regSwap_c_a();
    particles->regSwap_c_a();

    // Second RK step
    bssnSim->clearSrc();
    particles->addParticlesToBSSNSrc(bssnSim->fields);
    bssnSim->K2Calc();
    particles->RK2Step(bssnSim->fields);
    bssnSim->regSwap_c_a();
    particles->regSwap_c_a();

    // Third RK step
    bssnSim->clearSrc();
    particles->addParticlesToBSSNSrc(bssnSim->fields);
    bssnSim->K3Calc();
    particles->RK3Step(bssnSim->fields);
    bssnSim->regSwap_c_a();
    particles->regSwap_c_a();

    // Fourth RK step
    bssnSim->clearSrc();
    particles->addParticlesToBSSNSrc(bssnSim->fields);
    bssnSim->K4Calc();
    particles->RK4Step(bssnSim->fields);

    // Wrap up
    bssnSim->stepTerm();
    particles->stepTerm();
    // "current" data should be in the _p array.
  _timer["RK_steps"].stop();
}

void ParticleSim::runStep()
{
  runCommonStepTasks();

  initParticleStep();
  outputParticleStep();
  runParticleStep();
}

} /* namespace cosmo */
