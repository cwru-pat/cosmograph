#include "dust.h"

namespace cosmo
{

DustSim::DustSim()
{
  // just check to make sure we can use this class.
  if(USE_HARMONIC_ALPHA) {
    iodata->log("Error - not using synchronous gauge! You must use it for dust sims.");
    iodata->log("Please change this setting in cosmo_macros.h and recompile.");
    throw -1;
  }
}

void DustSim::init()
{
  _timer["init"].start();

  // initialize base class
  simInit();

  iodata->log("Running 'dust' type simulation.");
  staticSim = new Static();
  staticSim->init();
  
  iodata->log("Creating initial conditions.");
  setDustICs();

  _timer["init"].stop();
}

/**
 * @brief      Set conformally flat initial conditions for a w=0
 *  fluid in synchronous gauge.
 */
void DustSim::setDustICs()
{
  idx_t i, j, k;

  arr_t & DIFFr_a = *bssnSim->fields["DIFFr_a"];
  arr_t & DIFFphi_p = *bssnSim->fields["DIFFphi_p"];
  arr_t & DIFFphi_a = *bssnSim->fields["DIFFphi_a"];
  arr_t & DIFFphi_c = *bssnSim->fields["DIFFphi_c"];
  arr_t & DIFFphi_f = *bssnSim->fields["DIFFphi_f"];

  arr_t & DIFFD_a = *staticSim->fields["DIFFD_a"];

  auto & frw = bssnSim->frw;

  ICsData icd = cosmo_get_ICsData();
  iodata->log( "Generating ICs with peak at k = " + stringify(icd.peak_k) );
  iodata->log( "Generating ICs with peak amp. = " + stringify(icd.peak_amplitude) );

  // the conformal factor in front of metric is the solution to
  // d^2 exp(\phi) = -2*pi exp(5\phi) * \delta_rho
  // generate gaussian random field 1 + xi = exp(phi) (store xi in phi_p temporarily):
  set_gaussian_random_field(DIFFphi_p._array, fourier, &icd);

  // delta_rho = -lap(phi)/(1+xi)^5/2pi
  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k) {
    DIFFr_a[NP_INDEX(i,j,k)] = -0.5/PI/(
      pow(1.0 + DIFFphi_p[NP_INDEX(i,j,k)], 5.0)
    )*(
      double_derivative(i, j, k, 1, 1, DIFFphi_p)
      + double_derivative(i, j, k, 2, 2, DIFFphi_p)
      + double_derivative(i, j, k, 3, 3, DIFFphi_p)
    );
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
    idx_t idx = NP_INDEX(i,j,k);
    real_t rho_FRW = icd.rho_K_matter;
    real_t DIFFr = DIFFr_a[idx];
    real_t rho = rho_FRW + DIFFr;
    // phi_FRW = 0
    real_t DIFFphi = DIFFphi_a[idx];
    // phi = DIFFphi
    // DIFFK = 0

    DIFFD_a[idx] =
      rho_FRW*expm1(6.0*DIFFphi) + exp(6.0*DIFFphi)*DIFFr;

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
  iodata->log( "Average fluctuation density: " + stringify(average(DIFFD_a)) );
  iodata->log( "Std.dev fluctuation density: " + stringify(standard_deviation(DIFFD_a)) );
  if(min < 0.0)
  {
    iodata->log("Error: negative density in some regions.");
    throw -1;
  }

  #if USE_REFERENCE_FRW
    // Set values in reference FRW integrator
    real_t rho_FRW = icd.rho_K_matter;
    real_t K_frw = -sqrt(24.0*PI*rho_FRW);

    frw->set_phi(0.0);
    frw->set_K(K_frw);
    frw->addFluid(rho_FRW, 0.0 /* w=0 */);
  # else
    arr_t & DIFFK_p = *bssnSim->fields["DIFFK_p"];
    arr_t & DIFFK_a = *bssnSim->fields["DIFFK_a"];
    // add in FRW pieces to ICs
    // phi is unchanged
    // rho (D) and K get contribs
    // w=0 fluid only
    #pragma omp parallel for default(shared) private(i,j,k)
    LOOP3(i,j,k)
    {
      idx_t idx = NP_INDEX(i,j,k);
      real_t rho_FRW = icd.rho_K_matter;
      real_t D_FRW = rho_FRW; // on initial slice

      DIFFr_a[idx] += rho_FRW;

      DIFFK_a[idx] = -sqrt(24.0*PI*rho_FRW);
      DIFFK_p[idx] = -sqrt(24.0*PI*rho_FRW);

      DIFFD_a[idx] += D_FRW;
    }
  #endif
}

void DustSim::initDustStep()
{
  _timer["RK_steps"].start();
    bssnSim->stepInit();
    bssnSim->clearSrc();
    staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
  _timer["RK_steps"].stop();
}

void DustSim::outputDustStep()
{
  _timer["output"].start();
    prepBSSNOutput();
    io_bssn_fields_snapshot(iodata, step, bssnSim->fields);
    io_bssn_fields_powerdump(iodata, step, bssnSim->fields, fourier);
    io_bssn_dump_statistics(iodata, step, bssnSim->fields, bssnSim->frw);
    io_bssn_constraint_violation(iodata, step, bssnSim);
  _timer["output"].stop();
}

void DustSim::runDustStep()
{
  _timer["RK_steps"].start();
    // First RK step
    // source already set in initDustStep() (used for output)
    bssnSim->K1Calc();
    bssnSim->regSwap_c_a();

    // Second RK step
    bssnSim->clearSrc();
    staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
    bssnSim->K2Calc();
    bssnSim->regSwap_c_a();

    // Third RK step
    bssnSim->clearSrc();
    staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
    bssnSim->K3Calc();
    bssnSim->regSwap_c_a();

    // Fourth RK step
    bssnSim->clearSrc();
    staticSim->addBSSNSrc(bssnSim->fields, bssnSim->frw);
    bssnSim->K4Calc();

    // Wrap up
    bssnSim->stepTerm();
    // "current" data should be in the _p array.
  _timer["RK_steps"].stop();
}

void DustSim::runStep()
{
  runCommonStepTasks();

  initDustStep();
  outputDustStep();
  runDustStep();
}

} /* namespace cosmo */
