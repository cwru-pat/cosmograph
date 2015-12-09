#include "ICs.h"

namespace cosmo
{

ICsData cosmo_get_ICsData()
{
  ICsData icd = {0};

  icd.rho_K_matter = 3.0/PI/8.0; // matter density term from FRW equation

  real_t rho_K_lambda_frac = (real_t) stold(_config["rho_K_lambda_frac"]); // DE density
  icd.rho_K_lambda = rho_K_lambda_frac*icd.rho_K_matter;

  // power spectrum amplitude as a fraction of the density
  real_t peak_amplitude_frac = (real_t) stold(_config["peak_amplitude_frac"]); // fluctuation amplitude
  real_t peak_amplitude = peak_amplitude_frac*(1.0e-15); // scaling in arb. units

  real_t ic_spec_cut = (real_t) stold(_config["ic_spec_cut"]); // power spectrum cutoff parameter

  /* (peak scale in hubble units) * (to pixel scale) */
  icd.peak_k = (1.0/0.07)*H_LEN_FRAC;
  icd.peak_amplitude = peak_amplitude;
  icd.ic_spec_cut = ic_spec_cut; // cut spectrum off around p ~ ic_spec_cut

  icd.viol_amp = stold(_config["IC_viol_amp"]);

  return icd;
}

// analytic form of a power spectrum to use
// eg in LCDM, http://ned.ipac.caltech.edu/level5/Sept11/Norman/Norman2.html
real_t cosmo_power_spectrum(real_t k, ICsData *icd)
{
  real_t pre = icd->peak_amplitude;
  return pre/(1.0 + pow(fabs(k)/icd->peak_k, 4.0)/3.0)/pow(fabs(k)/icd->peak_k, 3.0);
}

// set a field to an arbitrary gaussian random field
void set_gaussian_random_field(real_t *field, Fourier *fourier, ICsData *icd)
{
  idx_t i, j, k;
  real_t px, py, pz, pmag;
  real_t scale;

  // populate "field" with random values
  std::random_device rd;
  std::mt19937 gen(9.0 /*rd()*/);
  std::normal_distribution<real_t> gaussian_distribution;
  std::uniform_real_distribution<double> angular_distribution(0.0, 2.0*PI);
   // calling these here before looping suppresses a warning (bug)
  gaussian_distribution(gen);
  angular_distribution(gen);

  // scale amplitudes in fourier space
  // don't expect to run at >512^3 anytime soon; loop over all momenta out to that resolution.
  // this won't work for a larger grid.
  idx_t NMAX = 512;
  for(i=0; i<NMAX; i++)
  {
    px = (real_t) (i<=NMAX/2 ? i : i-NMAX);
    for(j=0; j<NMAX; j++)
    {
      py = (real_t) (j<=NMAX/2 ? j : j-NMAX);
      for(k=0; k<NMAX/2+1; k++)
      {
        pz = (real_t) k;

        // generate the same random modes for all resolutions (up to NMAX)
        real_t rand_mag = gaussian_distribution(gen);
        real_t rand_phase = angular_distribution(gen);

        // only store momentum values for relevant bins
        if( fabs(px) < (real_t) NX/2+1 + 0.01 && fabs(py) < (real_t) NY/2+1 + 0.01 && fabs(pz) < (real_t) NZ/2+1 + 0.01 )
        {
          idx_t fft_index = FFT_NP_INDEX(
            px > -0.5 ? ROUND_2_IDXT(px) : NX + ROUND_2_IDXT(px),
            py > -0.5 ? ROUND_2_IDXT(py) : NY + ROUND_2_IDXT(py),
            pz > -0.5 ? ROUND_2_IDXT(pz) : NZ + ROUND_2_IDXT(pz)
          );

          pmag = sqrt(
            pw2(px * ( (real_t) N / (real_t) NX ) )
             + pw2(py * ( (real_t) N / (real_t) NY ))
             + pw2(pz * ( (real_t) N / (real_t) NZ ))
            );

          // Scale by power spectrum
          // don't want much power on scales smaller than ~3 pixels
          // Or scales p > 1/(3*dx)
          real_t cutoff = 1.0 / (
              1.0 + exp(10.0*(pmag - icd->ic_spec_cut))
          );
          scale = cutoff*sqrt(cosmo_power_spectrum(pmag, icd));

          (fourier->f_field)[fft_index][0] = scale*rand_mag*cos(rand_phase);
          (fourier->f_field)[fft_index][1] = scale*rand_mag*sin(rand_phase);
        }
      }
    }
  }

  // zero-mode (mean density)... set this to something later
  (fourier->f_field)[FFT_NP_INDEX(0,0,0)][0] = 0;
  (fourier->f_field)[FFT_NP_INDEX(0,0,0)][1] = 0;

  // FFT back; 'field' array should now be populated with a gaussian random
  // field and power spectrum given by cosmo_power_spectrum.
  fftw_execute_dft_c2r(fourier->p_c2r, fourier->f_field, field);

  return;
}

// doesn't specify monopole / expansion contribution
void set_conformal_ICs(
  std::map <std::string, arr_t *> & bssn_fields,
  Fourier *fourier, IOData *iod, FRW<real_t> *frw)
{
  idx_t i, j, k;
  real_t px, py, pz, p2;

  real_t * const DIFFphi_p = bssn_fields["DIFFphi_p"]->_array;
  real_t * const DIFFdustrho_p = bssn_fields["DIFFdustrho_p"]->_array;

  ICsData icd = cosmo_get_ICsData();
  LOG(iod->log, "Generating ICs with peak at k = " << icd.peak_k << "\n");
  LOG(iod->log, "Generating ICs with peak amp. = " << icd.peak_amplitude << "\n");

  // the conformal factor in front of metric is the solution to
  // d^2 exp(\phi) = -2*pi exp(5\phi) * \rho
  // generate gaussian random field xi = exp(phi) (use phi_p as a proxy):
  set_gaussian_random_field(DIFFphi_p, fourier, &icd);

  // phi = ln(xi)
  PARALLEL_LOOP3(i,j,k) {
    idx_t idx = NP_INDEX(i,j,k);
    DIFFphi_p[idx] = log1p(DIFFphi_p[idx]);
  }

  // rho = -lap(phi)/xi^5/2pi
  PARALLEL_LOOP3(i,j,k) {
    DIFFdustrho_p[INDEX(i,j,k)] = -0.5/PI*exp(-4.0*DIFFphi_p[INDEX(i,j,k)])*(
      double_derivative(i, j, k, 1, 1, bssn_fields["DIFFphi_p"])
      + double_derivative(i, j, k, 2, 2, bssn_fields["DIFFphi_p"])
      + double_derivative(i, j, k, 3, 3, bssn_fields["DIFFphi_p"])
      + pw2(derivative(i, j, k, 1, bssn_fields["DIFFphi_p"]))
      + pw2(derivative(i, j, k, 2, bssn_fields["DIFFphi_p"]))
      + pw2(derivative(i, j, k, 3, bssn_fields["DIFFphi_p"]))
    );
  }

  // Make sure min density value > 0
  // Set conserved density variable field
  real_t min = icd.rho_K_matter;
  real_t max = min;
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    real_t rho_FRW = icd.rho_K_matter;
    real_t DIFFdustrho = DIFFdustrho_p[idx];
    real_t rho = rho_FRW + DIFFdustrho;
    // phi_FRW = 0
    real_t DIFFphi = DIFFphi_p[idx];
    // phi = DIFFphi
    // DIFFK = 0

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
      LOG(iod->log, "Error: NaN energy density.\n");
      throw -1;
    }
  }

  LOG(iod->log, "Minimum fluid conservative conformal density: " << min << "\n");
  LOG(iod->log, "Maximum fluid conservative conformal density: " << max << "\n");
  LOG(iod->log, "Average fluctuation density: " << average(bssn_fields["DIFFdustrho_p"]) << "\n");
  LOG(iod->log, "Std.dev fluctuation density: " << standard_deviation(bssn_fields["DIFFdustrho_p"]) << "\n");
  if(min < 0.0) {
    LOG(iod->log, "Error: negative density in some regions.\n");
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
    // add in FRW pieces to ICs
    // phi is unchanged
    // rho (D) and K get contribs
    // w=0 fluid only
    PARALLEL_LOOP3(i,j,k)
    {
      idx_t idx = NP_INDEX(i,j,k);
      real_t rho_FRW = icd.rho_K_matter;
      real_t D_FRW = rho_FRW; // on initial slice

      DIFFdustrho_p[idx] += rho_FRW;
      bssn_fields["DIFFK_p"]->_array[idx] = -sqrt(24.0*PI*rho_FRW);
    }
  #endif

  PARALLEL_LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    bssn_fields["eta_p"]->_array[idx] = 2.0/3.0;
  }

}


void set_stability_test_ICs(
  std::map <std::string, arr_t *> & bssn_fields)
{
  idx_t i, j, k;

  std::random_device rd;
  std::mt19937 gen(7.0 /*rd()*/);
  std::uniform_real_distribution<real_t> dist((-1.0e-10)*50/NX*50/NX, (1.0e-10)*50/NX*50/NX);

  PARALLEL_LOOP3(i, j, k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    // default flat static vacuum spacetime.
    bssn_fields["DIFFgamma11_p"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma11_a"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma11_c"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma11_f"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma12_p"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma12_a"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma12_c"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma12_f"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma13_p"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma13_a"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma13_c"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma13_f"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma22_p"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma22_a"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma22_c"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma22_f"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma23_p"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma23_a"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma23_c"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma23_f"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma33_p"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma33_a"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma33_c"]->_array[idx] = dist(gen);
    bssn_fields["DIFFgamma33_f"]->_array[idx] = dist(gen);
    
    bssn_fields["DIFFphi_p"]->_array[idx] = dist(gen);
    bssn_fields["DIFFphi_a"]->_array[idx] = dist(gen);
    bssn_fields["DIFFphi_c"]->_array[idx] = dist(gen);
    bssn_fields["DIFFphi_f"]->_array[idx] = dist(gen);
    
    bssn_fields["A11_p"]->_array[idx]     = dist(gen);
    bssn_fields["A11_a"]->_array[idx]     = dist(gen);
    bssn_fields["A11_c"]->_array[idx]     = dist(gen);
    bssn_fields["A11_f"]->_array[idx]     = dist(gen);
    bssn_fields["A12_p"]->_array[idx]     = dist(gen);
    bssn_fields["A12_a"]->_array[idx]     = dist(gen);
    bssn_fields["A12_c"]->_array[idx]     = dist(gen);
    bssn_fields["A12_f"]->_array[idx]     = dist(gen);
    bssn_fields["A13_p"]->_array[idx]     = dist(gen);
    bssn_fields["A13_a"]->_array[idx]     = dist(gen);
    bssn_fields["A13_c"]->_array[idx]     = dist(gen);
    bssn_fields["A13_f"]->_array[idx]     = dist(gen);
    bssn_fields["A22_p"]->_array[idx]     = dist(gen);
    bssn_fields["A22_a"]->_array[idx]     = dist(gen);
    bssn_fields["A22_c"]->_array[idx]     = dist(gen);
    bssn_fields["A22_f"]->_array[idx]     = dist(gen);
    bssn_fields["A23_p"]->_array[idx]     = dist(gen);
    bssn_fields["A23_a"]->_array[idx]     = dist(gen);
    bssn_fields["A23_c"]->_array[idx]     = dist(gen);
    bssn_fields["A23_f"]->_array[idx]     = dist(gen);
    bssn_fields["A33_p"]->_array[idx]     = dist(gen);
    bssn_fields["A33_a"]->_array[idx]     = dist(gen);
    bssn_fields["A33_c"]->_array[idx]     = dist(gen);
    bssn_fields["A33_f"]->_array[idx]     = dist(gen);

    bssn_fields["DIFFK_p"]->_array[idx]   = dist(gen);
    bssn_fields["DIFFK_a"]->_array[idx]   = dist(gen);
    bssn_fields["DIFFK_c"]->_array[idx]   = dist(gen);
    bssn_fields["DIFFK_f"]->_array[idx]   = dist(gen);

    bssn_fields["Gamma1_p"]->_array[idx]  = dist(gen);
    bssn_fields["Gamma1_a"]->_array[idx]  = dist(gen);
    bssn_fields["Gamma1_c"]->_array[idx]  = dist(gen);
    bssn_fields["Gamma1_f"]->_array[idx]  = dist(gen);
    bssn_fields["Gamma2_p"]->_array[idx]  = dist(gen);
    bssn_fields["Gamma2_a"]->_array[idx]  = dist(gen);
    bssn_fields["Gamma2_c"]->_array[idx]  = dist(gen);
    bssn_fields["Gamma2_f"]->_array[idx]  = dist(gen);
    bssn_fields["Gamma3_p"]->_array[idx]  = dist(gen);
    bssn_fields["Gamma3_a"]->_array[idx]  = dist(gen);
    bssn_fields["Gamma3_c"]->_array[idx]  = dist(gen);
    bssn_fields["Gamma3_f"]->_array[idx]  = dist(gen);

    bssn_fields["DIFFalpha_p"]->_array[idx]  = 0.0;
    bssn_fields["DIFFalpha_a"]->_array[idx]  = 0.0;
    bssn_fields["DIFFalpha_c"]->_array[idx]  = 0.0;
    bssn_fields["DIFFalpha_f"]->_array[idx]  = 0.0;

    bssn_fields["DIFFdustrho_a"]->_array[NP_INDEX(i,j,k)] = 0.0 + 0.0*dist(gen);
  }

std::cout << "dist is: " << dist(gen) << "\n";

}

void set_linear_wave_ICs(
  std::map <std::string, arr_t *> & bssn_fields)
{
  idx_t i, j, k;

  PARALLEL_LOOP3(i,j,k)
  {
    bssn_fields["DIFFgamma22_p"]->_array[NP_INDEX(i,j,k)] = 1.0e-8*sin( 2.0*PI*((real_t) i)*dx );
    bssn_fields["DIFFgamma33_p"]->_array[NP_INDEX(i,j,k)] = -1.0e-8*sin( 2.0*PI*((real_t) i)*dx );
    bssn_fields["A22_p"]->_array[NP_INDEX(i,j,k)] = PI*1.0e-8*cos( 2.0*PI*((real_t) i)*dx );
    bssn_fields["A33_p"]->_array[NP_INDEX(i,j,k)] = -PI*1.0e-8*cos( 2.0*PI*((real_t) i)*dx );

    // FRW Background parameters
    // real_t rho = D_MATTER; // phi = 0

    // bssn_fields["r_a"][NP_INDEX(i,j,k)] = rho;

    // bssn_fields["K_p"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*rho);
    // bssn_fields["K_a"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*rho);
    // bssn_fields["K_f"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*rho);

  }

}

} // namespace cosmo
