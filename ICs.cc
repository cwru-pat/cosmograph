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
  real_t peak_amplitude = icd.rho_K_matter*peak_amplitude_frac*(1.0e-7)*(5.0/0.7); // scaling in arb. units

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
  real_t pre = icd->peak_amplitude*4.0/3.0;
  return pre*fabs(k)/icd->peak_k/(1.0 + pow(fabs(k)/icd->peak_k, 4.0)/3.0);
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
        if( fabs(px) < (real_t) N/2+1 + 0.01 && fabs(py) < (real_t) N/2+1 + 0.01 && fabs(pz) < (real_t) N/2+1 + 0.01 )
        {
          idx_t fft_index = FFT_NP_INDEX(
            px > -0.5 ? ROUND_2_IDXT(px) : N + ROUND_2_IDXT(px),
            py > -0.5 ? ROUND_2_IDXT(py) : N + ROUND_2_IDXT(py),
            pz > -0.5 ? ROUND_2_IDXT(pz) : N + ROUND_2_IDXT(pz)
          );

          pmag = sqrt(pw2(px) + pw2(py) + pw2(pz));

          // Scale by power spectrum
          // don't want much power on scales smaller than ~3 pixels
          // Or scales p > 1/(3*dx), or p > N/3
          real_t cutoff = 1.0 / (
              1.0 +
                exp(10.0*(fabs(px) - icd->ic_spec_cut))
                *exp(10.0*(fabs(py) - icd->ic_spec_cut))
                *exp(10.0*(fabs(pz) - icd->ic_spec_cut))
          );
          if(fabs(px)+0.01 > icd->ic_spec_cut || fabs(py)+0.01 > icd->ic_spec_cut || fabs(pz)+0.01 > icd->ic_spec_cut) {
            cutoff = 0.0;
          }
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
  std::map <std::string, real_t *> & bssn_fields,
  std::map <std::string, real_t *> & static_field,
  Fourier *fourier, IOData *iod)
{
  idx_t i, j, k;
  real_t px, py, pz, p2;
  ICsData icd = cosmo_get_ICsData();
  LOG(iod->log, "Generating ICs with peak at k = " << icd.peak_k << "\n");
  LOG(iod->log, "Generating ICs with peak amp. = " << icd.peak_amplitude << "\n");

  // the conformal factor in front of metric is the solution to
  // d^2 exp(\phi) = -2*pi exp(5\phi) * \rho
  // generate gaussian random field xi = exp(phi) (use phi_p as a proxy):
  set_gaussian_random_field(bssn_fields["phi_p"], fourier, &icd);
  // scale xi and grad^2 xi
  LOOP3(i,j,k) {
    bssn_fields["phi_p"][NP_INDEX(i,j,k)] += 1.0;
  }

  // rho = -lap(phi)/xi^5/2pi
  LOOP3(i,j,k) {
    bssn_fields["r_a"][NP_INDEX(i,j,k)] = -0.5/PI/(
      pow(bssn_fields["phi_p"][NP_INDEX(i,j,k)], 5.0)
    )*(
      double_derivative(i, j, k, 1, 1, bssn_fields["phi_p"])
      + double_derivative(i, j, k, 2, 2, bssn_fields["phi_p"])
      + double_derivative(i, j, k, 3, 3, bssn_fields["phi_p"])
    );
  }

  // phi = ln(xi)
  LOOP3(i,j,k) {
    bssn_fields["phi_a"][NP_INDEX(i,j,k)] = log(bssn_fields["phi_p"][NP_INDEX(i,j,k)]);
    bssn_fields["phi_f"][NP_INDEX(i,j,k)] = log(bssn_fields["phi_p"][NP_INDEX(i,j,k)]);
    bssn_fields["phi_p"][NP_INDEX(i,j,k)] = log(bssn_fields["phi_p"][NP_INDEX(i,j,k)]);
  }

  // Make sure min density value > 0
  // Set conserved density variable field
  real_t min = bssn_fields["r_a"][NP_INDEX(0,0,0)] + icd.rho_K_matter;
  real_t max = min;
  LOOP3(i,j,k)
  {
    // add bad ICs to K

    bssn_fields["K_p"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*(icd.rho_K_matter));
    bssn_fields["K_f"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*(icd.rho_K_matter));

    static_field["D_a"][NP_INDEX(i,j,k)] = icd.rho_K_matter + bssn_fields["r_a"][NP_INDEX(i,j,k)];
    static_field["D_a"][NP_INDEX(i,j,k)] *= exp(6.0*bssn_fields["phi_p"][NP_INDEX(i,j,k)]);

    if(static_field["D_a"][NP_INDEX(i,j,k)] < min)
    {
      min = static_field["D_a"][NP_INDEX(i,j,k)];
    }
    if(static_field["D_a"][NP_INDEX(i,j,k)] > max)
    {
      max = static_field["D_a"][NP_INDEX(i,j,k)];
    }
    if(static_field["D_a"][NP_INDEX(i,j,k)] != static_field["D_a"][NP_INDEX(i,j,k)])
    {
      LOG(iod->log, "Error: NaN energy density.\n");
      throw -1;
    }
  }

  LOG(iod->log, "Minimum fluid conservative conformal density: " << min << "\n");
  LOG(iod->log, "Maximum fluid conservative conformal density: " << max << "\n");
  LOG(iod->log, "Average fluid conservative conformal density: " << average(static_field["D_a"]) << "\n");
  LOG(iod->log, "Std.dev fluid conservative conformal density: " << standard_deviation(static_field["D_a"]) << "\n");
  if(min < 0.0) {
    LOG(iod->log, "Error: negative density in some regions.\n");
    throw -1;
  }

  // populate "field" with random values
  std::random_device rd;
  std::mt19937 gen(7.0 /*rd()*/);
  std::normal_distribution<real_t> gaussian_distribution;
  // calling these here before looping suppresses a warning (bug)
  gaussian_distribution(gen);

  // Add in cosmological constant
  real_t oldrho;
  LOOP3(i,j,k)
  {
    real_t viol = icd.viol_amp*gaussian_distribution(gen);

    oldrho = pw2(bssn_fields["K_p"][NP_INDEX(i,j,k)])/24.0/PI;
    bssn_fields["K_p"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*(icd.rho_K_lambda+oldrho)) + viol;
    bssn_fields["K_f"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*(icd.rho_K_lambda+oldrho)) + viol;
  }

}

} // namespace cosmo
