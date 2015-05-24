#include "ICs.h"

namespace cosmo
{

ICsData cosmo_get_ICsData()
{
  ICsData icd = {0};

  // H_LEN_FRAC is box side length in hubble units
  icd.rho_K_matter = 3.0/PI/8.0*pw2(H_LEN_FRAC/(N*dx)); // matter density term from FRW equation

  real_t rho_K_lambda_frac = (real_t) stold(_config["rho_K_lambda_frac"]); // DE density
  icd.rho_K_lambda = rho_K_lambda_frac*icd.rho_K_matter;

  // power spectrum amplitude as a fraction of the density
  real_t peak_amplitude_frac = (real_t) stold(_config["peak_amplitude_frac"]); // fluctuation amplitude
  real_t peak_amplitude = icd.rho_K_matter*peak_amplitude_frac*(1.0e-6); // scaling in arb. units

  real_t ic_spec_cut = (real_t) stold(_config["ic_spec_cut"]); // power spectrum cutoff parameter

  /* (peak scale in hubble units) * (to pixel scale) */
  icd.peak_k = (1.0/0.07)*H_LEN_FRAC; // need to scale this correctly?
  icd.peak_amplitude = peak_amplitude;
  icd.ic_spec_cut = ic_spec_cut; // cut spectrum off around p ~ ic_spec_cut
                                 // (max is p ~ sqrt(2.5)*N )

  return icd;
}

// analytic form of a power spectrum to use
// eg in LCDM, http://ned.ipac.caltech.edu/level5/Sept11/Norman/Norman2.html
real_t cosmo_power_spectrum(real_t k, ICsData *icd)
{
  real_t pre = icd->peak_amplitude*4.0/3.0;
  return pre*fabs(k)/(1.0 + pow(fabs(k)/icd->peak_k, 4.0)/3.0);
}

// set a field to an arbitrary gaussian random field
void set_gaussian_random_field(real_t *field, Fourier *fourier, ICsData *icd)
{
  idx_t i, j, k;
  real_t px, py, pz, pmag;
  real_t scale;

  // populate "field" with random values
  std::random_device rd;
  std::mt19937 gen(7.0 /*rd()*/);
  std::normal_distribution<real_t> gaussian_distribution;
  std::uniform_real_distribution<double> angular_distribution(0.0, 2.0*PI);
   // calling these here before looping suppresses a warning (bug)
  gaussian_distribution(gen);
  angular_distribution(gen);

  // scale amplitudes in fourier space
  for(i=0; i<N; i++)
  {
    px = (real_t) (i<=N/2 ? i : i-N);
    for(j=0; j<N; j++)
    {
      py = (real_t) (j<=N/2 ? j : j-N);
      for(k=0; k<N/2+1; k++)
      {
        pz = (real_t) k;
        pmag = sqrt(pw2(px) + pw2(py) + pw2(pz));

        // Scale by power spectrum
        // don't want much power on scales smaller than ~3 pixels
        // Or scales p > 1/(3*dx), or p > N/3
        real_t cutoff = 1.0/(1.0 + exp(10.0*(pmag - icd->ic_spec_cut)));
        scale = cutoff*sqrt(cosmo_power_spectrum(pmag, icd));

        real_t rand_mag = gaussian_distribution(gen);
        real_t rand_phase = angular_distribution(gen);

        (fourier->f_field)[FFT_NP_INDEX(i,j,k)][0] = scale*rand_mag*cos(rand_phase);
        (fourier->f_field)[FFT_NP_INDEX(i,j,k)][1] = scale*rand_mag*sin(rand_phase);

      }
      // run through more random numbers so ICs are similar
      // at different resolutions up to 256^3
      for(int x=N/2+1; x<256/2+1; x++) {
        gaussian_distribution(gen);
        angular_distribution(gen);
      }
    }
    // run through more random numbers so ICs are similar
    // at different resolutions up to 256^3
    for(int y=N; y<256; y++)
      for(int x=0; x<256/2+1; x++) {
        gaussian_distribution(gen);
        angular_distribution(gen);
      }
  }

  // zero-mode (mean density)... set this to something later
  (fourier->f_field[FFT_NP_INDEX(0,0,0)])[0] = 0;
  (fourier->f_field[FFT_NP_INDEX(0,0,0)])[1] = 0;

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
  real_t px, py, pz, p2i;
  ICsData icd = cosmo_get_ICsData();
  LOG(iod->log, "Generating ICs with peak at k = " << icd.peak_k << "\n");
  LOG(iod->log, "Generating ICs with peak amp. = " << icd.peak_amplitude << "\n");

  set_gaussian_random_field(static_field["D_a"], fourier, &icd);

  // the conformal factor in front of metric is the solution to
  // d^2 f = -2*pi f^5 * \rho
  // or 
  // d^2 f = -\rho_conformal
  // so assume -\rho_conformal was specified in UD_a, and that
  // we will need to convert it to \rho = f^-5 * \rho_conformal/(2pi).

  // FFT of conformal field
  fftw_execute_dft_r2c(fourier->p_r2c, static_field["D_a"], fourier->f_field);

  // scale amplitudes in fourier space;
  // solve for: \tilde{f} = \tilde{\rho_conformal}/k^2 (minus signs cancel)
  for(i=0; i<N; i++)
  {
    px = (real_t) (i<=N/2 ? i : i-N);
    for(j=0; j<N; j++)
    {
      py = (real_t) (j<=N/2 ? j : j-N);
      for(k=0; k<N/2+1; k++)
      {
        pz = (real_t) k;
        // Here we choose "1/k^2" such that the derivative stencil
        // applied later will agree with the metric solution we find.
        // For the usual second order laplacian stencil:
        // p2i = 1.0/4.0/(pw2(sin(PI*px/N)) + pw2(sin(PI*py/N)) + pw2(sin(PI*pz/N)));
        // for a sum of 4th-order order second derivative stencils,
        p2i = dx*dx*3.0/(
            16.0*( pw2(sin(PI*px/N)) + pw2(sin(PI*py/N)) + pw2(sin(PI*pz/N)) )
            - 1.0*( pw2(sin(2.0*PI*px/N)) + pw2(sin(2.0*PI*py/N)) + pw2(sin(2.0*PI*pz/N)) )
          );
        // account for fftw normalization here
        (fourier->f_field)[FFT_NP_INDEX(i,j,k)][0] *= p2i;
        (fourier->f_field)[FFT_NP_INDEX(i,j,k)][1] *= p2i;
      }
    }
  }
  // ignore the zero-mode (average) for now
  (fourier->f_field[FFT_NP_INDEX(0,0,0)])[0] = 0;
  (fourier->f_field[FFT_NP_INDEX(0,0,0)])[1] = 0;

  // temporarily use phi to store conformal factor; conformal factor is inverse fft of this
  // conformal factor array uses _p register initially, rather than _a
  fftw_execute_dft_c2r(fourier->p_c2r, fourier->f_field, bssn_fields["phi_p"]);
  // add in monopole contribution
  // scale density to correct "units"
  LOOP3(i,j,k)
  {
    bssn_fields["phi_p"][NP_INDEX(i,j,k)] /= POINTS; // account for FFTW transform normalization
    bssn_fields["phi_p"][NP_INDEX(i,j,k)] += 1.0;
    static_field["D_a"][NP_INDEX(i,j,k)] /= 2.0*PI*pow(bssn_fields["phi_p"][NP_INDEX(i,j,k)], 5);
  }
  // UD_a should now contain the physical density fluctuations
  LOG(iod->log, "Average physical fluctuation density (*16 pi) is: " << 16*PI*average(static_field["D_a"]) << "\n");

  // reconstruct BSSN conformal metric factor; \phi = log(\psi)
  LOOP3(i,j,k)
  {
    bssn_fields["phi_p"][NP_INDEX(i,j,k)] = log(bssn_fields["phi_p"][NP_INDEX(i,j,k)]);
    bssn_fields["phi_f"][NP_INDEX(i,j,k)] = bssn_fields["phi_p"][NP_INDEX(i,j,k)];
  }

  // check constraint: see if \grad^2 e^\phi = e^(5*\phi)*\rho
  real_t resid = 0.0;
  LOOP3(i,j,k)
  {
    // 4th-order stencil of e^\phi
    real_t grad2e = (
        - 1.0*(
          exp(bssn_fields["phi_p"][INDEX(i+2,j,k)]) + exp(bssn_fields["phi_p"][INDEX(i,j+2,k)]) + exp(bssn_fields["phi_p"][INDEX(i,j,k+2)])
          + exp(bssn_fields["phi_p"][INDEX(i-2,j,k)]) + exp(bssn_fields["phi_p"][INDEX(i,j-2,k)]) + exp(bssn_fields["phi_p"][INDEX(i,j,k-2)])
        )
        + 16.0*(
          exp(bssn_fields["phi_p"][INDEX(i+1,j,k)]) + exp(bssn_fields["phi_p"][INDEX(i,j+1,k)]) + exp(bssn_fields["phi_p"][INDEX(i,j,k+1)])
          + exp(bssn_fields["phi_p"][INDEX(i-1,j,k)]) + exp(bssn_fields["phi_p"][INDEX(i,j-1,k)]) + exp(bssn_fields["phi_p"][INDEX(i,j,k-1)])
        )
        - 90.0*exp(bssn_fields["phi_p"][NP_INDEX(i,j,k)])
      )/12.0/dx/dx;
    // conformal density
    real_t conformal_density = static_field["D_a"][NP_INDEX(i,j,k)]*exp(5.0*bssn_fields["phi_p"][NP_INDEX(i,j,k)]);

    // fractional constraint violation
    resid += fabs(grad2e + 2.0*PI*conformal_density) / sqrt(pw2(grad2e) + pw2(2.0*PI*conformal_density));
  }
  LOG(iod->log, "Average fractional constraint violation magnitude is: " << resid/POINTS << "\n");

  // Make sure min density value > 0
  real_t min = static_field["D_a"][NP_INDEX(0,0,0)] + icd.rho_K_matter;
  real_t max = min;
  resid = 0.0;

  // Scale density variable to conservative units
  LOOP3(i,j,k)
  {
    bssn_fields["K_p"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*(icd.rho_K_matter));
    bssn_fields["K_f"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*(icd.rho_K_matter));

    static_field["D_a"][NP_INDEX(i,j,k)] += icd.rho_K_matter;
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

  LOG(iod->log, "Minimum fluid 'density': " << min << "\n");
  LOG(iod->log, "Maximum fluid 'density': " << max << "\n");
  LOG(iod->log, "Average fluid 'density': " << average(static_field["D_a"]) << "\n");
  LOG(iod->log, "Std.dev fluid 'density': " << standard_deviation(static_field["D_a"]) << "\n");
  if(min < 0.0) {
    LOG(iod->log, "Error: negative density in some regions.\n");
    throw -1;
  }

  // Add in cosmological constant
  real_t oldrho;
  LOOP3(i,j,k)
  {
    oldrho = pw2(bssn_fields["K_p"][NP_INDEX(i,j,k)])/24.0/PI;
    bssn_fields["K_p"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*(icd.rho_K_lambda+oldrho));
    bssn_fields["K_f"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*(icd.rho_K_lambda+oldrho));
  }

}


void set_flat_dynamic_ICs(
  std::map <std::string, real_t *> & bssn_fields,
  std::map <std::string, real_t *> & hydro_fields,
  Fourier *fourier, IOData *iod)
{
  idx_t i, j, k;
  ICsData icd = cosmo_get_ICsData();

  // simple hamiltonian constraint; flat metric, \phi = 0, 
  set_gaussian_random_field(hydro_fields["UD_a"], fourier, &icd);
  // add in average
  LOOP3(i,j,k)
  {
    hydro_fields["UD_a"][NP_INDEX(i,j,k)] += icd.rho_K_matter;
  }

  // Make sure min density value > 0
  real_t min = hydro_fields["UD_a"][NP_INDEX(0,0,0)];
  real_t max = min;
  LOOP3(i,j,k)
  {
    if(hydro_fields["UD_a"][NP_INDEX(i,j,k)] < min)
    {
      min = hydro_fields["UD_a"][NP_INDEX(i,j,k)];
    }
    if(hydro_fields["UD_a"][NP_INDEX(i,j,k)] > max)
    {
      max = hydro_fields["UD_a"][NP_INDEX(i,j,k)];
    }
    if(hydro_fields["UD_a"][NP_INDEX(i,j,k)] != hydro_fields["UD_a"][NP_INDEX(i,j,k)])
    {
      LOG(iod->log, "Error: NaN energy density.\n");
      throw -1;
    }
  }

  LOG(iod->log, "Minimum fluid density: " << min << "\n");
  LOG(iod->log, "Maximum fluid density: " << max << "\n");
  LOG(iod->log, "Average fluid density: " << average(hydro_fields["UD_a"]) << "\n");
  LOG(iod->log, "Std.dev fluid density: " << standard_deviation(hydro_fields["UD_a"]) << "\n");
  if(min < 0.0) {
    LOG(iod->log, "Error: negative density in some regions.\n");
    throw -1;
  }

  // K = -sqrt(24pi*\rho)
  LOOP3(i,j,k)
  {
    bssn_fields["K_p"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*(icd.rho_K_lambda+hydro_fields["UD_a"][NP_INDEX(i,j,k)]));
    bssn_fields["K_f"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*(icd.rho_K_lambda+hydro_fields["UD_a"][NP_INDEX(i,j,k)]));
  }

  // Momentum constraint must now be satisfied
  // US_i = d_i K / 12 / pi
  LOOP3(i,j,k)
  {
    hydro_fields["US1_a"][NP_INDEX(i,j,k)] = (bssn_fields["K_p"][INDEX(i+1,j,k)] - bssn_fields["K_p"][INDEX(i-1,j,k)])/2.0/dx/12.0/PI;
    hydro_fields["US2_a"][NP_INDEX(i,j,k)] = (bssn_fields["K_p"][INDEX(i,j+1,k)] - bssn_fields["K_p"][INDEX(i,j-1,k)])/2.0/dx/12.0/PI;
    hydro_fields["US3_a"][NP_INDEX(i,j,k)] = (bssn_fields["K_p"][INDEX(i,j,k+1)] - bssn_fields["K_p"][INDEX(i,j,k-1)])/2.0/dx/12.0/PI;
  }

  // that's it...
  return;
}

void set_flat_static_ICs(
  std::map <std::string, real_t *> & bssn_fields,
  std::map <std::string, real_t *> & hydro_fields,
  Fourier *fourier, IOData *iod)
{
  idx_t i, j, k;

}

} // namespace cosmo
