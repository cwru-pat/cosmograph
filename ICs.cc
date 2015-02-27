#include "ICs.h"

namespace cosmo
{

// analytic form of a power spectrum to use
// eg in LCDM, http://ned.ipac.caltech.edu/level5/Sept11/Norman/Norman2.html
real_t cosmo_power_spectrum(real_t k, ICsData *icd)
{
  real_t pre = (1.0/N)*icd->peak_amplitude*4.0/3.0;
  return pre*fabs(k)/((real_t) N)/(1.0 + pow(fabs(k)/((real_t) N)/icd->peak_k, 4.0)/3.0);
}

// set a field to an arbitrary gaussian random field
void set_gaussian_random_field(real_t *field, Fourier *fourier, ICsData *icd)
{
  idx_t i, j, k;
  real_t px, py, pz, pmag;
  real_t scale;

  // populate "field" with random values
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<real_t> distribution;
  distribution(gen); // calling this here suppresses a warning (bug)
  LOOP3(i,j,k)
  {
    field[NP_INDEX(i,j,k)] = distribution(gen);
  }

  // FFT of random field
  fftw_execute_dft_r2c(fourier->p_r2c, field, fourier->f_field);

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
        // don't want much power on scales smaller than ~2 pixels
        // Or scales p > 1/(3*dx), or p > N/3
        real_t cutoff = 1.0/(1.0+exp(pmag-N/2.0));
        scale = cutoff*sqrt(cosmo_power_spectrum(pmag, icd));

        // fftw transform is unnormalized; account for an N^3 here
        (fourier->f_field)[FFT_NP_INDEX(i,j,k)][0] *= scale/((real_t) POINTS);
        (fourier->f_field)[FFT_NP_INDEX(i,j,k)][1] *= scale/((real_t) POINTS);
      }
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
  std::map <std::string, real_t *> & hydro_fields,
  Fourier *fourier,
  ICsData *icd,
  real_t rho_K_matter,
  real_t rho_K_lambda)
{
  idx_t i, j, k;
  real_t px, py, pz, p2i;

  set_gaussian_random_field(hydro_fields["UD_a"], fourier, icd);

  // working in conformal transverse-traceless decomposition,
  // the conformal factor in front of metric is the solution to
  // d^2 f = f^5 * \rho
  // or 
  // d^2 f = \rho_conformal
  // this assumes \rho_conformal*2pi was specified in UD_a, and that
  // we will need to convert it to \rho = f^-5 * \rho_conformal.

  // FFT of conformal field
  fftw_execute_dft_r2c(fourier->p_r2c, hydro_fields["UD_a"], fourier->f_field);

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
        // Here we choose the magnitude of k such that the derivative stencil
        // applied later will agree with the metric solution we find.
        p2i = 1.0/(pw2(2.0*sin(PI*px/N)) + pw2(2.0*sin(PI*py/N)) + pw2(2.0*sin(PI*pz/N)));
        // account for fftw normalization here
        (fourier->f_field)[FFT_NP_INDEX(i,j,k)][0] *= p2i/((real_t) POINTS);
        (fourier->f_field)[FFT_NP_INDEX(i,j,k)][1] *= p2i/((real_t) POINTS);
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
    bssn_fields["phi_p"][NP_INDEX(i,j,k)] += 1.0;
    hydro_fields["UD_a"][NP_INDEX(i,j,k)] /= 2.0*PI;
  }
  // UD_a should now contain the conformal density

  // reconstruct physical density, and store in UD
  LOOP3(i,j,k)
  {
    hydro_fields["UD_a"][NP_INDEX(i,j,k)] /= pow(bssn_fields["phi_p"][NP_INDEX(i,j,k)], 5);
  }

  // reconstruct BSSN conformal metric factor; \phi = log(\psi)
  LOOP3(i,j,k)
  {
    bssn_fields["phi_p"][NP_INDEX(i,j,k)] = log(bssn_fields["phi_p"][NP_INDEX(i,j,k)]);
    bssn_fields["phi_f"][NP_INDEX(i,j,k)] = bssn_fields["phi_p"][NP_INDEX(i,j,k)];
  }

  // check constraint: see if \grad^2 e^\phi = \bar{\rho} (conformal density fluctuations, stored in UD_a)
  real_t resid = 0.0;
  LOOP3(i,j,k)
  {
    real_t grad2e = exp(bssn_fields["phi_p"][INDEX(i+1,j,k)]) + exp(bssn_fields["phi_p"][INDEX(i,j+1,k)]) + exp(bssn_fields["phi_p"][INDEX(i,j,k+1)])
        + exp(bssn_fields["phi_p"][INDEX(i-1,j,k)]) + exp(bssn_fields["phi_p"][INDEX(i,j-1,k)]) + exp(bssn_fields["phi_p"][INDEX(i,j,k-1)])
        - 6.0*exp(bssn_fields["phi_p"][NP_INDEX(i,j,k)]);
    real_t cden = hydro_fields["UD_a"][NP_INDEX(i,j,k)]*exp(5.0*bssn_fields["phi_p"][NP_INDEX(i,j,k)]);
    resid += fabs(grad2e + 2.0*PI*cden); // should = 0
  }
  std::cout << "Total metric solution error residual is: " << resid << "\n";

  // Make sure min density value > 0
  real_t min = hydro_fields["UD_a"][NP_INDEX(0,0,0)] + rho_K_matter;
  real_t max = min;
  real_t oldrho;

  LOOP3(i,j,k)
  {
    hydro_fields["UD_a"][NP_INDEX(i,j,k)] += rho_K_matter;
    hydro_fields["UD_a"][NP_INDEX(i,j,k)] *= exp(6.0*bssn_fields["phi_p"][NP_INDEX(i,j,k)]);

    oldrho = pw2(bssn_fields["K_p"][NP_INDEX(i,j,k)])/24.0/PI;
    bssn_fields["K_p"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*(rho_K_matter+oldrho));
    bssn_fields["K_f"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*(rho_K_matter+oldrho));

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
      std::cout << "Error: NaN energy density.\n";
      throw -1;
    }
  }

  std::cout << "Minimum fluid 'density': " << min << "\n";
  std::cout << "Maximum fluid 'density': " << max << "\n";
  std::cout << "Average fluid 'density': " << average(hydro_fields["UD_a"]) << "\n";
  std::cout << "Std.dev fluid 'density': " << standard_deviation(hydro_fields["UD_a"]) << "\n";
  if(min < 0.0) {
    std::cout << "Error: negative density in some regions.\n";
    throw -1;
  }

  // Add in cosmological constant
  LOOP3(i,j,k)
  {
    oldrho = pw2(bssn_fields["K_p"][NP_INDEX(i,j,k)])/24.0/PI;
    bssn_fields["K_p"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*(rho_K_lambda+oldrho));
    bssn_fields["K_f"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*(rho_K_lambda+oldrho));
  }

}

void set_flat_ICs(
  std::map <std::string, real_t *> & bssn_fields,
  std::map <std::string, real_t *> & hydro_fields,
  Fourier *fourier,
  ICsData *icd,
  real_t rho_K_matter,
  real_t rho_K_lambda)
{
  idx_t i, j, k;

  // simple hamiltonian constraint; flat metric, \phi = 0, 
  set_gaussian_random_field(hydro_fields["UD_a"], fourier, icd);
  // add in average
  LOOP3(i,j,k)
  {
    hydro_fields["UD_a"][NP_INDEX(i,j,k)] += rho_K_matter;
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
      std::cout << "Error: NaN energy density.\n";
      throw -1;
    }
  }

  std::cout << "Minimum fluid density: " << min << "\n";
  std::cout << "Maximum fluid density: " << max << "\n";
  std::cout << "Average fluid density: " << average(hydro_fields["UD_a"]) << "\n";
  std::cout << "Std.dev fluid density: " << standard_deviation(hydro_fields["UD_a"]) << "\n";
  if(min < 0.0) {
    std::cout << "Error: negative density in some regions.\n";
    throw -1;
  }

  // K = -sqrt(24pi*\rho)
  LOOP3(i,j,k)
  {
    bssn_fields["K_p"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*(rho_K_lambda+hydro_fields["UD_a"][NP_INDEX(i,j,k)]));
    bssn_fields["K_f"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*(rho_K_lambda+hydro_fields["UD_a"][NP_INDEX(i,j,k)]));
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

} // namespace cosmo
