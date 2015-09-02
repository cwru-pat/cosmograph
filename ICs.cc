#include "ICs.h"

namespace cosmo
{

ICsData cosmo_get_ICsData()
{
  ICsData icd = {0};

  // H_LEN_FRAC is box side length in hubble units
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
        real_t cutoff = 1.0/( 1.0 + exp(10.0*(abs(px) - icd->ic_spec_cut))*exp(10.0*(abs(py) - icd->ic_spec_cut))*exp(10.0*(abs(pz) - icd->ic_spec_cut)) );
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
  // d^2 exp(\phi) = -2*pi exp(5\phi) * \rho
  // or
  // d^2 \phi = -2*pi exp(4\phi) * \rho - (d\phi)^2
  // or 
  // d^2 \phi = \rho_conformal
  // so given \rho_conformal = "D_a" (not actually the variable; just
  // a re-use of the field space for setting ICs), we will need to solve
  // for \phi and subsequently find the physical \rho
  // \rho = ("D_a" + (d\phi)^2) / (-2*pi) / exp(4\phi)

  // FFT of gaussian random field
  fftw_execute_dft_r2c(fourier->p_r2c, static_field["D_a"], fourier->f_field);

  // scale amplitudes in fourier space;
  // solve for: \tilde{\phi} = \tilde{-\rho_conformal}/k^2 (minus signs cancel)
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
        p2i = dx*dx*5040.0/4.0/(
            -  9.0*( pw2(sin(4.0*PI*px/N)) + pw2(sin(4.0*PI*py/N)) + pw2(sin(4.0*PI*pz/N)) )
            + 128.0*( pw2(sin(3.0*PI*px/N)) + pw2(sin(3.0*PI*py/N)) + pw2(sin(3.0*PI*pz/N)) )
            - 1008.0*( pw2(sin(2.0*PI*px/N)) + pw2(sin(2.0*PI*py/N)) + pw2(sin(2.0*PI*pz/N)) )
            + 8064.0*( pw2(sin(1.0*PI*px/N)) + pw2(sin(1.0*PI*py/N)) + pw2(sin(1.0*PI*pz/N)) )
          );
        // p2i = dx*dx*3.0/(
        //     16.0*( pw2(sin(PI*px/N)) + pw2(sin(PI*py/N)) + pw2(sin(PI*pz/N)) )
        //     - 1.0*( pw2(sin(2.0*PI*px/N)) + pw2(sin(2.0*PI*py/N)) + pw2(sin(2.0*PI*pz/N)) )
        //   );

        // account for fftw normalization here
        (fourier->f_field)[FFT_NP_INDEX(i,j,k)][0] *= p2i;
        (fourier->f_field)[FFT_NP_INDEX(i,j,k)][1] *= p2i;
      }
    }
  }
  // ignore the zero-mode (average) for now
  (fourier->f_field[FFT_NP_INDEX(0,0,0)])[0] = 0;
  (fourier->f_field[FFT_NP_INDEX(0,0,0)])[1] = 0;

  fftw_execute_dft_c2r(fourier->p_c2r, fourier->f_field, bssn_fields["phi_p"]);
  // account for FFTW transform normalization and set \phi_f
  LOOP3(i,j,k)
  {
    bssn_fields["phi_p"][NP_INDEX(i,j,k)] *= 1.0/POINTS;
    bssn_fields["phi_f"][NP_INDEX(i,j,k)] = bssn_fields["phi_p"][NP_INDEX(i,j,k)];
  }
  real_t phi_mean = average(bssn_fields["phi_p"]);
  LOG(iod->log, "Average conformal field value is: " << phi_mean << "\n");

  // reconstruct physical density and handle conformal factor normalization
  LOOP3(i,j,k)
  {
    static_field["D_a"][NP_INDEX(i,j,k)] = (
      // \rho_conformal + (d\phi)^2 => "D_a" + (d\phi)^2
      static_field["D_a"][NP_INDEX(i,j,k)] + (
        pw2(derivative(i, j, k, 1, bssn_fields["phi_p"]))
        + pw2(derivative(i, j, k, 2, bssn_fields["phi_p"]))
        + pw2(derivative(i, j, k, 3, bssn_fields["phi_p"]))
      )/pw2(12.0*dx)
    ) / (2.0*PI) / exp(4.0*bssn_fields["phi_p"][NP_INDEX(i,j,k)]);
  }
  // UD_a should now contain the physical density fluctuations
  real_t mean = average(static_field["D_a"]);
  LOG(iod->log, "Average physical fluctuation density (*16 pi) is: " << 16*PI*mean << "\n");
  LOG(iod->log, "St. Dev physical fluctuation density (*16 pi) is: " << standard_deviation(static_field["D_a"], mean) << "\n");

  // check constraint: see if \grad^2 e^\phi = e^(5*\phi)*\rho
  // (May differ at O(dx^2) from expression in terms of grad^2 \phi)
  real_t resid = 0.0;
  LOOP3(i,j,k)
  {
    // 4th-order stencil of e^\phi
    real_t grad2e = exp(bssn_fields["phi_p"][NP_INDEX(i,j,k)])*(
        // grad^2 phi
        (
          double_derivative(i, j, k, 1, 1, bssn_fields["phi_p"])
          + double_derivative(i, j, k, 2, 2, bssn_fields["phi_p"])
          + double_derivative(i, j, k, 3, 3, bssn_fields["phi_p"])
        )
        +
        // (d\phi)^2
        (
          pw2(derivative(i, j, k, 1, bssn_fields["phi_p"]))
          + pw2(derivative(i, j, k, 2, bssn_fields["phi_p"]))
          + pw2(derivative(i, j, k, 3, bssn_fields["phi_p"]))
        )
      );

    // conformal density
    // force density to agree with stencil used to calculate grad2e
    static_field["D_a"][NP_INDEX(i,j,k)] = -grad2e/2.0/PI/exp(5.0*bssn_fields["phi_p"][NP_INDEX(i,j,k)]);
    real_t density = static_field["D_a"][NP_INDEX(i,j,k)];

    // fractional constraint violation
    resid += fabs(grad2e + 2.0*PI*density*exp(5.0*bssn_fields["phi_p"][NP_INDEX(i,j,k)])) / sqrt(pw2(grad2e) + pw2(2.0*PI*density*exp(5.0*bssn_fields["phi_p"][NP_INDEX(i,j,k)])));
  }
  LOG(iod->log, "Average fractional constraint violation magnitude is: " << resid/POINTS << "\n");

  // Make sure min density value > 0
  real_t min = static_field["D_a"][NP_INDEX(0,0,0)] + icd.rho_K_matter;
  real_t max = min;
  resid = 0.0;

  // Scale density variable to conservative units
  LOOP3(i,j,k)
  {
    bssn_fields["phi_a"][NP_INDEX(i,j,k)] = bssn_fields["phi_p"][NP_INDEX(i,j,k)];

    // add bad ICs to K

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
