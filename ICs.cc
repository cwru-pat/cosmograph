#include "ICs.h"

namespace cosmo
{

// analytic form of a power spectrum to use
// eg in LCDM, http://ned.ipac.caltech.edu/level5/Sept11/Norman/Norman2.html
real_t cosmo_power_spectrum(real_t k, ICsData *icd)
{
  real_t pre = icd->peak_amplitude*4.0/3.0/icd->peak_k;
  return pre*abs(k)/(1.0 + pow(abs(k)/icd->peak_k, 4.0)/3.0);
}

// set a field to an arbitrary gaussian random field
void set_gaussian_random_field(real_t *field, Fourier *fourier, ICsData *icd)
{
  idx_t i,j,k;
  real_t px, py, pz, pmag;
  real_t scale;

  // populate "field" with random values
  std::default_random_engine generator;
  std::normal_distribution<real_t> distribution(0.0,1.0);
  LOOP3(i,j,k)
  {
    field[NP_INDEX(i,j,k)] = distribution(generator);
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

        scale = cosmo_power_spectrum(pmag, icd);

        (fourier->f_field)[FFT_NP_INDEX(i,j,k)][0] *= scale;
        (fourier->f_field)[FFT_NP_INDEX(i,j,k)][1] *= scale;
      }
    }
  }
  // zero-mode (mean density)... set this to something
  (fourier->f_field[FFT_NP_INDEX(0,0,0)])[0] = 0;
  (fourier->f_field[FFT_NP_INDEX(0,0,0)])[1] = 0;

  // FFT back; 'field' array should now be populated.
  fftw_execute_dft_c2r(fourier->p_c2r, fourier->f_field, field);
}

void set_physical_metric_and_density(
  std::map <std::string, real_t *> & bssn_fields,
  std::map <std::string, real_t *> & hydro_fields,
  Fourier *fourier)
{
  idx_t i, j, k;
  real_t px, py, pz, p2i;
  // working in conformal transverse-traceless decomposition,
  // the conformal factor in front of metric is the solution to
  // d^2 f = f^5 * \rho
  // or 
  // d^2 f = f^5 * \rho_conformal
  // this assumes \rho_conformal was specified in UD_a, and that
  // we will need to convert it to \rho = f^-5 * \rho_conformal.
  // It is possible to calculate f, we 

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
        p2i = 1/(pw2(2.0*sin(PI*px/N)) + pw2(2.0*sin(PI*py/N)) + pw2(2.0*sin(PI*pz/N)));
        (fourier->f_field)[FFT_NP_INDEX(i,j,k)][0] *= p2i;
        (fourier->f_field)[FFT_NP_INDEX(i,j,k)][1] *= p2i;
      }
    }
  }
  // ignore the zero-mode (average) for now
  (fourier->f_field[FFT_NP_INDEX(0,0,0)])[0] = 0;
  (fourier->f_field[FFT_NP_INDEX(0,0,0)])[1] = 0;

  // temporarily use phi as conformal factor; conformal factor is inverse fft of this
  fftw_execute_dft_c2r(fourier->p_c2r, fourier->f_field, bssn_fields["phi_a"]);
  
  // TODO: reconstruct physical metric
  // TODO: reconstruct physical density

}

} // namespace cosmo
