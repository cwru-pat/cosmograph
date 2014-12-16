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

void set_gaussian_random_field(void (*Pk)(real_t, real_t, real_t),
  real_t *field, ICsData *icd)
{
  idx_t i,j,k;
  real_t kmag;

  // populate "field" with random values
  std::default_random_engine generator;
  std::normal_distribution<real_t> distribution(0.0,1.0);
  LOOP3(i,j,k)
  {
    field[INDEX(i,j,k)] = distribution(generator);
  }

  // FFT of random field
  fftw_plan p = fftw_plan_r2r_3d(N, N, N,
      field, field,
      FFTW_R2HC, FFTW_R2HC, FFTW_R2HC, 
      FFTW_MEASURE
    );
  fftw_execute_r2r(p, field, field);

  // scale amplitudes of FFTs
  LOOP3(i,j,k)
  {
    kmag = sqrt(pw2(i) + pw2(j) + pw2(k%(N/2)))/((real_t) N);
    field[INDEX(i,j,k)] *= cosmo_power_spectrum(kmag, icd);
  }
}

} // namespace cosmo
