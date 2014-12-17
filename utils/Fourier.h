#ifndef COSMO_UTILS_FOURIER_H
#define COSMO_UTILS_FOURIER_H

#include <fftw3.h>

namespace cosmo
{

class Fourier
{
public:
  // FFT field
  fftw_complex *f_field;

  // plans for taking FFTs
  fftw_plan p_c2r;
  fftw_plan p_r2c;

  Fourier();
  ~Fourier();

  template<typename IT, typename RT>
  void Initialize(IT n, RT *field);
};


template<typename IT, typename RT>
void Fourier::Initialize(IT n, RT *field)
{
  //fftw_malloc
  f_field = (fftw_complex *) fftw_malloc(n*n*(n/2+1)
                                         *((long long) sizeof(fftw_complex)));

  // create plans
  p_r2c = fftw_plan_dft_r2c_3d(n, n, n,
                               field, f_field,
                               FFTW_MEASURE);
  p_c2r = fftw_plan_dft_c2r_3d(n, n, n,
                               f_field, field,
                               FFTW_MEASURE);
}


}

#endif
