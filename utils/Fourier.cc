
#include "Fourier.h"

namespace cosmo
{

Fourier::Fourier()
{
  // template for initialization; see Fourier::Initialize.
}

Fourier::~Fourier()
{
  // dealloc
#if USE_LONG_DOUBLES
  fftwl_free(f_field);
  fftwl_destroy_plan(p_r2c);
  fftwl_destroy_plan(p_c2r);
#else
  fftw_free(f_field);
  fftw_destroy_plan(p_r2c);
  fftw_destroy_plan(p_c2r);
#endif
  delete [] double_field;
}

}
