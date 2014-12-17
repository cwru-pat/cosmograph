
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
  fftw_free(f_field);
  fftw_destroy_plan(p_r2c);
  fftw_destroy_plan(p_c2r);
}

}
