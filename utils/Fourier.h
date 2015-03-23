#ifndef COSMO_UTILS_FOURIER_H
#define COSMO_UTILS_FOURIER_H

#include <fftw3.h>
#include <zlib.h>
#include <string>
#include <iostream>
#include "../cosmo_macros.h"

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

  template<typename RT, typename IOT>
  void powerDump(RT *in, IOT *iodata);
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

template<typename RT, typename IOT>
void Fourier::powerDump(RT *in, IOT *iodata)
{
  // Transform input array
  fftw_execute_dft_r2c(p_r2c, in, f_field);

  // average power over angles
  const int numbins = (int)(1.74*(N/2.0))+1; // Actual number of bins
  RT array_out[numbins];
  int numpoints[numbins]; // Number of points in each momentum bin
  RT p[numbins];
  RT f2[numbins]; // Values for each bin: Momentum, |F-k|^2, n_k

  double pmagnitude; // Total momentum (p) in units of lattice spacing, pmagnitude = Sqrt(px^2+py^2+pz^2).
                     // This also gives the bin index since bin spacing is set to equal lattice spacing.
  double fp2;
  int i, j, k, px, py, pz; // px, py, and pz are components of momentum in units of grid spacing

  // Initial magnitude of momentum in each bin
  for(i=0; i<numbins; i++) {
    f2[i] = 0.0;
    numpoints[i] = 0;
  }

  // Perform average over all angles here (~integral d\Omega).
  for(i=0; i<N; i++)
  {
    px = (i<=N/2 ? i : i-N);
    for(j=0; j<N; j++)
    {
      py = (j<=N/2 ? j : j-N);
      for(k=1; k<N/2; k++)
      {
        pz = k;
        pmagnitude = sqrt((RT) (pw2(px) + pw2(py) + pw2(pz)));
        fp2 = pw2(C_RE((f_field)[FFT_NP_INDEX(i,j,k)])) + pw2(C_IM((f_field)[FFT_NP_INDEX(i,j,k)]));
        numpoints[(int)pmagnitude] += 2;
        f2[(int)pmagnitude] += 2.*fp2;
      }

      pz = 0;
      k = 0;
      pmagnitude = sqrt((RT) (pw2(px) + pw2(py) + pw2( pz)));
      fp2 = pw2(C_RE((f_field)[FFT_NP_INDEX(i,j,k)])) + pw2(C_IM((f_field)[FFT_NP_INDEX(i,j,k)]));
      numpoints[(int)pmagnitude] += 1;
      f2[(int)pmagnitude] += fp2;

      pz = N/2;
      k = N/2;
      pmagnitude = sqrt((RT) (pw2(px) + pw2(py) + pw2(pz)));
      fp2 = pw2(C_RE((f_field)[FFT_NP_INDEX(i,j,k)])) + pw2(C_IM((f_field)[FFT_NP_INDEX(i,j,k)]));
      numpoints[(int)pmagnitude] += 1;
      f2[(int)pmagnitude] += fp2;
    }
  }

  for(i=0; i<numbins; i++)
  {
    // Converts sums to averages. (numpoints[i] should always be greater than zero.)
    if(numpoints[i] > 0)
    {
      array_out[i] = f2[i]/((double) numpoints[i]);
    }
    else
    {
      array_out[i] = 0.;
    }
  }

  // write data
  std::string filename = iodata->output_dir + "spec.dat.gz";
  char data[20];

  gzFile datafile = gzopen(filename.c_str(), "ab");
  if(datafile == Z_NULL) {
    printf("Error opening file: %s\n", filename.c_str());
    return;
  }

  for(i=0; i<numbins; i++)
  {
    // field values
    sprintf(data, "%g\t", array_out[i]);
    gzwrite(datafile, data, std::char_traits<char>::length(data));
  }
  gzwrite(datafile, "\n", std::char_traits<char>::length("\n")); 

  gzclose(datafile);
  return;
}


}

#endif
