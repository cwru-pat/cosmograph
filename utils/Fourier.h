#ifndef COSMO_UTILS_FOURIER_H
#define COSMO_UTILS_FOURIER_H

// FFT indexing
#define FFT_INDEX(i,j,k) ((NZ/2+1)*NY*((i+NX)%NX) + (NZ/2+1)*((j+NY)%NY) + ((k+NZ)%NZ))
// FFT indexing without periodicity
#define FFT_NP_INDEX(i,j,k) ((NZ/2+1)*NY*(i) + (NZ/2+1)*(j) + (k))

#include <fftw3.h>
#include <zlib.h>
#include <math.h>
#include <string>
#include <iostream>
#include "../cosmo_macros.h"
#include "../cosmo_types.h"
#include "../cosmo_globals.h"

namespace cosmo
{

class Fourier
{
public:
#if USE_LONG_DOUBLES
  typedef long double fft_rt;
  typedef fftwl_complex fft_ct;
#else
  typedef double fft_realt;
  typedef fftw_complex fft_ct;
#endif

  // FFT field
  fft_ct *f_field;
  fft_rt *double_field;

  // plans for taking FFTs
#if USE_LONG_DOUBLES
  fftwl_plan p_c2r;
  fftwl_plan p_r2c;
#else
  fftw_plan p_c2r;
  fftw_plan p_r2c;
#endif

  Fourier();
  ~Fourier();

  template<typename IT>
  void Initialize(IT nx, IT ny, IT nz);
  
  template<typename IT>
  void Initialize_1D(IT n);

#if USE_LONG_DOUBLES
  void execute_f_r2c() { fftwl_execute(p_r2c); }
  void execute_f_c2r() { fftwl_execute(p_c2r); }
#else
  void execute_f_r2c() { fftw_execute(p_r2c); }
  void execute_f_c2r() { fftw_execute(p_c2r); }
#endif

  template<typename RT, typename IOT>
  void powerDump(RT *in, IOT *iodata);

  template<typename IT, typename RT>
  void inverseLaplacian(RT *field);
};


/**
 * @brief Initialize a fourier class instance
 * @details Create fftw plans, allocate memory
 * 
 * @param nx points in x-direction
 * @param ny points in y-direction
 * @param nz points in z-direction
 */
template<typename IT>
void Fourier::Initialize(IT nx, IT ny, IT nz)
{
  // create plans
#if USE_LONG_DOUBLES
  //fftw_malloc
  f_field = (fft_ct *) fftwl_malloc(nx*ny*(nz/2+1)
                                         *((long long) sizeof(fft_ct)));
  double_field = new fft_rt[nx*ny*nz];

  p_r2c = fftwl_plan_dft_r2c_3d(nx, ny, nz,
                               double_field, f_field,
                               FFTW_MEASURE);
  p_c2r = fftwl_plan_dft_c2r_3d(nx, ny, nz,
                               f_field, double_field,
                               FFTW_MEASURE);
#else
  //fftw_malloc
  f_field = (fft_ct *) fftw_malloc(nx*ny*(nz/2+1)
                                         *((long long) sizeof(fft_ct)));
  double_field = new fft_rt[nx*ny*nz];

  p_r2c = fftw_plan_dft_r2c_3d(nx, ny, nz,
                               double_field, f_field,
                               FFTW_MEASURE);
  p_c2r = fftw_plan_dft_c2r_3d(nx, ny, nz,
                               f_field, double_field,
                               FFTW_MEASURE);
#endif
}

template<typename IT>
void Fourier::Initialize_1D(IT n)
{
  // create plans
#if USE_LONG_DOUBLES
  //fftw_malloc
  f_field = (fft_ct *) fftwl_malloc((n/2 +1)*((long long) sizeof(fft_ct)));
  double_field = new fft_rt[n];

  p_r2c = fftwl_plan_dft_r2c_1d(n,
                               double_field, f_field,FFTW_ESTIMATE
                               );
  p_c2r = fftwl_plan_dft_c2r_1d(n,
                               f_field, double_field,FFTW_ESTIMATE
			       );
#else
  //fftw_malloc
  f_field = (fft_ct *) fftw_malloc((n/2 +1)*((long long) sizeof(fft_ct)));
  double_field = new fft_rt[n];

  p_r2c = fftw_plan_dft_r2c_1d(n,
                               double_field, f_field,FFTW_ESTIMATE
                               );
  p_c2r = fftw_plan_dft_c2r_1d(n,
                               f_field, double_field,FFTW_ESTIMATE
             );
#endif
}


/**
 * @brief Compute a power spectrum and write to file
 * 
 * @param in Grid to compute the spectrum of
 * @param iodata reference to data structure containing output directory
 * information, iodata->dir()
 */
template<typename RT, typename IOT>
void Fourier::powerDump(RT *in, IOT *iodata)
{
  // Transform input array
  for(long int i=0; i<POINTS; ++i)
    double_field[i] = (fft_rt) in[i];

#if USE_LONG_DOUBLES
  fftwl_execute_dft_r2c(p_r2c, double_field, f_field);
#else
  fftw_execute_dft_r2c(p_r2c, double_field, f_field);
#endif

  // average power over angles
  const int numbins = (int) (sqrt(NX*NX + NY*NY + NZ*NZ)/2.0) + 1; // Actual number of bins
  RT * array_out = new RT[numbins];
  int * numpoints = new int[numbins]; // Number of points in each momentum bin
  RT * p = new RT[numbins];
  RT * f2 = new RT[numbins]; // Values for each bin: Momentum, |F-k|^2, n_k

  fft_rt pmagnitude; // Total momentum (p) in units of lattice spacing, pmagnitude = Sqrt(px^2+py^2+pz^2).
                     // This also gives the bin index since bin spacing is set to equal lattice spacing.
  fft_rt fp2;
  int i, j, k, px, py, pz; // px, py, and pz are components of momentum in units of grid spacing

  // Initial magnitude of momentum in each bin
  for(i=0; i<numbins; i++) {
    f2[i] = 0.0;
    numpoints[i] = 0;
  }

  // Perform average over all angles here (~integral d\Omega).
  for(i=0; i<NX; i++)
  {
    px = (i<=NX/2 ? i : i-NX);
    for(j=0; j<NY; j++)
    {
      py = (j<=NY/2 ? j : j-NY);
      for(k=1; k<NZ/2; k++)
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

      pz = NZ/2;
      k = NZ/2;
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
      array_out[i] = f2[i]/((fft_rt) numpoints[i]);
    }
    else
    {
      array_out[i] = 0.;
    }
  }

  // write data
  std::string filename = iodata->dir() + "spec.dat.gz";
  char data[20];

  gzFile datafile = gzopen(filename.c_str(), "ab");
  if(datafile == Z_NULL) {
    printf("Error opening file: %s\n", filename.c_str());
    return;
  }

  for(i=0; i<numbins; i++)
  {
    // field values
    sprintf(data, "%g\t", (fft_rt) array_out[i]);
    gzwrite(datafile, data, std::char_traits<char>::length(data));
  }
  gzwrite(datafile, "\n", std::char_traits<char>::length("\n")); 

  gzclose(datafile);
  delete[] array_out;
  delete[] numpoints;
  delete[] p;
  delete[] f2;

  return;
}


// compute inverse laplacian for input array
template<typename IT, typename RT>
void Fourier::inverseLaplacian(RT *field)
{
  IT i, j, k;
  RT px, py, pz, pmag;

  for(long int i=0; i<POINTS; ++i)
    double_field[i] = (fft_rt) field[i];

#if USE_LONG_DOUBLES
  fftwl_execute_dft_r2c(p_r2c, double_field, f_field);
#else
  fftw_execute_dft_r2c(p_r2c, double_field, f_field);
#endif

  for(i=0; i<NX; i++)
  {
    px = (RT) (i<=NX/2 ? i : i-NX);
    for(j=0; j<NY; j++)
    {
      py = (RT) (j<=NY/2 ? j : j-NY);
      for(k=0; k<NZ/2+1; k++)
      {
        pz = (RT) k;

        IT fft_index = FFT_NP_INDEX(i,j,k);

        pmag = sqrt( pw2(px) + pw2(py) + pw2(pz) )*2.0*PI/H_LEN_FRAC;

        f_field[fft_index][0] /= -pmag*pmag*POINTS;
        f_field[fft_index][1] /= -pmag*pmag*POINTS;
      }
    }
  }
  // zero mode?
  f_field[0][0] = 0;
  f_field[0][1] = 0;

#if USE_LONG_DOUBLES
  fftwl_execute_dft_c2r(p_c2r, f_field, double_field);
#else
  fftw_execute_dft_c2r(p_c2r, f_field, double_field);
#endif
  

  for(long int i=0; i<POINTS; ++i)
    field[i] = (RT) double_field[i];
}

} // namespace cosmo

#endif
