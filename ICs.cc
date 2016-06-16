#include "ICs.h"

namespace cosmo
{

ICsData cosmo_get_ICsData()
{
  ICsData icd = {0};

  icd.rho_K_matter = 3.0/PI/8.0; // matter density term from FRW equation

  real_t rho_K_lambda_frac = (real_t) stold(_config["rho_K_lambda_frac"]); // DE density
  icd.rho_K_lambda = rho_K_lambda_frac*icd.rho_K_matter;

  // power spectrum amplitude as a fraction of the density
  real_t peak_amplitude_frac = (real_t) stold(_config["peak_amplitude_frac"]); // fluctuation amplitude
  real_t peak_amplitude = peak_amplitude_frac*(1.0e-15); // scaling in arb. units

  real_t ic_spec_cut = (real_t) stold(_config["ic_spec_cut"]); // power spectrum cutoff parameter

  /* (peak scale in hubble units) * (to pixel scale) */
  icd.peak_k = (1.0/0.07)*H_LEN_FRAC;
  icd.peak_amplitude = peak_amplitude;
  icd.ic_spec_cut = ic_spec_cut; // cut spectrum off around p ~ ic_spec_cut

  icd.viol_amp = stold(_config["IC_viol_amp"]);

  return icd;
}

/**
 * @brief A cosmologically-motivated power spectrum for the conformal factor
 * @details Return the power spectrum amplitude for a particular wavenumber.
 *  The power spectrum describing the density goes roughly as 
 *    P(k) ~ k^{ 1} for small k
 *    P(k) ~ k^{-3} for large k
 *  (see, eg, http://ned.ipac.caltech.edu/level5/Sept11/Norman/Norman2.html)
 *  An analytic function with this form is the function
 *    P(k) ~ k / ( 1 + k^4 )
 *  The conformal factor obeys k^2 \phi ~ \rho; so \phi ~ \rho / k^2.
 *  The power spectrum < \phi \phi > then goes as 1/k^4 < \rho \rho >.
 *  Thus, the power spectrum for the conformal factor scales as
 *    P(k) ~ 1 / k^3 / ( 1 + k^4 )
 *  Which after additional scalings, this function returns.
 * 
 * @param k wavenumber of interest
 * @param icd initial condition data structre
 * 
 * @return power spectrum amplitude
 */
real_t cosmo_power_spectrum(real_t k, ICsData *icd)
{
  real_t pre = icd->peak_amplitude;
  return pre/(1.0 + pow(fabs(k)/icd->peak_k, 4.0)/3.0)/pow(fabs(k)/icd->peak_k, 3.0);
}

// set a field to an arbitrary gaussian random field

void set_gaussian_random_field(real_t *field, Fourier *fourier, ICsData *icd)
{
  idx_t i, j, k;
  real_t px, py, pz, pmag;
  real_t scale;

  // populate "field" with random values
  std::random_device rd;
  std::mt19937 gen(9.0 /*rd()*/);
  std::normal_distribution<real_t> gaussian_distribution;
  std::uniform_real_distribution<double> angular_distribution(0.0, 2.0*PI);
   // calling these here before looping suppresses a warning (bug)
  gaussian_distribution(gen);
  angular_distribution(gen);

  // scale amplitudes in fourier space
  // don't expect to run at >512^3 anytime soon; loop over all momenta out to that resolution.
  // this won't work for a larger grid.
  idx_t NMAX = 512;
  for(i=0; i<NMAX; i++)
  {
    px = (real_t) (i<=NMAX/2 ? i : i-NMAX);
    for(j=0; j<NMAX; j++)
    {
      py = (real_t) (j<=NMAX/2 ? j : j-NMAX);
      for(k=0; k<NMAX/2+1; k++)
      {
        pz = (real_t) k;

        // generate the same random modes for all resolutions (up to NMAX)
        real_t rand_mag = gaussian_distribution(gen);
        real_t rand_phase = angular_distribution(gen);

        // only store momentum values for relevant bins
        if( fabs(px) < (real_t) NX/2+1 + 0.01 && fabs(py) < (real_t) NY/2+1 + 0.01 && fabs(pz) < (real_t) NZ/2+1 + 0.01 )
        {
          idx_t fft_index = FFT_NP_INDEX(
            px > -0.5 ? ROUND_2_IDXT(px) : NX + ROUND_2_IDXT(px),
            py > -0.5 ? ROUND_2_IDXT(py) : NY + ROUND_2_IDXT(py),
            pz > -0.5 ? ROUND_2_IDXT(pz) : NZ + ROUND_2_IDXT(pz)
          );

          pmag = sqrt(
            pw2(px * ( (real_t) N / (real_t) NX ) )
             + pw2(py * ( (real_t) N / (real_t) NY ))
             + pw2(pz * ( (real_t) N / (real_t) NZ ))
            );

          // Scale by power spectrum
          // don't want much power on scales smaller than ~3 pixels
          // Or scales p > 1/(3*dx)
          real_t cutoff = 1.0 / (
              1.0 + exp(10.0*(pmag - icd->ic_spec_cut))
          );
          scale = cutoff*sqrt(cosmo_power_spectrum(pmag, icd));

          (fourier->f_field)[fft_index][0] = scale*rand_mag*cos(rand_phase);
          (fourier->f_field)[fft_index][1] = scale*rand_mag*sin(rand_phase);
        }
      }
    }
  }

  // zero-mode (mean density)... set this to something later
  (fourier->f_field)[FFT_NP_INDEX(0,0,0)][0] = 0;
  (fourier->f_field)[FFT_NP_INDEX(0,0,0)][1] = 0;

  // FFT back; 'field' array should now be populated with a gaussian random
  // field and power spectrum given by cosmo_power_spectrum.
  fftw_execute_dft_c2r(fourier->p_c2r, fourier->f_field, field);

  return;
}

/**
 * @brief Initialize vector of rays
 * @details Initialize a vector of rays with given rays, given velocity
 *  directions (ray velocity magnitude is normalized by the raytrace class)
 * 
 * @param rays reference to RayTrace class rays should belong to
 * @param n_rays number of rays
 */
void init_ray_vector(std::vector<RayTrace<real_t, idx_t> *> * rays, idx_t n_rays)
{
  idx_t i = 0;
  RaytraceData<real_t> rd = {0};

  std::mt19937 gen(7.0);
  std::uniform_real_distribution<real_t> dist(0.0, 1.0);
  dist(gen);
  for(i=0; i<n_rays; i++)
  {
    // rays distributed (uniformly) around a sphere
    real_t theta = 2.0*PI*dist(gen);
    real_t U = 2.0*dist(gen) - 1.0;
    // sphere radius R = L/2
    real_t R = N*dx/2.0;
    // inside some over/underdensity (hopefully)
    real_t X0, Y0, Z0;
    if(i<n_rays/3)
    {
      X0 = 0.382813*N*dx;
      Y0 = 0.476563*N*dx;
      Z0 = 0.351563*N*dx;
    }
    else if(i<2*n_rays/3)
    {
      X0 = 0.460938*N*dx;
      Y0 = 0.554688*N*dx;
      Z0 = 0.492188*N*dx;
    }
    else
    {
      X0 = 0.0429688*N*dx;
      Y0 = 0.0360625*N*dx;
      Z0 = 0.0390625*N*dx;
    }
    // ray position starts at an observer (integrating back in time...)
    rd.x[0] = X0;
    rd.x[1] = Y0;
    rd.x[2] = Z0;
    // velocity directed inwards
    // normalization of V is enforced by raytrace class
    rd.V[0] = -1.0*std::sqrt(1.0-U*U)*std::cos(theta);
    rd.V[1] = -1.0*std::sqrt(1.0-U*U)*std::sin(theta);
    rd.V[2] = -1.0*U;

    // energy, angle in arb. units
    rd.E = 1.0;
    rd.Phi = 1.0;

    RayTrace<real_t, idx_t> * ray;
    ray = new RayTrace<real_t, idx_t> (-std::fabs(dt), dx, rd);
    rays->push_back( ray );
  }
}


/** -- Need to fix up below --- **/

/**
 * @brief      AwA stability test
 *
 * @param      bssn_fields   The bssn fields
 * @param      static_field  The static field
 */
void set_stability_test_ICs(map_t & bssn_fields, map_t & static_field)
{
  idx_t i, j, k;

  std::random_device rd;
  std::mt19937 gen(7.0 /*rd()*/);
  std::uniform_real_distribution<real_t> dist((-1.0e-10)*50/NX*50/NX, (1.0e-10)*50/NX*50/NX);

  LOOP3(i, j, k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    // default flat static vacuum spacetime.
    bssn_fields["DIFFgamma11_p"][idx] = dist(gen);
    bssn_fields["DIFFgamma11_a"][idx] = dist(gen);
    bssn_fields["DIFFgamma11_c"][idx] = dist(gen);
    bssn_fields["DIFFgamma11_f"][idx] = dist(gen);
    bssn_fields["DIFFgamma12_p"][idx] = dist(gen);
    bssn_fields["DIFFgamma12_a"][idx] = dist(gen);
    bssn_fields["DIFFgamma12_c"][idx] = dist(gen);
    bssn_fields["DIFFgamma12_f"][idx] = dist(gen);
    bssn_fields["DIFFgamma13_p"][idx] = dist(gen);
    bssn_fields["DIFFgamma13_a"][idx] = dist(gen);
    bssn_fields["DIFFgamma13_c"][idx] = dist(gen);
    bssn_fields["DIFFgamma13_f"][idx] = dist(gen);
    bssn_fields["DIFFgamma22_p"][idx] = dist(gen);
    bssn_fields["DIFFgamma22_a"][idx] = dist(gen);
    bssn_fields["DIFFgamma22_c"][idx] = dist(gen);
    bssn_fields["DIFFgamma22_f"][idx] = dist(gen);
    bssn_fields["DIFFgamma23_p"][idx] = dist(gen);
    bssn_fields["DIFFgamma23_a"][idx] = dist(gen);
    bssn_fields["DIFFgamma23_c"][idx] = dist(gen);
    bssn_fields["DIFFgamma23_f"][idx] = dist(gen);
    bssn_fields["DIFFgamma33_p"][idx] = dist(gen);
    bssn_fields["DIFFgamma33_a"][idx] = dist(gen);
    bssn_fields["DIFFgamma33_c"][idx] = dist(gen);
    bssn_fields["DIFFgamma33_f"][idx] = dist(gen);
    
    bssn_fields["DIFFphi_p"][idx] = dist(gen);
    bssn_fields["DIFFphi_a"][idx] = dist(gen);
    bssn_fields["DIFFphi_c"][idx] = dist(gen);
    bssn_fields["DIFFphi_f"][idx] = dist(gen);
    
    bssn_fields["A11_p"][idx]     = dist(gen);
    bssn_fields["A11_a"][idx]     = dist(gen);
    bssn_fields["A11_c"][idx]     = dist(gen);
    bssn_fields["A11_f"][idx]     = dist(gen);
    bssn_fields["A12_p"][idx]     = dist(gen);
    bssn_fields["A12_a"][idx]     = dist(gen);
    bssn_fields["A12_c"][idx]     = dist(gen);
    bssn_fields["A12_f"][idx]     = dist(gen);
    bssn_fields["A13_p"][idx]     = dist(gen);
    bssn_fields["A13_a"][idx]     = dist(gen);
    bssn_fields["A13_c"][idx]     = dist(gen);
    bssn_fields["A13_f"][idx]     = dist(gen);
    bssn_fields["A22_p"][idx]     = dist(gen);
    bssn_fields["A22_a"][idx]     = dist(gen);
    bssn_fields["A22_c"][idx]     = dist(gen);
    bssn_fields["A22_f"][idx]     = dist(gen);
    bssn_fields["A23_p"][idx]     = dist(gen);
    bssn_fields["A23_a"][idx]     = dist(gen);
    bssn_fields["A23_c"][idx]     = dist(gen);
    bssn_fields["A23_f"][idx]     = dist(gen);
    bssn_fields["A33_p"][idx]     = dist(gen);
    bssn_fields["A33_a"][idx]     = dist(gen);
    bssn_fields["A33_c"][idx]     = dist(gen);
    bssn_fields["A33_f"][idx]     = dist(gen);

    bssn_fields["DIFFK_p"][idx]   = dist(gen);
    bssn_fields["DIFFK_a"][idx]   = dist(gen);
    bssn_fields["DIFFK_c"][idx]   = dist(gen);
    bssn_fields["DIFFK_f"][idx]   = dist(gen);

    bssn_fields["Gamma1_p"][idx]  = dist(gen);
    bssn_fields["Gamma1_a"][idx]  = dist(gen);
    bssn_fields["Gamma1_c"][idx]  = dist(gen);
    bssn_fields["Gamma1_f"][idx]  = dist(gen);
    bssn_fields["Gamma2_p"][idx]  = dist(gen);
    bssn_fields["Gamma2_a"][idx]  = dist(gen);
    bssn_fields["Gamma2_c"][idx]  = dist(gen);
    bssn_fields["Gamma2_f"][idx]  = dist(gen);
    bssn_fields["Gamma3_p"][idx]  = dist(gen);
    bssn_fields["Gamma3_a"][idx]  = dist(gen);
    bssn_fields["Gamma3_c"][idx]  = dist(gen);
    bssn_fields["Gamma3_f"][idx]  = dist(gen);

    bssn_fields["DIFFalpha_p"][idx]  = 0.0;
    bssn_fields["DIFFalpha_a"][idx]  = 0.0;
    bssn_fields["DIFFalpha_c"][idx]  = 0.0;
    bssn_fields["DIFFalpha_f"][idx]  = 0.0;

    static_field["DIFFD_a"][NP_INDEX(i,j,k)] = 0.0 + 0.0*dist(gen);
  }
}

void set_linear_wave_ICs(map_t & bssn_fields)
{
  idx_t i, j, k;

  LOOP3(i,j,k)
  {
    bssn_fields["DIFFgamma22_p"][NP_INDEX(i,j,k)] = 1.0e-8*sin( 2.0*PI*((real_t) i)*dx );
    bssn_fields["DIFFgamma22_a"][NP_INDEX(i,j,k)] = bssn_fields["DIFFgamma22_p"][NP_INDEX(i,j,k)];
    bssn_fields["DIFFgamma22_f"][NP_INDEX(i,j,k)] = bssn_fields["DIFFgamma22_p"][NP_INDEX(i,j,k)];

    bssn_fields["DIFFgamma33_p"][NP_INDEX(i,j,k)] = -1.0e-8*sin( 2.0*PI*((real_t) i)*dx );
    bssn_fields["DIFFgamma33_a"][NP_INDEX(i,j,k)] = bssn_fields["DIFFgamma33_p"][NP_INDEX(i,j,k)];
    bssn_fields["DIFFgamma33_f"][NP_INDEX(i,j,k)] = bssn_fields["DIFFgamma33_p"][NP_INDEX(i,j,k)];

    bssn_fields["A22_p"][NP_INDEX(i,j,k)] = PI*1.0e-8*cos( 2.0*PI*((real_t) i)*dx );
    bssn_fields["A22_a"][NP_INDEX(i,j,k)] = bssn_fields["A22_p"][NP_INDEX(i,j,k)];
    bssn_fields["A22_f"][NP_INDEX(i,j,k)] = bssn_fields["A22_p"][NP_INDEX(i,j,k)];

    bssn_fields["A33_p"][NP_INDEX(i,j,k)] = -PI*1.0e-8*cos( 2.0*PI*((real_t) i)*dx );
    bssn_fields["A33_a"][NP_INDEX(i,j,k)] = bssn_fields["A33_p"][NP_INDEX(i,j,k)];
    bssn_fields["A33_f"][NP_INDEX(i,j,k)] = bssn_fields["A33_p"][NP_INDEX(i,j,k)];

    // FRW Background parameters
    // real_t rho = D_MATTER; // phi = 0

    // bssn_fields["r_a"][NP_INDEX(i,j,k)] = rho;

    // bssn_fields["K_p"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*rho);
    // bssn_fields["K_a"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*rho);
    // bssn_fields["K_f"][NP_INDEX(i,j,k)] = -sqrt(24.0*PI*rho);
  }
}

} // namespace cosmo
