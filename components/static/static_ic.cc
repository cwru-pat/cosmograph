#include "static_ic.h"
#include "../../cosmo_includes.h"
#include "../../cosmo_types.h"
#include "../../cosmo_globals.h"
#include "../../ICs/ICs.h"
#include "../../utils/math.h"

#include <boost/math/special_functions/spherical_harmonic.hpp>

#include <utility>
#include <random>

namespace cosmo
{

typedef std::pair<real_t, real_t> complex_t;
#define m_idx(l, m) ((l)+(m))


/**
 * @brief Gaussian random field ICs
 * 
 * General steps:
 *  1) Generate the Newtonian gauge variables:
 *     phi = psi, \dot{phi} = 0, \delta_u
 *  2) Compute the synchronous gauge variables
 *  3) Fix \phi, \bar{gamma}_{ij}, K
 *  4) Solve for V, U in CTT decomposition -> obtain \bar{A}_{ij}^{\rm NL}
 *  5) Use the Hamiltonian constraint residual to solve for the density field
 */
void static_ic_set_random(BSSN * bssn, Static * stat, Lambda * lambda,
  Fourier * fourier, IOData * iodata)
{
  idx_t i, j, k;

  // Background cosmology, a_FRW = 1
  real_t rho_FRW = 3.0/PI/8.0;
  real_t Omega_L = std::stod(_config("Omega_L", "1.0e-6"));
  real_t p0 = std::stod(_config("p0", "7.0"));
  real_t P = H_LEN_FRAC*H_LEN_FRAC*1.0e-15*std::stod(_config("P", "1.0"));
  real_t p_cut = std::stod(_config("p_cut", "1.0"));
  real_t rho_L = Omega_L * rho_FRW;
  lambda->setLambda(rho_L);

  // final metric fields
  arr_t & DIFFphi_p = *bssn->fields["DIFFphi_p"];
  arr_t & DIFFK_p = *bssn->fields["DIFFK_p"];
  arr_t & A11_p = *bssn->fields["A11_p"];
  arr_t & A12_p = *bssn->fields["A12_p"];
  arr_t & A13_p = *bssn->fields["A13_p"];
  arr_t & A22_p = *bssn->fields["A22_p"];
  arr_t & A23_p = *bssn->fields["A23_p"];
  arr_t & A33_p = *bssn->fields["A33_p"];
  arr_t & DIFFD_a = *stat->fields["DIFFD_a"];

  // Extra / auxiliary fields for intermediate computations
  // (re-use _c register, it gets overwritten later anyways)
  arr_t & phi_N = *bssn->fields["DIFFphi_c"];
  arr_t & lap_phi_N = *bssn->fields["DIFFK_c"];
  arr_t & invlape6pd1K = *bssn->fields["DIFFgamma11_c"];
  arr_t & invlape6pd2K = *bssn->fields["DIFFgamma12_c"];
  arr_t & invlape6pd3K = *bssn->fields["DIFFgamma13_c"];
  arr_t & W1 = *bssn->fields["A22_c"];
  arr_t & W2 = *bssn->fields["A23_c"];
  arr_t & W3 = *bssn->fields["A33_c"];

  // 1.a) Synchronous-gauge Newtonian potential:
  set_gaussian_random_Phi_N(phi_N, fourier, P, p0, p_cut);

  // 1.b) Preliminary metric variables: phi, K
# pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    
    lap_phi_N[idx] = laplacian(i, j, k, phi_N);

    DIFFphi_p[idx] = log1p(-10.0/3.0*phi_N[idx])/4.0;
    DIFFK_p[idx] = -3.0*(1.0 + 2.0*phi_N[idx]) + 2.0/3.0*lap_phi_N[idx];
  }

  // 1.c) Vector contribution to the extrinsic curvature, W^i
  // compute invlapK
# pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k) {
    idx_t idx = NP_INDEX(i,j,k);
    real_t e6p = exp(6.0*DIFFphi_p[idx]);
    invlape6pd1K[idx] = e6p*derivative(i,j,k,1,DIFFK_p);
    invlape6pd2K[idx] = e6p*derivative(i,j,k,2,DIFFK_p);
    invlape6pd3K[idx] = e6p*derivative(i,j,k,3,DIFFK_p);
  }
  fourier->inverseLaplacian <idx_t, real_t> (invlape6pd1K._array);
  fourier->inverseLaplacian <idx_t, real_t> (invlape6pd2K._array);
  fourier->inverseLaplacian <idx_t, real_t> (invlape6pd3K._array);
  // compute W^i = W_i
# pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    W1[idx] = -derivative(i,j,k,1,phi_N)/3.0 + invlape6pd1K[idx]/2.0;
    W2[idx] = -derivative(i,j,k,2,phi_N)/3.0 + invlape6pd2K[idx]/2.0;
    W3[idx] = -derivative(i,j,k,3,phi_N)/3.0 + invlape6pd3K[idx]/2.0;
  }

  // 1.d) Aij components
# pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    // these are the CTT conformal Aij
    real_t CTT2BSSNAij = exp(-6.0*DIFFphi_p[idx]);
    real_t DkWk = derivative(i,j,k,1,W1) + derivative(i,j,k,2,W2) + derivative(i,j,k,3,W3);
    A11_p[idx] = CTT2BSSNAij * ( derivative(i,j,k,1,W1) + derivative(i,j,k,1,W1) + 2.0/3.0*double_derivative(i,j,k,1,1,phi_N)
      - 2.0/3.0*DkWk - 2.0/9.0*lap_phi_N[idx] );
    A12_p[idx] = CTT2BSSNAij * ( derivative(i,j,k,1,W2) + derivative(i,j,k,2,W1) + 2.0/3.0*double_derivative(i,j,k,1,2,phi_N) );
    A13_p[idx] = CTT2BSSNAij * ( derivative(i,j,k,1,W3) + derivative(i,j,k,3,W1) + 2.0/3.0*double_derivative(i,j,k,1,3,phi_N) );
    A22_p[idx] = CTT2BSSNAij * ( derivative(i,j,k,2,W2) + derivative(i,j,k,2,W2) + 2.0/3.0*double_derivative(i,j,k,2,2,phi_N)
      - 2.0/3.0*DkWk - 2.0/9.0*lap_phi_N[idx] );
    A23_p[idx] = CTT2BSSNAij * ( derivative(i,j,k,2,W3) + derivative(i,j,k,3,W2) + 2.0/3.0*double_derivative(i,j,k,2,3,phi_N) );
    A33_p[idx] = CTT2BSSNAij * ( derivative(i,j,k,3,W3) + derivative(i,j,k,3,W3) + 2.0/3.0*double_derivative(i,j,k,3,3,phi_N)
      - 2.0/3.0*DkWk - 2.0/9.0*lap_phi_N[idx] );
  }

  // 1.e) Density field
# pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    real_t AijAij = A11_p[idx]*A11_p[idx] + A22_p[idx]*A22_p[idx] + A33_p[idx]*A33_p[idx]
      + 2.0*(A12_p[idx]*A12_p[idx] + A13_p[idx]*A13_p[idx] + A23_p[idx]*A23_p[idx]);
    real_t lap_e_p = exp(DIFFphi_p[idx]) * ( laplacian(i,j,k,DIFFphi_p)
      + pw2(derivative(i,j,k,1,DIFFphi_p)) + pw2(derivative(i,j,k,2,DIFFphi_p)) + pw2(derivative(i,j,k,3,DIFFphi_p)) );
    real_t rho_ADM = 1.0/2.0/PI * ( DIFFK_p[idx]*DIFFK_p[idx]/12.0 - AijAij/8.0 - exp(-5.0*DIFFphi_p[idx])*lap_e_p );
    real_t rho_0 = rho_ADM - rho_L;

    if(rho_0 < 0.0)
    {
      iodata->log("Error: negative density in some regions.");
      throw -1;
    }

    DIFFD_a[idx] = exp(6.0*DIFFphi_p[idx])*rho_0;
  }

  iodata->log( "Average fluctuation density: " + stringify(average(DIFFD_a)) );
  iodata->log( "Std.dev fluctuation density: " + stringify(standard_deviation(DIFFD_a)) );
}


void static_ic_set_sinusoid_3d(BSSN * bssn, Static * stat, Lambda * lambda, Fourier * fourier,
  IOData * iodata)
{

  idx_t i, j, k;

  arr_t & DIFFr_a = *bssn->fields["DIFFr_a"];
  arr_t & DIFFphi_p = *bssn->fields["DIFFphi_p"];
  arr_t & DIFFD_a = *stat->fields["DIFFD_a"];

  real_t rho_FRW = 3.0/PI/8.0;

  real_t Omega_L = std::stod(_config("Omega_L", "0.0"));
  real_t rho_m = (1.0 - Omega_L) * rho_FRW;
  real_t rho_L = Omega_L * rho_FRW;
  lambda->setLambda(rho_L);

  real_t A = H_LEN_FRAC*H_LEN_FRAC*std::stod(_config("peak_amplitude_frac", "0.001"));

  // the conformal factor in front of metric is the solution to
  // d^2 exp(\phi) = -2*pi exp(5\phi) * \delta_rho
  // generate random mode in \phi
  // delta_rho = -(lap e^\phi)/e^(4\phi)/2pi
  real_t phix = std::stod(_config("phix", "0.0"));
  real_t phiy = phix, phiz = phix;

  // grid values
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    real_t x_frac = ((real_t) i / (real_t) NX),
      y_frac = ((real_t) j / (real_t) NY), z_frac = ((real_t) k / (real_t) NZ);
    
    real_t phi = A*
      (sin(2.0*PI*x_frac + phix) + sin(2.0*PI*y_frac + phiy) + sin(2.0*PI*z_frac + phiz));

    real_t DIFFrho = -exp(-4.0*phi)/PI/2.0*(
      pw2(A*2.0*PI / H_LEN_FRAC * cos(2.0*PI*x_frac + phix))
      + pw2(A*2.0*PI / H_LEN_FRAC * cos(2.0*PI*y_frac + phiy))
      + pw2(A*2.0*PI / H_LEN_FRAC * cos(2.0*PI*z_frac + phiz))
      -  ( A * pw2(2.0*PI / H_LEN_FRAC) * sin(2.0*PI*x_frac + phix)
           + A * pw2(2.0*PI / H_LEN_FRAC) * sin(2.0*PI*y_frac + phiy)
           + A * pw2(2.0*PI / H_LEN_FRAC) * sin(2.0*PI*z_frac + phiz))
    );

    // These aren't difference vars
    DIFFphi_p[NP_INDEX(i,j,k)] = phi;
    DIFFr_a[idx] = DIFFrho;

// // debugging: throw away field
// DIFFD_a[idx] = phi;
  }
// std::cout << std::setprecision(17);
// std::cout << "field[0] = " << DIFFD_a[0] << "; ";
// fourier->inverseLaplacian <idx_t, real_t> (DIFFD_a._array);
// std::cout << "lap/lap field = " << laplacian(0,0,0,DIFFD_a) << "\n";

  // Make sure min density value > 0
  // Set conserved density variable field
  real_t min = rho_m;
  real_t max = min;
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    real_t DIFFr = DIFFr_a[idx];
    real_t rho = rho_m + DIFFr;
    // phi_FRW = 0
    real_t DIFFphi = DIFFphi_p[idx];
    // phi = DIFFphi
    // DIFFK = 0

    DIFFD_a[idx] =
      rho_m*expm1(6.0*DIFFphi) + exp(6.0*DIFFphi)*DIFFr;

    if(rho < min)
    {
      min = rho;
    }
    if(rho > max)
    {
      max = rho;
    }
    if(rho != rho)
    {
      iodata->log("Error: NaN energy density.");
      throw -1;
    }
  }

  iodata->log( "Minimum fluid density: " + stringify(min) );
  iodata->log( "Maximum fluid density: " + stringify(max) );
  iodata->log( "Average fluctuation density: " + stringify(average(DIFFD_a)) );
  iodata->log( "Std.dev fluctuation density: " + stringify(standard_deviation(DIFFD_a)) );
  if(min < 0.0)
  {
    iodata->log("Error: negative density in some regions.");
    throw -1;
  }

# if USE_REFERENCE_FRW
  // Set values in reference FRW integrator
  auto & frw = bssn->frw;
  real_t K_FRW = -sqrt(24.0*PI*rho_FRW);
  frw->set_phi(0.0);
  frw->set_K(K_FRW);
  frw->addFluid(rho_m, 0.0 /* w=0 */);
# else
  arr_t & DIFFK_p = *bssn->fields["DIFFK_p"];
  arr_t & DIFFK_a = *bssn->fields["DIFFK_a"];
  // add in FRW pieces to ICs
  // phi is unchanged
  // rho (D) and K get contribs
  // w=0 fluid only
# pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    real_t D_FRW = rho_m; // on initial slice

    DIFFr_a[idx] += rho_m;

    DIFFK_a[idx] = -sqrt(24.0*PI*rho_FRW);
    DIFFK_p[idx] = -sqrt(24.0*PI*rho_FRW);

    DIFFD_a[idx] += D_FRW;
  }
# endif

}


/**
 * @brief Sinusoidal mode ICs
 */
void static_ic_set_sinusoid(BSSN * bssn, Static * stat, Lambda * lambda, Fourier * fourier,
  IOData * iodata)
{
  idx_t i, j, k;

  arr_t & DIFFr_a = *bssn->fields["DIFFr_a"];
  arr_t & DIFFphi_p = *bssn->fields["DIFFphi_p"];
  arr_t & DIFFD_a = *stat->fields["DIFFD_a"];

  real_t rho_FRW = 3.0/PI/8.0;

  real_t Omega_L = std::stod(_config("Omega_L", "0.0"));
  real_t rho_m = (1.0 - Omega_L) * rho_FRW;
  real_t rho_L = Omega_L * rho_FRW;
  lambda->setLambda(rho_L);

  real_t A = H_LEN_FRAC*H_LEN_FRAC*std::stod(_config("peak_amplitude_frac", "0.001"));

  // the conformal factor in front of metric is the solution to
  // d^2 exp(\phi) = -2*pi exp(5\phi) * \delta_rho
  // generate random mode in \phi
  // delta_rho = -(lap e^\phi)/e^(4\phi)/2pi
  real_t phix = std::stod(_config("phix", "0.0"));
  real_t twopi_L = 2.0*PI/H_LEN_FRAC;
  real_t pw2_twopi_L = twopi_L*twopi_L;
  // grid values
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    real_t x = ((real_t) i / (real_t) NX);
    real_t phi = A*sin(2.0*PI*x + phix);
    real_t DIFFrho = -exp(-4.0*phi)/PI/2.0*(
        pw2(twopi_L*A*cos(2.0*PI*x + phix))
        - pw2_twopi_L*A*sin(2.0*PI*x + phix)
      );

    // These aren't difference vars
    DIFFphi_p[NP_INDEX(i,j,k)] = phi;
    DIFFr_a[idx] = DIFFrho;

  }

  // Make sure min density value > 0
  // Set conserved density variable field
  real_t min = rho_m;
  real_t max = min;
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    real_t DIFFr = DIFFr_a[idx];
    real_t rho = rho_m + DIFFr;
    // phi_FRW = 0
    real_t DIFFphi = DIFFphi_p[idx];
    // phi = DIFFphi
    // DIFFK = 0

    DIFFD_a[idx] =
      rho_m*expm1(6.0*DIFFphi) + exp(6.0*DIFFphi)*DIFFr;

    if(rho < min)
    {
      min = rho;
    }
    if(rho > max)
    {
      max = rho;
    }
    if(rho != rho)
    {
      iodata->log("Error: NaN energy density.");
      throw -1;
    }
  }

  iodata->log( "Minimum fluid density: " + stringify(min) );
  iodata->log( "Maximum fluid density: " + stringify(max) );
  iodata->log( "Average fluctuation density: " + stringify(average(DIFFD_a)) );
  iodata->log( "Std.dev fluctuation density: " + stringify(standard_deviation(DIFFD_a)) );
  if(min < 0.0)
  {
    iodata->log("Error: negative density in some regions.");
    throw -1;
  }

# if USE_REFERENCE_FRW
  // Set values in reference FRW integrator
  auto & frw = bssn->frw;
  real_t K_FRW = -sqrt(24.0*PI*rho_FRW);
  frw->set_phi(0.0);
  frw->set_K(K_FRW);
  frw->addFluid(rho_m, 0.0 /* w=0 */);
# else
  arr_t & DIFFK_p = *bssn->fields["DIFFK_p"];
  arr_t & DIFFK_a = *bssn->fields["DIFFK_a"];
  // add in FRW pieces to ICs
  // phi is unchanged
  // rho (D) and K get contribs
  // w=0 fluid only
# pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    real_t D_FRW = rho_m; // on initial slice

    DIFFr_a[idx] += rho_m;

    DIFFK_a[idx] = -sqrt(24.0*PI*rho_FRW);
    DIFFK_p[idx] = -sqrt(24.0*PI*rho_FRW);

    DIFFD_a[idx] += D_FRW;
  }
# endif
}


void static_ic_set_semianalytic(
  BSSN * bssn, Static * stat, Lambda * lambda, Fourier * fourier,
  IOData * iodata)
{
  idx_t i, j, k;

  arr_t & DIFFr_a = *bssn->fields["DIFFr_a"];
  arr_t & DIFFphi_p = *bssn->fields["DIFFphi_p"];
  arr_t & DIFFD_a = *stat->fields["DIFFD_a"];
  arr_t & DIFFK_p = *bssn->fields["DIFFK_p"];

  real_t rho_FRW = 3.0/PI/8.0;
  real_t K_FRW = -3.0;

  real_t Omega_L = std::stod(_config("Omega_L", "0.0"));
  real_t rho_m = (1.0 - Omega_L) * rho_FRW;
  real_t rho_L = Omega_L * rho_FRW;
  lambda->setLambda(rho_L);
  real_t L = H_LEN_FRAC;
  real_t A = std::stod(_config("peak_amplitude", "0.001"))*0.026699*L*L;
  // grid values
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    real_t x = ((real_t) i) * dx;
    DIFFphi_p[idx] = std::log1p( A*std::sin(2.0*PI*x/L) );
    DIFFr_a[idx] = rho_m + 2.0*A*PI*std::sin((2*PI*x)/L)
      / ( L*L*std::pow(1 + A*std::sin((2*PI*x)/L), 5) );
    DIFFK_p[idx] = K_FRW;
    DIFFD_a[idx] = exp(6.0*DIFFphi_p[idx])*DIFFr_a[idx];
  }

  // Make sure min density value > 0
  // Set conserved density variable field
  real_t min = rho_m;
  real_t max = min;
  LOOP3(i,j,k)
  {
    real_t rho = DIFFr_a[NP_INDEX(i,j,k)];
    if(rho < min)
    {
      min = rho;
    }
    if(rho > max)
    {
      max = rho;
    }
    if(rho != rho)
    {
      iodata->log("Error: NaN energy density.");
      throw -1;
    }
  }

  iodata->log( "Minimum fluid density: " + stringify(min) );
  iodata->log( "Maximum fluid density: " + stringify(max) );
  iodata->log( "Average fluctuation density: " + stringify(average(DIFFD_a)) );
  iodata->log( "Std.dev fluctuation density: " + stringify(standard_deviation(DIFFD_a)) );
  if(min < 0.0)
  {
    iodata->log("Error: negative density in some regions.");
    throw -1;
  }
}



/**
 * @brief Spherical "shell" of perturbations around an observer
 */
void static_ic_set_sphere(BSSN * bssn, Static * stat, IOData * iodata)
{
  idx_t i, j, k;

  arr_t & DIFFr_a = *bssn->fields["DIFFr_a"];
  arr_t & DIFFphi_p = *bssn->fields["DIFFphi_p"];
  arr_t & DIFFphi_a = *bssn->fields["DIFFphi_a"];
  arr_t & DIFFphi_f = *bssn->fields["DIFFphi_f"];

  arr_t & DIFFD_a = *stat->fields["DIFFD_a"];

  // shell amplitude
  const real_t A = stod(_config("shell_amplitude", "1e-5"));
  // Shell described by only one fixed l:
  const idx_t l = stoi(_config("shell_angular_scale_l", "1"));
  iodata->log( "Generating ICs with shell angular scale of l = " + stringify(l) );
  iodata->log( "Generating ICs with peak amp. = " + stringify(A) );

  // spherical shell of perturbations in phi0field

  // place shell around center of box
  real_t x0 = (NX-0.001)*dx/2.0;
  real_t y0 = (NY-0.001)*dx/2.0;
  real_t z0 = (NZ-0.001)*dx/2.0;

  // place spherical "shell" of fluctuations at r = NX/4, 1/2-way between observer and boundary 
  real_t r_shell = NX*dx / 4.0;
  // shell width
  real_t shell_width = NX*dx / 40.0;

  // Angular fluctuations in shell described by spherical harmonic coeffs, a_lm's,
  complex_t * alms = new complex_t[m_idx(l,l)+1];
  const real_t seed = stod(_config("mt19937_seed", "7"));
  std::mt19937 gen(seed);
  std::normal_distribution<> normal_dist(0.0, 1.0);
  std::uniform_real_distribution<> uniform_dist(0.0, 2.0*PI);

  std::cout << "normal_dist(gen) = " << normal_dist(gen) << ", uniform_dist(gen) = " << uniform_dist(gen) << "\n";

  // zero mode:
  alms[m_idx(l, 0)].first = normal_dist(gen);
  alms[m_idx(l, 0)].second = 0;
  // positive modes:
  for(int m = 1; m <= l; m++)
  {
    real_t phase = uniform_dist(gen);
    real_t amp = normal_dist(gen);
    alms[m_idx(l,m)].first = amp*std::cos(phase);
    alms[m_idx(l,m)].second = amp*std::sin(phase);
  }
  // negative modes:
  for(int m = -l; m <= -1; m++)
  {
    real_t Condon_Shortley_phase = std::abs(m) % 2 ? -1.0 : 1.0; // 0 (false) if m even, 1 (true) if odd
    alms[m_idx(l,m)].first = Condon_Shortley_phase*alms[m_idx(l,std::abs(m))].first;
    alms[m_idx(l,m)].second = -Condon_Shortley_phase*alms[m_idx(l,std::abs(m))].second;
  }

  for(int m = -l; m <= l; m++)
  {
    std::cout << "Amp. of a_{" << l << "," << m << "} = " << alms[m_idx(l,m)].first
      << " + " << alms[m_idx(l,m)].second << "i\n";
  }

  LOOP3(i,j,k) {
    idx_t idx = NP_INDEX(i,j,k);

    real_t x = i*dx;
    real_t y = j*dx;
    real_t z = k*dx;

    real_t r = sqrt( pw2(x - x0) + pw2(y - y0) + pw2(z - z0) );
    real_t theta = acos((z - z0) / r); // "theta" (polar angle)
    real_t phi = atan2(y - y0, x - x0); // "phi" (azimuthal angle)

    real_t DIFFphi_r = 0.0;
    real_t DIFFphi_i = 0.0;
    for(int m=-l; m<=l; m++)
    {
      real_t Y_r = boost::math::spherical_harmonic_r(l, m, theta, phi);
      real_t Y_i = boost::math::spherical_harmonic_i(l, m, theta, phi);

      DIFFphi_r += alms[m_idx(l,m)].first*Y_r - alms[m_idx(l,m)].second*Y_i;
      DIFFphi_i += alms[m_idx(l,m)].first*Y_i + alms[m_idx(l,m)].second*Y_r;
    }
    if(std::abs(DIFFphi_i) > 1e-6)
    {
      iodata->log("Significant non-zero imaginary component of solution exists!");
      throw -1;
    }

    // gaussian profile shell of fluctuations
    real_t U_r = A*std::exp( -pw2((r - r_shell)/2.0/shell_width) );
    // cosine profile
    // real_t U_r = (r < r_shell-shell_width || r > r_shell+shell_width ) ? 0 : A*(1+std::cos(PI*(r-r_shell)/shell_width));

    DIFFphi_p[idx] = U_r*DIFFphi_r;


    if(i==NX/2 && j==NY/2 && k==NZ/2)
      std::cout << "At a grid point close to the middle, r = " << r
                << ", r_shell = " << r_shell
                << ", shell_width = " << shell_width
                << ", U_r = " << U_r
                << ", DIFFphi_p = " << DIFFphi_p[idx] << "\n";
  }
  // cleanup
  delete [] alms;

  // delta_rho = -lap(phi)/(1+xi)^5/2pi
# pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k) {
    DIFFr_a[NP_INDEX(i,j,k)] = -0.5/PI/(
      pow(1.0 + DIFFphi_p[NP_INDEX(i,j,k)], 5.0)
    )*(
      double_derivative(i, j, k, 1, 1, DIFFphi_p)
      + double_derivative(i, j, k, 2, 2, DIFFphi_p)
      + double_derivative(i, j, k, 3, 3, DIFFphi_p)
    );
  }

  // phi = ln(xi)
# pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k) {
    idx_t idx = NP_INDEX(i,j,k);
    DIFFphi_a[idx] = log1p(DIFFphi_p[idx]);
    DIFFphi_f[idx] = log1p(DIFFphi_p[idx]);
    DIFFphi_p[idx] = log1p(DIFFphi_p[idx]);
  }

  // Make sure min density value > 0
  // Set conserved density variable field
  real_t min = 3.0/PI/8.0;
  real_t max = min;
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    real_t rho_FRW = 3.0/PI/8.0;
    real_t DIFFr = DIFFr_a[idx];
    real_t rho = rho_FRW + DIFFr;
    // phi_FRW = 0
    real_t DIFFphi = DIFFphi_a[idx];
    // phi = DIFFphi
    // DIFFK = 0

    DIFFD_a[idx] =
      rho_FRW*expm1(6.0*DIFFphi) + exp(6.0*DIFFphi)*DIFFr;

    if(rho < min)
    {
      min = rho;
    }
    if(rho > max)
    {
      max = rho;
    }
    if(rho != rho)
    {
      iodata->log("Error: NaN energy density.");
      throw -1;
    }
  }

  iodata->log( "Minimum fluid density: " + stringify(min) );
  iodata->log( "Maximum fluid density: " + stringify(max) );
  iodata->log( "Average fluctuation density: " + stringify(average(DIFFD_a)) );
  iodata->log( "Std.dev fluctuation density: " + stringify(standard_deviation(DIFFD_a)) );
  if(min < 0.0)
  {
    iodata->log("Error: negative density in some regions.");
    throw -1;
  }

# if USE_REFERENCE_FRW
  // Set values in reference FRW integrator
  real_t rho_FRW = 3.0/PI/8.0;
  real_t K_frw = -sqrt(24.0*PI*rho_FRW);

  auto & frw = bssn->frw;
  frw->set_phi(0.0);
  frw->set_K(K_frw);
  frw->addFluid(rho_FRW, 0.0 /* w=0 */);
# else
  arr_t & DIFFK_p = *bssn->fields["DIFFK_p"];
  arr_t & DIFFK_a = *bssn->fields["DIFFK_a"];
  // add in FRW pieces to ICs
  // phi is unchanged
  // rho (D) and K get contribs
  // w=0 fluid only
# pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    real_t rho_FRW = 3.0/PI/8.0;
    real_t D_FRW = rho_FRW; // on initial slice

    DIFFr_a[idx] += rho_FRW;

    DIFFK_a[idx] = -sqrt(24.0*PI*rho_FRW);
    DIFFK_p[idx] = -sqrt(24.0*PI*rho_FRW);

    DIFFD_a[idx] += D_FRW;
  }
# endif
}

}
