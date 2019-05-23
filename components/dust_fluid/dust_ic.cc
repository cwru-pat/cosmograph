#include "dust_ic.h"
#include "../../cosmo_includes.h"
#include "../../cosmo_types.h"
#include "../../cosmo_globals.h"
#include "../../utils/Fourier.h"
#include "../../utils/math.h"
#include "../../ICs/ICs.h"

namespace cosmo
{

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
void dust_ic_set_random( BSSN * bssn, Dust * dust, Lambda * lambda,
  Fourier * fourier, IOData * iodata )
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

  // initial metric fields
  arr_t & DIFFphi_p = *bssn->fields["DIFFphi_p"];
  arr_t & DIFFK_p = *bssn->fields["DIFFK_p"];
  arr_t & A11_p = *bssn->fields["A11_p"];
  arr_t & A12_p = *bssn->fields["A12_p"];
  arr_t & A13_p = *bssn->fields["A13_p"];
  arr_t & A22_p = *bssn->fields["A22_p"];
  arr_t & A23_p = *bssn->fields["A23_p"];
  arr_t & A33_p = *bssn->fields["A33_p"];
  arr_t & D_p = dust->D._array_p;

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

    D_p[idx] = exp(6.0*DIFFphi_p[idx])*rho_0;
  }

  iodata->log( "Average fluctuation density: " + stringify(average(D_p)) );
  iodata->log( "Std.dev fluctuation density: " + stringify(standard_deviation(D_p)) );
}

} // namespace cosmo
