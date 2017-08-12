#include "bardeen.h"
#include "../../utils/math.h"

namespace cosmo
{

void Bardeen::setPotentials()
{
  idx_t i, j, k;

  // compute conformal factor, time-derivatives (assumes dust universe)
  arr_t & DIFFphi_a = *bssn->fields["DIFFphi_a"];
  real_t a = exp( 2.0*( bssn->frw->get_phi() + conformal_average(DIFFphi_a, DIFFphi_a, bssn->frw->get_phi()) ) );
  real_t dadt = 1.0/std::sqrt(a); // reliant upon on dust universe
  real_t d2adt2 = -1.0/2.0/a/a; // reliant upon on dust universe

  // construct h_ij components, time derivatives
# pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    BSSNData bd = {0};
    bssn->set_bd_values(i, j, k, &bd); // TODO: remove redundant computations here?

    real_t e4phi = std::exp(4.0*bd.phi);
    h11[idx] = e4phi*bd.gamma11 - a*a;
    h22[idx] = e4phi*bd.gamma22 - a*a;
    h33[idx] = e4phi*bd.gamma33 - a*a;
    h12[idx] = e4phi*bd.gamma12;
    h13[idx] = e4phi*bd.gamma13;
    h23[idx] = e4phi*bd.gamma23;

    // dt gamma_ij = -2.0*K_ij
    // K_{ij} = e^{4\phi} * A_{ij} + gamma_{ij}*K/3
    dt_h11[idx] = -2.0*e4phi*(bd.A11 + bd.gamma11*bd.K/3.0) - 2.0*a*dadt;
    dt_h22[idx] = -2.0*e4phi*(bd.A22 + bd.gamma22*bd.K/3.0) - 2.0*a*dadt;
    dt_h33[idx] = -2.0*e4phi*(bd.A33 + bd.gamma33*bd.K/3.0) - 2.0*a*dadt;
    dt_h12[idx] = -2.0*e4phi*(bd.A12 + bd.gamma12*bd.K/3.0);
    dt_h13[idx] = -2.0*e4phi*(bd.A13 + bd.gamma13*bd.K/3.0);
    dt_h23[idx] = -2.0*e4phi*(bd.A23 + bd.gamma23*bd.K/3.0);

    // dt K_{ij}
    real_t ev_phi = bssn->ev_DIFFK(&bd) + -bd.K_FRW/6.0;
    real_t ev_K = bssn->ev_DIFFK(&bd) + bd.K_FRW*bd.K_FRW/3.0 + 4.0*PI*bd.rho_FRW;
    d2t_h11[idx] = -2.0*( 4.0*(ev_phi)*e4phi*(bd.A11 + bd.gamma11*bd.K/3.0)
                    + e4phi*( bssn->ev_A11(&bd) + bssn->ev_DIFFgamma11(&bd)*bd.K/3.0 + bd.gamma11*ev_K/3.0 )
                    + dadt*dadt + a*d2adt2 );
    d2t_h22[idx] = -2.0*( 4.0*(ev_phi)*e4phi*(bd.A22 + bd.gamma22*bd.K/3.0)
                    + e4phi*( bssn->ev_A22(&bd) + bssn->ev_DIFFgamma22(&bd)*bd.K/3.0 + bd.gamma22*ev_K/3.0 )
                    + dadt*dadt + a*d2adt2 );
    d2t_h33[idx] = -2.0*( 4.0*(ev_phi)*e4phi*(bd.A33 + bd.gamma33*bd.K/3.0)
                    + e4phi*( bssn->ev_A33(&bd) + bssn->ev_DIFFgamma33(&bd)*bd.K/3.0 + bd.gamma33*ev_K/3.0 )
                    + dadt*dadt + a*d2adt2 );
    d2t_h12[idx] = -2.0*( 4.0*(ev_phi)*e4phi*(bd.A12 + bd.gamma12*bd.K/3.0)
                    + e4phi*( bssn->ev_A12(&bd) + bssn->ev_DIFFgamma12(&bd)*bd.K/3.0 + bd.gamma12*ev_K/3.0 ) );
    d2t_h13[idx] = -2.0*( 4.0*(ev_phi)*e4phi*(bd.A13 + bd.gamma13*bd.K/3.0)
                    + e4phi*( bssn->ev_A13(&bd) + bssn->ev_DIFFgamma13(&bd)*bd.K/3.0 + bd.gamma13*ev_K/3.0 ) );
    d2t_h23[idx] = -2.0*( 4.0*(ev_phi)*e4phi*(bd.A23 + bd.gamma23*bd.K/3.0)
                    + e4phi*( bssn->ev_A23(&bd) + bssn->ev_DIFFgamma23(&bd)*bd.K/3.0 + bd.gamma23*ev_K/3.0 ) );
  }

  // construct A (and its time derivatives) in increments:
  // (A.1) construct d_i d_j h_{ij}
# pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    A[idx] = (
      double_derivative(i, j, k, 1, 1, h11) + double_derivative(i, j, k, 2, 2, h22) + double_derivative(i, j, k, 3, 3, h33)
      + 2.0*(double_derivative(i, j, k, 1, 2, h12) + double_derivative(i, j, k, 1, 3, h13) + double_derivative(i, j, k, 2, 3, h23))
    );
    dt_A[idx] = (
      double_derivative(i, j, k, 1, 1, dt_h11) + double_derivative(i, j, k, 2, 2, dt_h22) + double_derivative(i, j, k, 3, 3, dt_h33)
      + 2.0*(double_derivative(i, j, k, 1, 2, dt_h12) + double_derivative(i, j, k, 1, 3, dt_h13) + double_derivative(i, j, k, 2, 3, dt_h23))
    );
    d2t_A[idx] = (
      double_derivative(i, j, k, 1, 1, d2t_h11) + double_derivative(i, j, k, 2, 2, d2t_h22) + double_derivative(i, j, k, 3, 3, d2t_h33)
      + 2.0*(double_derivative(i, j, k, 1, 2, d2t_h12) + double_derivative(i, j, k, 1, 3, d2t_h13) + double_derivative(i, j, k, 2, 3, d2t_h23))
    );
  }
  // (A.2) compute inverse laplacian of (A.1)
  fourier->inverseLaplacian <idx_t, real_t> (A._array);
  fourier->inverseLaplacian <idx_t, real_t> (dt_A._array);
  fourier->inverseLaplacian <idx_t, real_t> (d2t_A._array);
  // (A.3) subtract from trace and /(2a^2)
# pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    
    real_t h_tr = h11[idx] + h22[idx] + h33[idx];
    real_t dt_h_tr = dt_h11[idx] + dt_h22[idx] + dt_h33[idx];
    real_t d2t_h_tr = d2t_h11[idx] + d2t_h22[idx] + d2t_h33[idx];

    // defn's for semantic clarity
    real_t djdk_d2_hjk = A[idx];
    real_t dt_djdk_d2_hjk = dt_A[idx];
    real_t d2t_djdk_d2_hjk = d2t_A[idx];
    
    // final values
    A[idx] = ( h_tr - djdk_d2_hjk )/a/a/2.0;
    dt_A[idx] = -2.0*A[idx]*dadt/a + ( dt_h_tr - dt_djdk_d2_hjk )/a/a/2.0;
    d2t_A[idx] =  -2.0*A[idx]*( pw2(dadt/a) + d2adt2/a ) - 4.0*dadt/a*dt_A[idx]
                  + ( d2t_h_tr - d2t_djdk_d2_hjk )/a/a/2.0;
  }
  // fix monopoles: <h_tr> / a^2 = <3A> + <d^2 B> = <3A> (periodic spacetime)
  //             => <A> = <h_tr> / 3 / a^2
  real_t avg_A = average(A);
  real_t avg_dt_A = average(dt_A);
  real_t avg_d2t_A = average(d2t_A);
  real_t avg_h_tr = 0, avg_dt_h_tr = 0, avg_d2t_h_tr = 0;
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    avg_h_tr += h11[idx] + h22[idx] + h33[idx];
    avg_dt_h_tr += dt_h11[idx] + dt_h22[idx] + dt_h33[idx];
    avg_d2t_h_tr += d2t_h11[idx] + d2t_h22[idx] + d2t_h33[idx];
  }
  avg_h_tr /= POINTS;
  avg_dt_h_tr /= POINTS;
  avg_d2t_h_tr /= POINTS;

  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    A[idx] += avg_h_tr/3.0/a/a - avg_A;
    dt_A[idx] += avg_dt_h_tr/3.0/a/a - 2.0*avg_h_tr*dadt/3.0/a/a/a - avg_dt_A;
    d2t_A[idx] += avg_d2t_h_tr/3.0/a/a - 2.0*dadt*avg_dt_h_tr/3.0/a/a/a
     + 6.0*avg_h_tr*dadt*dadt/3.0/a/a/a/a - 2.0*avg_h_tr*d2adt2/3.0/a/a/a - 2.0*avg_dt_h_tr*dadt/3.0/a/a/a
     - avg_d2t_A;
  }

  // construct B (and its time derivatives) in increments:
  // (B.1) trace - 3A
# pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    
    real_t h_tr = h11[idx] + h22[idx] + h33[idx];
    real_t dt_h_tr = dt_h11[idx] + dt_h22[idx] + dt_h33[idx];
    real_t d2t_h_tr = d2t_h11[idx] + d2t_h22[idx] + d2t_h33[idx];
    
    B[idx] = ( h_tr/a/a - 3.0*A[idx] );
    dt_B[idx] = ( dt_h_tr/a/a - 2.0/a/a/a*dadt*h_tr - 3.0*dt_A[idx] );
    d2t_B[idx] = ( d2t_h_tr/a/a - 4.0/a/a/a*dadt*dt_h_tr
      - 2.0/a/a/a*d2adt2*h_tr + 6.0/a/a/a/a*dadt*dadt*h_tr - 3.0*d2t_A[idx] );
  }

  // inverse laplacian of
  fourier->inverseLaplacian <idx_t, real_t> (B._array);
  fourier->inverseLaplacian <idx_t, real_t> (dt_B._array);
  fourier->inverseLaplacian <idx_t, real_t> (d2t_B._array);

  // scalar metric fields obtained... get Bardeen potentials:
# pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    Phi[idx] = -a/2.0*( 2.0*dadt*dt_B[idx] + a*d2t_B[idx] );
    Psi[idx] = -1.0/2.0*A[idx] + a*dadt*dt_B[idx]/2.0;
  }

}

}
