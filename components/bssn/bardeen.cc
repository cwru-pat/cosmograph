#include "bardeen.h"
#include "../../utils/math.h"

namespace cosmo
{

real_t abs_dder(idx_t i, idx_t j, idx_t k, arr_t & field)
{
  return std::abs(double_derivative(i,j,k,1,1,field))
   + std::abs(double_derivative(i,j,k,1,2,field))
   + std::abs(double_derivative(i,j,k,1,3,field))
   + std::abs(double_derivative(i,j,k,2,2,field))
   + std::abs(double_derivative(i,j,k,2,3,field))
   + std::abs(double_derivative(i,j,k,3,3,field));
}

/**
 * @brief      Compute Bardeen & vector potentials.
 * Assumes no reference solultion is used (this should be checked in the
 * constructor).
 */
void Bardeen::setPotentials(real_t elapsed_sim_time)
{
  idx_t i, j, k;

  if( bssn->frw->get_K() != 0 || bssn->frw->get_phi() != 0)
  {
    std::cout << "Bardeen variables incompatible with use of reference metric."
     << " SVT and Bardeen fields will not be computed.";
    return;
  }

  // compute conformal factor, time-derivatives (assumes dust universe)
  arr_t & DIFFphi_a = *bssn->fields["DIFFphi_a"];
  arr_t & DIFFalpha_a = *bssn->fields["DIFFalpha_a"];
  arr_t & DIFFK_a = *bssn->fields["DIFFK_a"];
#if USE_Z4c_DAMPING
  arr_t & theta_a = *bssn->fields["theta_a"];
#endif
  
  // TODO: set these differently? Explore different choices?
  real_t alpha_avg = 1.0 + conformal_average(DIFFalpha_a, DIFFphi_a, 0.0);
  real_t K_avg = conformal_average(DIFFK_a, DIFFphi_a, 0.0);
  real_t phi_avg = conformal_average(DIFFphi_a, DIFFphi_a, 0.0);

#if USE_Z4c_DAMPING
  K_avg += conformal_average(theta_a, DIFFphi_a, 0.0);
#endif
  
  real_t a, dadt, H, d2adt2;
  if(use_matter_scale_factor)
  {
    real_t tI = 2.0/3.0; // t_I = 2/(3H) = 2/3 in units where H_I = 1 (sim units)
    real_t t = tI + elapsed_sim_time;
    a = std::pow(t/tI, 2.0/3.0);
    H = 2.0/3.0/t;
    dadt = H*a;
    d2adt2 = a*(H*H - 2.0/3.0/t/t);
  }
  else
  {
    // define a = e^(2<\phi>)
    a = exp( 2.0*phi_avg );
    // a' ~ e^(2<\phi>)' ~ a*2\phi' ~ -1/3*a*<alpha>*<K>
    dadt = -1.0/3.0*a*alpha_avg*K_avg;
    H = dadt/a;
    // a'' ~ H*a' + a*2\phi''
    // Set after \phi'' is computed,
    // Don't use d2adt2 this until then. (valgrind can catch.)
  }

  // Store time-derivative of BSSN metric components for later differentiation
#pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    BSSNData bd = {0};
    bssn->set_bd_values(i, j, k, &bd);

    // stores d/dt \bar{\gamma}_ij
    dt_g11[idx] = bssn->ev_DIFFgamma11(&bd);
    dt_g12[idx] = bssn->ev_DIFFgamma12(&bd);
    dt_g13[idx] = bssn->ev_DIFFgamma13(&bd);
    dt_g22[idx] = bssn->ev_DIFFgamma22(&bd);
    dt_g23[idx] = bssn->ev_DIFFgamma23(&bd);
    dt_g33[idx] = bssn->ev_DIFFgamma33(&bd);

    // stores d/dt \beta^i
#if USE_BSSN_SHIFT
    dt_beta1[idx] = bssn->ev_beta1(&bd);
    dt_beta2[idx] = bssn->ev_beta2(&bd);
    dt_beta3[idx] = bssn->ev_beta3(&bd);
#else
    dt_beta1[idx] = 0.0;
    dt_beta2[idx] = 0.0;
    dt_beta3[idx] = 0.0;
#endif

    // stores d/dt phi
    dt_phi[idx] = bssn->ev_DIFFphi(&bd);
  }


  // Construct d2t \bar{\gamma}_ij, d2t phi
#pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    BSSNData bd = {0};
    bssn->set_bd_values(i, j, k, &bd);

    real_t dkbetak = bd.d1beta1 + bd.d2beta2 + bd.d3beta3;
    real_t dkdtbetak = derivative(i,j,k,1,dt_beta1) + derivative(i,j,k,2,dt_beta2) + derivative(i,j,k,3,dt_beta3);
    real_t dtalpha = bssn->ev_DIFFalpha(&bd);

    // stores d^2/dt^2 \bar{\gamma}_ij, incl. Macro for calc.
#define dt_g21 dt_g12
#define dt_g31 dt_g13
#define dt_g32 dt_g23
#define D2T_g(I,J) -2.0*dtalpha*bd.A##I##J - 2.0*bd.alpha*bssn->ev_A##I##J(&bd) \
      + dt_beta1[idx]*bd.d1g##I##J + dt_beta2[idx]*bd.d2g##I##J + dt_beta3[idx]*bd.d3g##I##J \
      + bd.beta1*derivative(i,j,k,1,dt_g##I##J) + bd.beta2*derivative(i,j,k,2,dt_g##I##J) + bd.beta3*derivative(i,j,k,3,dt_g##I##J) \
      + dt_g##I##1[idx]*bd.d##J##beta1 + dt_g##I##2[idx]*bd.d##J##beta2 + dt_g##I##3[idx]*bd.d##J##beta3 \
      + dt_g##J##1[idx]*bd.d##I##beta1 + dt_g##J##2[idx]*bd.d##I##beta2 + dt_g##J##3[idx]*bd.d##I##beta3 \
      + bd.gamma##I##1*derivative(i,j,k,J,dt_beta1) + bd.gamma##I##2*derivative(i,j,k,J,dt_beta2) + bd.gamma##I##3*derivative(i,j,k,J,dt_beta3) \
      + bd.gamma##J##1*derivative(i,j,k,I,dt_beta1) + bd.gamma##J##2*derivative(i,j,k,I,dt_beta2) + bd.gamma##J##3*derivative(i,j,k,I,dt_beta3) \
      - 2.0/3.0*( dt_g##I##J[idx]*dkbetak + bd.gamma##I##J*dkdtbetak )

    d2t_g11[idx] = D2T_g(1, 1);
    d2t_g12[idx] = D2T_g(1, 2);
    d2t_g13[idx] = D2T_g(1, 3);
    d2t_g22[idx] = D2T_g(2, 2);
    d2t_g23[idx] = D2T_g(2, 3);
    d2t_g33[idx] = D2T_g(3, 3);

    // stores d^2/dt^2 phi
    d2t_phi[idx] = -1.0/6.0*( dtalpha*(bd.K + 2.0 * bd.theta) + bd.alpha*bssn->ev_DIFFK(&bd) )
      + dt_beta1[idx]*bd.d1phi + dt_beta2[idx]*bd.d2phi + dt_beta3[idx]*bd.d3phi
      + bd.beta1*derivative(i,j,k,1,dt_phi) + bd.beta2*derivative(i,j,k,2,dt_phi) + bd.beta3*derivative(i,j,k,3,dt_phi)
      + 1.0/6.0*dkdtbetak;
  }
  // set second derivative of a
  // a'' ~ H*a' + a*2\phi''
  if(!use_matter_scale_factor)
    d2adt2 = H*dadt + 2.0*a*conformal_average(d2t_phi, DIFFphi_a, 0.0);

  // construct h_ij, h_0i components, time derivatives
#pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    BSSNData bd = {0};
    bssn->set_bd_values(i, j, k, &bd); // TODO: remove redundant computations here?

    real_t e4phi = std::exp(4.0*bd.phi);

    // h_ij
    h11[idx] = e4phi*bd.gamma11 - a*a;
    h22[idx] = e4phi*bd.gamma22 - a*a;
    h33[idx] = e4phi*bd.gamma33 - a*a;
    h12[idx] = e4phi*bd.gamma12;
    h13[idx] = e4phi*bd.gamma13;
    h23[idx] = e4phi*bd.gamma23;

    // dt h_ij
#define DT_h(I,J) e4phi*(4.0*dt_phi[idx]*bd.gamma##I##J + dt_g##I##J[idx])
    dt_h11[idx] = DT_h(1,1) - 2.0*a*dadt;
    dt_h22[idx] = DT_h(2,2) - 2.0*a*dadt;
    dt_h33[idx] = DT_h(3,3) - 2.0*a*dadt;
    dt_h12[idx] = DT_h(1,2);
    dt_h13[idx] = DT_h(1,3);
    dt_h23[idx] = DT_h(2,3);

    // d2t h_ij
#define D2T_h(I,J) e4phi*(std::pow(4.0*dt_phi[idx],2)*bd.gamma##I##J \
        + 8.0*dt_phi[idx]*dt_g##I##J[idx] + 4.0*d2t_phi[idx]*bd.gamma##I##J \
        + d2t_g##I##J[idx])
    d2t_h11[idx] = D2T_h(1,1) - 2.0*(dadt*dadt + a*d2adt2);
    d2t_h22[idx] = D2T_h(2,2) - 2.0*(dadt*dadt + a*d2adt2);
    d2t_h33[idx] = D2T_h(3,3) - 2.0*(dadt*dadt + a*d2adt2);
    d2t_h12[idx] = D2T_h(1,2);
    d2t_h13[idx] = D2T_h(1,3);
    d2t_h23[idx] = D2T_h(2,3);

    // h0i
#define h0I(I) e4phi*(bd.gamma##I##1*bd.beta1 + bd.gamma##I##2*bd.beta2 + bd.gamma##I##3*bd.beta3)
    h01[idx] = h0I(1);
    h02[idx] = h0I(2);
    h03[idx] = h0I(3);

    // d/dt h0i
#define DT_h0I(I) 4.0*dt_phi[idx]*(h0I(I)) + e4phi*( \
        dt_g##I##1[idx]*bd.beta1 + dt_g##I##2[idx]*bd.beta2 + dt_g##I##3[idx]*bd.beta3 \
        + bd.gamma##I##1*dt_beta1[idx] + bd.gamma##I##2*dt_beta2[idx] + bd.gamma##I##3*dt_beta3[idx])
    h01[idx] = DT_h0I(1);
    h02[idx] = DT_h0I(2);
    h03[idx] = DT_h0I(3);

    // "E" SVT scalar
    E[idx] = bd.alpha*bd.alpha - 1.0
      + bd.gamma11*bd.beta1*bd.beta1 + bd.gamma22*bd.beta2*bd.beta2 + bd.gamma33*bd.beta3*bd.beta3
      + 2.0*(bd.gamma12*bd.beta1*bd.beta2 + bd.gamma13*bd.beta1*bd.beta3 + bd.gamma23*bd.beta2*bd.beta3 );
  }


  // construct A (and its time derivatives) in increments:
  // (A.1) construct d_i d_j h_{ij}
#pragma omp parallel for default(shared) private(i, j, k)
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
#pragma omp parallel for default(shared) private(i, j, k)
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
    dt_A[idx] = -2.0*A[idx]*H + ( dt_h_tr - dt_djdk_d2_hjk )/a/a/2.0;
    d2t_A[idx] =  -2.0*A[idx]*( pw2(H) + d2adt2/a ) - 4.0*H*dt_A[idx]
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
#pragma omp parallel for default(shared) private(i, j, k)
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


  // Construct F, dt_F
#pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    // d^2 F
    F[idx] = ( derivative(i,j,k,1,h01) + derivative(i,j,k,2,h02) + derivative(i,j,k,3,h03) ) / a;
    // d^2 F'
    dt_F[idx] = ( derivative(i,j,k,1,dt_h01) + derivative(i,j,k,2,dt_h02) + derivative(i,j,k,3,dt_h03) ) / a
      - H*F[idx];
  }
  fourier->inverseLaplacian <idx_t, real_t> (F._array);
  fourier->inverseLaplacian <idx_t, real_t> (dt_F._array);


  // All scalar metric fields have been obtained, so compute Bardeen potentials:
#pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    Psi[idx] = -1.0/2.0*( A[idx] + 2.0*H*(a*F[idx] - 0.5*a*a*dt_B[idx]) );
    Phi[idx] = 1.0/2.0*( E[idx]
      + 2.0*(dadt*F[idx] + a*dt_F[idx])
      - 2.0*a*dadt*dt_B[idx] - a*a*d2t_B[idx]
    );
  }


  // Construct vector potentials.
#pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    // G components
    G1[idx] = derivative(i,j,k,1,F) - h01[idx]/a;
    G2[idx] = derivative(i,j,k,2,F) - h02[idx]/a;
    G3[idx] = derivative(i,j,k,3,F) - h03[idx]/a;

    // d^2 C_i
    C1[idx] = ( derivative(i,j,k,2,h12) + derivative(i,j,k,3,h13)
      - derivative(i,j,k,1,h22) - derivative(i,j,k,1,h33) )/a/a 
      + 2.0*derivative(i,j,k,1,A);
    C2[idx] = ( derivative(i,j,k,1,h12) + derivative(i,j,k,3,h23)
      - derivative(i,j,k,2,h11) - derivative(i,j,k,2,h33) )/a/a 
      + 2.0*derivative(i,j,k,2,A);
    C3[idx] = ( derivative(i,j,k,1,h13) + derivative(i,j,k,2,h23)
      - derivative(i,j,k,3,h11) - derivative(i,j,k,3,h22) )/a/a 
      + 2.0*derivative(i,j,k,3,A);
  }
  fourier->inverseLaplacian <idx_t, real_t> (C1._array);
  fourier->inverseLaplacian <idx_t, real_t> (C2._array);
  fourier->inverseLaplacian <idx_t, real_t> (C3._array);


  // Construct tensor potentials.
#pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

#define DIJ(I,J) h##I##J[idx]/a/a - (I==J?1.0:0.0)*A[idx] - double_derivative(i,j,k,I,J,B) \
          - derivative(i,j,k,I,C##J) - derivative(i,j,k,J,C##I)
    D11[idx] = DIJ(1,1);
    D12[idx] = DIJ(1,2);
    D13[idx] = DIJ(1,3);
    D22[idx] = DIJ(2,2);
    D23[idx] = DIJ(2,3);
    D33[idx] = DIJ(3,3);
  }

  // linearized constraint violation
#pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    
    lin_viol[idx] = E[idx] + A[idx] - a*a*d2t_B[idx] - 3.0*a*dadt*dt_B[idx]
      + 2.0*a*dt_F[idx] + 4.0*dadt*F[idx];

    lin_viol_der_mag[idx] = abs_dder(i,j,k,E) + abs_dder(i,j,k,A)
      + a*a*abs_dder(i,j,k,d2t_B) + std::abs(3.0*a*dadt)*abs_dder(i,j,k,dt_B)
      + 2.0*a*abs_dder(i,j,k,dt_F) + std::abs(4.0*dadt)*abs_dder(i,j,k,F);
  }
#pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    lin_viol_der[idx] = abs_dder(i,j,k,lin_viol);
  }

  real_t mean_viol = average(lin_viol_der);
  real_t viol_scale = average(lin_viol_der_mag);
  real_t max_viol = max(lin_viol_der);
  real_t std_viol = standard_deviation(lin_viol_der, mean_viol);

  viols[0] = mean_viol / viol_scale;
  viols[1] = std_viol / viol_scale;
  viols[2] = max_viol / viol_scale;
  viols[3] = max_viol;
  viols[4] = a / exp( 2.0*phi_avg );
  viols[5] = dadt / (-1.0/3.0*a*alpha_avg*K_avg);
  viols[6] = d2adt2 / ( H*dadt + 2.0*a*conformal_average(d2t_phi, DIFFphi_a, 0.0) );

  // TODO: Does ( G - a*dt_C ) ~ 1/a^2 ? (vector modes)
}


void Bardeen::getSVTViolations(real_t * viols_copyto)
{
  for(int i=0; i<7; ++i)
    viols_copyto[i] = viols[i];
}


} // namespace
