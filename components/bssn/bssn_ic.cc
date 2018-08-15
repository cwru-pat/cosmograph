#include "bssn_ic.h"
#include "../../cosmo_types.h"
#include "../../cosmo_globals.h"

namespace cosmo
{

/**
 * @brief      Set flat vacuum spacetime initial conditions,
 * plus small "noise" fluctuations.
 *
 * @param      bssn  BSSN instance
 */
void bssn_ic_awa_stability(BSSN * bssn, real_t A)
{
  idx_t i, j, k;

  std::random_device rd;
  std::mt19937 gen(7.0 /*rd()*/);
  std::uniform_real_distribution<real_t> dist(
      -A*50/NX*50/NX, A*50/NX*50/NX
    );

  arr_t & DIFFgamma11_p = *bssn->fields["DIFFgamma11_p"];
  arr_t & DIFFgamma12_p = *bssn->fields["DIFFgamma12_p"];
  arr_t & DIFFgamma13_p = *bssn->fields["DIFFgamma13_p"];
  arr_t & DIFFgamma22_p = *bssn->fields["DIFFgamma22_p"];
  arr_t & DIFFgamma23_p = *bssn->fields["DIFFgamma23_p"];
  arr_t & DIFFgamma33_p = *bssn->fields["DIFFgamma33_p"];
  arr_t & DIFFphi_p = *bssn->fields["DIFFphi_p"];
  arr_t & A11_p = *bssn->fields["A11_p"];
  arr_t & A12_p = *bssn->fields["A12_p"];
  arr_t & A13_p = *bssn->fields["A13_p"];
  arr_t & A22_p = *bssn->fields["A22_p"];
  arr_t & A23_p = *bssn->fields["A23_p"];
  arr_t & A33_p = *bssn->fields["A33_p"];
  arr_t & DIFFK_p = *bssn->fields["DIFFK_p"];
  arr_t & Gamma1_p = *bssn->fields["Gamma1_p"];
  arr_t & Gamma2_p = *bssn->fields["Gamma2_p"];
  arr_t & Gamma3_p = *bssn->fields["Gamma3_p"];

  LOOP3(i, j, k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    DIFFgamma11_p[idx] = dist(gen);
    DIFFgamma12_p[idx] = dist(gen);
    DIFFgamma13_p[idx] = dist(gen);
    DIFFgamma22_p[idx] = dist(gen);
    DIFFgamma23_p[idx] = dist(gen);
    DIFFgamma33_p[idx] = dist(gen);
    DIFFphi_p[idx] = dist(gen);
    A11_p[idx]     = dist(gen);
    A12_p[idx]     = dist(gen);
    A13_p[idx]     = dist(gen);
    A22_p[idx]     = dist(gen);
    A23_p[idx]     = dist(gen);
    A33_p[idx]     = dist(gen);
    DIFFK_p[idx]   = dist(gen);
    Gamma1_p[idx]  = dist(gen);
    Gamma2_p[idx]  = dist(gen);
    Gamma3_p[idx]  = dist(gen);
  }
}

/**
 * @brief Initialize with a "linear" wave propagating in the x-direction.
 * @details AwA test.
 */
void bssn_ic_awa_linear_wave(BSSN * bssn)
{
  bssn_ic_awa_linear_wave(bssn, 1.0e-8, 1);
}

/**
 * @brief Initialize with a "linear" wave propagating in the 'dir'-direction.
 * @details AwA test.
 */
void bssn_ic_awa_linear_wave(BSSN * bssn, real_t A, int dir)
{
  idx_t i, j, k;

  arr_t & DIFFgamma11_p = *bssn->fields["DIFFgamma11_p"];
  arr_t & DIFFgamma22_p = *bssn->fields["DIFFgamma22_p"];
  arr_t & DIFFgamma33_p = *bssn->fields["DIFFgamma33_p"];

  arr_t & A11_p = *bssn->fields["A11_p"];
  arr_t & A22_p = *bssn->fields["A22_p"];
  arr_t & A33_p = *bssn->fields["A33_p"];

  LOOP3(i,j,k)
  {
    // position of dir
    real_t w = 0;
    switch(dir)
    {
      case 1 :
        w = ((real_t) i)*dx;
        DIFFgamma22_p[NP_INDEX(i,j,k)] = A*sin( 2.0*PI*w );
        DIFFgamma33_p[NP_INDEX(i,j,k)] = -A*sin( 2.0*PI*w );
        A22_p[NP_INDEX(i,j,k)] = PI*A*cos( 2.0*PI*w );
        A33_p[NP_INDEX(i,j,k)] = -PI*A*cos( 2.0*PI*w );
        break;
      case 2 :
        w = ((real_t) j)*dx;
        DIFFgamma11_p[NP_INDEX(i,j,k)] = A*sin( 2.0*PI*w );
        DIFFgamma33_p[NP_INDEX(i,j,k)] = -A*sin( 2.0*PI*w );
        A11_p[NP_INDEX(i,j,k)] = PI*A*cos( 2.0*PI*w );
        A33_p[NP_INDEX(i,j,k)] = -PI*A*cos( 2.0*PI*w );
        break;
      case 3 :
        w = ((real_t) k)*dx;
        DIFFgamma11_p[NP_INDEX(i,j,k)] = A*sin( 2.0*PI*w );
        DIFFgamma22_p[NP_INDEX(i,j,k)] = -A*sin( 2.0*PI*w );
        A11_p[NP_INDEX(i,j,k)] = PI*A*cos( 2.0*PI*w );
        A22_p[NP_INDEX(i,j,k)] = -PI*A*cos( 2.0*PI*w );
        break;
    }
  }
}

/**
 * @brief Initialize with a "linear" wave propagating in the x-direction.
 * @details AwA test.
 */
void bssn_ic_awa_diagonal_linear_wave(BSSN * bssn)
{
  // TODO
}

/**
 * @brief Initialize with a "linear" wave propagating in the x-direction,
 * on top of a deSitter spacetime.
 * @details Solution should behave, at linear order, as a damped wave equation.
 * Should be run as a vacuum sim (todo: clean up)
 */
void bssn_ic_awa_linear_wave_desitter(BSSN * bssn)
{
  bssn_ic_awa_linear_wave(bssn, 1.0e-8, 1);

  real_t rho = 1.0e-3;
  auto & frw = bssn->frw;

  real_t K = -sqrt(24.0*PI*rho);

  frw->set_phi(0.0);
  frw->set_K(K);
  frw->addFluid(rho, -1.0 /* w=-1 */);

  // idx_t i, j, k;

  // arr_t & DIFFr_a = *bssn->fields["DIFFr_a"];
  // arr_t & DIFFS_a = *bssn->fields["DIFFS_a"];
  // arr_t & DIFFK_p = *bssn->fields["DIFFK_p"];

  // LOOP3(i,j,k)
  // {
  //   // FRW Background parameters
  //   real_t rho = 1.0e-3;

  //   DIFFr_a[NP_INDEX(i,j,k)] = rho; // constant density
  //                                   // not evolved for a vacuum sim
  //   DIFFS_a[NP_INDEX(i,j,k)] = -3.0*rho;
  //   DIFFK_p[NP_INDEX(i,j,k)] = -sqrt(24.0*PI*rho);
  // }
}

/**
 * @brief      AwA Gauge Wave Test
 */
void bssn_ic_awa_gauge_wave(BSSN * bssn)
{
  bssn_ic_awa_gauge_wave(bssn, 1);
}

/**
 * @brief      AwA Gauge Wave Test, setting wave propagation direction
 */
void bssn_ic_awa_gauge_wave(BSSN * bssn, int dir)
{
  idx_t i, j, k;

  real_t A = 0.5;

  arr_t & DIFFphi_p = *bssn->fields["DIFFphi_p"];
  arr_t & DIFFK_p = *bssn->fields["DIFFK_p"];
  arr_t & DIFFalpha_p = *bssn->fields["DIFFalpha_p"];

  arr_t & DIFFgamma11_p = *bssn->fields["DIFFgamma11_p"];
  arr_t & DIFFgamma22_p = *bssn->fields["DIFFgamma22_p"];
  arr_t & DIFFgamma33_p = *bssn->fields["DIFFgamma33_p"];
  
  arr_t & A11_p = *bssn->fields["A11_p"];
  arr_t & A22_p = *bssn->fields["A22_p"];
  arr_t & A33_p = *bssn->fields["A33_p"];

  arr_t & Gamma1_p = *bssn->fields["Gamma1_p"];
  arr_t & Gamma2_p = *bssn->fields["Gamma2_p"];
  arr_t & Gamma3_p = *bssn->fields["Gamma3_p"];

  LOOP3(i,j,k)
  {
    // position of dir
    real_t w = 0;
    switch(dir)
    {
      case 1 :
        w = ((real_t) i)*dx;
        break;
      case 2 :
        w = ((real_t) j)*dx;
        break;
      case 3 :
        w = ((real_t) k)*dx;
        break;
    }

    // H(t = 0)
    real_t H = A*sin( 2.0*PI*w );
    real_t dtH = -2.0*PI*A*cos( 2.0*PI*w );
    real_t Kxx = dtH/(2.0*sqrt(1.0-H));
    real_t K = Kxx / (1.0 - H);

    DIFFphi_p[NP_INDEX(i,j,k)] = log1p(-H)/12.0;
    DIFFalpha_p[NP_INDEX(i,j,k)] = sqrt(1.0-H) - 1.0;
    DIFFK_p[NP_INDEX(i,j,k)] = K;

    DIFFgamma11_p[NP_INDEX(i,j,k)] = pow(1.0-H, -1.0/3.0) - 1.0;
    DIFFgamma22_p[NP_INDEX(i,j,k)] = pow(1.0-H, -1.0/3.0) - 1.0;
    DIFFgamma33_p[NP_INDEX(i,j,k)] = pow(1.0-H, -1.0/3.0) - 1.0;

    A11_p[NP_INDEX(i,j,k)] = - K/3.0*pow(1.0-H, -1.0/3.0);
    A22_p[NP_INDEX(i,j,k)] = - K/3.0*pow(1.0-H, -1.0/3.0);
    A33_p[NP_INDEX(i,j,k)] = - K/3.0*pow(1.0-H, -1.0/3.0);

    switch(dir)
    {
      case 1 :
        DIFFgamma11_p[NP_INDEX(i,j,k)] = pow(1.0-H, 2.0/3.0) - 1.0;
        A11_p[NP_INDEX(i,j,k)] = Kxx*pow(1.0-H, -1.0/3.0) - K/3.0*pow(1.0-H, 2.0/3.0);
        Gamma1_p[NP_INDEX(i,j,k)] = 2.0/3.0*dtH*pow(1.0-H, -5.0/3.0);
        break;
      case 2 :
        DIFFgamma22_p[NP_INDEX(i,j,k)] = pow(1.0-H, 2.0/3.0) - 1.0;
        A22_p[NP_INDEX(i,j,k)] = Kxx*pow(1.0-H, -1.0/3.0) - K/3.0*pow(1.0-H, 2.0/3.0);
        Gamma2_p[NP_INDEX(i,j,k)] = 2.0/3.0*dtH*pow(1.0-H, -5.0/3.0);
        break;
      case 3 :
        DIFFgamma33_p[NP_INDEX(i,j,k)] = pow(1.0-H, 2.0/3.0) - 1.0;
        A33_p[NP_INDEX(i,j,k)] = Kxx*pow(1.0-H, -1.0/3.0) - K/3.0*pow(1.0-H, 2.0/3.0);
        Gamma3_p[NP_INDEX(i,j,k)] = 2.0/3.0*dtH*pow(1.0-H, -5.0/3.0);
        break;
      default:
        throw -1;
        break;
    }

  }
}

void bssn_ic_awa_shifted_gauge_wave(BSSN * bssn)
{
  bssn_ic_awa_shifted_gauge_wave(bssn, 1);
}

/**
 * @brief      AwA Shifted Gauge Wave Test
 */
void bssn_ic_awa_shifted_gauge_wave(BSSN * bssn, int dir)
{
# if ! USE_BSSN_SHIFT
  std::cerr << "USE_BSSN_SHIFT must be enabled for shifted gauge wave! (cmake using -DCOSMO_USE_BSSN_SHIFT=1)" << std::endl;
  throw -1;
# endif

  idx_t i, j, k;

  real_t A = 0.5;

  arr_t & DIFFphi_p = *bssn->fields["DIFFphi_p"];
  arr_t & DIFFK_p = *bssn->fields["DIFFK_p"];
  arr_t & DIFFalpha_p = *bssn->fields["DIFFalpha_p"];

  arr_t & DIFFgamma11_p = *bssn->fields["DIFFgamma11_p"];
  arr_t & DIFFgamma22_p = *bssn->fields["DIFFgamma22_p"];
  arr_t & DIFFgamma33_p = *bssn->fields["DIFFgamma33_p"];
  
  arr_t & A11_p = *bssn->fields["A11_p"];
  arr_t & A22_p = *bssn->fields["A22_p"];
  arr_t & A33_p = *bssn->fields["A33_p"];

  arr_t & Gamma1_p = *bssn->fields["Gamma1_p"];
  arr_t & Gamma2_p = *bssn->fields["Gamma2_p"];
  arr_t & Gamma3_p = *bssn->fields["Gamma3_p"];

  arr_t & beta1_p = *bssn->fields["beta1_p"];
  arr_t & beta2_p = *bssn->fields["beta2_p"];
  arr_t & beta3_p = *bssn->fields["beta3_p"];

  LOOP3(i,j,k)
  {
    // position of dir
    real_t w = 0;
    switch(dir)
    {
      case 1 :
        w = ((real_t) i)*dx;
        break;
      case 2 :
        w = ((real_t) j)*dx;
        break;
      case 3 :
        w = ((real_t) k)*dx;
        break;
    }

    // H(t = 0)
    real_t H = A*sin( 2.0*PI*w );
    real_t dtH = -2.0*PI*A*cos( 2.0*PI*w );
    real_t Kxx = dtH/(2.0*sqrt(1.0+H));
    real_t K = Kxx / (1.0 + H);

    DIFFphi_p[NP_INDEX(i,j,k)] = log1p(H)/12.0;
    DIFFalpha_p[NP_INDEX(i,j,k)] = sqrt(1.0/(1.0+H)) - 1.0;
    DIFFK_p[NP_INDEX(i,j,k)] = K;

    DIFFgamma11_p[NP_INDEX(i,j,k)] = pow(1.0+H, -1.0/3.0) - 1.0;
    DIFFgamma22_p[NP_INDEX(i,j,k)] = pow(1.0+H, -1.0/3.0) - 1.0;
    DIFFgamma33_p[NP_INDEX(i,j,k)] = pow(1.0+H, -1.0/3.0) - 1.0;

    A11_p[NP_INDEX(i,j,k)] = - K/3.0*pow(1.0+H, -1.0/3.0);
    A22_p[NP_INDEX(i,j,k)] = - K/3.0*pow(1.0+H, -1.0/3.0);
    A33_p[NP_INDEX(i,j,k)] = - K/3.0*pow(1.0+H, -1.0/3.0);

    switch(dir)
    {
      case 1 :
        Gamma1_p[NP_INDEX(i,j,k)] = -2.0/3.0*dtH*pow(1.0+H, -5.0/3.0);
        beta1_p[NP_INDEX(i,j,k)] = -H/(1.0+H);
        DIFFgamma11_p[NP_INDEX(i,j,k)] = pow(1.0+H, 2.0/3.0) - 1.0;
        A11_p[NP_INDEX(i,j,k)] = Kxx*pow(1.0+H, -1.0/3.0) - K/3.0*pow(1.0+H, 2.0/3.0);
        break;
      case 2 :
        Gamma2_p[NP_INDEX(i,j,k)] = -2.0/3.0*dtH*pow(1.0+H, -5.0/3.0);
        beta2_p[NP_INDEX(i,j,k)] = -H/(1.0+H);
        DIFFgamma22_p[NP_INDEX(i,j,k)] = pow(1.0+H, 2.0/3.0) - 1.0;
        A22_p[NP_INDEX(i,j,k)] = Kxx*pow(1.0+H, -1.0/3.0) - K/3.0*pow(1.0+H, 2.0/3.0);
        break;
      case 3 :
        Gamma3_p[NP_INDEX(i,j,k)] = -2.0/3.0*dtH*pow(1.0+H, -5.0/3.0);
        beta3_p[NP_INDEX(i,j,k)] = -H/(1.0+H);
        DIFFgamma33_p[NP_INDEX(i,j,k)] = pow(1.0+H, 2.0/3.0) - 1.0;
        A33_p[NP_INDEX(i,j,k)] = Kxx*pow(1.0+H, -1.0/3.0) - K/3.0*pow(1.0+H, 2.0/3.0);
        break;
      default:
        throw -1;
        break;
    }
  }
}

/**
 * @brief      AwA Shifted Gauge Wave Test
 */
void bssn_ic_kasner(BSSN * bssn, real_t px)
{
  idx_t i, j, k;

  // gamma_ij = diag( t^(2px), t^(2py), t^(2pz) )
  // pick px,
  // px + py + pz = 1
  // px^2 + py^2 + pz^2 = 1
  // => py = ( 1 - px - sqrt( 1 + 2*px - 3px^2 ) )/2
  //    pz = ( 1 - px + sqrt( 1 + 2*px - 3px^2 ) )/2
  // det(gamma) = t^2
  // Let t_i = 1 => DIFFphi = 0, DIFFgamma_ij
  // 
  // K = - 1/(2 det(gamma)) d/dt det(gamma)
  //   = -1/(2t) * 2
  // => K = -1
  // 
  // K_ij = -1/2 * d_t gamma_ij
  //      = -1/2 * diag( 2px*t^(2px-1), 2py*t^(2py-1), 2pz*t^(2pz-1) )
  // 
  // Check: K = gamma^ij K_ij = -1 
  // 
  // A_ij = K_ij - K gamma_ij / 3
  //      = -1/2 * diag( 2px*t^(2px-1), 2py*t^(2py-1), 2pz*t^(2pz-1) ) - (-1/3)*diag( t^(2px), t^(2py), t^(2pz) )
  //      = -diag( px, py, pz ) + diag( 1, 1, 1 )/3
  //      = diag( 1/3-px, 1/3-py, 1/3-pz )

  real_t py = ( 1.0 - px - std::sqrt( 1.0 + 2.0*px - 3.0*px*px ) )/2.0;
  real_t pz = ( 1.0 - px + std::sqrt( 1.0 + 2.0*px - 3.0*px*px ) )/2.0;

  arr_t & DIFFK_p = *bssn->fields["DIFFK_p"];  
  arr_t & A11_p = *bssn->fields["A11_p"];
  arr_t & A22_p = *bssn->fields["A22_p"];
  arr_t & A33_p = *bssn->fields["A33_p"];

  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    DIFFK_p[idx] = -1.0;
    A11_p[idx] = 1.0/3.0 - px;
    A22_p[idx] = 1.0/3.0 - py;
    A33_p[idx] = 1.0/3.0 - pz;
  }
}


} // namespace cosmo
