#include "scalar_ic.h"
#include "../../cosmo_includes.h"
#include "../../cosmo_types.h"
#include "../../cosmo_globals.h"
#include "../../utils/Fourier.h"
#include "../../utils/math.h"
#if USE_MULTIGRID
#include "../elliptic_solver/multigrid.h"
#endif

namespace cosmo
{

/**
 * @brief Initialize a small-amplitude wave on a flat metric
 *  (use a stationary gaussian wave packet for now)
 *  NB: nonzero constraint violation
 */
void scalar_ic_set_wave(BSSN * bssn, Scalar * scalar)
{
  // BSSN is already initialized to flat, just initialize scalar fields
  arr_t & phi = scalar->phi._array_p; // gaussian 
  arr_t & psi1 = scalar->psi1._array_p; // derivative of phi in x-dir
  arr_t & psi2 = scalar->psi3._array_p; // derivative of phi in y-dir
  arr_t & psi3 = scalar->psi2._array_p; // derivative of phi in z-dir

  arr_t & K_p = *bssn->fields["DIFFK_p"]; // extrinsic curvature
  arr_t & K_a = *bssn->fields["DIFFK_a"]; // extrinsic curvature

  // iterators
  idx_t i, j, k;
  // gaussian parameters
  real_t amplitude = 1.0e-10;
  real_t sigx = (real_t) NX / 10.0;

  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    phi[INDEX(i,j,k)] = amplitude*exp(
        -pw2(i - ((real_t) NX-1)/2.0)/pw2(sigx)
      );
  }

  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    psi1[INDEX(i,j,k)] = derivative(i, j, k, 1, phi);
    psi2[INDEX(i,j,k)] = derivative(i, j, k, 2, phi);
    psi3[INDEX(i,j,k)] = derivative(i, j, k, 3, phi);
  }

  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    K_a[INDEX(i,j,k)] = K_p[INDEX(i,j,k)] = -std::sqrt(24*PI*scalar->V(phi[INDEX(i,j,k)]));
  }

  return;
}

/**
 * @brief Initialize scalar field to zero: assumes potential with a minimum,
 * constant value here, so a cosmological constant.
 */
void scalar_ic_set_Lambda(BSSN * bssn, Scalar * scalar)
{
  idx_t i, j, k;

  arr_t & K_p = *bssn->fields["DIFFK_p"]; // extrinsic curvature
  arr_t & K_a = *bssn->fields["DIFFK_a"]; // extrinsic curvature

  real_t Lambda = scalar->V(0);
  LOOP3(i,j,k)
  {
    K_a[INDEX(i,j,k)] = K_p[INDEX(i,j,k)] = -std::sqrt(24 * PI * Lambda);
  }
}

/**
 * @brief semi-analytic test initial conditions: a sinusoidal fluctuation in a
 * particular direction & corresponding metric
 */
void scalar_ic_set_semianalytic_test(BSSN * bssn, Scalar * scalar)
{
  idx_t i, j, k;

  arr_t & phi_p = *bssn->fields["DIFFphi_p"];
  arr_t & phi_a = *bssn->fields["DIFFphi_a"];

  //chosing the analytic solution first
  #pragma omp parallel for
  LOOP3(i, j, k)
  {
    phi_p[INDEX(i,j,k)] = 1.0 + 0.01 * std::sin(4.0 * PI *( (real_t)i / NX - 0.125));
    phi_a[INDEX(i,j,k)] = 1.0 + 0.01 * std::sin(4.0 * PI *( (real_t)i / NX - 0.125));
  }

  arr_t & phi = scalar->phi._array_p; // field
  arr_t & psi1 = scalar->psi1._array_p; // derivative of phi in x-dir
  arr_t & psi2 = scalar->psi3._array_p; // derivative of phi in y-dir
  arr_t & psi3 = scalar->psi2._array_p; // derivative of phi in z-dir
  
  arr_t & K_p = *bssn->fields["DIFFK_p"]; // extrinsic curvature
  arr_t & K_a = *bssn->fields["DIFFK_a"]; // extrinsic curvature

  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    K_a[INDEX(i,j,k)] = K_p[INDEX(i,j,k)] = 9.763423957014197;
  }
  real_t Lambda = 1.0;
  
  real_t * temp, * der_bak;
  temp = new real_t[NX], der_bak = new real_t[NX];
  real_t lap_dif = 0;
  for(i = 0; i < NX; i++)
  {
    if(2.0 * i <= NX)
      temp[i] = std::sqrt(std::fabs(
      (pw2(4.0 * PI) * 0.01  *  std::sin(4.0 * PI *( (real_t)i / NX - 0.125) )
       +(-2.0 * PI * Lambda + pw2(K_a[INDEX(i,0,0)])/12.0 )*
       std::pow(phi_a[INDEX(i,0,0)], 5.0) )
      /(PI * phi_a[INDEX(i,0,0)] ) )) ;
    else
      temp[i]= -  std::sqrt(std::fabs(
      (pw2(4.0 * PI) * 0.01  *  std::sin(4.0 * PI *( (real_t)i / NX - 0.125) )
       +(-2.0 * PI * Lambda + pw2(K_a[INDEX(i,0,0)])/12.0 )*
       std::pow(phi_a[INDEX(i,0,0)], 5.0) )
      /(PI * phi_a[INDEX(i,0,0)] ) )) ;
    der_bak[i] = temp[i];
    for(j = 0; j < NY; j++)
      for(k = 0; k < NZ; k++)
  phi[INDEX(i,j,k)] = temp[i];
    lap_dif = std::max(lap_dif, std::fabs(double_derivative(i,0,0,1,1,phi_p)
            +pw2(4.0 * PI) * 0.01  *  std::sin(4.0 * PI *( (real_t)i / NX - 0.125) ) ) );
  }

  std::cout << "Difference between exact lap and discrete lap is: "<<lap_dif<<"\n";
  Fourier * fourier;
  fourier = new Fourier();

  fourier->Initialize_1D(NX, temp);
  fourier->execute_f_r2c(0);

  for(i = 0; i < NX/2 +1; i++)
  {
    if(i > 0)
    {
      std::swap(fourier->f_field[i][0], fourier->f_field[i][1]);
      fourier->f_field[i][0] = dx * fourier->f_field[i][0] /
          ( 2.0 * PI * (real_t) i / NX);
      fourier->f_field[i][1] = - dx * fourier->f_field[i][1] /
          ( 2.0 * PI * (real_t) i / NX);
    }
    else
    {
      fourier->f_field[i][0] = fourier->f_field[i][1] = 0;
    }
  }
  fourier->execute_f_c2r(0);

  LOOP3(i,j,k)
  {
    phi[INDEX(i,j,k)] = temp[i]/NX;
    phi_p[INDEX(i,j,k)] =(real_t) std::log(phi_p[INDEX(i,j,k)]);
    phi_a[INDEX(i,j,k)] =(real_t) std::log(phi_a[INDEX(i,j,k)]);    
  }

  real_t max_deviation =- 1e100;

  for(i = 0; i < NX; i++)
  {
    max_deviation = std::max(max_deviation,std::fabs( derivative(i,0,0,1,phi) - der_bak[i]));
  }
  std::cout<<"The maximum deviation for solving scalar field under analytic solution of phi under odx8 is: "<<max_deviation;
  std::cout<<"\n";

  // initialize psi according to values in phi
  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    psi1[INDEX(i,j,k)] = derivative(i, j, k, 1, phi);
    psi2[INDEX(i,j,k)] = derivative(i, j, k, 2, phi);
    psi3[INDEX(i,j,k)] = derivative(i, j, k, 3, phi);
  }
}

#if USE_MULTIGRID
/**
 * @brief Use the multigrid solver to solve for metric factors given
 * a particular scalar field implementation.
 */
void scalar_ic_set_multigrid(BSSN * bssn, Scalar * scalar)
{
  idx_t i, j, k;

  // Choose a configuration for the scalar fields first:
  arr_t & phi = scalar->phi._array_p; // field
  arr_t & psi1 = scalar->psi1._array_p; // derivative of phi in x-dir
  arr_t & psi2 = scalar->psi2._array_p; // derivative of phi in y-dir
  arr_t & psi3 = scalar->psi3._array_p; // derivative of phi in z-dir

  std::random_device rd;
  std::mt19937 gen(7.0 /*rd()*/);
  std::uniform_real_distribution<real_t> dist(0, 2.0*PI);

  // cutoff @ "ic_spec_cut"; maybe initialize this field
  // according to some power spectrum?
  real_t n_max = std::stoi(_config["n_max"]);
  real_t phi_0 = std::stod(_config["phi_0"]);
  real_t delta = std::stod(_config["delta_phi"]);
  
  // background value
  LOOP3(i,j,k)
    phi[INDEX(i,j,k)] = phi_0;

  // sum over different modes
  for(int n = -n_max; n <= n_max; ++n)
  {
    if(n != 0)
    {
      // random phases
      real_t x_phase = dist(gen),
             y_phase = dist(gen),
             z_phase = dist(gen);

      #pragma omp parallel for default(shared) private(i,j,k)
      LOOP3(i,j,k)
      {
        // some sinusoidal modes
        phi[INDEX(i,j,k)] += delta*(
                cos(2.0*PI*((real_t) n/NX)*i + x_phase )
                 + cos(2.0*PI*((real_t) n/NY)*j + y_phase )
                 + cos(2.0*PI*((real_t) n/NZ)*k + z_phase )
            );
      }
    }
  }

  // initialize psi according to values in phi
  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    psi1[INDEX(i,j,k)] = derivative(i, j, k, 1, phi);
    psi2[INDEX(i,j,k)] = derivative(i, j, k, 2, phi);
    psi3[INDEX(i,j,k)] = derivative(i, j, k, 3, phi);
  }

  // PI is zero for now

  // compute background/average K
  real_t K_src = 0;
  LOOP3(i, j, k)
  {
    K_src += 3.0*scalar->V(phi[INDEX(i,j,k)])
      +3.0/2.0*(
        pw2(psi1[INDEX(i,j,k)]) + pw2(psi2[INDEX(i,j,k)]) + pw2(psi3[INDEX(i,j,k)])
      );
  }
  K_src = -std::sqrt(K_src/NX/NY/NZ);

  arr_t & K_p = *bssn->fields["DIFFK_p"]; // extrinsic curvature
  arr_t & K_a = *bssn->fields["DIFFK_a"]; // extrinsic curvature
  arr_t & K0 = *bssn->fields["K0_a"]; // initial extrinsic curvature

  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    idx_t idx = INDEX(i,j,k);
    K0[idx] = K_a[idx] = K_p[idx] = K_src;
  }

  // solve for BSSN fields using multigrid class:
  real_t relaxation_tolerance = std::stod(_config["relaxation_tolerance"]);
  FASMultigrid multigrid (N, N*dx, 4, relaxation_tolerance);

  idx_t u_exp[2] = { 1, 5 };
  multigrid.build_rho(2, u_exp);
  LOOP3(i, j, k)
  {
    real_t value = PI*(
        pw2(psi1[INDEX(i,j,k)]) + pw2(psi2[INDEX(i,j,k)]) + pw2(psi3[INDEX(i,j,k)])
      );
    multigrid.setPolySrcAtPt(i, j, k, 0, value);

    value = 2.0* PI*scalar->V(phi[INDEX(i,j,k)]) - K_a[INDEX(i,j,k)]*K_a[INDEX(i,j,k)]/12.0;
    multigrid.setPolySrcAtPt(i, j, k, 1, value);
  }
  multigrid.initializeRhoHeirarchy();

  if(std::stoi(_config["debug_multigrid"]))
  {
    multigrid.printSourceStrip(0, 5);
    multigrid.printSourceStrip(1, 5);
  }

  // set initial guess and solve using 3 V-cycles
  multigrid.setTrialSolution(0);
  multigrid.VCycles(std::stoi(_config["num_v_cycles"]));

  if(std::stoi(_config["debug_multigrid"]))
  {
    multigrid.printSolutionStrip(1);
    multigrid.printSolutionStrip(5);
  }

  // set bssn phi using multigrid solution
  arr_t & phi_p = *bssn->fields["DIFFphi_p"];
  arr_t & phi_a = *bssn->fields["DIFFphi_a"];

  REAL_T * u = multigrid.getSolution();
  #pragma omp parallel for
  LOOP3(i, j, k)
  {
    phi_p[INDEX(i,j,k)] = (real_t) std::log(std::abs(u[INDEX(i,j,k)]));
    phi_a[INDEX(i,j,k)] = (real_t) std::log(std::abs(u[INDEX(i,j,k)]));
  }

  return;
}
#endif // if USE_MULTIGRID

} // namespace cosmo
