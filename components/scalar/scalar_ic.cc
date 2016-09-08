#include "scalar_ic.h"
#include "../../cosmo_includes.h"
#include "../../cosmo_types.h"
#include "../../cosmo_globals.h"
#include "../../utils/Fourier.h"
#include "../../utils/math.h"
#if USE_MULTIGRID
#include "../elliptic_solver/full_multigrid.h"
#endif

namespace cosmo
{

/**
 * @brief Initialize a small-amplitude wave on a flat metric
 *  (use a stationary Gaussian wave packet for now)
 *  NB: nonzero constraint violation
 */
void scalar_ic_set_wave(BSSN * bssn, Scalar * scalar)
{
  // BSSN is already initialized to flat, just initialize scalar fields
  arr_t & phi = scalar->phi._array_p; // Gaussian 
  arr_t & psi1 = scalar->psi1._array_p; // derivative of phi in x-dir
  arr_t & psi2 = scalar->psi3._array_p; // derivative of phi in y-dir
  arr_t & psi3 = scalar->psi2._array_p; // derivative of phi in z-dir

  arr_t & K_p = *bssn->fields["DIFFK_p"]; // extrinsic curvature
  arr_t & K_a = *bssn->fields["DIFFK_a"]; // extrinsic curvature

  // iterators
  idx_t i, j, k;
  // Gaussian parameters
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
void scalar_ic_set_semianalytic_test(BSSN * bssn, Scalar * scalar,
  IOData * iodata)
{
  idx_t i, j, k;

  arr_t & phi_p = *bssn->fields["DIFFphi_p"];
  arr_t & phi_a = *bssn->fields["DIFFphi_a"];

  // choosing the analytic solution first
  #pragma omp parallel for
  LOOP3(i, j, k)
  {
    phi_p[INDEX(i,j,k)] = 1.0 + 0.01 * std::sin(4.0 * PI *( (real_t)i / NX - 0.125));
    phi_a[INDEX(i,j,k)] = 1.0 + 0.01 * std::sin(4.0 * PI *( (real_t)i / NX - 0.125));
  }

  arr_t & phi = scalar->phi._array_p; // field
  arr_t & psi1 = scalar->psi1._array_p; // derivative of phi in x-dir
  arr_t & psi2 = scalar->psi2._array_p; // derivative of phi in y-dir
  arr_t & psi3 = scalar->psi3._array_p; // derivative of phi in z-dir
  
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

  iodata->log("The difference between exact lap and discrete lap is: " + stringify(lap_dif));
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
    max_deviation = std::max(max_deviation, std::fabs( derivative(i,0,0,1,phi) - der_bak[i]));
  }
  iodata->log("The maximum deviation of the numerical and analytic solution of phi using odx"
    + stringify(STENCIL_ORDER) + " stencils is: " + stringify(max_deviation));

  // initialize psi according to values in phi
  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    psi1[INDEX(i,j,k)] = derivative(i, j, k, 1, phi);
    psi2[INDEX(i,j,k)] = derivative(i, j, k, 2, phi);
    psi3[INDEX(i,j,k)] = derivative(i, j, k, 3, phi);
  }
}

void scalar_ic_set_full_equations(BSSN * bssn, Scalar * scalar, IOData * iodata)
{
  idx_t i, j, k;

    // Choose a configuration for the scalar fields first:
  arr_t & phi = scalar->phi._array_p; // field
  arr_t & psi1 = scalar->psi1._array_p; // derivative of phi in x-dir
  arr_t & psi2 = scalar->psi2._array_p; // derivative of phi in y-dir
  arr_t & psi3 = scalar->psi3._array_p; // derivative of phi in z-dir

  arr_t & Pi = scalar->Pi._array_p;

  real_t dt_phi = 0.1;
  
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
    Pi[INDEX(i,j,k)] = -dt_phi;
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

  arr_t & phi_p = *bssn->fields["DIFFphi_p"];
  arr_t & phi_a = *bssn->fields["DIFFphi_a"];


  arr_t & A11_p = *bssn->fields["A11_p"];
  arr_t & A11_a = *bssn->fields["A11_a"];

  arr_t & A12_p = *bssn->fields["A12_p"];
  arr_t & A12_a = *bssn->fields["A12_a"];

  arr_t & A13_p = *bssn->fields["A13_p"];
  arr_t & A13_a = *bssn->fields["A13_a"];

  arr_t & A22_p = *bssn->fields["A22_p"];
  arr_t & A22_a = *bssn->fields["A22_a"];

  arr_t & A23_p = *bssn->fields["A23_p"];
  arr_t & A23_a = *bssn->fields["A23_a"];

  arr_t & A33_p = *bssn->fields["A33_p"];
  arr_t & A33_a = *bssn->fields["A33_a"];
  
  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    idx_t idx = INDEX(i,j,k);
    K0[idx] = K_a[idx] = K_p[idx] = K_src;
  }


  
  arr_t * X = new arr_t [4];

  for(i = 1; i < 4; i++)
    X[i].init(NX, NY, NZ);

  X[0].nx = phi_p.nx;
  X[0].ny = phi_p.ny;
  X[0].nz = phi_p.nz;
  X[0].pts = phi_p.pts;
  X[0]._array = phi_p._array;
  
  idx_t molecule_n[4] = {18, 5, 5, 5};

  real_t relaxation_tolerance = std::stod(_config["relaxation_tolerance"]);

  atom atom_tmp = {0};

  FASMultigrid multigrid(X, 4, molecule_n, 4, 5, relaxation_tolerance);

  /*Starting adding all the terms in equations*******************************/


  multigrid.eqns[0][0].init(1, 1);
  //adding laplacian term
  atom_tmp.type =  multigrid.atom_type::lap;
  atom_tmp.u_id = 0;
  multigrid.eqns[0][0].add_atom(atom_tmp);

  //adding \Psi^-7 term
  for(i = 1; i <= 3; i++)
    for(j = 1; j <= 3; j++)
    {
      atom_tmp.type = i+1;
      atom_tmp.u_id = j;
      if(i == j)
        multigrid.eqns[0][3*(i-1)+j].init(3, 0.5 - 1.0/6.0);
      else
        multigrid.eqns[0][3*(i-1)+j].init(3, 0.25);
      multigrid.eqns[0][3*(i-1)+j].add_atom(atom_tmp);

      multigrid.eqns[0][3*(i-1)+j].add_atom(atom_tmp);

      atom_tmp.type = 1;
      atom_tmp.u_id = 0;
      atom_tmp.value = -7;
      multigrid.eqns[0][3*(i-1)+j].add_atom(atom_tmp);
    }
  for(i = 1; i <= 3; i++)
    for(j = i+1; j <= 3; j++)
    {
      atom_tmp.type = i+1;
      atom_tmp.u_id = j;
      multigrid.eqns[0][(9 + i+j - 2)].init(3, 0.5);
      multigrid.eqns[0][(9 + i+j - 2)].add_atom(atom_tmp);
      
      atom_tmp.type = j+1;
      atom_tmp.u_id = i;
      multigrid.eqns[0][(9 + i+j - 2)].add_atom(atom_tmp);

      atom_tmp.type = 1;
      atom_tmp.u_id = 0;
      atom_tmp.value = -7;
      multigrid.eqns[0][(9 + i+j - 2)].add_atom(atom_tmp);
    }

  for(i = 1; i <= 3; i++)
    for(j = i+1; j <= 3; j++)
    {
      multigrid.eqns[0][(12 + i+j - 2)].init(3, -1.0/3.0);
      atom_tmp.type = i+1;
      atom_tmp.u_id = i;
      multigrid.eqns[0][(12+i+j-2)].add_atom(atom_tmp);

      atom_tmp.type = j+1;
      atom_tmp.u_id = j;
      multigrid.eqns[0][(12+i+j-2)].add_atom(atom_tmp);

      atom_tmp.type = 1;
      atom_tmp.u_id = 0;
      atom_tmp.value = -7;
      multigrid.eqns[0][(12 + i+j - 2)].add_atom(atom_tmp);
    }

  
  multigrid.eqns[0][16].init(2, 1.0);
  multigrid.eqns[0][17].init(2, 1.0);

  atom_tmp.type = multigrid.atom_type::const_f;
  multigrid.eqns[0][16].add_atom(atom_tmp);
  multigrid.eqns[0][17].add_atom(atom_tmp);
  
  atom_tmp.type = multigrid.atom_type::poly;
  atom_tmp.u_id = 0;
  atom_tmp.value = 5;
  multigrid.eqns[0][16].add_atom(atom_tmp);

  atom_tmp.value = multigrid.atom_type::poly;
  multigrid.eqns[0][17].add_atom(atom_tmp);
  

  
  
  

  real_t avg1 = 0.0, avg5 = 0.0;
  LOOP3(i, j, k)
  {
    real_t value = PI*(
        pw2(psi1[INDEX(i,j,k)]) + pw2(psi2[INDEX(i,j,k)]) + pw2(psi3[INDEX(i,j,k)])
      );
    avg1 += value;
    multigrid.setPolySrcAtPt(0, 17, i, j, k, value); //set value for term 1

    value = 2.0* PI*scalar->V(phi[INDEX(i,j,k)]) - K_a[INDEX(i,j,k)]*K_a[INDEX(i,j,k)]/12.0 + PI * dt_phi * dt_phi;
    multigrid.setPolySrcAtPt(0, 16, i, j, k, value); //set value for term 2
    avg5 += value;
  }
  avg1 = avg1/N/N/N;
  avg5 = avg5/N/N/N;

  /*Finishing adding Hamiltonian constraint**************************************************/

  multigrid.eqns[1][0].init(1, 1.0);
  multigrid.eqns[1][1].init(1, 1.0/3.0);
  multigrid.eqns[1][2].init(1, 1.0/3.0);
  multigrid.eqns[1][3].init(1, 1.0/3.0);
  multigrid.eqns[1][4].init(2, 1.0);

  multigrid.eqns[2][0].init(1, 1.0);
  multigrid.eqns[2][1].init(1, 1.0/3.0);
  multigrid.eqns[2][2].init(1, 1.0/3.0);
  multigrid.eqns[2][3].init(1, 1.0/3.0);
  multigrid.eqns[2][4].init(2, 1.0);

  multigrid.eqns[3][0].init(1, 1.0);
  multigrid.eqns[3][1].init(1, 1.0/3.0);
  multigrid.eqns[3][2].init(1, 1.0/3.0);
  multigrid.eqns[3][3].init(1, 1.0/3.0);
  multigrid.eqns[3][4].init(2, 1.0);

  //adding terms to eqn 1
  atom_tmp.type = multigrid.atom_type::lap;
  atom_tmp.u_id = 1;
  multigrid.eqns[1][0].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::der11;
  atom_tmp.u_id = 1;
  multigrid.eqns[1][1].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::der12;
  atom_tmp.u_id = 2;
  multigrid.eqns[1][2].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::der13;
  atom_tmp.u_id = 3;
  multigrid.eqns[1][3].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::const_f;
  multigrid.eqns[1][4].add_atom(atom_tmp);
  atom_tmp.type = multigrid.atom_type::poly;
  atom_tmp.u_id = 0;
  atom_tmp.value = 6;
  multigrid.eqns[1][4].add_atom(atom_tmp);

  
  //adding terms to eqn 2
  atom_tmp.type = multigrid.atom_type::lap;
  atom_tmp.u_id = 2;
  multigrid.eqns[2][0].add_atom(atom_tmp);


  atom_tmp.type = multigrid.atom_type::der12;
  atom_tmp.u_id = 1;
  multigrid.eqns[2][1].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::der22;
  atom_tmp.u_id = 2;
  multigrid.eqns[2][2].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::der23;
  atom_tmp.u_id = 3;
  multigrid.eqns[2][3].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::const_f;
  multigrid.eqns[2][4].add_atom(atom_tmp);
  atom_tmp.type = multigrid.atom_type::poly;
  atom_tmp.u_id = 0;
  atom_tmp.value = 6;
  multigrid.eqns[2][4].add_atom(atom_tmp);
 
  //adding terms to eqn 3
  atom_tmp.type = multigrid.atom_type::lap;
  atom_tmp.u_id = 3;
  multigrid.eqns[3][0].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::der13;
  atom_tmp.u_id = 1;
  multigrid.eqns[3][1].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::der23;
  atom_tmp.u_id = 2;
  multigrid.eqns[3][2].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::der33;
  atom_tmp.u_id = 3;
  multigrid.eqns[3][3].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::const_f;
  multigrid.eqns[3][4].add_atom(atom_tmp);
  atom_tmp.type = multigrid.atom_type::poly;
  atom_tmp.u_id = 0;
  atom_tmp.value = 6;
  multigrid.eqns[3][4].add_atom(atom_tmp);

  LOOP3(i,j,k)
  {
    idx_t idx = INDEX(i,j,k);
    real_t value = 8.0 * PI * psi1[idx] * dt_phi;
    multigrid.setPolySrcAtPt(1, 4, i, j, k, value);
    value = 8.0 * PI * psi2[idx] * dt_phi;
    multigrid.setPolySrcAtPt(2, 4, i, j, k, value);
    value = 8.0 * PI * psi3[idx] * dt_phi;
    multigrid.setPolySrcAtPt(3, 4, i, j, k, value);
  }

  //finishing adding atoms to each term!!!!!!!!!!!!!!!

  multigrid.initializeRhoHeirarchy();


  
  iodata->log("The average value of coefficient of the first order term is: " + stringify(avg1));
  iodata->log("The average value of coefficient of the fifth order term is: " + stringify(avg5));
  iodata->log("The suggested initial value of multigrid solver is: " + stringify(std::pow(-avg1/avg5,1.0/4.0)));
  iodata->log("The estimated value of gradient energy/potential is:" + stringify(-avg5/PI/2.0/scalar->V(1)));
  iodata->log("The ratio of H_LEN_FRAC/H_0^-1 is: " + stringify(dx*N/(3.0/K_src)));

  
  LOOP3(i, j, k)
  {
    idx_t idx = INDEX(i,j,k);
    X[0][idx] = std::pow(-avg1/avg5,1.0/4.0);
    X[1][idx] = X[2][idx] = X[3][idx] = 0.0;
  }
  //  std::cout<<std::max(psi2.max(),psi3.max())<<"\n";
  multigrid.VCycles(std::stoi(_config["num_v_cycles"]));

  LOOP3(i,j,k)
  {
    idx_t idx = INDEX(i, j, k);

    real_t temp = 0.0;

    for(idx_t kk = 1; kk <=3; kk++)
      temp += derivative(i, j, k, kk, X[kk]);
    
    A11_p[idx] = A11_a[idx] = std::pow(phi_p[idx], -6.0) * (derivative(i, j, k, 1, X[1]) + derivative(i, j, k, 1, X[1]) - 2.0 * temp / 3.0);

    A12_p[idx] = A12_a[idx] = std::pow(phi_p[idx], -6.0) * (derivative(i, j, k, 1, X[2]) + derivative(i, j, k, 2, X[1]));

    A13_p[idx] = A13_a[idx] = std::pow(phi_p[idx], -6.0) * (derivative(i, j, k, 1, X[3]) + derivative(i, j, k, 3, X[1]));

    A22_p[idx] = A22_a[idx] = std::pow(phi_p[idx], -6.0) * (derivative(i, j, k, 2, X[2]) + derivative(i, j, k, 2, X[2]) - 2.0 * temp / 3.0);
    
    A23_p[idx] = A23_a[idx] = std::pow(phi_p[idx], -6.0) * (derivative(i, j, k, 2, X[3]) + derivative(i, j, k, 3, X[2]));

    A33_p[idx] = A33_a[idx] = std::pow(phi_p[idx], -6.0) * (derivative(i, j, k, 3, X[3]) + derivative(i, j, k, 3, X[3]) - 2.0 * temp / 3.0);

    
    phi_p[idx] = std::log(fabs(phi_p[idx]));
    phi_a[idx] = phi_p[idx];
  }
  
}
  
  
void scalar_ic_set_Bowen_York(BSSN * bssn, Scalar * scalar, IOData * iodata)
{
  idx_t i, j, k;
  arr_t * X = new arr_t [3]; //creating 3 variables which can be converted to A^{ij}
  //
  for(i=0; i < 3; i++)
    X[i].init(NX, NY, NZ);

  idx_t molecule_n[3] = {4, 4, 4}; //three terms for each equation

  real_t relaxation_tolerance = std::stod(_config["relaxation_tolerance"]);
  
  FASMultigrid multigrid(X, 3, molecule_n, 4, 5, relaxation_tolerance);
  
  atom atom_tmp = {0};

  multigrid.eqns[0][0].init(1, 1.0);
  multigrid.eqns[0][1].init(1, 1.0/3.0);
  multigrid.eqns[0][2].init(1, 1.0/3.0);
  multigrid.eqns[0][3].init(1, 1.0/3.0);

  multigrid.eqns[1][0].init(1, 1.0);
  multigrid.eqns[1][1].init(1, 1.0/3.0);
  multigrid.eqns[1][2].init(1, 1.0/3.0);
  multigrid.eqns[1][3].init(1, 1.0/3.0);

  multigrid.eqns[2][0].init(1, 1.0);
  multigrid.eqns[2][1].init(1, 1.0/3.0);
  multigrid.eqns[2][2].init(1, 1.0/3.0);
  multigrid.eqns[2][3].init(1, 1.0/3.0);

  //adding terms to eqn 1
  atom_tmp.type = multigrid.atom_type::lap;
  atom_tmp.u_id = 0;
  multigrid.eqns[0][0].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::der11;
  atom_tmp.u_id = 0;
  multigrid.eqns[0][1].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::der12;
  atom_tmp.u_id = 1;
  multigrid.eqns[0][2].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::der13;
  atom_tmp.u_id = 2;
  multigrid.eqns[0][3].add_atom(atom_tmp);

  //adding terms to eqn 2
  atom_tmp.type = multigrid.atom_type::lap;
  atom_tmp.u_id = 1;
  multigrid.eqns[1][0].add_atom(atom_tmp);


  atom_tmp.type = multigrid.atom_type::der12;
  atom_tmp.u_id = 0;
  multigrid.eqns[1][1].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::der22;
  atom_tmp.u_id = 1;
  multigrid.eqns[1][2].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::der23;
  atom_tmp.u_id = 2;
  multigrid.eqns[1][3].add_atom(atom_tmp);
 
  //adding terms to eqn 3
  atom_tmp.type = multigrid.atom_type::lap;
  atom_tmp.u_id = 2;
  multigrid.eqns[2][0].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::der13;
  atom_tmp.u_id = 0;
  multigrid.eqns[2][1].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::der23;
  atom_tmp.u_id = 1;
  multigrid.eqns[2][2].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::der33;
  atom_tmp.u_id = 2;
  multigrid.eqns[2][3].add_atom(atom_tmp);
 
  //finishing adding atoms to each term!!!!!!!!!!!!!!!

  std::random_device rd;
  std::mt19937 gen(7.0 /*rd()*/);
  std::uniform_real_distribution<real_t> dist(0, 2.0*PI);

  
  LOOP3(i,j,k)
  {

    idx_t idx = INDEX(i, j, k);
    idx_t n = 1;
    real_t x_phase = dist(gen),
             y_phase = dist(gen),
             z_phase = dist(gen);
    X[0][idx] = 1 + 0.01 * (
                cos(2.0*PI*((real_t) n/NX)*i + x_phase )
                 + cos(2.0*PI*((real_t) n/NY)*j + y_phase )
                + cos(2.0*PI*((real_t) n/NZ)*k + z_phase ));

    x_phase = dist(gen);
    y_phase = dist(gen);
    z_phase = dist(gen);

    X[1][idx] = 1 + 0.01 * (
                cos(2.0*PI*((real_t) n/NX)*i + x_phase )
                 + cos(2.0*PI*((real_t) n/NY)*j + y_phase )
                + cos(2.0*PI*((real_t) n/NZ)*k + z_phase ));

    x_phase = dist(gen);
    y_phase = dist(gen);
    z_phase = dist(gen);
    X[2][idx] = 1 + 0.01 * (
                cos(2.0*PI*((real_t) n/NX)*i + x_phase )
                 + cos(2.0*PI*((real_t) n/NY)*j + y_phase )
                + cos(2.0*PI*((real_t) n/NZ)*k + z_phase ));
  }

    multigrid.VCycles(std::stoi(_config["num_v_cycles"]));
}
  
#if USE_MULTIGRID
/**
 * @brief Use the multigrid solver to solve for metric factors given
 * a particular scalar field implementation.
 */
void scalar_ic_set_multigrid(BSSN * bssn, Scalar * scalar, IOData * iodata)
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

  iodata->log("K_0 = " + stringify(K_src) + ", H_0 = "
      + stringify(-K_src/3.0) + ", and k/H_0 = "
      + stringify(2.0*PI/(N*dx)/(-K_src/3.0))
    );

   // solve for BSSN fields using multigrid class:

  
  idx_t molecule_n[] = {3};
  arr_t & phi_p = *bssn->fields["DIFFphi_p"];
  arr_t & phi_a = *bssn->fields["DIFFphi_a"];

  arr_t * phi_ini = & phi_p;


  
  FASMultigrid multigrid(phi_ini, 1, molecule_n, 4, 5, relaxation_tolerance);

  
  
  atom atom_tmp = {0};


  //initializing equations
  multigrid.eqns[0][0].init(1, 1);
  multigrid.eqns[0][1].init(2, 1);
  multigrid.eqns[0][2].init(2, 1);
  
  //adding terms to eqn
  //add first laplacian term
  atom_tmp.type = multigrid.atom_type::lap;
  atom_tmp.u_id = 0;
  multigrid.eqns[0][0].add_atom(atom_tmp);

  //add second term
  atom_tmp.type = multigrid.atom_type::const_f;
  multigrid.eqns[0][1].add_atom(atom_tmp);
  
  atom_tmp.type = multigrid.atom_type::poly;
  atom_tmp.u_id = 0;
  atom_tmp.value = 1;
  multigrid.eqns[0][1].add_atom(atom_tmp);

  //add third term
  atom_tmp.type = multigrid.atom_type::const_f;
  multigrid.eqns[0][2].add_atom(atom_tmp);

  atom_tmp.type = multigrid.atom_type::poly;
  atom_tmp.u_id = 0;
  atom_tmp.value = 5;
  multigrid.eqns[0][2].add_atom(atom_tmp);

  
  
  real_t avg1 = 0.0, avg5 = 0.0;
  LOOP3(i, j, k)
  {
    real_t value = PI*(
        pw2(psi1[INDEX(i,j,k)]) + pw2(psi2[INDEX(i,j,k)]) + pw2(psi3[INDEX(i,j,k)])
      );
    avg1 += value;
    multigrid.setPolySrcAtPt(0, 1, i, j, k, value); //set value for term 1

    value = 2.0* PI*scalar->V(phi[INDEX(i,j,k)]) - K_a[INDEX(i,j,k)]*K_a[INDEX(i,j,k)]/12.0;
    multigrid.setPolySrcAtPt(0, 2, i, j, k, value); //set value for term 2
    avg5 += value;

  }
  
  avg1 = avg1/N/N/N;
  avg5 = avg5/N/N/N;
  multigrid.initializeRhoHeirarchy();
  iodata->log("The average value of coefficient of the first order term is: " + stringify(avg1));
  iodata->log("The average value of coefficient of the fifth order term is: " + stringify(avg5));
  iodata->log("The suggested initial value of multigrid solver is: " + stringify(std::pow(-avg1/avg5,1.0/4.0)));
  iodata->log("The estimated value of gradient energy/potential is:" + stringify(-avg5/PI/2.0/scalar->V(1)));
  iodata->log("The ratio of H_LEN_FRAC/H_0^-1 is: " + stringify(dx*N/(3.0/K_src)));

  LOOP3(i, j, k)
  {
    idx_t idx = INDEX(i,j,k);
    phi_ini[0][idx] = std::pow(-avg1/avg5,1.0/4.0);
  }

  multigrid.VCycles(std::stoi(_config["num_v_cycles"]));

  LOOP3(i,j,k)
  {
    idx_t idx = INDEX(i, j, k);
    phi_p[idx] = std::log(fabs(phi_p[idx]));
    phi_a[idx] = phi_p[idx];
  }
  
  
  
  return;
}
#endif // if USE_MULTIGRID

} // namespace cosmo
