#include "sheets_ic.h"
#include "../../cosmo_includes.h"
#include "../../cosmo_types.h"
#include "../../cosmo_globals.h"
#include "../../utils/Fourier.h"
#include "../../utils/math.h"


namespace cosmo
{

// Assigning initial sheets samplings configuration
// by calculating \rho(x) based on initial choice of \phi(x).
// detail is in the note
void sheets_ic_sinusoid_3d(
  BSSN *bssnSim, Sheet *sheetSim, IOData * iodata, real_t & tot_mass)
{
  iodata->log("Setting sinusoidal ICs.");
  idx_t i, j, k;

  // conformal factor
  arr_t & DIFFphi_p = *bssnSim->fields["DIFFphi_p"];
  // DIFFK is initially zero
  arr_t & DIFFK_p = *bssnSim->fields["DIFFK_p"];
  // matter sources
  arr_t & DIFFr_a = *bssnSim->fields["DIFFr_a"];

  arr_t & Dx = sheetSim->Dx._array_p;
  arr_t & Dy = sheetSim->Dy._array_p;
  arr_t & Dz = sheetSim->Dz._array_p;
  
  real_t A = sheetSim->lx*sheetSim->lx*std::stod(_config("peak_amplitude", "0.0001"));
  iodata->log( "Generating ICs with peak amp. = " + stringify(A) );

  real_t rho_FRW = 3.0/PI/8.0;
  real_t K_FRW = -sqrt(24.0*PI*rho_FRW);
  iodata->log( "FRW density is " + stringify(rho_FRW) + ", total mass in simulation volume is "
    + stringify(rho_FRW*sheetSim->lx*sheetSim->ly*sheetSim->lz));

  // the conformal factor in front of metric is the solution to
  // d^2 exp(\phi) = -2*pi exp(5\phi) * \delta_rho
  // generate random mode in \phi
  // delta_rho = -(lap e^\phi)/e^(5\phi)/2pi
  real_t phix = 0, phiy = 0, phiz = 0;
  real_t twopi_L = 2.0*PI / sheetSim->lx;
  real_t pw2_twopi_L = twopi_L*twopi_L;
  
  // grid values
  //#pragma omp parallel for
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    real_t x_frac = ((real_t) i / (real_t) NX),
      y_frac = ((real_t) j / (real_t) NY), z_frac = ((real_t) k / (real_t) NZ);
    
    real_t phi = A*
      (sin(2.0*PI*x_frac + phix) + sin(2.0*PI*y_frac + phiy) + sin(2.0*PI*z_frac + phiz));

    real_t rho = rho_FRW -exp(-4.0*phi)/PI/2.0*(
      pw2(A*2.0*PI / sheetSim->lx * cos(2.0*PI*x_frac + phix))
      + pw2(A*2.0*PI / sheetSim->ly * cos(2.0*PI*y_frac + phiy))
      + pw2(A*2.0*PI / sheetSim->lz * cos(2.0*PI*z_frac + phiz))
      -  ( A * pw2(2.0*PI / sheetSim->lx) * sin(2.0*PI*x_frac + phix)
           + A * pw2(2.0*PI / sheetSim->ly) * sin(2.0*PI*y_frac + phiy)
           + A * pw2(2.0*PI / sheetSim->lz) * sin(2.0*PI*z_frac + phiz))
    );

    // These aren't difference vars
    DIFFphi_p[NP_INDEX(i,j,k)] = phi;
    DIFFK_p[idx] = K_FRW;
    DIFFr_a[idx] = rho;
  }
  
  idx_t integration_points_per_dx = std::stod(_config("integration_points_per_dx", "1000"));
  idx_t integration_points_per_dy = integration_points_per_dx;
  idx_t integration_points_per_dz = integration_points_per_dx;
  idx_t integration_points_x = NX * integration_points_per_dx;
  idx_t integration_points_y = NY * integration_points_per_dx;
  idx_t integration_points_z = NZ * integration_points_per_dx;
  std::cout << "Setting initial conditions using " << integration_points_x << " integration_points" << std::endl;
  real_t integration_interval_x = sheetSim->lx / integration_points_x;
  real_t integration_interval_y = sheetSim->ly / integration_points_y;
  real_t integration_interval_z = sheetSim->lz / integration_points_z;

  tot_mass = 0;

  real_t tot_vol = 0;
  // compute total mass in simulation, mass per tracer particle

#pragma omp parallel for collapse(2) reduction(+:tot_mass, tot_vol)
  for(idx_t i=0; i<integration_points_x; ++i)
    for(idx_t j=0; j<integration_points_y; ++j)
      for(idx_t k=0; k<integration_points_z; ++k)
      {
        real_t x_frac = ((real_t) i / (real_t) integration_points_x),
          y_frac = ((real_t) j / (real_t) integration_points_y),
          z_frac = ((real_t) k / (real_t) integration_points_z);
    
        real_t phi = A*
          (sin(2.0*PI*x_frac + phix) + sin(2.0*PI*y_frac + phiy) + sin(2.0*PI*z_frac + phiz));

        real_t rho = rho_FRW -exp(-4.0*phi)/PI/2.0*(
          pw2(A*2.0*PI / sheetSim->lx * cos(2.0*PI*x_frac + phix))
          + pw2(A*2.0*PI / sheetSim->ly * cos(2.0*PI*y_frac + phiy))
          + pw2(A*2.0*PI / sheetSim->lz * cos(2.0*PI*z_frac + phiz))
          -  ( A * pw2(2.0*PI / sheetSim->lx) * sin(2.0*PI*x_frac + phix)
               + A * pw2(2.0*PI / sheetSim->ly) * sin(2.0*PI*y_frac + phiy)
               + A * pw2(2.0*PI / sheetSim->lz) * sin(2.0*PI*z_frac + phiz))
        );

        real_t rootdetg = std::exp(6.0*phi);

        if(rho < 0)
        {
          std::cout<<"Density "<<rho<<" is smaller than zero at "
                   <<i<<" "<<j<<" "<<k<<"\n";
          throw(-1);
        }
        
        tot_mass += rho * rootdetg *
          integration_interval_x * integration_interval_y * integration_interval_z;
        
        tot_vol += rootdetg * integration_interval_x * integration_interval_y * integration_interval_z;
      }

  // density on the sheets
  // NOT mass in each voxel
  real_t rho_s = tot_mass / sheetSim->lx / sheetSim->ly / sheetSim->lz;

  std::cout << "Total mass and total volum are " << tot_mass<<" "<<tot_vol<<".\n";

  // s1(x,y,z), s2(y,z), s3(z)
  arr_t s1(NX, NY, NZ) , s2(NX, NY, NZ), s3(NX, NY, NZ);

  // just for test
  //  arr_t g(NX, NY , NZ), f(NX, NY, NZ);
  
  real_t cur_z = 0;
    
  for(int k = 0; k < integration_points_z; k++)
  {
    real_t cur_xy = 0;
    #pragma omp parallel for collapse(2) reduction(+:cur_xy)
    for(int i = 0; i < integration_points_x; i++)
      for(int j = 0; j < integration_points_y; j++)
      {
        real_t x_frac = (real_t)i/(real_t) integration_points_x;
        real_t y_frac = (real_t)j/(real_t) integration_points_y;

        real_t z_frac = ((real_t) k / (real_t) integration_points_z);
    
        real_t phi = A*
          (sin(2.0*PI*x_frac + phix) + sin(2.0*PI*y_frac + phiy) + sin(2.0*PI*z_frac + phiz));

        real_t rho = rho_FRW -exp(-4.0*phi)/PI/2.0*(
          pw2(A*2.0*PI / sheetSim->lx * cos(2.0*PI*x_frac + phix))
          + pw2(A*2.0*PI / sheetSim->ly * cos(2.0*PI*y_frac + phiy))
          + pw2(A*2.0*PI / sheetSim->lz * cos(2.0*PI*z_frac + phiz))
          -  ( A * pw2(2.0*PI / sheetSim->lx) * sin(2.0*PI*x_frac + phix)
               + A * pw2(2.0*PI / sheetSim->ly) * sin(2.0*PI*y_frac + phiy)
               + A * pw2(2.0*PI / sheetSim->lz) * sin(2.0*PI*z_frac + phiz))
        );

        cur_xy += rho * integration_interval_x * integration_interval_y;
      }

    cur_z += cur_xy * integration_interval_z / rho_s / sheetSim->lx / sheetSim->ly
      - integration_interval_z;

    if( (k%integration_points_per_dz == 0))
    {
#pragma omp parallel for collapse(2)
      for(int i = 0; i < NX; i++)
        for(int j = 0; j < NY; j++)
        {
          s3[NP_INDEX(i, j, k / integration_points_per_dz)] = cur_z;
        }
    }
  }

  for(int k = 0; k < NZ; k++)
  {
    real_t cur_xy_num = 0;
    real_t cur_xy_den = 0;
    
    #pragma omp parallel for collapse(2) reduction(+:cur_xy_den)
    for(int j = 0; j < integration_points_y; j++)
    {
      for(int i = 0; i < integration_points_x; i++)
      {
        real_t x_frac = (real_t)i/(real_t) integration_points_x;
        real_t y_frac = (real_t)j/(real_t) integration_points_y;

        real_t z_frac = ((real_t) k / (real_t) NZ);
    
        real_t phi = A*
          (sin(2.0*PI*x_frac + phix) + sin(2.0*PI*y_frac + phiy) + sin(2.0*PI*z_frac + phiz));

        real_t rho = rho_FRW -exp(-4.0*phi)/PI/2.0*(
          pw2(A*2.0*PI / sheetSim->lx * cos(2.0*PI*x_frac + phix))
          + pw2(A*2.0*PI / sheetSim->ly * cos(2.0*PI*y_frac + phiy))
          + pw2(A*2.0*PI / sheetSim->lz * cos(2.0*PI*z_frac + phiz))
          -  ( A * pw2(2.0*PI / sheetSim->lx) * sin(2.0*PI*x_frac + phix)
               + A * pw2(2.0*PI / sheetSim->ly) * sin(2.0*PI*y_frac + phiy)
               + A * pw2(2.0*PI / sheetSim->lz) * sin(2.0*PI*z_frac + phiz))
        );

        cur_xy_den += rho * integration_interval_x * integration_interval_y;
      }
    }

    for(int j = 0; j < integration_points_y; j++)
    {
      real_t cur_x = 0;
      for(int i = 0; i < integration_points_x; i++)
      {
        real_t x_frac = (real_t)i/(real_t) integration_points_x;
        real_t y_frac = (real_t)j/(real_t) integration_points_y;

        real_t z_frac = ((real_t) k / (real_t) NZ);
    
        real_t phi = A*
          (sin(2.0*PI*x_frac + phix) + sin(2.0*PI*y_frac + phiy) + sin(2.0*PI*z_frac + phiz));

        real_t rho = rho_FRW -exp(-4.0*phi)/PI/2.0*(
          pw2(A*2.0*PI / sheetSim->lx * cos(2.0*PI*x_frac + phix))
          + pw2(A*2.0*PI / sheetSim->ly * cos(2.0*PI*y_frac + phiy))
          + pw2(A*2.0*PI / sheetSim->lz * cos(2.0*PI*z_frac + phiz))
          -  ( A * pw2(2.0*PI / sheetSim->lx) * sin(2.0*PI*x_frac + phix)
               + A * pw2(2.0*PI / sheetSim->ly) * sin(2.0*PI*y_frac + phiy)
               + A * pw2(2.0*PI / sheetSim->lz) * sin(2.0*PI*z_frac + phiz))
        );

        cur_x += rho * integration_interval_x;
      }

      cur_xy_num +=  cur_x * integration_interval_y * sheetSim->ly / cur_xy_den - integration_interval_y;
      
      if(  j%integration_points_per_dy == 0)
      {
        for(int i = 0; i < NX; i++)
          s2[NP_INDEX(i, j/integration_points_per_dy, k)] = cur_xy_num;
      }
      
    }
  }

#pragma omp parallel for collapse(2)
  for(int k = 0; k < NZ; k++)
  {
    for(int j = 0; j < NY; j++)
    {
      real_t cur_x_num = 0;
      real_t cur_x_den = 0;
    
      for(int i = 0; i < integration_points_x; i++)
      {
        real_t x_frac = (real_t)i/(real_t) integration_points_x;
        real_t y_frac = (real_t)j/(real_t) NY;

        real_t z_frac = ((real_t) k / (real_t) NZ);
    
        real_t phi = A*
          (sin(2.0*PI*x_frac + phix) + sin(2.0*PI*y_frac + phiy) + sin(2.0*PI*z_frac + phiz));

        real_t rho = rho_FRW -exp(-4.0*phi)/PI/2.0*(
          pw2(A*2.0*PI / sheetSim->lx * cos(2.0*PI*x_frac + phix))
          + pw2(A*2.0*PI / sheetSim->ly * cos(2.0*PI*y_frac + phiy))
          + pw2(A*2.0*PI / sheetSim->lz * cos(2.0*PI*z_frac + phiz))
          -  ( A * pw2(2.0*PI / sheetSim->lx) * sin(2.0*PI*x_frac + phix)
               + A * pw2(2.0*PI / sheetSim->ly) * sin(2.0*PI*y_frac + phiy)
               + A * pw2(2.0*PI / sheetSim->lz) * sin(2.0*PI*z_frac + phiz))
        );

        cur_x_den += rho * integration_interval_x;
      }

      for(int i = 0; i < integration_points_x; i++)
      {
        real_t x_frac = (real_t)i/(real_t) integration_points_x;
        real_t y_frac = (real_t)j/(real_t) NY;

        real_t z_frac = ((real_t) k / (real_t) NZ);
    
        real_t phi = A*
          (sin(2.0*PI*x_frac + phix) + sin(2.0*PI*y_frac + phiy) + sin(2.0*PI*z_frac + phiz));

        real_t rho = rho_FRW -exp(-4.0*phi)/PI/2.0*(
          pw2(A*2.0*PI / sheetSim->lx * cos(2.0*PI*x_frac + phix))
          + pw2(A*2.0*PI / sheetSim->ly * cos(2.0*PI*y_frac + phiy))
          + pw2(A*2.0*PI / sheetSim->lz * cos(2.0*PI*z_frac + phiz))
          -  ( A * pw2(2.0*PI / sheetSim->lx) * sin(2.0*PI*x_frac + phix)
               + A * pw2(2.0*PI / sheetSim->ly) * sin(2.0*PI*y_frac + phiy)
               + A * pw2(2.0*PI / sheetSim->lz) * sin(2.0*PI*z_frac + phiz))
        );
        
        cur_x_num += rho * integration_interval_x  * sheetSim->lx / cur_x_den - integration_interval_x;

        if( i%integration_points_per_dx == 0)
        {
          s1[NP_INDEX(i/integration_points_per_dx, j, k)] =
            cur_x_num;
        }
      }
    }
  }

  real_t max_err = 0, max_err_i = 0;

  #pragma omp parallel for collapse(2) reduction(max:max_err)
  for(int i = 0 ; i < NX; i ++)
    for(int j = 0; j < NY; j++)
      for(int k = 0; k < NZ; k++)
      {
        real_t err = fabs(DIFFr_a[NP_INDEX(i, j, k)] / rho_s
                          - ((derivative(i, j, k, 1, s1) + 1.0)
                             * (derivative(0, j, k, 2, s2) + 1.0)
                             * (derivative(0, 0, k, 3, s3) + 1.0)));
        if(err > max_err) max_err_i = i;
        max_err = std::max(max_err,err);
      }
  
  std::cout<<"Max error to Eq. 10 is "<<max_err<<" at "<<max_err_i<<"\n";

  const int b_search_interation_limit = 20;
  real_t max_inverse_deviation = 0;

  // reverse to get s3
#pragma omp parallel for
  for(int k =0; k < sheetSim->ns3; k++)
  {
    real_t cur_s3 = ((real_t)k / sheetSim->ns3) * sheetSim->lz;

    real_t z = cur_s3;

    real_t dz_cur = 0;

    real_t dz_lower = -z, dz_upper = (real_t)sheetSim->lz - z;
        
    int iter_cnt = 0;
    real_t x_idx, y_idx, z_idx;
    while(iter_cnt <= b_search_interation_limit)
    {
      dz_cur = (dz_lower + dz_upper) / 2.0;

      z_idx = (z + dz_cur) / dx;

      if(z+dz_cur + s3.getTriCubicInterpolatedValue(0,0,z_idx) < cur_s3)
      {
        dz_lower = dz_cur;
      }
      else
      {
        dz_upper = dz_cur;
      }
      iter_cnt++;
    }

    for(int i = 0; i < sheetSim->ns1; i++)
      for(int j = 0; j < sheetSim->ns2; j++)
        Dz(i, j, k) = dz_cur;

  }

  // reverse to get s2
#pragma omp parallel for collapse(2)
  for(int j =0; j < sheetSim->ns2; j++)
  {
    for(int k =0; k < sheetSim->ns3; k++)
    {
      real_t cur_s2 = ((real_t)j / sheetSim->ns2) * sheetSim->ly;
      real_t cur_s3 = ((real_t)k / sheetSim->ns3) * sheetSim->lz;

      real_t z = cur_s3 + Dz(0, 0, k);
      real_t y = cur_s2;

      real_t dy_cur = 0;

      real_t dy_lower = -y, dy_upper = (real_t)sheetSim->ly - y;
        
      int iter_cnt = 0;
      real_t x_idx, y_idx, z_idx;
      while(iter_cnt <= b_search_interation_limit)
      {
        dy_cur = (dy_lower + dy_upper) / 2.0;

        y_idx = (y + dy_cur) / dx;
        z_idx = z / dx;

        if(y+dy_cur + s2.getTriCubicInterpolatedValue(0,y_idx,z_idx) < cur_s2)
        {
          dy_lower = dy_cur;
        }
        else
        {
          dy_upper = dy_cur;
        }
        iter_cnt++;
      }

      for(int i = 0; i < sheetSim->ns1; i++)
          Dy(i, j, k) = dy_cur;

    }
  }

  // reverset to get s1
#pragma omp parallel for collapse(2) reduction(max:max_inverse_deviation)
  for(int i =0; i < sheetSim->ns1; i++)
  {
    for(int j = 0; j < sheetSim->ns2; j++)
    {
      for(int k =0; k < sheetSim->ns3; k ++)
      {
        real_t cur_s1 = ((real_t)i / sheetSim->ns1) * sheetSim->lx;
        real_t cur_s2 = ((real_t)j / sheetSim->ns2) * sheetSim->ly;
        real_t cur_s3 = ((real_t)k / sheetSim->ns3) * sheetSim->lz;

        real_t x = cur_s1;
        real_t y = cur_s2 + Dy(0, j, k);
        real_t z = cur_s3 + Dz(0 , 0, k);

        real_t dx_cur = 0;

        real_t dx_lower = -x, dx_upper = (real_t)sheetSim->lx - x;
        
        int iter_cnt = 0;
        real_t x_idx, y_idx, z_idx;
        while(iter_cnt <= b_search_interation_limit)
        {
          dx_cur = (dx_lower + dx_upper) / 2.0;
          x_idx = (x + dx_cur) / dx;
          y_idx = y / dx;
          z_idx = z / dx;

          if(x+dx_cur + s1.getTriCubicInterpolatedValue(x_idx,y_idx,z_idx) < cur_s1)
          {
            dx_lower = dx_cur;
          }
          else
          {
            dx_upper = dx_cur;
          }
          iter_cnt++;
        }
        Dx(i, j, k) = dx_cur;
        max_inverse_deviation = std::max(max_inverse_deviation,
                                         fabs(x+dx_cur+s1.getTriCubicInterpolatedValue(x_idx,y_idx,z_idx) - cur_s1));
        max_inverse_deviation = std::max(max_inverse_deviation,
                                         fabs(y+s2.getTriCubicInterpolatedValue(x_idx,y_idx,z_idx) - cur_s2));
        max_inverse_deviation = std::max(max_inverse_deviation,
                                         fabs(z+s3.getTriCubicInterpolatedValue(x_idx,y_idx,z_idx) - cur_s3));

      }
    }
  }
  

  std::cout<<"The inversion of function brings an deviation of "<<max_inverse_deviation<<"\n";
    
  real_t max_dev = 0;
  int maxi, maxj, maxk;
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    real_t x_frac = ((real_t) i / (real_t) NX),
      y_frac = ((real_t) j / (real_t) NY), z_frac = ((real_t) k / (real_t) NZ);
    
    real_t phi = A*
      (sin(2.0*PI*x_frac + phix) + sin(2.0*PI*y_frac + phiy) + sin(2.0*PI*z_frac + phiz));

    real_t rho = rho_FRW -exp(-4.0*phi)/PI/2.0*(
      pw2(A*2.0*PI / sheetSim->lx * cos(2.0*PI*x_frac + phix))
      + pw2(A*2.0*PI / sheetSim->ly * cos(2.0*PI*y_frac + phiy))
      + pw2(A*2.0*PI / sheetSim->lz * cos(2.0*PI*z_frac + phiz))
      -  ( A * pw2(2.0*PI / sheetSim->lx) * sin(2.0*PI*x_frac + phix)
           + A * pw2(2.0*PI / sheetSim->ly) * sin(2.0*PI*y_frac + phiy)
           + A * pw2(2.0*PI / sheetSim->lz) * sin(2.0*PI*z_frac + phiz))
    );

    // These aren't difference vars
    max_dev = std::max(max_dev, fabs(DIFFr_a[NP_INDEX(i,j,k)] - rho));
  }
  std::cout<<"max dev is "<<max_dev<<" at "<<maxi<<" "<<maxj<<" "<<maxk<<"\n";
}

void sheets_ic_sinusoid(
  BSSN *bssnSim, Sheet *sheetSim, Lambda * lambda, IOData * iodata, real_t & tot_mass)
{
  iodata->log("Setting sinusoidal ICs.");
  idx_t i, j, k;

  // conformal factor
  arr_t & DIFFphi_p = *bssnSim->fields["DIFFphi_p"];
  // DIFFK is initially zero
  arr_t & DIFFK_p = *bssnSim->fields["DIFFK_p"];
  // matter sources
  arr_t & DIFFr_a = *bssnSim->fields["DIFFr_a"];

  arr_t & Dx = sheetSim->Dx._array_p;
  
  real_t A = sheetSim->lx*sheetSim->lx*std::stod(_config("peak_amplitude", "0.0001"));
  iodata->log( "Generating ICs with peak amp. = " + stringify(A) );

  real_t K_FRW = -3.0;
  real_t rho_FRW = 3.0/PI/8.0;
  
  real_t Omega_L = std::stod(_config("Omega_L", "0.0"));
  real_t rho_m = (1.0 - Omega_L) * rho_FRW;
  real_t rho_L = Omega_L * rho_FRW;
  lambda->setLambda(rho_L);

  iodata->log( "Total density is " + stringify(rho_FRW) + ". Matter density is "
    + stringify(rho_m) + ", total matter mass in simulation volume is "
    + stringify(rho_m*sheetSim->lx*sheetSim->ly*sheetSim->lz));

  // the conformal factor in front of metric is the solution to
  // d^2 exp(\phi) = -2*pi exp(5\phi) * \delta_rho
  // generate random mode in \phi
  // delta_rho = -(lap e^\phi)/e^(4\phi)/2pi
  real_t phix = 0;
  real_t twopi_L = 2.0*PI / sheetSim->lx;
  real_t pw2_twopi_L = twopi_L*twopi_L;
  
  // grid values
  //#pragma omp parallel for
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    real_t x_frac = ((real_t) i / (real_t) NX);
    real_t phi = A*sin(2.0*PI*x_frac + phix);
    real_t rho = rho_m + -exp(-4.0*phi)/PI/2.0*(
      pw2(twopi_L*A*cos(2.0*PI*x_frac + phix))
      - pw2_twopi_L*A*sin(2.0*PI*x_frac + phix)
    );

    // These aren't difference vars
    DIFFphi_p[NP_INDEX(i,j,k)] = phi;
    DIFFK_p[idx] = K_FRW;
    DIFFr_a[idx] = rho;
  }

  real_t integration_points = sheetSim->ns1 * std::stod(_config("integration_points_per_dx", "1000"));
  std::cout << "Setting initial conditions using " << integration_points << " integration_points" << std::endl;
  real_t integration_interval = sheetSim->lx / integration_points;
  tot_mass = 0;

  // compute total mass in simulation, mass per tracer particle
  for(idx_t i=0; i<integration_points; ++i)
  {
    real_t x_frac = i/(real_t) integration_points;

    real_t phi = A*sin(2.0*PI*x_frac + phix);
    real_t rho = rho_m + -exp(-4.0*phi)/PI/2.0*(
      pw2(twopi_L*A*cos(2.0*PI*x_frac + phix))
      - pw2_twopi_L*A*sin(2.0*PI*x_frac + phix)
    );

    real_t rootdetg = std::exp(6.0*phi);

    tot_mass += rho * rootdetg * integration_interval * sheetSim->ly * sheetSim->lz;
  }

  real_t mass_per_tracer = tot_mass / (real_t) (sheetSim->ns1);

  std::cout << "Total mass and mass_per_tracer are " << tot_mass
    << ", " << mass_per_tracer << ".\n";


  // Cumulatively integrate density, deposit particles when integral reaches a particle mass
  idx_t cur_s1 = 1; // (boundary condition: Dx(x=0) = 0, so start positioning s1=1 particle
  real_t cur_mass = 0;
  for(i=0; i<=integration_points; ++i)
  {
    real_t x_frac = i/(real_t) integration_points;
    real_t x = sheetSim->lx * x_frac;

    real_t phi = A*sin(2.0*PI*x_frac + phix);
    real_t rho = rho_m + -exp(-4.0*phi)/PI/2.0*(
      pw2(twopi_L*A*cos(2.0*PI*x_frac + phix))
      - pw2_twopi_L*A*sin(2.0*PI*x_frac + phix)
    );


    real_t rootdetg = std::exp(6.0*phi);
    cur_mass += rho * rootdetg * integration_interval * sheetSim->ly * sheetSim->lz;

    if(cur_mass >= mass_per_tracer)
    {
      for(j=0; j<sheetSim->ns2; ++j)
        for(k=0; k<sheetSim->ns3; ++k)
        {
          Dx(cur_s1, j, k) = x - sheetSim->_S1IDXtoX0(cur_s1)
           + integration_interval * (cur_mass - mass_per_tracer) / mass_per_tracer;
        }
      cur_mass = cur_mass - mass_per_tracer;
      cur_s1++;
      if(cur_s1 == sheetSim->ns1) break;
    }
  }

  if(cur_s1 < sheetSim->ns1 - 1)
  {
    std::cout<<"Error in setting initial distribution!\n";
    throw(-1);
  }
}


// Assigning initial sheets samplings configuration
// by calculating \rho(x) based on initial choice of \phi(x).
// detail is in the note
void sheets_ic_sinusoid_3d_diffusion(
  BSSN *bssnSim, Sheet *sheetSim, IOData * iodata, real_t & tot_mass)
{
  iodata->log("Setting sinusoidal 3D ICs with diffusion method");
  idx_t i, j, k;


    // conformal factor
  arr_t & DIFFphi_p = *bssnSim->fields["DIFFphi_p"];
  // DIFFK is initially zero
  arr_t & DIFFK_p = *bssnSim->fields["DIFFK_p"];
  // matter sources
  arr_t & DIFFr_a = *bssnSim->fields["DIFFr_a"];

  arr_t & Dx_p = sheetSim->Dx._array_p;
  arr_t & Dy_p = sheetSim->Dy._array_p;
  arr_t & Dz_p = sheetSim->Dz._array_p;

  arr_t & Dx_a = sheetSim->Dx._array_a;
  arr_t & Dy_a = sheetSim->Dy._array_a;
  arr_t & Dz_a = sheetSim->Dz._array_a;
  
  real_t A = sheetSim->lx*sheetSim->lx*std::stod(_config("peak_amplitude", "0.0001"));

  // setting iteration stuff
  real_t precision_goal = sheetSim->lx * sheetSim->lx * std::stod(_config("precision_goal", "1e-5"));
  real_t damping_coef = sheetSim->lx * sheetSim->lx * std::stod(_config("damping_coef", "0.1"));
  
  iodata->log( "Generating ICs with peak amp. = " + stringify(A) );

  real_t rho_FRW = 3.0/PI/8.0;
  real_t K_FRW = -sqrt(24.0*PI*rho_FRW);
  iodata->log( "FRW density is " + stringify(rho_FRW) + ", total mass in simulation volume is "
    + stringify(rho_FRW*sheetSim->lx*sheetSim->ly*sheetSim->lz));

  // the conformal factor in front of metric is the solution to
  // d^2 exp(\phi) = -2*pi exp(5\phi) * \delta_rho
  // generate random mode in \phi
  // delta_rho = -(lap e^\phi)/e^(5\phi)/2pi
  real_t phix = 0, phiy = 0, phiz = 0;
  real_t twopi_L = 2.0*PI / sheetSim->lx;
  real_t pw2_twopi_L = twopi_L*twopi_L;

  // defining difference of \rho
  arr_t rho(NX, NY, NZ), rho_err(NX, NY, NZ);
  // variables to store gradiant of delta \rho
  arr_t drhodx(NX, NY, NZ), drhody(NX, NY, NZ), drhodz(NX, NY, NZ);

  arr_t fourier_temp(NX, NY, NZ);

  arr_t d1phi(NX, NY, NZ), d2phi(NX, NY, NZ), d3phi(NX, NY, NZ);
  
  // grid values
  //#pragma omp parallel for
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    real_t x_frac = ((real_t) i / (real_t) NX),
      y_frac = ((real_t) j / (real_t) NY), z_frac = ((real_t) k / (real_t) NZ);
    
    real_t phi = A*
      (sin(2.0*PI*x_frac + phix) + sin(2.0*PI*y_frac + phiy) + sin(2.0*PI*z_frac + phiz));

    real_t rho_c = rho_FRW -exp(-4.0*phi)/PI/2.0*(
      pw2(A*2.0*PI / sheetSim->lx * cos(2.0*PI*x_frac + phix))
      + pw2(A*2.0*PI / sheetSim->ly * cos(2.0*PI*y_frac + phiy))
      + pw2(A*2.0*PI / sheetSim->lz * cos(2.0*PI*z_frac + phiz))
      -  ( A * pw2(2.0*PI / sheetSim->lx) * sin(2.0*PI*x_frac + phix)
           + A * pw2(2.0*PI / sheetSim->ly) * sin(2.0*PI*y_frac + phiy)
           + A * pw2(2.0*PI / sheetSim->lz) * sin(2.0*PI*z_frac + phiz))
    );

    // These aren't difference vars
    DIFFphi_p[NP_INDEX(i,j,k)] = phi;
    DIFFK_p[idx] = K_FRW;
    // set target \rho to rho temperarily
    rho[idx] = rho_c;

    //    real_t rootdetg = std::exp(6.0*phi);
    
  }
  

  idx_t integration_points_per_dx = std::stod(_config("integration_points_per_dx", "1000"));
  idx_t integration_points_per_dy = integration_points_per_dx;
  idx_t integration_points_per_dz = integration_points_per_dx;
  idx_t integration_points_x = NX * integration_points_per_dx;
  idx_t integration_points_y = NY * integration_points_per_dx;
  idx_t integration_points_z = NZ * integration_points_per_dx;
  std::cout << "Setting initial conditions using " << integration_points_x << " integration_points" << std::endl;
  real_t integration_interval_x = sheetSim->lx / integration_points_x;
  real_t integration_interval_y = sheetSim->ly / integration_points_y;
  real_t integration_interval_z = sheetSim->lz / integration_points_z;

  tot_mass = 0;

  real_t tot_vol = 0;
  // compute total mass in simulation, mass per tracer particle

#pragma omp parallel for collapse(2) reduction(+:tot_mass, tot_vol)
  for(idx_t i=0; i<integration_points_x; ++i)
    for(idx_t j=0; j<integration_points_y; ++j)
      for(idx_t k=0; k<integration_points_z; ++k)
      {
        real_t x_frac = ((real_t) i / (real_t) integration_points_x),
          y_frac = ((real_t) j / (real_t) integration_points_y),
          z_frac = ((real_t) k / (real_t) integration_points_z);
    
        real_t phi = A*
          (sin(2.0*PI*x_frac + phix) + sin(2.0*PI*y_frac + phiy) + sin(2.0*PI*z_frac + phiz));

        real_t rho_c = rho_FRW -exp(-4.0*phi)/PI/2.0*(
          pw2(A*2.0*PI / sheetSim->lx * cos(2.0*PI*x_frac + phix))
          + pw2(A*2.0*PI / sheetSim->ly * cos(2.0*PI*y_frac + phiy))
          + pw2(A*2.0*PI / sheetSim->lz * cos(2.0*PI*z_frac + phiz))
          -  ( A * pw2(2.0*PI / sheetSim->lx) * sin(2.0*PI*x_frac + phix)
               + A * pw2(2.0*PI / sheetSim->ly) * sin(2.0*PI*y_frac + phiy)
               + A * pw2(2.0*PI / sheetSim->lz) * sin(2.0*PI*z_frac + phiz))
        );

        real_t rootdetg = std::exp(6.0*phi);

        if(rho_c < 0)
        {
          std::cout<<"Density "<<rho_c<<" is smaller than zero at "
                   <<i<<" "<<j<<" "<<k<<"\n";
          throw(-1);
        }
        
        tot_mass += rho_c * rootdetg * 
          integration_interval_x * integration_interval_y * integration_interval_z;
        
        tot_vol += rootdetg * integration_interval_x * integration_interval_y * integration_interval_z;
      }

  std::cout<<"Total mass is "<<tot_mass<<"\n";

  // setting proper initial guess

#pragma omp parallel for collapse(2)
  for(idx_t i = 0; i < NX; i++)
    for(idx_t j = 0; j < NY; j++)
      for(idx_t k = 0; k < NZ; k++)
      {
        idx_t idx = NP_INDEX(i,j,k);

        real_t x_frac = ((real_t) i / (real_t) NX),
          y_frac = ((real_t) j / (real_t) NY), z_frac = ((real_t) k / (real_t) NZ);
    
        real_t phi = A*
          (sin(2.0*PI*x_frac + phix) + sin(2.0*PI*y_frac + phiy) + sin(2.0*PI*z_frac + phiz));

        real_t rho_c = rho_FRW -exp(-4.0*phi)/PI/2.0*(
          pw2(A*2.0*PI / sheetSim->lx * cos(2.0*PI*x_frac + phix))
          + pw2(A*2.0*PI / sheetSim->ly * cos(2.0*PI*y_frac + phiy))
          + pw2(A*2.0*PI / sheetSim->lz * cos(2.0*PI*z_frac + phiz))
          -  ( A * pw2(2.0*PI / sheetSim->lx) * sin(2.0*PI*x_frac + phix)
               + A * pw2(2.0*PI / sheetSim->ly) * sin(2.0*PI*y_frac + phiy)
               + A * pw2(2.0*PI / sheetSim->lz) * sin(2.0*PI*z_frac + phiz))
        );
        real_t rootdetg = std::exp(6.0*phi);
        fourier_temp[idx] = sheetSim->lx * sheetSim->ly * sheetSim->lz * rho_c * rootdetg / tot_mass - 1.0;
      }

    // trying to set better initial guess, but not working very well 

  Fourier * fourier;
  fourier = new Fourier();
  fourier->Initialize(NX, NY, NZ,
    bssnSim->fields["DIFFphi_a"]->_array);

  fourier->inverseLaplacian <idx_t, real_t> (fourier_temp._array);

  // generating arrays of derivative of phi

  #pragma omp parallel for collapse(2)
  for(idx_t i = 0; i < NX; i++)
    for(idx_t j = 0; j < NY; j++)
      for(idx_t k = 0; k < NZ; k++)
      {
        idx_t idx = NP_INDEX(i,j,k);
        d1phi[idx] = derivative(i, j, k, 1, fourier_temp);
        d2phi[idx] = derivative(i, j, k, 2, fourier_temp);
        d3phi[idx] = derivative(i, j, k, 3, fourier_temp);
      }

  // staring inverse process
  const int b_search_interation_limit = 20;
  real_t max_inverse_deviation = 0;

  // reverse to get s3
#pragma omp parallel for
  for(int k =0; k < sheetSim->ns3; k++)
  {
    real_t cur_s3 = ((real_t)k / sheetSim->ns3) * sheetSim->lz;

    real_t z = cur_s3;

    real_t dz_cur = 0;

    real_t dz_lower = -z, dz_upper = (real_t)sheetSim->lz - z;
        
    int iter_cnt = 0;
    real_t x_idx, y_idx, z_idx;
    while(iter_cnt <= b_search_interation_limit)
    {
      dz_cur = (dz_lower + dz_upper) / 2.0;

      z_idx = (z + dz_cur) / dx;

      if(z+dz_cur + d3phi.getTriCubicInterpolatedValue(0,0,z_idx) < cur_s3)
      {
        dz_lower = dz_cur;
      }
      else
      {
        dz_upper = dz_cur;
      }
      iter_cnt++;
    }

    for(int i = 0; i < sheetSim->ns1; i++)
      for(int j = 0; j < sheetSim->ns2; j++)
        Dz_p(i, j, k) = dz_cur;

  }

  // reverse to get s2
#pragma omp parallel for collapse(2)
  for(int j =0; j < sheetSim->ns2; j++)
  {
    for(int k =0; k < sheetSim->ns3; k++)
    {
      real_t cur_s2 = ((real_t)j / sheetSim->ns2) * sheetSim->ly;
      real_t cur_s3 = ((real_t)k / sheetSim->ns3) * sheetSim->lz;

      real_t z = cur_s3 + Dz_p(0, 0, k);
      real_t y = cur_s2;

      real_t dy_cur = 0;

      real_t dy_lower = -y, dy_upper = (real_t)sheetSim->ly - y;
        
      int iter_cnt = 0;
      real_t x_idx, y_idx, z_idx;
      while(iter_cnt <= b_search_interation_limit)
      {
        dy_cur = (dy_lower + dy_upper) / 2.0;

        y_idx = (y + dy_cur) / dx;
        z_idx = z / dx;

        if(y+dy_cur + d2phi.getTriCubicInterpolatedValue(0,y_idx,z_idx) < cur_s2)
        {
          dy_lower = dy_cur;
        }
        else
        {
          dy_upper = dy_cur;
        }
        iter_cnt++;
      }

      for(int i = 0; i < sheetSim->ns1; i++)
          Dy_p(i, j, k) = dy_cur;

    }
  }

  // reverset to get s1
#pragma omp parallel for collapse(2) reduction(max:max_inverse_deviation)
  for(int i =0; i < sheetSim->ns1; i++)
  {
    for(int j = 0; j < sheetSim->ns2; j++)
    {
      for(int k =0; k < sheetSim->ns3; k ++)
      {
        real_t cur_s1 = ((real_t)i / sheetSim->ns1) * sheetSim->lx;
        real_t cur_s2 = ((real_t)j / sheetSim->ns2) * sheetSim->ly;
        real_t cur_s3 = ((real_t)k / sheetSim->ns3) * sheetSim->lz;

        real_t x = cur_s1;
        real_t y = cur_s2 + Dy_p(0, j, k);
        real_t z = cur_s3 + Dz_p(0 , 0, k);

        real_t dx_cur = 0;

        real_t dx_lower = -x, dx_upper = (real_t)sheetSim->lx - x;
        
        int iter_cnt = 0;
        real_t x_idx, y_idx, z_idx;
        while(iter_cnt <= b_search_interation_limit)
        {
          dx_cur = (dx_lower + dx_upper) / 2.0;
          x_idx = (x + dx_cur) / dx;
          y_idx = y / dx;
          z_idx = z / dx;

          if(x+dx_cur + d1phi.getTriCubicInterpolatedValue(x_idx,y_idx,z_idx) < cur_s1)
          {
            dx_lower = dx_cur;
          }
          else
          {
            dx_upper = dx_cur;
          }
          iter_cnt++;
        }
        Dx_p(i, j, k) = dx_cur;
        max_inverse_deviation = std::max(max_inverse_deviation,
                                         fabs(x+dx_cur+d1phi.getTriCubicInterpolatedValue(x_idx,y_idx,z_idx) - cur_s1));
        max_inverse_deviation = std::max(max_inverse_deviation,
                                         fabs(y+d2phi.getTriCubicInterpolatedValue(x_idx,y_idx,z_idx) - cur_s2));
        max_inverse_deviation = std::max(max_inverse_deviation,
                                         fabs(z+d3phi.getTriCubicInterpolatedValue(x_idx,y_idx,z_idx) - cur_s3));

      }
    }
  }
  

  std::cout<<"The inversion of function brings an deviation of "<<max_inverse_deviation<<"\n";



  // setting DIFFr_a as array to store differences
  // difference equals to rho for now!

  real_t max_err = 1e99;
  real_t previous_err = 1e100;
  idx_t ns1 = Dx_a.nx;
  idx_t ns2 = Dx_a.ny;
  idx_t ns3 = Dx_a.nz;

  idx_t iter_cnt = 0;
  
  // doing iteration
  // stop when max_err increase 
  while(max_err <= previous_err)
  {
    previous_err = max_err;
    if(iter_cnt >= 2000)
    {
      std::cout<<"Jump out since iteration counts are too many!\n";
      break;
    }
    Dx_a = Dx_p;
    Dy_a = Dy_p;
    Dz_a = Dz_p;

    #pragma omp parallel for
    for(idx_t i = 0; i < rho.pts; i++)
      DIFFr_a._array[i] = 0;
    
    // depositing current rho to DIFFr_a
    sheetSim->addBSSNSource(bssnSim, tot_mass);

    // calcuating \delta \rho
#pragma omp parallel for
    for(idx_t i = 0; i < rho.pts; i++)
      rho_err._array[i] = rho._array[i] - DIFFr_a._array[i];

    max_err = rho_err.abs_max();

    iter_cnt++;
    
    std::cout<<"Itering step "<<iter_cnt
             <<" "<<"with delta rho max as "<<max_err<<"\n";


    // generating gradiant to rho_err
#pragma omp parallel for collapse(2)
    for(idx_t i = 0; i < NX; i++)
      for(idx_t j = 0; j < NY; j++)
        for(idx_t k = 0; k < NZ; k++)
        {
          idx_t idx = NP_INDEX(i,j,k);
          drhodx[idx] = derivative(i, j, k, 1, rho_err);
          drhody[idx] = derivative(i, j, k, 2, rho_err);
          drhodz[idx] = derivative(i, j, k, 3, rho_err);
        }
      
    
    
      // doing updating 
#pragma omp parallel for collapse(2)
  for(idx_t s1=0; s1<ns1; ++s1)
    for(idx_t s2=0; s2<ns2; ++s2)
      for(idx_t s3=0; s3<ns3; ++s3)
      {
        real_t x_pt = Dx_a(s1,s2,s3) + sheetSim->_S1IDXtoX0(s1);
        real_t y_pt = Dy_a(s1,s2,s3) + sheetSim->_S2IDXtoY0(s2);
        real_t z_pt = Dz_a(s1,s2,s3) + sheetSim->_S3IDXtoZ0(s3);

        real_t x_idx = x_pt/sheetSim->dx;
        real_t y_idx = y_pt/sheetSim->dy;
        real_t z_idx = z_pt/sheetSim->dz;        


        
        // should do interpolation at position (Dx, Dy, Dz)
        Dx_p(s1, s2, s3) += damping_coef * drhodx.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        Dy_p(s1, s2, s3) += damping_coef * drhody.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        Dz_p(s1, s2, s3) += damping_coef * drhodz.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
      }

  }

  // dump the result
  io_dump_3dslice(iodata, Dx_p, "Dx");
  io_dump_3dslice(iodata, Dy_p, "Dy");
  io_dump_3dslice(iodata, Dz_p, "Dz");  
  
}

void sheets_ic_sinusoid_1d_diffusion(
  BSSN *bssnSim, Sheet *sheetSim, IOData * iodata, real_t & tot_mass)
{
  iodata->log("Setting sinusoidal 3D ICs with diffusion method");
  idx_t i, j, k;


    // conformal factor
  arr_t & DIFFphi_p = *bssnSim->fields["DIFFphi_p"];
  // DIFFK is initially zero
  arr_t & DIFFK_p = *bssnSim->fields["DIFFK_p"];
  // matter sources
  arr_t & DIFFr_a = *bssnSim->fields["DIFFr_a"];

  arr_t & Dx_p = sheetSim->Dx._array_p;
  arr_t & Dx_a = sheetSim->Dx._array_a;
  
  real_t A = sheetSim->lx*sheetSim->lx*std::stod(_config("peak_amplitude", "0.0001"));

  // setting iteration stuff
  real_t precision_goal = sheetSim->lx * sheetSim->lx * std::stod(_config("precision_goal", "1e-5"));
  real_t damping_coef = sheetSim->lx * sheetSim->lx * std::stod(_config("damping_coef", "0.1"));
  
  iodata->log( "Generating ICs with peak amp. = " + stringify(A) );

  real_t rho_FRW = 3.0/PI/8.0;
  real_t K_FRW = -sqrt(24.0*PI*rho_FRW);
  iodata->log( "FRW density is " + stringify(rho_FRW) + ", total mass in simulation volume is "
    + stringify(rho_FRW*sheetSim->lx*sheetSim->ly*sheetSim->lz));

  // the conformal factor in front of metric is the solution to
  // d^2 exp(\phi) = -2*pi exp(5\phi) * \delta_rho
  // generate random mode in \phi
  // delta_rho = -(lap e^\phi)/e^(5\phi)/2pi
  real_t phix = 0;
  real_t twopi_L = 2.0*PI / sheetSim->lx;
  real_t pw2_twopi_L = twopi_L*twopi_L;

  // defining difference of \rho
  arr_t rho(NX, NY, NZ), rho_err(NX, NY, NZ);
  // variables to store gradiant of delta \rho
  arr_t drhodx(NX, NY, NZ);
  
  // grid values
  //#pragma omp parallel for
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    real_t x_frac = ((real_t) i / (real_t) NX);
    
    real_t phi = A*(sin(2.0*PI*x_frac + phix));

    real_t rho_c = rho_FRW + -exp(-4.0*phi)/PI/2.0*(
      pw2(twopi_L*A*cos(2.0*PI*x_frac + phix))
      - pw2_twopi_L*A*sin(2.0*PI*x_frac + phix)
    );
    // These aren't difference vars
    DIFFphi_p[NP_INDEX(i,j,k)] = phi;
    DIFFK_p[idx] = K_FRW;
    // set target \rho to rho temperarily
    rho[idx] = rho_c;
  }
  

  idx_t integration_points = std::stod(_config("integration_points_per_dx", "1000"));
  std::cout << "Setting initial conditions using " << integration_points<< " integration_points" << std::endl;
  real_t integration_interval = sheetSim->lx / integration_points;

  tot_mass = 0;

    // compute total mass in simulation, mass per tracer particle
  for(idx_t i=0; i<integration_points; ++i)
  {
    real_t x_frac = i/(real_t) integration_points;

    real_t phi = A*sin(2.0*PI*x_frac + phix);
    real_t rho_c = rho_FRW + -exp(-4.0*phi)/PI/2.0*(
      pw2(twopi_L*A*cos(2.0*PI*x_frac + phix))
      - pw2_twopi_L*A*sin(2.0*PI*x_frac + phix)
    );

    real_t rootdetg = std::exp(6.0*phi);

    tot_mass += rho_c * rootdetg * integration_interval * sheetSim->ly * sheetSim->lz;
  }

  real_t mass_per_tracer = tot_mass / (real_t) (sheetSim->ns1);

  std::cout << "Total mass and mass_per_tracer are " << tot_mass
    << ", " << mass_per_tracer << ".\n";

  // setting DIFFr_a as array to store differences
  // difference equals to rho for now!

  real_t max_err = 1e100;

  idx_t ns1 = Dx_a.nx;
  idx_t ns2 = Dx_a.ny;
  idx_t ns3 = Dx_a.nz;

  idx_t iter_cnt = 0;

  
  // doing iteration
  while(max_err >= precision_goal)
  {
    if(iter_cnt >= 5000)
      break;
    
    Dx_a = Dx_p;


    
    #pragma omp parallel for
    for(idx_t i = 0; i < rho.pts; i++)
      DIFFr_a._array[i] = 0;
    
    // depositing current rho to DIFFr_a
    sheetSim->addBSSNSource(bssnSim, tot_mass);

    
    // calcuating \delta \rho
#pragma omp parallel for
    for(idx_t i = 0; i < rho.pts; i++)
      rho_err._array[i] = rho._array[i] - DIFFr_a._array[i];
    
    max_err = rho_err.abs_max();

    iter_cnt++;
    
    std::cout<<"Itering step "<<iter_cnt
             <<" "<<"with delta rho max as "<<max_err<<"\n";


    // generating gradiant to rho_err
#pragma omp parallel for collapse(2)
    for(idx_t i = 0; i < NX; i++)
      for(idx_t j = 0; j < NY; j++)
        for(idx_t k = 0; k < NZ; k++)
        {
          idx_t idx = NP_INDEX(i,j,k);
          drhodx[idx] = derivative(i, j, k, 1, rho_err);
        }
      
    
      // doing updating 
#pragma omp parallel for collapse(2)
  for(idx_t s1=0; s1<ns1; ++s1)
    for(idx_t s2=0; s2<ns2; ++s2)
      for(idx_t s3=0; s3<ns3; ++s3)
      {
        real_t x_pt = Dx_a(s1,s2,s3) + sheetSim->_S1IDXtoX0(s1);

        real_t x_idx = x_pt/sheetSim->dx;

        
        // should do interpolation at position (Dx, Dy, Dz)
        Dx_p(s1, s2, s3) += damping_coef * drhodx.getTriCubicInterpolatedValue(x_idx, 0, 0);
      }

  }

  
}


} // namespace cosmo
