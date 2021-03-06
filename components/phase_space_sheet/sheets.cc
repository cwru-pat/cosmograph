#include "sheets.h"
#include "../../utils/math.h"
#include "../../cosmo_includes.h"
#include "../../cosmo_globals.h"
#include "sheets_macros.h"

namespace cosmo
{

Sheet::Sheet():
  ns1(std::stoi(_config["ns1"])),
  ns2(std::stoi(_config["ns2"])),
  ns3(std::stoi(_config["ns3"])),
  Dx(ns1, ns2, ns3, dt),
  Dy(ns1, ns2, ns3, dt),
  Dz(ns1, ns2, ns3, dt),
  vx(ns1, ns2, ns3, dt),
  vy(ns1, ns2, ns3, dt),
  vz(ns1, ns2, ns3, dt),
  d1alpha_a(NX, NY, NZ),
  d2alpha_a(NX, NY, NZ),
  d3alpha_a(NX, NY, NZ),
  d1gammai11_a(NX, NY, NZ),
  d1gammai22_a(NX, NY, NZ),
  d1gammai33_a(NX, NY, NZ),
  d1gammai12_a(NX, NY, NZ),
  d1gammai13_a(NX, NY, NZ),
  d1gammai23_a(NX, NY, NZ),
  d2gammai11_a(NX, NY, NZ),
  d2gammai22_a(NX, NY, NZ),
  d2gammai33_a(NX, NY, NZ),
  d2gammai12_a(NX, NY, NZ),
  d2gammai13_a(NX, NY, NZ),
  d2gammai23_a(NX, NY, NZ),
  d3gammai11_a(NX, NY, NZ),
  d3gammai22_a(NX, NY, NZ),
  d3gammai33_a(NX, NY, NZ),
  d3gammai12_a(NX, NY, NZ),
  d3gammai13_a(NX, NY, NZ),
  d3gammai23_a(NX, NY, NZ),
  d1beta1_a(NX, NY, NZ),
  d1beta2_a(NX, NY, NZ),
  d1beta3_a(NX, NY, NZ),
  d2beta1_a(NX, NY, NZ),
  d2beta2_a(NX, NY, NZ),
  d2beta3_a(NX, NY, NZ),
  d3beta1_a(NX, NY, NZ),
  d3beta2_a(NX, NY, NZ),
  d3beta3_a(NX, NY, NZ)
{
  dx = H_LEN_FRAC / (real_t) COSMO_N;
  dy = dx;
  dz = dx;

  lx = NX * dx;
  ly = NY * dx;
  lz = NZ * dx;

  std::cout << "Initiaizing sheet class with lx,ly,lz = " << lx << "," << ly << "," << lz
    << ", dx,dy,dz = " << dx << "," << dy << "," << dz << std::endl;

  carrier_count_scheme = static_cast<carrierCountScheme> (std::stoi(_config("carrier_count_scheme","1")));
  deposit = static_cast<depositScheme> (std::stoi(_config("deposit_scheme","1")));

  follow_null_geodesics = !!std::stoi(_config("follow_null_geodesics", "0"));
  rescale_sheet = std::stod(_config("rescale_sheet", "1.0"));
  if(rescale_sheet == 1.0) rescale_sheet = 0.0;
  ray_bundle_epsilon = std::stod(_config("ray_bundle_epsilon","1.0")) / (real_t POINTS);
  det_g_obs = 0.0;

  carriers_per_dx = std::stoi(_config("carriers_per_dx","1"));
  carriers_per_dy = std::stoi(_config("carriers_per_dy","1"));
  carriers_per_dz = std::stoi(_config("carriers_per_dz","1"));
}


Sheet::~Sheet()
{
  // Anything to do?
}

void Sheet::setDt(real_t dt)
{
  // set dt for all RK4 register fields
  Dx.setDt(dt);
  Dy.setDt(dt);
  Dz.setDt(dt);
  vx.setDt(dt);
  vy.setDt(dt);
  vz.setDt(dt);
}

void Sheet::_MassDeposit(real_t weight, real_t x_idx, real_t y_idx,
  real_t z_idx, arr_t &rho)
{
  switch (deposit)
  {
    case PCS:
      _PCSDeposit(weight, x_idx, y_idx, z_idx, rho);
      break;
    case CINT:
      _CINTDeposit(weight, x_idx, y_idx, z_idx, rho);
      break;
    case CIC:
    default:
      _CICDeposit(weight, x_idx, y_idx, z_idx, rho);
      break;
  }
}

void Sheet::_CICDeposit(real_t weight, real_t x_idx, real_t y_idx,
  real_t z_idx, arr_t &rho)
{
  idx_t ix, iy, iz,
        ixp, iyp, izp;

  real_t x_f, y_f, z_f,
         x_h, y_h, z_h;

  // gridpoint index "left" of x_idx
  ix = (x_idx < 0 ? (idx_t) x_idx - 1 : (idx_t) x_idx );
  ixp = ix + 1;
  x_f = x_idx - (real_t) ix;
  x_h = 1.0 - x_f;

  iy = (y_idx < 0 ? (idx_t) y_idx - 1 : (idx_t) y_idx );
  iyp = iy + 1;
  y_f = y_idx - (real_t) iy;
  y_h = 1.0 - y_f;

  iz = (z_idx < 0 ? (idx_t) z_idx - 1 : (idx_t) z_idx );
  izp = iz + 1;
  z_f = z_idx - (real_t) iz;
  z_h = 1.0 - z_f;

#pragma omp atomic    
  rho(ix, iy, iz) += x_h*y_h*z_h*weight;
#pragma omp atomic
  rho(ix, iy, izp) += x_h*y_h*z_f*weight;
#pragma omp atomic
  rho(ix, iyp, iz) += x_h*y_f*z_h*weight;
#pragma omp atomic
  rho(ix, iyp, izp) += x_h*y_f*z_f*weight;

#pragma omp atomic
  rho(ixp, iy, iz) += x_f*y_h*z_h*weight;
#pragma omp atomic
  rho(ixp, iy, izp) += x_f*y_h*z_f*weight;
#pragma omp atomic
  rho(ixp, iyp, iz) += x_f*y_f*z_h*weight;
#pragma omp atomic
  rho(ixp, iyp, izp) += x_f*y_f*z_f*weight;

}

/**
 * @brief      Piecewise cubic spline deposition
 */
void Sheet::_PCSDeposit(real_t weight, real_t x_idx, real_t y_idx,
  real_t z_idx, arr_t &rho)
{

  if(ns2 == 1 && NY == 1 && ns3 == 1 && NZ == 1)
  {
    idx_t ix = (x_idx < 0 ? (idx_t) x_idx - 1 : (idx_t) x_idx );

    real_t norm = 0.0;
    for(idx_t i=-1; i<=2; ++i)
    {
      real_t s = std::abs( ix+i-x_idx );
      
      if(s<1.0)
      {
        norm += (4.0 - 6.0*s*s + 3.0*std::pow(std::abs(s), 3))/6.0;
      }
      else if(s<2.0 && s>=1.0)
      {
        norm += std::pow(2.0 - std::abs(s), 3)/6.0;
      }
    }

    real_t pcs;
    for(idx_t i=-1; i<=2; ++i)
    {
      real_t s = std::abs( ix+i-x_idx );
          
      if(s<1.0)
      {
        pcs = (4.0 - 6.0*s*s + 3.0*std::pow(std::abs(s), 3))/6.0;
#pragma omp atomic
        rho(ix+i, 0, 0) += pcs*weight/norm;
      }
      else if(s<2.0 && s>=1.0)
      {
        pcs = std::pow(2.0 - std::abs(s), 3)/6.0;
#pragma omp atomic
        rho(ix+i, 0, 0) += pcs*weight/norm;
      }
    }

  }
  else
  {
    idx_t ix = (x_idx < 0 ? (idx_t) x_idx - 1 : (idx_t) x_idx );
    idx_t iy = (y_idx < 0 ? (idx_t) y_idx - 1 : (idx_t) y_idx );
    idx_t iz = (z_idx < 0 ? (idx_t) z_idx - 1 : (idx_t) z_idx );

    real_t norm = 0.0;
    for(idx_t i=-1; i<=2; ++i)
      for(idx_t j=-1; j<=2; ++j)
        for(idx_t k=-1; k<=2; ++k)
        {
          real_t s = std::sqrt( std::pow(ix+i-x_idx, 2)
            + std::pow(iy+j-y_idx, 2)
            + std::pow(iz+k-z_idx, 2) );
          
          if(s<1.0)
          {
            norm += (4.0 - 6.0*s*s + 3.0*std::pow(std::abs(s), 3))/6.0;
          }
          else if(s<2.0 && s>=1.0)
          {
            norm += std::pow(2.0 - std::abs(s), 3)/6.0;
          }
        }

    real_t pcs;
    for(idx_t i=-1; i<=2; ++i)
      for(idx_t j=-1; j<=2; ++j)
        for(idx_t k=-1; k<=2; ++k)
        {
          real_t s = std::sqrt( std::pow(ix+i-x_idx, 2)
            + std::pow(iy+j-y_idx, 2)
            + std::pow(iz+k-z_idx, 2) );
          
          if(s<1.0)
          {
            pcs = (4.0 - 6.0*s*s + 3.0*std::pow(std::abs(s), 3))/6.0;
#pragma omp atomic
            rho(ix+i, iy+j, iz+k) += pcs*weight/norm;
          }
          else if(s<2.0 && s>=1.0)
          {
            pcs = std::pow(2.0 - std::abs(s), 3)/6.0;
#pragma omp atomic
            rho(ix+i, iy+j, iz+k) += pcs*weight/norm;
          }
        }
  }

}


/**
 * @brief      Deposition based on cubic interoplation (CINT)
 */
void Sheet::_CINTDeposit(real_t weight, real_t x_idx, real_t y_idx,
  real_t z_idx, arr_t &rho)
{
  // gridpoint index "left" of x_idx
  idx_t ix = (x_idx < 0 ? (idx_t) x_idx - 1 : (idx_t) x_idx );
  real_t x_f = x_idx - (real_t) ix;

  idx_t iy = (y_idx < 0 ? (idx_t) y_idx - 1 : (idx_t) y_idx );
  real_t y_f = y_idx - (real_t) iy;

  idx_t iz = (z_idx < 0 ? (idx_t) z_idx - 1 : (idx_t) z_idx );
  real_t z_f = z_idx - (real_t) iz;

  real_t x_wts[4] = { x_f*x_f*(2.0-x_f)-x_f, x_f*x_f*(3.0*x_f-5.0)+2.0,
         x_f*x_f*(4.0-3.0*x_f)+x_f, x_f*x_f*(x_f-1.0) };
  real_t y_wts[4] = { y_f*y_f*(2.0-y_f)-y_f, y_f*y_f*(3.0*y_f-5.0)+2.0,
         y_f*y_f*(4.0-3.0*y_f)+y_f, y_f*y_f*(y_f-1.0) };
  real_t z_wts[4] = { z_f*z_f*(2.0-z_f)-z_f, z_f*z_f*(3.0*z_f-5.0)+2.0,
         z_f*z_f*(4.0-3.0*z_f)+z_f, z_f*z_f*(z_f-1.0) };

  for(idx_t i=-1; i<=2; ++i)
    for(idx_t j=-1; j<=2; ++j)
      for(idx_t k=-1; k<=2; ++k)
      {        
#pragma omp atomic
          rho(ix+i, iy+j, iz+k) += weight*x_wts[i+1]*y_wts[j+1]*z_wts[k+1]/2.0/2.0/2.0;
      }

}


/**
* Get Min/max x/y/z coordinates at voxel corners
*/
real_t Sheet::_getXRangeInSVoxel(register_t & DX, idx_t s1_idx, idx_t s2_idx,
  idx_t s3_idx, real_t X0_lower, real_t X0_upper)
{
  real_t X_vals_on_bounds[8] = {
    DX(s1_idx,s2_idx,s3_idx) + X0_lower,
    DX(s1_idx+1,s2_idx,s3_idx) + X0_upper,
    DX(s1_idx,s2_idx+1,s3_idx) + X0_lower,
    DX(s1_idx,s2_idx,s3_idx+1) + X0_lower,
    DX(s1_idx+1,s2_idx+1,s3_idx) + X0_upper,
    DX(s1_idx+1,s2_idx,s3_idx+1) + X0_upper,
    DX(s1_idx,s2_idx+1,s3_idx+1) + X0_lower,
    DX(s1_idx+1,s2_idx+1,s3_idx+1) + X0_upper
  };

  real_t X_min = std::numeric_limits<real_t>::max();
  real_t X_max = std::numeric_limits<real_t>::min();
  for(int i=0; i<8; ++i)
  {
    if(X_vals_on_bounds[i] > X_max) X_max = X_vals_on_bounds[i];
    if(X_vals_on_bounds[i] < X_min) X_min = X_vals_on_bounds[i];
  }

  return X_max - X_min;
}


/**
 * Compute conribution to rho(x) from data in a phase-space
 * Sheet voxel and add to rho(x) grid. Do so via 1501.01959 mass
 * deposition scheme.
 * TODO: improve; consider higher-order or analytic versions of this
 */
void Sheet::addBSSNSource(BSSN *bssn, real_t tot_mass)
{
  _timer["_pushsheetToStressTensor"].start();

  arr_t & DIFFr_a = *bssn->fields["DIFFr_a"];
  arr_t & DIFFS_a = *bssn->fields["DIFFS_a"];
  arr_t & S1_a = *bssn->fields["S1_a"];
  arr_t & S2_a = *bssn->fields["S2_a"];
  arr_t & S3_a = *bssn->fields["S3_a"];
  arr_t & STF11_a = *bssn->fields["STF11_a"];
  arr_t & STF12_a = *bssn->fields["STF12_a"];
  arr_t & STF13_a = *bssn->fields["STF13_a"];
  arr_t & STF22_a = *bssn->fields["STF22_a"];
  arr_t & STF23_a = *bssn->fields["STF23_a"];
  arr_t & STF33_a = *bssn->fields["STF33_a"];

  arr_t & DIFFphi_a = *bssn->fields["DIFFphi_a"];
  arr_t & DIFFgamma11_a = *bssn->fields["DIFFgamma11_a"];
  arr_t & DIFFgamma22_a = *bssn->fields["DIFFgamma22_a"];
  arr_t & DIFFgamma33_a = *bssn->fields["DIFFgamma33_a"];

  arr_t & DIFFgamma12_a = *bssn->fields["DIFFgamma12_a"];
  arr_t & DIFFgamma13_a = *bssn->fields["DIFFgamma13_a"];
  arr_t & DIFFgamma23_a = *bssn->fields["DIFFgamma23_a"];
  
  idx_t num_x_carriers, num_y_carriers, num_z_carriers;

  //    int i, j, k;
# pragma omp parallel for collapse(2)
  for(idx_t s1=0; s1<ns1; ++s1)
    for(idx_t s2=0; s2<ns2; ++s2)
      for(idx_t s3=0; s3<ns3; ++s3)
      {
        switch(carrier_count_scheme)
        {
        case per_dx:
          num_x_carriers = carriers_per_dx == 0 ? 1 : carriers_per_dx*(
            (idx_t) (0.5 + _getXRangeInSVoxel(Dx, s1, s2, s3, _S1IDXtoX0(s1), _S1IDXtoX0(s1+1)) / dx ) );
          num_y_carriers = carriers_per_dy == 0 ? 1 : carriers_per_dy*(
            (idx_t) (0.5 + _getXRangeInSVoxel(Dy, s1, s2, s3, _S2IDXtoY0(s2), _S2IDXtoY0(s2+1)) / dy ) );
          num_z_carriers = carriers_per_dz == 0 ? 1 : carriers_per_dz*(
            (idx_t) (0.5 + _getXRangeInSVoxel(Dz, s1, s2, s3, _S3IDXtoZ0(s3), _S3IDXtoZ0(s3+1)) / dz ) );
          break;
        case per_ds:
        default:
          num_x_carriers = carriers_per_dx;
          num_y_carriers = carriers_per_dy;
          num_z_carriers = carriers_per_dz;
          break;
        }


        if(num_x_carriers <= 0) num_x_carriers = 1;
        if(num_y_carriers <= 0) num_y_carriers = 1;
        if(num_z_carriers <= 0) num_z_carriers = 1;
        
        idx_t num_carriers = num_x_carriers*num_y_carriers*num_z_carriers;

        idx_t i, j, k;

        real_t f_Dx[64], f_Dy[64], f_Dz[64];
        real_t a_Dx[64], a_Dy[64], a_Dz[64];

        if(num_x_carriers != 1 || num_y_carriers != 1 || num_z_carriers != 1)
        {
        
          for(i=0; i<4; ++i)
            for(j=0; j<4; ++j)
              for(k=0; k<4; ++k)
              {
                f_Dx[i*16 + j*4 + k] = Dx(s1-1+i, s2-1+j, s3-1+k);
                f_Dy[i*16 + j*4 + k] = Dy(s1-1+i, s2-1+j, s3-1+k);
                f_Dz[i*16 + j*4 + k] = Dz(s1-1+i, s2-1+j, s3-1+k);
              }
          compute_tricubic_coeffs(a_Dx, f_Dx);
          compute_tricubic_coeffs(a_Dy, f_Dy);
          compute_tricubic_coeffs(a_Dz, f_Dz);
        }
        // distribute mass from all carriers
        for(i=0; i<num_x_carriers; ++i)
          for(j=0; j<num_y_carriers; ++j)
            for(k=0; k<num_z_carriers; ++k)
            {

              real_t carrier_s1d = (real_t) i / (real_t) num_x_carriers;
              real_t carrier_s2d = (real_t) j / (real_t) num_y_carriers;
              real_t carrier_s3d = (real_t) k / (real_t) num_z_carriers;

              real_t carrier_s1 = s1 + carrier_s1d;
              real_t carrier_s2 = s2 + carrier_s2d;
              real_t carrier_s3 = s3 + carrier_s3d;

              real_t carrier_x_idx, carrier_y_idx, carrier_z_idx;
              
              if(num_x_carriers == 1 && num_y_carriers == 1 && num_z_carriers == 1)
              {
                // special case if only one carrier
                carrier_x_idx = ( _S1IDXtoX0(s1) + Dx(s1, s2, s3) ) / dx;
                carrier_y_idx = ( _S2IDXtoY0(s2) + Dy(s1, s2, s3) ) / dy;
                carrier_z_idx = ( _S3IDXtoZ0(s3) + Dz(s1, s2, s3) ) / dz;
              }
              else
              {
                carrier_x_idx = ( _S1IDXtoX0(carrier_s1)
                                         + evaluate_interpolation(a_Dx, carrier_s1d, carrier_s2d, carrier_s3d) ) / dx;
                carrier_y_idx = ( _S2IDXtoY0(carrier_s2)
                                         + evaluate_interpolation(a_Dy, carrier_s1d, carrier_s2d, carrier_s3d) ) / dy;
                carrier_z_idx = ( _S3IDXtoZ0(carrier_s3)
                                         + evaluate_interpolation(a_Dz, carrier_s1d, carrier_s2d, carrier_s3d) ) / dz;
              }

              real_t u1 = vx._array_a.getTriCubicInterpolatedValue(carrier_s1, carrier_s2, carrier_s3);
              real_t u2 = vy._array_a.getTriCubicInterpolatedValue(carrier_s1, carrier_s2, carrier_s3);
              real_t u3 = vz._array_a.getTriCubicInterpolatedValue(carrier_s1, carrier_s2, carrier_s3);

              real_t phi = DIFFphi_a.getTriCubicInterpolatedValue(
                carrier_x_idx, carrier_y_idx, carrier_z_idx);

              real_t DIFFgamma11 = DIFFgamma11_a.getTriCubicInterpolatedValue(
                carrier_x_idx, carrier_y_idx, carrier_z_idx); 
              real_t DIFFgamma22 = DIFFgamma22_a.getTriCubicInterpolatedValue(
                carrier_x_idx, carrier_y_idx, carrier_z_idx) ; 
              real_t DIFFgamma33 = DIFFgamma33_a.getTriCubicInterpolatedValue(
                carrier_x_idx, carrier_y_idx, carrier_z_idx) ; 

              real_t DIFFgamma12 = DIFFgamma12_a.getTriCubicInterpolatedValue(
                carrier_x_idx, carrier_y_idx, carrier_z_idx); 
              real_t DIFFgamma23 = DIFFgamma23_a.getTriCubicInterpolatedValue(
                carrier_x_idx, carrier_y_idx, carrier_z_idx); 
              real_t DIFFgamma13 = DIFFgamma13_a.getTriCubicInterpolatedValue(
                carrier_x_idx, carrier_y_idx, carrier_z_idx); 

              real_t gammai11 = std::exp(-4.0*phi)*(1.0 + DIFFgamma22 + DIFFgamma33 - pw2(DIFFgamma23) + DIFFgamma22*DIFFgamma33);
              real_t gammai22 = std::exp(-4.0*phi)*(1.0 + DIFFgamma11 + DIFFgamma33 - pw2(DIFFgamma13) + DIFFgamma11*DIFFgamma33);
              real_t gammai33 = std::exp(-4.0*phi)*(1.0 + DIFFgamma11 + DIFFgamma22 - pw2(DIFFgamma12) + DIFFgamma11*DIFFgamma22);
              real_t gammai12 = std::exp(-4.0*phi)*(DIFFgamma13*DIFFgamma23 - DIFFgamma12*(1.0 + DIFFgamma33));
              real_t gammai13 = std::exp(-4.0*phi)*(DIFFgamma12*DIFFgamma23 - DIFFgamma13*(1.0 + DIFFgamma22));
              real_t gammai23 = std::exp(-4.0*phi)*(DIFFgamma12*DIFFgamma13 - DIFFgamma23*(1.0 + DIFFgamma11));

              real_t rootdetg = std::exp(6.0*phi);
        
              real_t W = 0.0;
              if(follow_null_geodesics)
              {
                W = std::sqrt(
                  gammai11*u1*u1 + gammai22*u2*u2 + gammai33*u3*u3
                  + 2.0*( gammai12*u1*u2 + gammai13*u1*u3 + gammai23*u2*u3 )
                );
              }
              else
              {
                W = std::sqrt( 1.0 +
                  gammai11*u1*u1 + gammai22*u2*u2 + gammai33*u3*u3
                  + 2.0*( gammai12*u1*u2 + gammai13*u1*u3 + gammai23*u2*u3 )
                );
              }
              
              real_t mass = tot_mass / (real_t) num_carriers / ns1/ns2/ns3;
              
              real_t MnA = mass / W / dx / dy / dz / rootdetg;

              real_t rho = MnA*W*W;
              real_t S = rho - MnA;
              real_t S1 = MnA*W*u1;
              real_t S2 = MnA*W*u2;
              real_t S3 = MnA*W*u3;

              // STF is not trace free yet
              real_t STF11 = (MnA*u1*u1);
              real_t STF12 = (MnA*u1*u2);
              real_t STF13 = (MnA*u1*u3);
              real_t STF22 = (MnA*u2*u2);
              real_t STF23 = (MnA*u2*u3);
              real_t STF33 = (MnA*u3*u3);

              _MassDeposit(rho, carrier_x_idx, carrier_y_idx, carrier_z_idx, DIFFr_a);

              _MassDeposit(S, carrier_x_idx, carrier_y_idx, carrier_z_idx, DIFFS_a);
              
              _MassDeposit(S1, carrier_x_idx, carrier_y_idx, carrier_z_idx, S1_a);
              _MassDeposit(S2, carrier_x_idx, carrier_y_idx, carrier_z_idx, S2_a);
              _MassDeposit(S3, carrier_x_idx, carrier_y_idx, carrier_z_idx, S3_a);

              _MassDeposit(STF11, carrier_x_idx, carrier_y_idx, carrier_z_idx, STF11_a);
              _MassDeposit(STF22, carrier_x_idx, carrier_y_idx, carrier_z_idx, STF22_a);
              _MassDeposit(STF33, carrier_x_idx, carrier_y_idx, carrier_z_idx, STF33_a);
              _MassDeposit(STF12, carrier_x_idx, carrier_y_idx, carrier_z_idx, STF12_a);
              _MassDeposit(STF13, carrier_x_idx, carrier_y_idx, carrier_z_idx, STF13_a);
              _MassDeposit(STF23, carrier_x_idx, carrier_y_idx, carrier_z_idx, STF23_a);
            }
        
      }

  idx_t i, j, k;
# pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i, j, k)
  {
    BSSNData bd = {0};
    bssn->set_bd_values(i, j, k, &bd);
    bssn->enforceTFSIJ(&bd);
  }


  _timer["_pushsheetToStressTensor"].stop();
}

void Sheet::rescaleFieldPerturbations(arr_t & field, real_t multiplier)
{
  idx_t i, j, k;
  real_t avg = average(field);
  LOOP3(i,j,k)
  {
    field[NP_INDEX(i,j,k)] = avg + multiplier*( field[NP_INDEX(i,j,k)] - avg );
  }
}

void Sheet::rescaleVelocityPerturbations(arr_t & ux, arr_t & uy, arr_t & uz, real_t multiplier)
{

  // velocities are perturbations around \vec{u} = \hat{r}
  // Can scale coordinates relative to average magnitude
  real_t avg_u = 0.0;
  idx_t s2 = 0, s3 = 0;
# pragma omp parallel for reduction(+:avg_u)
  for(idx_t s1=0; s1<ns1; ++s1)
  {
    avg_u += std::sqrt( pw2(ux(s1,s2,s3)) + pw2(uy(s1,s2,s3)) + pw2(uz(s1,s2,s3)) );
  }
  avg_u /= ns1;

  std::ifstream vecFile(_config["healpix_vecs_file"]);
  idx_t r = 0;
  while (!vecFile.eof() && r < ns1)
  {
    idx_t s1 = r;
    // "background" position
    real_t ux_bar, uy_bar, uz_bar;
    vecFile >> ux_bar;
    vecFile >> uy_bar;
    vecFile >> uz_bar;
    real_t mag_u = std::sqrt( pw2(ux_bar) + pw2(uy_bar) + pw2(uz_bar) );
    ux_bar *= avg_u/mag_u;
    uy_bar *= avg_u/mag_u;
    uz_bar *= avg_u/mag_u;

    for(idx_t s2=0; s2<3; ++s2)
    {
      idx_t s3 = 0.0;
      ux(s1,s2,s3) = ux_bar + multiplier*( ux(s1,s2,s3) - ux_bar );
      uy(s1,s2,s3) = uy_bar + multiplier*( uy(s1,s2,s3) - uy_bar );
      uz(s1,s2,s3) = uz_bar + multiplier*( uz(s1,s2,s3) - uz_bar );
    }
    r++;
  }
  vecFile.close();

}


void Sheet::rescalePositionPerturbations(arr_t & dx, arr_t & dy, arr_t & dz, real_t multiplier)
{
  // coordinate positions are "perturbations" around \vec{x} = \vec{r}
  // Can scale coordinates relative to average magnitude
  real_t avg_r = 0.0;
  idx_t s2 = 0, s3 = 0;
# pragma omp parallel for reduction(+:avg_r)
  for(idx_t s1=0; s1<ns1; ++s1)
  {
    real_t x_pt = dx(s1,s2,s3) + _S1IDXtoX0(s1);
    real_t y_pt = dy(s1,s2,s3) + _S2IDXtoY0(s2);
    real_t z_pt = dz(s1,s2,s3) + _S3IDXtoZ0(s3);
    avg_r += std::sqrt( pw2(x_pt) + pw2(y_pt) + pw2(z_pt) );
  }
  avg_r /= ns1;

  std::ifstream vecFile(_config["healpix_vecs_file"]);
  idx_t r = 0;
  while (!vecFile.eof() && r < ns1)
  {
    idx_t s1 = r;
    // "background" position
    real_t x_bar, y_bar, z_bar;
    vecFile >> x_bar;
    vecFile >> y_bar;
    vecFile >> z_bar;
    real_t mag_r = -std::sqrt( pw2(x_bar) + pw2(y_bar) + pw2(z_bar) );
    x_bar *= avg_r/mag_r;
    y_bar *= avg_r/mag_r;
    z_bar *= avg_r/mag_r;

    for(idx_t s2=0; s2<3; ++s2)
    {
      idx_t s3 = 0.0;
      Dx(s1,s2,s3) = x_bar + multiplier*( dx(s1,s2,s3) + _S1IDXtoX0(s1) - x_bar ) - _S1IDXtoX0(s1);
      Dy(s1,s2,s3) = y_bar + multiplier*( dy(s1,s2,s3) + _S2IDXtoY0(s2) - y_bar ) - _S2IDXtoY0(s2);
      Dz(s1,s2,s3) = z_bar + multiplier*( dz(s1,s2,s3) + _S3IDXtoZ0(s3) - z_bar ) - _S3IDXtoZ0(s3);
    }
    r++;
  }
  vecFile.close();

}


void Sheet::rescaleAllFieldPerturbations(BSSN *bssn, real_t multiplier)
{
  // rescale metric perturbations that get used in the sheet evolution
  rescaleFieldPerturbations(*bssn->fields["DIFFphi_a"], multiplier);
  rescaleFieldPerturbations(*bssn->fields["DIFFgamma11_a"], multiplier);
  rescaleFieldPerturbations(*bssn->fields["DIFFgamma22_a"], multiplier);
  rescaleFieldPerturbations(*bssn->fields["DIFFgamma33_a"], multiplier);
  rescaleFieldPerturbations(*bssn->fields["DIFFgamma12_a"], multiplier);
  rescaleFieldPerturbations(*bssn->fields["DIFFgamma13_a"], multiplier);
  rescaleFieldPerturbations(*bssn->fields["DIFFgamma23_a"], multiplier);
  rescaleFieldPerturbations(*bssn->fields["DIFFalpha_a"], multiplier);
#if USE_BSSN_SHIFT
  rescaleFieldPerturbations(*bssn->fields["beta1_a"], multiplier);
  rescaleFieldPerturbations(*bssn->fields["beta2_a"], multiplier);
  rescaleFieldPerturbations(*bssn->fields["beta3_a"], multiplier);
#endif

  rescalePositionPerturbations(Dx._array_a, Dy._array_a, Dz._array_a, multiplier);
  rescalePositionPerturbations(Dx._array_c, Dy._array_c, Dz._array_c, multiplier);

  rescaleVelocityPerturbations(vx._array_a, vy._array_a, vz._array_a, multiplier);
  rescaleVelocityPerturbations(vx._array_c, vy._array_c, vz._array_c, multiplier);
}


void Sheet::RKStep(BSSN *bssn)
{
  if(rescale_sheet) rescaleAllFieldPerturbations(bssn, rescale_sheet);

  idx_t i, j, k;

  arr_t & DIFFphi_a = *bssn->fields["DIFFphi_a"];
  arr_t & DIFFgamma11_a = *bssn->fields["DIFFgamma11_a"];
  arr_t & DIFFgamma22_a = *bssn->fields["DIFFgamma22_a"];
  arr_t & DIFFgamma33_a = *bssn->fields["DIFFgamma33_a"];

  arr_t & DIFFgamma12_a = *bssn->fields["DIFFgamma12_a"];
  arr_t & DIFFgamma13_a = *bssn->fields["DIFFgamma13_a"];
  arr_t & DIFFgamma23_a = *bssn->fields["DIFFgamma23_a"];

  arr_t & DIFFalpha_a = *bssn->fields["DIFFalpha_a"];

#if USE_BSSN_SHIFT
  arr_t & beta1_a = *bssn->fields["beta1_a"];
  arr_t & beta2_a = *bssn->fields["beta2_a"];
  arr_t & beta3_a = *bssn->fields["beta3_a"];
#endif


# pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i, j, k)
  {
#if USE_BSSN_SHIFT
    d1beta1_a(i, j, k) = derivative(i, j, k, 1, beta1_a);
    d1beta2_a(i, j, k) = derivative(i, j, k, 1, beta2_a);
    d1beta3_a(i, j, k) = derivative(i, j, k, 1, beta3_a);
    d2beta1_a(i, j, k) = derivative(i, j, k, 2, beta1_a);
    d2beta2_a(i, j, k) = derivative(i, j, k, 2, beta2_a);
    d2beta3_a(i, j, k) = derivative(i, j, k, 2, beta3_a);
    d3beta1_a(i, j, k) = derivative(i, j, k, 3, beta1_a);
    d3beta2_a(i, j, k) = derivative(i, j, k, 3, beta2_a);
    d3beta3_a(i, j, k) = derivative(i, j, k, 3, beta3_a);
#endif

    d1alpha_a(i, j, k) = derivative(i, j, k, 1, DIFFalpha_a);
    d2alpha_a(i, j, k) = derivative(i, j, k, 2, DIFFalpha_a);
    d3alpha_a(i, j, k) = derivative(i, j, k, 3, DIFFalpha_a);


    real_t phi = DIFFphi_a(i, j, k);

    real_t DIFFgamma11 = DIFFgamma11_a(i, j, k);
    real_t DIFFgamma22 = DIFFgamma22_a(i, j, k);
    real_t DIFFgamma33 = DIFFgamma33_a(i, j, k);

    real_t DIFFgamma12 = DIFFgamma12_a(i, j, k);
    real_t DIFFgamma23 = DIFFgamma23_a(i, j, k);
    real_t DIFFgamma13 = DIFFgamma13_a(i, j, k);
    
    real_t gammai11 = std::exp(-4.0*phi)*(1.0 + DIFFgamma22 + DIFFgamma33 - pw2(DIFFgamma23) + DIFFgamma22*DIFFgamma33);
    real_t gammai22 = std::exp(-4.0*phi)*(1.0 + DIFFgamma11 + DIFFgamma33 - pw2(DIFFgamma13) + DIFFgamma11*DIFFgamma33);
    real_t gammai33 = std::exp(-4.0*phi)*(1.0 + DIFFgamma11 + DIFFgamma22 - pw2(DIFFgamma12) + DIFFgamma11*DIFFgamma22);
    real_t gammai12 = std::exp(-4.0*phi)*(DIFFgamma13*DIFFgamma23 - DIFFgamma12*(1.0 + DIFFgamma33));
    real_t gammai13 = std::exp(-4.0*phi)*(DIFFgamma12*DIFFgamma23 - DIFFgamma13*(1.0 + DIFFgamma22));
    real_t gammai23 = std::exp(-4.0*phi)*(DIFFgamma12*DIFFgamma13 - DIFFgamma23*(1.0 + DIFFgamma11));

    SET_GAMMAI_DER(1);
    SET_GAMMAI_DER(2);
    SET_GAMMAI_DER(3);
  }

#pragma omp parallel for default(shared) private(i, j, k)
  for(i=0; i<ns1; ++i)
    for(j=0; j<ns2; ++j)
      for(k=0; k<ns3; ++k)
      {
        real_t x_pt = Dx._a(i,j,k) + _S1IDXtoX0(i);
        real_t y_pt = Dy._a(i,j,k) + _S2IDXtoY0(j);
        real_t z_pt = Dz._a(i,j,k) + _S3IDXtoZ0(k);

        real_t x_idx = x_pt/dx;
        real_t y_idx = y_pt/dy;
        real_t z_idx = z_pt/dz;
 
        real_t u1 = vx._a(i, j, k);
        real_t u2 = vy._a(i, j, k);
        real_t u3 = vz._a(i, j, k);

        
        real_t phi = DIFFphi_a.getTriCubicInterpolatedValue(
          x_idx, y_idx, z_idx);

        real_t DIFFalpha = DIFFalpha_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);

#if USE_BSSN_SHIFT
        real_t beta1 = beta1_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t beta2 = beta2_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t beta3 = beta3_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
#else
        real_t beta1 = 0, beta2 = 0, beta3 = 0;
#endif
        real_t DIFFgamma11 = DIFFgamma11_a.getTriCubicInterpolatedValue(
          x_idx, y_idx, z_idx);
        real_t DIFFgamma22 = DIFFgamma22_a.getTriCubicInterpolatedValue(
          x_idx, y_idx, z_idx);
        real_t DIFFgamma33 = DIFFgamma33_a.getTriCubicInterpolatedValue(
          x_idx, y_idx, z_idx);

        real_t DIFFgamma12 = DIFFgamma12_a.getTriCubicInterpolatedValue(
          x_idx, y_idx, z_idx); 
        real_t DIFFgamma23 = DIFFgamma23_a.getTriCubicInterpolatedValue(
          x_idx, y_idx, z_idx); 
        real_t DIFFgamma13 = DIFFgamma13_a.getTriCubicInterpolatedValue(
          x_idx, y_idx, z_idx); 

        real_t gammai11 = std::exp(-4.0*phi)*(1.0 + DIFFgamma22 + DIFFgamma33 - pw2(DIFFgamma23) + DIFFgamma22*DIFFgamma33);
        real_t gammai22 = std::exp(-4.0*phi)*(1.0 + DIFFgamma11 + DIFFgamma33 - pw2(DIFFgamma13) + DIFFgamma11*DIFFgamma33);
        real_t gammai33 = std::exp(-4.0*phi)*(1.0 + DIFFgamma11 + DIFFgamma22 - pw2(DIFFgamma12) + DIFFgamma11*DIFFgamma22);
        real_t gammai12 = std::exp(-4.0*phi)*(DIFFgamma13*DIFFgamma23 - DIFFgamma12*(1.0 + DIFFgamma33));
        real_t gammai13 = std::exp(-4.0*phi)*(DIFFgamma12*DIFFgamma23 - DIFFgamma13*(1.0 + DIFFgamma22));
        real_t gammai23 = std::exp(-4.0*phi)*(DIFFgamma12*DIFFgamma13 - DIFFgamma23*(1.0 + DIFFgamma11));

        real_t d1alpha = d1alpha_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d2alpha = d2alpha_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d3alpha = d3alpha_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);

#if USE_BSSN_SHIFT
        real_t d1beta1 = d1beta1_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d1beta2 = d1beta2_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d1beta3 = d1beta3_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);

        real_t d2beta1 = d2beta1_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d2beta2 = d2beta2_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d2beta3 = d2beta3_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);

        real_t d3beta1 = d3beta1_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d3beta2 = d3beta2_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d3beta3 = d3beta3_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
#else
        real_t d1beta1 = 0, d1beta2 = 0, d1beta3 = 0;
        real_t d2beta1 = 0, d2beta2 = 0, d2beta3 = 0;
        real_t d3beta1 = 0, d3beta2 = 0, d3beta3 = 0;
#endif

        real_t d1gammai11 = d1gammai11_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d1gammai22 = d1gammai22_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d1gammai33 = d1gammai33_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d1gammai12 = d1gammai12_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d1gammai13 = d1gammai13_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d1gammai23 = d1gammai23_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);

// optimize this for 1d
#if NY == 1 && NZ == 1
        real_t d2gammai11 = 0.0;
        real_t d2gammai22 = 0.0;
        real_t d2gammai33 = 0.0;
        real_t d2gammai12 = 0.0;
        real_t d2gammai13 = 0.0;
        real_t d2gammai23 = 0.0;

        real_t d3gammai11 = 0.0;
        real_t d3gammai22 = 0.0;
        real_t d3gammai33 = 0.0;
        real_t d3gammai12 = 0.0;
        real_t d3gammai13 = 0.0;
        real_t d3gammai23 = 0.0;
#else
        real_t d2gammai11 = d2gammai11_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d2gammai22 = d2gammai22_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d2gammai33 = d2gammai33_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d2gammai12 = d2gammai12_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d2gammai13 = d2gammai13_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d2gammai23 = d2gammai23_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);

        real_t d3gammai11 = d3gammai11_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d3gammai22 = d3gammai22_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d3gammai33 = d3gammai33_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d3gammai12 = d3gammai12_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d3gammai13 = d3gammai13_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
        real_t d3gammai23 = d3gammai23_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
#endif

        real_t W = 0.0;
        if(follow_null_geodesics)
        {
          W = std::sqrt(
            gammai11 * u1 * u1 + gammai22 * u2 * u2 + gammai33 * u3 * u3
            + 2.0 * (gammai12 * u1 * u2 + gammai13 * u1 * u3 + gammai23 * u2 * u3 )
          );
        }
        else
        {
          W = std::sqrt( 1.0 +
            gammai11 * u1 * u1 + gammai22 * u2 * u2 + gammai33 * u3 * u3
            + 2.0 * (gammai12 * u1 * u2 + gammai13 * u1 * u3 + gammai23 * u2 * u3 )
          );
        }

        real_t U0 = W/(DIFFalpha + 1.0);
        
        Dx._c(i,j,k) = (gammai11 * u1 + gammai12 * u2 + gammai13 * u3) / U0 - beta1;
        Dy._c(i,j,k) = (gammai12 * u1 + gammai22 * u2 + gammai23 * u3) / U0 - beta2;
        Dz._c(i,j,k) = (gammai13 * u1 + gammai23 * u2 + gammai33 * u3) / U0 - beta3;

        vx._c(i,j,k) = -1.0*W*d1alpha + u1*d1beta1 + u2*d1beta2 + u3*d1beta3
          -0.5 / U0 * (
            d1gammai11 * u1 * u1 + d1gammai22 * u2 * u2 + d1gammai33 * u3 * u3
            + 2.0 * (d1gammai12 * u1 * u2 + d1gammai13 * u1 * u3 + d1gammai23 * u2 * u3)
          );

        vy._c(i,j,k) = -1.0*W*d2alpha + u1*d2beta1 + u2*d2beta2 + u3*d2beta3
          -0.5 / U0 * (
            d2gammai11 * u1 * u1 + d2gammai22 * u2 * u2 + d2gammai33 * u3 * u3
            + 2.0 * (d2gammai12 * u1 * u2 + d2gammai13 * u1 * u3 + d2gammai23 * u2 * u3)
          );

        vz._c(i,j,k) = -1.0*W*d3alpha + u1*d3beta1 + u2*d3beta2 + u3*d3beta3
          -0.5 / U0 * (
            d3gammai11 * u1 * u1 + d3gammai22 * u2 * u2 + d3gammai33 * u3 * u3
            + 2.0 * (d3gammai12 * u1 * u2 + d3gammai13 * u1 * u3 + d3gammai23 * u2 * u3)
          );
      }

  if(rescale_sheet) rescaleAllFieldPerturbations(bssn, 1.0/rescale_sheet);
}

void Sheet::stepInit()
{
  Dx.stepInit();
  Dy.stepInit();
  Dz.stepInit();
  vx.stepInit();
  vy.stepInit();
  vz.stepInit();
}

void Sheet::K1Finalize()
{
  Dx.K1Finalize();
  Dy.K1Finalize();
  Dz.K1Finalize();
  vx.K1Finalize();
  vy.K1Finalize();
  vz.K1Finalize();
}

void Sheet::K2Finalize()
{
  Dx.K2Finalize();
  Dy.K2Finalize();
  Dz.K2Finalize();
  vx.K2Finalize();
  vy.K2Finalize();
  vz.K2Finalize();
}

void Sheet::K3Finalize()
{
  Dx.K3Finalize();
  Dy.K3Finalize();
  Dz.K3Finalize();
  vx.K3Finalize();
  vy.K3Finalize();
  vz.K3Finalize();
}

void Sheet::K4Finalize()
{
  Dx.K4Finalize();
  Dy.K4Finalize();
  Dz.K4Finalize();
  vx.K4Finalize();
  vy.K4Finalize();
  vz.K4Finalize();
}

std::vector<real_t> Sheet::getgammaIJ(idx_t s1, idx_t s2, idx_t s3, BSSN *bssnSim)
{
  std::vector<real_t> gammaIJ (7);

  arr_t & DIFFphi_p = *bssnSim->fields["DIFFphi_p"];
  arr_t & DIFFgamma11_p = *bssnSim->fields["DIFFgamma11_p"];
  arr_t & DIFFgamma22_p = *bssnSim->fields["DIFFgamma22_p"];
  arr_t & DIFFgamma33_p = *bssnSim->fields["DIFFgamma33_p"];
  arr_t & DIFFgamma12_p = *bssnSim->fields["DIFFgamma12_p"];
  arr_t & DIFFgamma13_p = *bssnSim->fields["DIFFgamma13_p"];
  arr_t & DIFFgamma23_p = *bssnSim->fields["DIFFgamma23_p"];

  real_t x_pt = Dx._p(s1,s2,s3) + _S1IDXtoX0(s1);
  real_t y_pt = Dy._p(s1,s2,s3) + _S2IDXtoY0(s2);
  real_t z_pt = Dz._p(s1,s2,s3) + _S3IDXtoZ0(s3);

  real_t x_idx = x_pt/dx;
  real_t y_idx = y_pt/dy;
  real_t z_idx = z_pt/dz;

  real_t phi = DIFFphi_p.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
  real_t detgamma = std::exp(12.0*phi);
  real_t gamma11 = std::exp(4.0*phi)*(1.0+DIFFgamma11_p.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx));
  real_t gamma22 = std::exp(4.0*phi)*(1.0+DIFFgamma22_p.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx));
  real_t gamma33 = std::exp(4.0*phi)*(1.0+DIFFgamma33_p.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx));
  real_t gamma12 = std::exp(4.0*phi)*DIFFgamma12_p.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
  real_t gamma13 = std::exp(4.0*phi)*DIFFgamma23_p.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
  real_t gamma23 = std::exp(4.0*phi)*DIFFgamma13_p.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);

  gammaIJ[0] = gamma11;
  gammaIJ[1] = gamma22;
  gammaIJ[2] = gamma33;
  gammaIJ[3] = gamma12;
  gammaIJ[4] = gamma13;
  gammaIJ[5] = gamma23;
  gammaIJ[6] = detgamma;

  return gammaIJ;
}

std::vector<real_t> Sheet::getgammaiIJ(idx_t s1, idx_t s2, idx_t s3, BSSN *bssnSim)
{
  std::vector<real_t> gammaiIJ (6);

  arr_t & DIFFphi_p = *bssnSim->fields["DIFFphi_p"];
  arr_t & DIFFgamma11_p = *bssnSim->fields["DIFFgamma11_p"];
  arr_t & DIFFgamma22_p = *bssnSim->fields["DIFFgamma22_p"];
  arr_t & DIFFgamma33_p = *bssnSim->fields["DIFFgamma33_p"];
  arr_t & DIFFgamma12_p = *bssnSim->fields["DIFFgamma12_p"];
  arr_t & DIFFgamma13_p = *bssnSim->fields["DIFFgamma13_p"];
  arr_t & DIFFgamma23_p = *bssnSim->fields["DIFFgamma23_p"];

  real_t x_pt = Dx._p(s1,s2,s3) + _S1IDXtoX0(s1);
  real_t y_pt = Dy._p(s1,s2,s3) + _S2IDXtoY0(s2);
  real_t z_pt = Dz._p(s1,s2,s3) + _S3IDXtoZ0(s3);

  real_t x_idx = x_pt/dx;
  real_t y_idx = y_pt/dy;
  real_t z_idx = z_pt/dz;

  real_t phi = DIFFphi_p.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
  real_t DIFFgamma11 = DIFFgamma11_p.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
  real_t DIFFgamma22 = DIFFgamma22_p.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
  real_t DIFFgamma33 = DIFFgamma33_p.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
  real_t DIFFgamma12 = DIFFgamma12_p.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
  real_t DIFFgamma23 = DIFFgamma23_p.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
  real_t DIFFgamma13 = DIFFgamma13_p.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
  real_t gammai11 = std::exp(-4.0*phi)*(1.0 + DIFFgamma22 + DIFFgamma33 - pw2(DIFFgamma23) + DIFFgamma22*DIFFgamma33);
  real_t gammai22 = std::exp(-4.0*phi)*(1.0 + DIFFgamma11 + DIFFgamma33 - pw2(DIFFgamma13) + DIFFgamma11*DIFFgamma33);
  real_t gammai33 = std::exp(-4.0*phi)*(1.0 + DIFFgamma11 + DIFFgamma22 - pw2(DIFFgamma12) + DIFFgamma11*DIFFgamma22);
  real_t gammai12 = std::exp(-4.0*phi)*(DIFFgamma13*DIFFgamma23 - DIFFgamma12*(1.0 + DIFFgamma33));
  real_t gammai13 = std::exp(-4.0*phi)*(DIFFgamma12*DIFFgamma23 - DIFFgamma13*(1.0 + DIFFgamma22));
  real_t gammai23 = std::exp(-4.0*phi)*(DIFFgamma12*DIFFgamma13 - DIFFgamma23*(1.0 + DIFFgamma11));

  gammaiIJ[0] = gammai11;
  gammaiIJ[1] = gammai22;
  gammaiIJ[2] = gammai33;
  gammaiIJ[3] = gammai12;
  gammaiIJ[4] = gammai13;
  gammaiIJ[5] = gammai23;

  return gammaiIJ;
}

// dot product of vectors (not ok for 4-vectors)
real_t dot_cov_spatial_vectors(real_t * v1, real_t * v2, std::vector<real_t> gammaiIJ)
{
  real_t gammai11 = gammaiIJ[0];
  real_t gammai22 = gammaiIJ[1];
  real_t gammai33 = gammaiIJ[2];
  real_t gammai12 = gammaiIJ[3];
  real_t gammai13 = gammaiIJ[4];
  real_t gammai23 = gammaiIJ[5];

  return (
      gammai11*v1[0]*v2[0] + gammai22*v1[1]*v2[1] + gammai33*v1[2]*v2[2]
      + gammai12*v1[0]*v2[1] + gammai13*v1[0]*v2[2] + gammai23*v1[1]*v2[2]
      + gammai12*v1[1]*v2[0] + gammai13*v1[2]*v2[1] + gammai23*v1[2]*v2[1]
    );
}
real_t mag_cov_spatial_vector(real_t * v1, std::vector<real_t> gammaiIJ)
{
  return std::sqrt( dot_cov_spatial_vectors(v1, v1, gammaiIJ) );
}
real_t dot_cont_spatial_vectors(real_t * v1, real_t * v2, std::vector<real_t> gammaIJ)
{
  real_t gamma11 = gammaIJ[0];
  real_t gamma22 = gammaIJ[1];
  real_t gamma33 = gammaIJ[2];
  real_t gamma12 = gammaIJ[3];
  real_t gamma13 = gammaIJ[4];
  real_t gamma23 = gammaIJ[5];

  return (
      gamma11*v1[0]*v2[0] + gamma22*v1[1]*v2[1] + gamma33*v1[2]*v2[2]
      + gamma12*v1[0]*v2[1] + gamma13*v1[0]*v2[2] + gamma23*v1[1]*v2[2]
      + gamma12*v1[1]*v2[0] + gamma13*v1[2]*v2[1] + gamma23*v1[2]*v2[1]
    );
}
real_t mag_cont_spatial_vector(real_t * v1, std::vector<real_t> gammaIJ)
{
  return std::sqrt( dot_cont_spatial_vectors(v1, v1, gammaIJ) );
}
real_t dot_cov_cont_spatial_vectors(real_t * v1, real_t * v2)
{
  return (
      v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    );
}

real_t dot_4_vectors_vvg(real_t v1[4], real_t v2[4], real_t g[4][4])
{
  real_t mag = 0.0;
  for(idx_t i=0; i<4; ++i)
    for(idx_t j=0; j<4; ++j)
      mag += g[i][j]*v1[i]*v2[j];
  return mag;
}

real_t dot_4_vectors_vv(real_t v1[4], real_t v2[4])
{
  real_t mag = 0.0;
  for(idx_t i=0; i<4; ++i)
    mag += v1[i]*v2[i];
  return mag;
}

/**
 * @brief      Function to get metric/data for a particular particle;
 * Currently only ok for zero shift (maybe ok in general?); returns data from _p register
 */
std::vector<real_t> Sheet::getRayDataAtS(idx_t s, BSSN *bssnSim, Lambda * lambda)
{
  std::vector<real_t> sheet_data (9);
  idx_t s1 = s, s2 = 0, s3 = 0;

  arr_t & DIFFr_a = *bssnSim->fields["DIFFr_a"];
  arr_t & S1_a = *bssnSim->fields["S1_a"];
  arr_t & S2_a = *bssnSim->fields["S2_a"];
  arr_t & S3_a = *bssnSim->fields["S3_a"];
  arr_t & DIFFalpha_p = *bssnSim->fields["DIFFalpha_p"];

  std::vector<real_t> gammaiIJ = getgammaiIJ(s1, s2, s3, bssnSim);
  real_t gammai11 = gammaiIJ[0];
  real_t gammai22 = gammaiIJ[1];
  real_t gammai33 = gammaiIJ[2];
  real_t gammai12 = gammaiIJ[3];
  real_t gammai13 = gammaiIJ[4];
  real_t gammai23 = gammaiIJ[5];

  std::vector<real_t> gammaIJ = getgammaIJ(s1, s2, s3, bssnSim);
  real_t gamma11 = gammaIJ[0];
  real_t gamma22 = gammaIJ[1];
  real_t gamma33 = gammaIJ[2];
  real_t gamma12 = gammaIJ[3];
  real_t gamma13 = gammaIJ[4];
  real_t gamma23 = gammaIJ[5];

  real_t x_pt = Dx._p(s1,s2,s3) + _S1IDXtoX0(s1);
  real_t y_pt = Dy._p(s1,s2,s3) + _S2IDXtoY0(s2);
  real_t z_pt = Dz._p(s1,s2,s3) + _S3IDXtoZ0(s3);
  real_t k1 = vx._p(s1,s2,s3);
  real_t k2 = vy._p(s1,s2,s3);
  real_t k3 = vz._p(s1,s2,s3);
  
  real_t x_idx = x_pt/dx;
  real_t y_idx = y_pt/dy;
  real_t z_idx = z_pt/dz;

  real_t rho_ADM = DIFFr_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
  real_t alpha = DIFFalpha_p.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
  real_t S1 = S1_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
  real_t S2 = S2_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
  real_t S3 = S3_a.getTriCubicInterpolatedValue(x_idx, y_idx, z_idx);
  
  real_t W_k = std::sqrt(
    gammai11*k1*k1 + gammai22*k2*k2 + gammai33*k3*k3
    + 2.0*( gammai12*k1*k2 + gammai13*k1*k3 + gammai23*k2*k3 )
  );

  real_t rho_ADM_m = rho_ADM - lambda->getLambda(); // matter only density
  real_t W_u = 1.0/std::sqrt( 1.0 - 
    ( gammai11*S1*S1 + gammai22*S2*S2 + gammai33*S3*S3 + 2.0*( gammai12*S1*S2 + gammai13*S1*S3 + gammai23*S2*S3 ) ) / rho_ADM_m/rho_ADM_m
  ); // fluid W-factor
  real_t rho_m = rho_ADM_m / W_u / W_u; // density in matter rest frame

  // observed photon frequency
  real_t u1 = S1/rho_m/W_u, u2 = S2/rho_m/W_u, u3 = S3/rho_m/W_u;
  real_t w = W_u*W_k - (
      gammai11*k1*u1 + gammai12*k1*u2 + gammai13*k1*u3
      + gammai12*k2*u1 + gammai22*k2*u2 + gammai23*k2*u3
      + gammai13*k3*u1 + gammai23*k3*u2 + gammai33*k3*u3
    );

  // angular diameter for a coordinate observer at rest

  // covariant k_0, u_0
  real_t k0 = -alpha*W_k;
  real_t k[4] = { k0, k1, k2, k3 };
  real_t u0 = -alpha*W_u;
  real_t u[4] = { u0, u1, u2, u3 };
  // direction vector d
  real_t d[4];
  for(idx_t i=0; i<4; ++i)
    d[i] = k[i]/w - u[i];
  // 4-metric (no shift...)
  real_t g[4][4] = {0};
  g[0][0] = -alpha*alpha;
  g[1][1] = gamma11; g[2][1] = gamma21; g[3][1] = gamma31;
  g[1][2] = gamma12; g[2][2] = gamma22; g[3][2] = gamma32;
  g[1][3] = gamma13; g[2][3] = gamma23; g[3][3] = gamma33;
  // screen projection matrix
  real_t Smn[4][4] = {0};
  for(idx_t i=0; i<4; ++i)
    for(idx_t j=0; j<4; ++j)
      Smn[i][j] = g[i][j] + u[i]*u[j] - d[i]*d[j];

  // get separation vectors; time separation is zero
  real_t sep1[4] = { 0, Dx._p(s1,1,0) - Dx._p(s1,0,0),
    Dy._p(s1,1,0) + _S2IDXtoY0((idx_t) 1) - Dy._p(s1,0,0), Dz._p(s1,1,0) - Dz._p(s1,0,0) };
  real_t sep2[4] = { 0, Dx._p(s1,2,0) - Dx._p(s1,0,0),
    Dy._p(s1,2,0) + _S2IDXtoY0((idx_t) 2) - Dy._p(s1,0,0), Dz._p(s1,2,0) - Dz._p(s1,0,0) };

  // Get screen vectors
  // Use separation vectors as guesses (dont worry about co/contra, just a guess), and orthonormalize
  real_t V1[4], v1[4];
  real_t V2[4], v2[4];
  for(int i=0; i<4; ++i)
  {
    V1[i] = sep1[i];
    V2[i] = sep2[i];
  }
  // Project & normalize V1
  real_t v1norm = std::sqrt( dot_4_vectors_vvg(V1, V1, Smn) );
  for(int i=0; i<4; i++)
    v1[i] = ( Smn[i][0]*V1[0] + Smn[i][1]*V1[1] + Smn[i][2]*V1[2] + Smn[i][3]*V1[3] )/v1norm;
  // project and normalize V2
  real_t v1V2 = 0.0;
  for(int i=0; i<4; i++)
    v1V2 += v1[i]*V2[i];
  real_t SVV = dot_4_vectors_vvg(V2, V2, Smn);
  real_t v2norm = std::sqrt( std::abs( SVV - v1V2*v1V2 ) );
  for(int i=0; i<4; i++)
    v2[i] = ( Smn[i][0]*V2[0] + Smn[i][1]*V2[1] + Smn[i][2]*V2[2] + Smn[i][3]*V2[3] - v1[i]*v1V2 ) / v2norm;

  real_t D11 = dot_4_vectors_vv(sep1, v1);
  real_t D12 = dot_4_vectors_vv(sep1, v2);
  real_t D21 = dot_4_vectors_vv(sep2, v1);
  real_t D22 = dot_4_vectors_vv(sep2, v2);
  real_t DA = std::sqrt(std::abs(D11*D22 - D12*D21)) / std::atan(ray_bundle_epsilon);

  sheet_data[0] = x_pt;
  sheet_data[1] = y_pt;
  sheet_data[2] = z_pt;

  sheet_data[3] = u1;
  sheet_data[4] = u2;
  sheet_data[5] = u3;

  sheet_data[6] = w;
  sheet_data[7] = rho_m;

  sheet_data[8] = DA;

  return sheet_data;
}


} // namespace
