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

  carrier_count_scheme = static_cast<carrierCountScheme> (std::stoi(_config["carrier_count_scheme"]));
  deposit = static_cast<depositScheme> (std::stoi(_config["deposit_scheme"]));

  carriers_per_dx = std::stod(_config["carriers_per_dx"]);
  carriers_per_dy = std::stod(_config["carriers_per_dy"]);
  carriers_per_dz = std::stod(_config["carriers_per_dz"]);
}


Sheet::~Sheet()
{
  // Anything to do?
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
  #pragma omp parallel for collapse(2)        
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
        
              real_t W = std::sqrt( 1.0 +
                gammai11*u1*u1 + gammai22*u2*u2 + gammai33*u3*u3
                + 2.0*( gammai12*u1*u2 + gammai13*u1*u3 + gammai23*u2*u3 )
              );
              
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


void Sheet::RKStep(BSSN *bssn)
{
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
          
        real_t W = std::sqrt( 1.0 +
                              gammai11 * u1 * u1 + gammai22 * u2 * u2 + gammai33 * u3 * u3
                              + 2.0 * (gammai12 * u1 * u2 + gammai13 * u1 * u3 + gammai23 * u2 * u3 )
        );

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

} // namespace
