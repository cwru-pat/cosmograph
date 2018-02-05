#include "sheets_ic.h"
#include "../../cosmo_includes.h"
#include "../../cosmo_types.h"
#include "../../cosmo_globals.h"
#include "../../utils/Fourier.h"
#include "../../utils/math.h"


namespace cosmo
{

void sheets_ic_sinusoid(
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
  
  real_t A = std::stod(_config("peak_amplitude", "0.0001"));
  iodata->log( "Generating ICs with peak amp. = " + stringify(A) );

  real_t rho_FRW = 3.0/PI/8.0;
  real_t K_FRW = -sqrt(24.0*PI*rho_FRW);

  // the conformal factor in front of metric is the solution to
  // d^2 exp(\phi) = -2*pi exp(5\phi) * \delta_rho
  // generate random mode in \phi
  // delta_rho = -(lap e^\phi)/e^(4\phi)/2pi
  real_t phix = 2.77;
  real_t twopi_L = 2.0*PI/H_LEN_FRAC;
  real_t pw2_twopi_L = twopi_L*twopi_L;
  // grid values
  for(i=0; i<NX; ++i)
    for(j=0; j<NY; ++j)
      for(k=0; k<NZ; ++k)
      {
        idx_t idx = NP_INDEX(i,j,k);

        real_t x = ((real_t) i / (real_t) NX);
        real_t phi = A*sin(2.0*PI*x + phix);
        real_t rho = rho_FRW + -exp(-4.0*phi)/PI/2.0*(
          pw2(twopi_L*A*cos(2.0*PI*x + phix))
          - pw2_twopi_L*A*sin(2.0*PI*x + phix)
        );

        // These aren't difference vars
        DIFFphi_p[NP_INDEX(i,j,k)] = phi;
        DIFFK_p[idx] = K_FRW;
        DIFFr_a[idx] = rho;
      }

  // particle values
  // parallelizing may break this, be careful
  //  idx_t particles_per_dx = std::stoi(_config("particles_per_dx", "1"));
  real_t integration_interval = std::stod(_config("integration_interval", "0.01"));


  int cur_s1 = 1;
  real_t cur_mass = 0;

  tot_mass = 0;

  for(real_t cur_x = 0; cur_x < H_LEN_FRAC; cur_x += integration_interval)  
  {
    real_t x = cur_x;

    real_t phi = A*sin(2.0*PI*x + phix);
    real_t rho = rho_FRW + -exp(-4.0*phi)/PI/2.0*(
      pw2(twopi_L*A*cos(2.0*PI*x + phix))
      - pw2_twopi_L*A*sin(2.0*PI*x + phix)
    );

    real_t rootdetg = std::exp(6.0*phi);

    tot_mass += rho * rootdetg * integration_interval;
  }

  real_t mass_per_voxel = tot_mass / (double)Dx.nx;

  for(real_t cur_x = 0; cur_x < H_LEN_FRAC; cur_x += integration_interval)
  {
    real_t x = cur_x;
    real_t phi = A*sin(2.0*PI*x + phix);
    real_t rho = rho_FRW + -exp(-4.0*phi)/PI/2.0*(
      pw2(twopi_L*A*cos(2.0*PI*x + phix))
      - pw2_twopi_L*A*sin(2.0*PI*x + phix)
    );

    real_t rootdetg = std::exp(6.0*phi);
    
    cur_mass += rootdetg * integration_interval * rho;

    if(cur_mass > mass_per_voxel)
    {
      for(j=0; j<Dx.ny; ++j)
        for(k=0; k<Dx.nz; ++k)
        {
          Dx(cur_s1, j, k) = integration_interval * (cur_mass - mass_per_voxel)
            / mass_per_voxel + cur_x;
        }
      cur_mass = cur_mass - mass_per_voxel;
      cur_s1++;
    }
  }
  if(cur_s1 < Dx.nx - 1)
  {
    for(j=0; j<NY; ++j)
      for(k=0; k<NZ; ++k)
      {
        Dx(Dx.nx-1, j, k) = H_LEN_FRAC;
      }
  }
  


  // iodata->log( "Minimum fluid density: " + stringify(min) );
  // iodata->log( "Maximum fluid density: " + stringify(max) );
  // iodata->log( "Average fluctuation density: " + stringify(average(DIFFr_a)) );
  // iodata->log( "Std.dev fluctuation density: " + stringify(standard_deviation(DIFFr_a)) );
  // if(min < 0.0)
  // {
  //   iodata->log("Error: negative density in some regions.");
  //   throw -1;
  // }

  
}


}
