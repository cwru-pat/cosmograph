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
  return;
}

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
  
  real_t A = sheetSim->lx*sheetSim->lx*std::stod(_config("peak_amplitude", "0.0001"));
  iodata->log( "Generating ICs with peak amp. = " + stringify(A) );

  real_t rho_FRW = 3.0/PI/8.0;
  real_t K_FRW = -sqrt(24.0*PI*rho_FRW);
  iodata->log( "FRW density is " + stringify(rho_FRW) + ", total mass in simulation volume is "
    + stringify(rho_FRW*sheetSim->lx*sheetSim->ly*sheetSim->lz));

  // the conformal factor in front of metric is the solution to
  // d^2 exp(\phi) = -2*pi exp(5\phi) * \delta_rho
  // generate random mode in \phi
  // delta_rho = -(lap e^\phi)/e^(4\phi)/2pi
  real_t phix = 0;
  real_t twopi_L = 2.0*PI / sheetSim->lx;
  real_t pw2_twopi_L = twopi_L*twopi_L;
  
  // grid values
#pragma omp parallel for
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    real_t x_frac = ((real_t) i / (real_t) NX);
    real_t phi = A*sin(2.0*PI*x_frac + phix);
    real_t rho = rho_FRW + -exp(-4.0*phi)/PI/2.0*(
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
    real_t rho = rho_FRW + -exp(-4.0*phi)/PI/2.0*(
      pw2(twopi_L*A*cos(2.0*PI*x_frac + phix))
      - pw2_twopi_L*A*sin(2.0*PI*x_frac + phix)
    );


    tot_mass += rho  * integration_interval * sheetSim->ly * sheetSim->lz;
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
    real_t rho = rho_FRW + -exp(-4.0*phi)/PI/2.0*(
      pw2(twopi_L*A*cos(2.0*PI*x_frac + phix))
      - pw2_twopi_L*A*sin(2.0*PI*x_frac + phix)
    );

    
    cur_mass += rho *  integration_interval * sheetSim->ly * sheetSim->lz;

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

  
} // namespace cosmo
