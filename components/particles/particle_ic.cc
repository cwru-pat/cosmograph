#include "particle_ic.h"
#include "../../cosmo_includes.h"
#include "../../cosmo_types.h"
#include "../../cosmo_globals.h"
#include "../../ICs/ICs.h"
#include "../../utils/math.h"

namespace cosmo
{

/**
 * @brief Initialize particles from gaussian random field data
 * @details Initialize particles from gaussian random field data, so particle
 *  masses are that needed to recreate the corresponding density field.
 */
void particle_ic_set_random(BSSN * bssnSim, Particles * particles, Fourier * fourier,
  IOData * iodata)
{
  idx_t i, j, k;
  ICsData icd = cosmo_get_ICsData();
  real_t rho_FRW = icd.rho_K_matter;
  arr_t & DIFFr = *bssnSim->fields["DIFFr_a"];
  arr_t & DIFFphi_p = *bssnSim->fields["DIFFphi_p"];
  arr_t & DIFFphi_a = *bssnSim->fields["DIFFphi_a"];
  arr_t & DIFFphi_f = *bssnSim->fields["DIFFphi_f"];
  iodata->log( "Generating ICs with peak at k = " + stringify(icd.peak_k) );
  iodata->log( "Generating ICs with peak amp. = " + stringify(icd.peak_amplitude) );

  // the conformal factor in front of metric is the solution to
  // d^2 exp(\phi) = -2*pi exp(5\phi) * \rho
  // generate gaussian random field 1 + xi = exp(phi) (use phi_p as a proxy):
  set_gaussian_random_field(DIFFphi_p, fourier, &icd);

  // rho = -lap(phi)/xi^5/2pi
# pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i, j, k)
  {
    DIFFr[NP_INDEX(i,j,k)] = rho_FRW - 0.5/PI/(
      pow(1.0 + DIFFphi_p[NP_INDEX(i,j,k)], 5.0)
    )*(
      double_derivative(i, j, k, 1, 1, DIFFphi_p)
      + double_derivative(i, j, k, 2, 2, DIFFphi_p)
      + double_derivative(i, j, k, 3, 3, DIFFphi_p)
    );
  }

  LOOP3(i, j, k)
  {
#   if 0
    // 1 particle per gridpoint
    real_t rho = DIFFr[INDEX(i,j,k)];
    real_t rootdetg = std::exp(6.0*DIFFphi_p[INDEX(i,j,k)]);
    real_t mass = rho*dx*dx*dx*rootdetg;
    Particle<real_t> particle = {0};
    particle.X[0] = i*dx;
    particle.X[1] = j*dx;
    particle.X[2] = k*dx;
    particle.M = mass;
    particles->addParticle( particle );
#   elif 1
    // 8 particles per gridpoint
    // linear interpolation to get particles 1/2-way between
    // neighboring gridpoints

    // TODO: doc this or something
    real_t m[2][2][2];
    for(int fx=0; fx<=1; ++fx)
      for(int fy=0; fy<=1; ++fy)
        for(int fz=0; fz<=1; ++fz)
        {
          real_t rho = DIFFr[INDEX(i+fx,j+fy,k+fz)];
          real_t rootdetg = std::exp(6.0*DIFFphi_p[INDEX(i+fx,j+fy,k+fz)]);
          m[fx][fy][fz] = rho*dx*dx*dx*rootdetg / 8.0; // 8 particles per gridpoint
        }

    for(int fx=0; fx<=1; ++fx)
      for(int fy=0; fy<=1; ++fy)
        for(int fz=0; fz<=1; ++fz)
        {
          // interpolated mass
          real_t mass = (
              m[0][0][0] + fx*fy*fz*m[1][1][1]
              + fx*m[1][0][0] + fy*m[0][1][0] + fz*m[0][0][1]
              + fx*fy*m[1][1][0] + fy*fz*m[0][1][1] + fx*fz*m[1][0][1]
              ) / std::pow(2.0, fx+fy+fz);

          Particle<real_t> particle = {0};
          particle.X[0] = (i+fx/2.0)*dx;
          particle.X[1] = (j+fy/2.0)*dx;
          particle.X[2] = (k+fz/2.0)*dx;
          particle.M = mass;
          particles->addParticle( particle );
        }
#   else
    // 3^3 particles per gridpoint
    // linear interpolation to get particles 1/2-way between
    // neighboring gridpoints

    // TODO: doc this or something
    real_t m[2][2][2];
    for(int fx=0; fx<=1; ++fx)
      for(int fy=0; fy<=1; ++fy)
        for(int fz=0; fz<=1; ++fz)
        {
          real_t rho = DIFFr[INDEX(i+fx,j+fy,k+fz)];
          real_t rootdetg = std::exp(6.0*DIFFphi_p[INDEX(i+fx,j+fy,k+fz)]);
          m[fx][fy][fz] = rho*dx*dx*dx*rootdetg / 27.0; // 27 particles per gridpoint
        }
    for(int fx=0; fx<=2; ++fx)
      for(int fy=0; fy<=2; ++fy)
        for(int fz=0; fz<=2; ++fz)
        {
          // interpolated mass
          real_t mass = (
              (3.0-fx)*(3.0-fy)*(3.0-fz)*m[0][0][0] + fx*fy*fz*m[1][1][1]
              + fx*(3.0-fy)*(3.0-fz)*m[1][0][0] + (3.0-fx)*fy*(3.0-fz)*m[0][1][0] + (3.0-fx)*(3.0-fy)*fz*m[0][0][1]
              + fx*fy*(3.0-fz)*m[1][1][0] + (3.0-fx)*fy*fz*m[0][1][1] + fx*(3.0-fy)*fz*m[1][0][1]
            ) / 27.0;
          Particle<real_t> particle = {0};
          particle.X[0] = (i+fx/3.0)*dx;
          particle.X[1] = (j+fy/3.0)*dx;
          particle.X[2] = (k+fz/3.0)*dx;
          particle.M = mass;
          particles->addParticle( particle );
        }
#   endif
  }

  // phi = ln(xi)
  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k) {
    idx_t idx = NP_INDEX(i,j,k);
    DIFFphi_a[idx] = log1p(DIFFphi_p[idx]);
    DIFFphi_f[idx] = log1p(DIFFphi_p[idx]);
    DIFFphi_p[idx] = log1p(DIFFphi_p[idx]);
  }

  // Make sure min density value > 0
  // Set conserved density variable field
  real_t min = icd.rho_K_matter;
  real_t max = min;
  LOOP3(i,j,k)
  {
    real_t rho = DIFFr[NP_INDEX(i,j,k)];

    if(rho < min)
    {
      min = rho;
    }
    if(rho > max)
    {
      max = rho;
    }
    if(rho != rho)
    {
      iodata->log("Error: NaN energy density.");
      throw -1;
    }
  }

  iodata->log( "Minimum fluid density: " + stringify(min) );
  iodata->log( "Maximum fluid density: " + stringify(max) );
  iodata->log( "Average fluctuation density: " + stringify(average(DIFFr)) );
  iodata->log( "Std.dev fluctuation density: " + stringify(standard_deviation(DIFFr)) );
  if(min < 0.0)
  {
    iodata->log( "Error: negative density in some regions.");
    throw -1;
  }

  arr_t & DIFFK_p = *bssnSim->fields["DIFFK_p"];
  arr_t & DIFFK_a = *bssnSim->fields["DIFFK_a"];
  #pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    DIFFK_a[idx] = -sqrt(24.0*PI*rho_FRW);
    DIFFK_p[idx] = -sqrt(24.0*PI*rho_FRW);
  }
}

/**
 * @brief Initialize particles per sinusoidal mode
 */
void particle_ic_set_sinusoid(BSSN * bssnSim, Particles * particles, IOData * iodata)
{
  iodata->log("Setting sinusoidal ICs.");
  idx_t i, j, k;

  // conformal factor
  arr_t & DIFFphi_p = *bssnSim->fields["DIFFphi_p"];
  // DIFFK is initially zero
  arr_t & DIFFK_p = *bssnSim->fields["DIFFK_p"];
  // matter sources
  arr_t & DIFFr_a = *bssnSim->fields["DIFFr_a"];

  real_t A = H_LEN_FRAC*H_LEN_FRAC*std::stod(_config("peak_amplitude_frac", "0.0001"));
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
  for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
      for(k=0; k<Nz; ++k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    real_t x = ((real_t) i / (real_t) Nx);
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
  idx_t particles_per_dx = std::stoi(_config("particles_per_dx", "1"));
  iodata->log("Particles per dx: " + stringify(particles_per_dx));
  for(i=0; i<Nx*particles_per_dx; ++i)
    for(j=0; j<Ny; ++j)
      for(k=0; k<Nz; ++k)
  {
    real_t x = ((real_t) i / (real_t) Nx / (real_t) particles_per_dx);

    real_t phi = A*sin(2.0*PI*x + phix);

    // \rho at a few/adjacent points
    real_t rho = rho_FRW + -exp(-4.0*phi)/PI/2.0*(
        pw2(twopi_L*A*cos(2.0*PI*x + phix))
        - pw2_twopi_L*A*sin(2.0*PI*x + phix)
      );
    real_t xp = x + dx/particles_per_dx;
    real_t xm = x - dx/particles_per_dx;
    real_t rhop = rho_FRW + -exp(-4.0*phi)/PI/2.0*(
        pw2(twopi_L*A*cos(2.0*PI*xp + phix))
        - pw2_twopi_L*A*sin(2.0*PI*xp + phix)
      );
    real_t rhom = rho_FRW + -exp(-4.0*phi)/PI/2.0*(
        pw2(twopi_L*A*cos(2.0*PI*xm + phix))
        - pw2_twopi_L*A*sin(2.0*PI*xm + phix)
      );
    // deconvolution to get mass
    // mass is distributed across nearby points;
    // try to counteract this
    // TODO: tune this?
    real_t stren = std::stod(_config("deconvolution_strength", "1.0"));
    rho = -stren*rhop + (1.0+2.0*stren)*rho - stren*rhom;


    real_t rootdetg = std::exp(6.0*phi);
    Particle<real_t> particle = {0};
    particle.X[0] = ((real_t) i)/((real_t) particles_per_dx)*dx;
    particle.X[1] = j*dx;
    particle.X[2] = k*dx;
    particle.M = rho*(dx/particles_per_dx)*dx*dx*rootdetg;
    particles->addParticle( particle );
  }


  // find min/max density,
  // Make sure min density value > 0
  real_t min = rho_FRW;
  real_t max = min;
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    real_t rho = rho_FRW + DIFFr_a[idx];

    if(rho < min)
    {
      min = rho;
    }
    if(rho > max)
    {
      max = rho;
    }
    if(rho != rho)
    {
      iodata->log("Error: NaN energy density.");
      throw -1;
    }
  }

  iodata->log( "Minimum fluid density: " + stringify(min) );
  iodata->log( "Maximum fluid density: " + stringify(max) );
  iodata->log( "Average fluctuation density: " + stringify(average(DIFFr_a)) );
  iodata->log( "Std.dev fluctuation density: " + stringify(standard_deviation(DIFFr_a)) );
  if(min < 0.0)
  {
    iodata->log("Error: negative density in some regions.");
    throw -1;
  }
}

void particle_ic_set_sinusoid_to_compare(BSSN * bssnSim, Particles * particles, IOData * iodata)
{
  iodata->log("Setting sinusoidal ICs.");
  idx_t i, j, k;

  // conformal factor
  arr_t & DIFFphi_p = *bssnSim->fields["DIFFphi_p"];
  // DIFFK is initially zero
  arr_t & DIFFK_p = *bssnSim->fields["DIFFK_p"];
  // matter sources
  arr_t & DIFFr_a = *bssnSim->fields["DIFFr_a"];

  real_t A = H_LEN_FRAC*H_LEN_FRAC*std::stod(_config("peak_amplitude_frac", "0.0001"));
  iodata->log( "Generating ICs with peak amp. = " + stringify(A) );

  real_t rho_FRW = 3.0/PI/8.0;
  real_t K_FRW = -sqrt(24.0*PI*rho_FRW);

  // the conformal factor in front of metric is the solution to
  // d^2 exp(\phi) = -2*pi exp(5\phi) * \delta_rho
  // generate random mode in \phi
  // delta_rho = -(lap e^\phi)/e^(4\phi)/2pi
  real_t phix = 0;
  real_t twopi_L = 2.0*PI/H_LEN_FRAC;
  real_t pw2_twopi_L = twopi_L*twopi_L;
  // grid values
  for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
      for(k=0; k<Nz; ++k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    real_t x = ((real_t) i / (real_t) Nx);
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
  idx_t particles_per_dx = std::stoi(_config("particles_per_dx", "1"));
  iodata->log("Particles per dx: " + stringify(particles_per_dx));
  for(i=0; i<Nx*particles_per_dx; ++i)
    for(j=0; j<Ny; ++j)
      for(k=0; k<Nz; ++k)
  {
    real_t x = ((real_t) i / (real_t) Nx / (real_t) particles_per_dx);

    real_t phi = A*sin(2.0*PI*x + phix);

    // \rho at a few/adjacent points
    real_t rho = rho_FRW + -exp(-4.0*phi)/PI/2.0*(
        pw2(twopi_L*A*cos(2.0*PI*x + phix))
        - pw2_twopi_L*A*sin(2.0*PI*x + phix)
      );
    real_t xp = x + dx/particles_per_dx;
    real_t xm = x - dx/particles_per_dx;
    real_t rhop = rho_FRW + -exp(-4.0*phi)/PI/2.0*(
        pw2(twopi_L*A*cos(2.0*PI*xp + phix))
        - pw2_twopi_L*A*sin(2.0*PI*xp + phix)
      );
    real_t rhom = rho_FRW + -exp(-4.0*phi)/PI/2.0*(
        pw2(twopi_L*A*cos(2.0*PI*xm + phix))
        - pw2_twopi_L*A*sin(2.0*PI*xm + phix)
      );
    // deconvolution to get mass
    // mass is distributed across nearby points;
    // try to counteract this
    // TODO: tune this?
    real_t stren = std::stod(_config("deconvolution_strength", "1.0"));
    rho = -stren*rhop + (1.0+2.0*stren)*rho - stren*rhom;


    real_t rootdetg = std::exp(6.0*phi);
    Particle<real_t> particle = {0};
    particle.X[0] = ((real_t) i)/((real_t) particles_per_dx)*dx;
    particle.X[1] = j*dx;
    particle.X[2] = k*dx;
    particle.M = rho*(dx/particles_per_dx)*dx*dx*rootdetg;
    particles->addParticle( particle );
  }


  // find min/max density,
  // Make sure min density value > 0
  real_t min = rho_FRW;
  real_t max = min;
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    real_t rho = rho_FRW + DIFFr_a[idx];

    if(rho < min)
    {
      min = rho;
    }
    if(rho > max)
    {
      max = rho;
    }
    if(rho != rho)
    {
      iodata->log("Error: NaN energy density.");
      throw -1;
    }
  }

  iodata->log( "Minimum fluid density: " + stringify(min) );
  iodata->log( "Maximum fluid density: " + stringify(max) );
  iodata->log( "Average fluctuation density: " + stringify(average(DIFFr_a)) );
  iodata->log( "Std.dev fluctuation density: " + stringify(standard_deviation(DIFFr_a)) );
  if(min < 0.0)
  {
    iodata->log("Error: negative density in some regions.");
    throw -1;
  }
}

  
/**
 * @brief Initialize particles per vector mode ID
 */
void particle_ic_set_vectorpert(BSSN * bssnSim, Particles * particles,
  IOData * iodata)
{
  iodata->log("Setting ICs based on a vector mode fluctuation.");
  // assumes functions vary in the y-direction
  idx_t i, j, k;

  // K is initially K_FRW
  arr_t & DIFFK_p = *bssnSim->fields["DIFFK_p"];
  // A_xy contains non-FRW fluctuation
  arr_t & A12_p = *bssnSim->fields["A12_p"];
  // beta1 (x-shift) to counter fluid velocity
  arr_t & beta1_p = *bssnSim->fields["beta1_p"];

  real_t b = std::stod(_config("peak_amplitude", "0.01"));
  iodata->log( "Generating ICs with b = " + stringify(b) );
  real_t use_initial_shift = std::stoi(_config("use_initial_shift", "1"));
  if(use_initial_shift)
  {
    iodata->log( "Using initial shift." );
  }
  else
  {
    iodata->log( "Not using initial shift." );
  }

  real_t rho_FRW = 3.0/PI/8.0;
  real_t K_FRW = -sqrt(24.0*PI*rho_FRW);
  real_t L = H_LEN_FRAC;
  real_t B = b*L*L;
  real_t phase = 2.0*PI*0.5*dx/L;

  // grid values
  for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
      for(k=0; k<Nz; ++k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    real_t y = j*dx;

    DIFFK_p[idx] = K_FRW;

    real_t Axy = B*PI/L*std::cos(2.0*PI*y/L + phase);
    A12_p[idx] = Axy;

    // more stable when using shift?
    if(use_initial_shift)
    {
      beta1_p[idx] = B*std::sin(2.0*PI*y/L + phase);
    }
  }

  // particle values
  idx_t particles_per_dy = std::stoi(_config("particles_per_dy", "1"));
  iodata->log("Particles per dx: " + stringify(particles_per_dy));


  for(i=0; i<Nx; ++i)
    for(j=0; j<Ny; ++j)
      for(k=0; k<Nz; ++k)
        for(int p=0; p<particles_per_dy; ++p)
  {
    real_t p_frac = dx/particles_per_dy;

    // particle position + adjacent particle positions
    real_t ys[3] = { j*dx + (p-1)*p_frac, j*dx + p*p_frac, j*dx + (p+1)*p_frac };

    // determined via Hamiltonian constraint:
    real_t Axy, Sx;
    real_t rho, Ux, W, M; // variables to construct (some intermediate)
    real_t MWs[3], MUs[3]; // variables to deconvolve
    for(int s=0; s<3; ++s)
    {
      Axy = B*PI/L*std::cos(2.0*PI*ys[s]/L + phase);
      Sx = -B*PI/L/L/4.0*std::sin(2.0*PI*ys[s]/L + phase);

      rho = (K_FRW*K_FRW/12.0 - 2.0*Axy*Axy/8.0)/2.0/PI;
      Ux = Sx/std::sqrt(rho*rho - Sx*Sx);
      W = std::sqrt(1.0 + Ux*Ux);

      // Particle mass
      M = rho*(dx/particles_per_dy)*dx*dx/W;

      MWs[s] = M*W;
      MUs[s] = M*Ux;
    }

    // de-convolved variables
    real_t stren = std::stod(_config("deconvolution_strength", "1.0"));
    real_t MW = -stren*MWs[0] + (1.0+2.0*stren)*MWs[1] - stren*MWs[2];
    real_t MU = -stren*MUs[0] + (1.0+2.0*stren)*MUs[1] - stren*MUs[2];

    Ux = (MU > 0 ? 1.0 : -1.0) / std::sqrt( pw2(MW/MU) - 1.0 );
    M = MU/Ux;

    Particle<real_t> particle = {0};
    particle.X[0] = i*dx;
    particle.X[1] = ys[1];
    particle.X[2] = k*dx;
    particle.U[0] = Ux;
    particle.M = M;

    particles->addParticle( particle );
  }

}

}
