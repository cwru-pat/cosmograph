#include "particles.h"
#include "particles_macros.h"
#include "../../utils/math.h"
#include "../../utils/Timer.h"
#include "../../IO/io.h"
#include "../../cosmo_globals.h"

namespace cosmo
{

Particles::Particles()
{
  particles = new particle_vec();
}

Particles::~Particles()
{
  // anything to do?
}

particle_vec * Particles::getParticleVec()
{
  return particles;
}

/**
 * @brief Create random particles
 * @details Add particles with randomized positions, zero velocity,
 * and unit mass to internal particles vector.
 * 
 * @param n_particles number of particles to create
 */
void Particles::init(idx_t n_particles)
{
  idx_t i = 0;

  std::mt19937 gen(7.0);
  std::uniform_real_distribution<real_t> dist(0.0, 1.0);
  dist(gen);

  for(i=0; i<n_particles; i++)
  {
    Particle<real_t> particle = {0};

    // Randomized position
    particle.X[0] = 0.0*dist(gen)*COSMO_N*dx;
    particle.X[1] = 0.0*dist(gen)*COSMO_N*dx;
    particle.X[2] = 0.0*dist(gen)*COSMO_N*dx;
    // Mass in units TBD
    particle.M = 1.0*dx*dx*dx;

    particle.U[0] = dist(gen)/3.0;
    particle.U[1] = dist(gen)/13.0;
    particle.U[2] = dist(gen)/2.0;

    ParticleRegister<real_t> particle_register;
    particle_register.p_p = particle;
    particle_register.p_a = particle;
    particle_register.p_c = particle;
    particle_register.p_f = particle;

    particles->push_back( particle_register );
  }
}

/**
 * @brief Add an individual particle
 * @details Add an individual particle to the internal list of particles
 * 
 * @param particle Particle struct
 */
void Particles::addParticle(Particle<real_t> particle)
{
  ParticleRegister<real_t> particle_register;
  particle_register.p_p = particle;
  particle_register.p_a = particle;
  particle_register.p_c = particle;
  particle_register.p_f = particle;
  particles->push_back( particle_register );
  return;
}

/**
 * @brief Compute non-integer index at a point.
 */
real_t Particles::getFractionalIndex(real_t x)
{
  return real_t_mod(x, COSMO_N*dx)/dx;
}

/**
 * @brief Returns the index "below" a point (floor of getFractionalIndex)
 * 
 * @return index associated with position
 */
idx_t Particles::getIndexBelow(real_t x)
{
  return (idx_t) getFractionalIndex(x); // floor to get 
}

/**
 * @brief Returns the index nearest a point (rounds getFractionalIndex)
 */
idx_t Particles::getNearestIndex(real_t x)
{
  return (idx_t) ( getFractionalIndex(x) + 0.5 );
}

/**
 * @brief Sets a fractional distance between grid points
 * @details Stores the fractional distance between grid points of a given
 * point. Eg, for dx = 0.5, the value 1.23 is 46/100 the way between
 * gridpoints at 1.0 and 1.5, so 0.46 is computed. Following this, supplying
 * the point (1.23, 4.56, 7.89) as an argument for X would cause x_d to be
 * set to (0.46, 0.12, 0.78)
 * 
 * @param X true position
 * @param x_d fractional position (passed by reference; values are mutated)
 */
void Particles::setX_d(real_t X[3], real_t x_d[3])
{
  for(idx_t i=0; i<3; i++)
  {
    real_t frac_idx = getFractionalIndex(X[i]);
    real_t idx_below = (real_t) getIndexBelow(X[i]);
    x_d[i] = frac_idx - idx_below;
  }
}

/**
 * @brief Linearly interpolate metric quantities from corners of a "box"
 * @details Linearly interpolate metric quantities from corners of a "box"
 * 
 * @param corner_pp_in A 2x2x2 array containing metric values at corners
 * of the "box"
 * @param x_d fractional position between corners inside the "box" to
 * interpolate values at (see setX_d)
 * @return [description]
 */
ParticleMetricPrimitives<real_t> Particles::interpolatePrimitivesFromCorners(
  ParticleMetricPrimitives<real_t> corner_pp_in[2][2][2], real_t x_d[3])
{
  ParticleMetricPrimitives<real_t> pp = {0};

  pp.rootdetg = PARTICLES_INTERPOLATION(rootdetg);
  pp.alpha = PARTICLES_INTERPOLATION(alpha);

  for(idx_t i=0; i<3; i++)
  {
    pp.beta[i] = PARTICLES_INTERPOLATION(beta[i]);
    pp.dalpha[i] = PARTICLES_INTERPOLATION(dalpha[i]);
    for(idx_t j=0; j<3; j++)
    {
      pp.dbeta[i][j] = PARTICLES_INTERPOLATION(dbeta[i][j]);
    }
  }

  for(idx_t i=0; i<6; i++)
  {
    pp.gi[i] = PARTICLES_INTERPOLATION(gi[i]);
    for(idx_t j=0; j<3; j++)
    {
      pp.dgi[j][i] = PARTICLES_INTERPOLATION(dgi[j][i]);
    }
  }

  return pp;
}

/**
 * @brief Only set some metric components
 * @details Only set some metric components:
 *  gamma_ij, rootdetg, and alpha
 *  see also: getInterpolatedPrimitivesIncomplete
 * 
 * @param corner_pp_in 2x2x2 array of metric primitives at corners (incomplete)
 * @return metric primitives (incomplete)
 */
ParticleMetricPrimitives<real_t> Particles::interpolatePrimitivesFromCornersIncomplete(
  ParticleMetricPrimitives<real_t> corner_pp_in[2][2][2], real_t x_d[3])
{
  ParticleMetricPrimitives<real_t> pp = {0};
  pp.rootdetg = PARTICLES_INTERPOLATION(rootdetg);
  pp.alpha = PARTICLES_INTERPOLATION(alpha);
  for(idx_t i=0; i<6; i++)
  {
    pp.gi[i] = PARTICLES_INTERPOLATION(gi[i]);
  }
  return pp;
}

real_t Particles::linearInterpolation(
  real_t C000, real_t C001, real_t C010, real_t C011, /* values at "C"orners, C_{x,y,z} */
  real_t C100, real_t C101, real_t C110, real_t C111, /* "binary" order */
  real_t x_d[3] /* normalized position within cube */
  )
{
  real_t C00 = C000*(1.0 - x_d[0]) + C100*x_d[0];
  real_t C01 = C001*(1.0 - x_d[0]) + C101*x_d[0];
  real_t C10 = C010*(1.0 - x_d[0]) + C110*x_d[0];
  real_t C11 = C011*(1.0 - x_d[0]) + C111*x_d[0];

  real_t C0 = C00*(1.0 - x_d[1]) + C10*x_d[1];
  real_t C1 = C01*(1.0 - x_d[1]) + C11*x_d[1];

  // return "C"
  return C0*(1.0 - x_d[2]) + C1*x_d[2];
}

/**
 * @brief      Interpolate some metric primitives required for geodesic integration
 * @details    Interpolate some metric primitives (eg, only some variables in a
 *  ParticleMetricPrimitives struct: alpha, phi, gamma_ij) required for geodesic
 *  integration near a particular particle. See also: getInterpolatedPrimitives
 *
 * @param      p            Instance of particle to compute interpolated metric
 *  quantities near (Particle type).
 * @param      bssn_fields  Map from bssn class to fields.
 *
 * @return     Returns ParticleMetricPrimitives with values near the particle
 */
ParticleMetricPrimitives<real_t> Particles::getInterpolatedPrimitivesIncomplete(Particle<real_t> * p,
  map_t & bssn_fields)
{
  arr_t & DIFFalpha_a = *bssn_fields["DIFFalpha_a"];
  arr_t & DIFFphi_a = *bssn_fields["DIFFphi_a"];
  arr_t & DIFFgamma11_a = *bssn_fields["DIFFgamma11_a"];
  arr_t & DIFFgamma12_a = *bssn_fields["DIFFgamma12_a"];
  arr_t & DIFFgamma13_a = *bssn_fields["DIFFgamma13_a"];
  arr_t & DIFFgamma22_a = *bssn_fields["DIFFgamma22_a"];
  arr_t & DIFFgamma23_a = *bssn_fields["DIFFgamma23_a"];
  arr_t & DIFFgamma33_a = *bssn_fields["DIFFgamma33_a"];

  // Linear interpolant
  ParticleMetricPrimitives<real_t> corner_pp[2][2][2];
  for(idx_t i=0; i<2; i++)
    for(idx_t j=0; j<2; j++)
      for(idx_t k=0; k<2; k++)
      {
        idx_t x_idx = getIndexBelow(p->X[0]) + i;
        idx_t y_idx = getIndexBelow(p->X[1]) + j;
        idx_t z_idx = getIndexBelow(p->X[2]) + k;
        idx_t idx = INDEX(x_idx, y_idx, z_idx);

        ParticleMetricPrimitives<real_t> pp = {0};

        pp.rootdetg = std::exp(6.0*DIFFphi_a[idx]);
        pp.alpha = DIFFalpha_a[idx];

        pp.gi[aIDX(1,1)] = std::exp(-4.0*DIFFphi_a[idx])*(1.0 + DIFFgamma22_a[idx] + DIFFgamma33_a[idx] - pw2(DIFFgamma23_a[idx]) + DIFFgamma22_a[idx]*DIFFgamma33_a[idx]);
        pp.gi[aIDX(2,2)] = std::exp(-4.0*DIFFphi_a[idx])*(1.0 + DIFFgamma11_a[idx] + DIFFgamma33_a[idx] - pw2(DIFFgamma13_a[idx]) + DIFFgamma11_a[idx]*DIFFgamma33_a[idx]);
        pp.gi[aIDX(3,3)] = std::exp(-4.0*DIFFphi_a[idx])*(1.0 + DIFFgamma11_a[idx] + DIFFgamma22_a[idx] - pw2(DIFFgamma12_a[idx]) + DIFFgamma11_a[idx]*DIFFgamma22_a[idx]);
        pp.gi[aIDX(1,2)] = std::exp(-4.0*DIFFphi_a[idx])*(DIFFgamma13_a[idx]*DIFFgamma23_a[idx] - DIFFgamma12_a[idx]*(1.0 + DIFFgamma33_a[idx]));
        pp.gi[aIDX(1,3)] = std::exp(-4.0*DIFFphi_a[idx])*(DIFFgamma12_a[idx]*DIFFgamma23_a[idx] - DIFFgamma13_a[idx]*(1.0 + DIFFgamma22_a[idx]));
        pp.gi[aIDX(2,3)] = std::exp(-4.0*DIFFphi_a[idx])*(DIFFgamma12_a[idx]*DIFFgamma13_a[idx] - DIFFgamma23_a[idx]*(1.0 + DIFFgamma11_a[idx]));

        corner_pp[i][j][k] = pp;
      }
  /* end i, j, k loops */

  real_t x_d[3];
  setX_d(p->X, x_d);
  return interpolatePrimitivesFromCornersIncomplete(corner_pp, x_d);
}

/**
 * @brief      Interpolate metric primitives required for geodesic integration
 * @details    Interpolate metric primitives (eg, variables in a
 *  ParticleMetricPrimitives struct) required for geodesic integration near a
 *  particular particle.
 *
 * @param      p            Instance of particle to compute interpolated metric
 *  quantities near (Particle type).
 * @param      bssn_fields  Map from bssn class to fields.
 *
 * @return     Returns ParticleMetricPrimitives with values near the particle
 */
ParticleMetricPrimitives<real_t> Particles::getInterpolatedPrimitives(Particle<real_t> * p,
  map_t & bssn_fields)
{
  arr_t & DIFFalpha_a = *bssn_fields["DIFFalpha_a"];
  
# if USE_BSSN_SHIFT
    arr_t & beta1_a = *bssn_fields["beta1_a"];
    arr_t & beta2_a = *bssn_fields["beta2_a"];
    arr_t & beta3_a = *bssn_fields["beta3_a"];
# endif

  arr_t & DIFFphi_a = *bssn_fields["DIFFphi_a"];

  arr_t & DIFFgamma11_a = *bssn_fields["DIFFgamma11_a"];
  arr_t & DIFFgamma12_a = *bssn_fields["DIFFgamma12_a"];
  arr_t & DIFFgamma13_a = *bssn_fields["DIFFgamma13_a"];
  arr_t & DIFFgamma22_a = *bssn_fields["DIFFgamma22_a"];
  arr_t & DIFFgamma23_a = *bssn_fields["DIFFgamma23_a"];
  arr_t & DIFFgamma33_a = *bssn_fields["DIFFgamma33_a"];

  // Linear interpolant
  ParticleMetricPrimitives<real_t> corner_pp[2][2][2];
  for(idx_t i=0; i<2; i++)
    for(idx_t j=0; j<2; j++)
      for(idx_t k=0; k<2; k++)
      {
        idx_t x_idx = getIndexBelow(p->X[0]) + i;
        idx_t y_idx = getIndexBelow(p->X[1]) + j;
        idx_t z_idx = getIndexBelow(p->X[2]) + k;
        idx_t idx = INDEX(x_idx, y_idx, z_idx);

        ParticleMetricPrimitives<real_t> pp;

        pp.rootdetg = std::exp(6.0*DIFFphi_a[idx]);
        pp.alpha = DIFFalpha_a[idx];

#       if USE_BSSN_SHIFT
          pp.beta[0] = beta1_a[idx];
          pp.beta[1] = beta2_a[idx];
          pp.beta[2] = beta3_a[idx];
#       else
          pp.beta[0] = 0;
          pp.beta[1] = 0;
          pp.beta[2] = 0;
#       endif

        pp.gi[aIDX(1,1)] = std::exp(-4.0*DIFFphi_a[idx])*(1.0 + DIFFgamma22_a[idx] + DIFFgamma33_a[idx] - pw2(DIFFgamma23_a[idx]) + DIFFgamma22_a[idx]*DIFFgamma33_a[idx]);
        pp.gi[aIDX(2,2)] = std::exp(-4.0*DIFFphi_a[idx])*(1.0 + DIFFgamma11_a[idx] + DIFFgamma33_a[idx] - pw2(DIFFgamma13_a[idx]) + DIFFgamma11_a[idx]*DIFFgamma33_a[idx]);
        pp.gi[aIDX(3,3)] = std::exp(-4.0*DIFFphi_a[idx])*(1.0 + DIFFgamma11_a[idx] + DIFFgamma22_a[idx] - pw2(DIFFgamma12_a[idx]) + DIFFgamma11_a[idx]*DIFFgamma22_a[idx]);
        pp.gi[aIDX(1,2)] = std::exp(-4.0*DIFFphi_a[idx])*(DIFFgamma13_a[idx]*DIFFgamma23_a[idx] - DIFFgamma12_a[idx]*(1.0 + DIFFgamma33_a[idx]));
        pp.gi[aIDX(1,3)] = std::exp(-4.0*DIFFphi_a[idx])*(DIFFgamma12_a[idx]*DIFFgamma23_a[idx] - DIFFgamma13_a[idx]*(1.0 + DIFFgamma22_a[idx]));
        pp.gi[aIDX(2,3)] = std::exp(-4.0*DIFFphi_a[idx])*(DIFFgamma12_a[idx]*DIFFgamma13_a[idx] - DIFFgamma23_a[idx]*(1.0 + DIFFgamma11_a[idx]));

        for(idx_t a=0; a<3; a++)
        {
          pp.dalpha[a] = DER(DIFFalpha_a);

#         if USE_BSSN_SHIFT
            pp.dbeta[a][0] = DER(beta1_a);
            pp.dbeta[a][1] = DER(beta2_a);
            pp.dbeta[a][2] = DER(beta3_a);
#         else
            pp.dbeta[a][0] = 0;
            pp.dbeta[a][1] = 0;
            pp.dbeta[a][2] = 0;
#         endif

          pp.dgi[a][aIDX(1,1)] = -4.0*DER(DIFFphi_a)*pp.gi[aIDX(1,1)]
            + std::exp(-4.0*DIFFphi_a[idx])*(DER(DIFFgamma22_a) + DER(DIFFgamma33_a) - 2.0*DIFFgamma23_a[idx]*DER(DIFFgamma23_a) + DER(DIFFgamma22_a)*DIFFgamma33_a[idx] + DIFFgamma22_a[idx]*DER(DIFFgamma33_a));
          pp.dgi[a][aIDX(2,2)] = -4.0*DER(DIFFphi_a)*pp.gi[aIDX(2,2)]
            + std::exp(-4.0*DIFFphi_a[idx])*(DER(DIFFgamma11_a) + DER(DIFFgamma33_a) - 2.0*DIFFgamma13_a[idx]*DER(DIFFgamma13_a) + DER(DIFFgamma11_a)*DIFFgamma33_a[idx] + DIFFgamma11_a[idx]*DER(DIFFgamma33_a));
          pp.dgi[a][aIDX(3,3)] = -4.0*DER(DIFFphi_a)*pp.gi[aIDX(3,3)]
            + std::exp(-4.0*DIFFphi_a[idx])*(DER(DIFFgamma11_a) + DER(DIFFgamma22_a) - 2.0*DIFFgamma12_a[idx]*DER(DIFFgamma12_a) + DER(DIFFgamma11_a)*DIFFgamma22_a[idx] + DIFFgamma11_a[idx]*DER(DIFFgamma22_a));
          pp.dgi[a][aIDX(1,2)] = -4.0*DER(DIFFphi_a)*pp.gi[aIDX(1,2)]
            + std::exp(-4.0*DIFFphi_a[idx])*(DER(DIFFgamma13_a)*DIFFgamma23_a[idx] + DIFFgamma13_a[idx]*DER(DIFFgamma23_a) - DER(DIFFgamma12_a)*(1.0 + DIFFgamma33_a[idx]) - DIFFgamma12_a[idx]*DER(DIFFgamma33_a));
          pp.dgi[a][aIDX(1,3)] = -4.0*DER(DIFFphi_a)*pp.gi[aIDX(1,3)]
            + std::exp(-4.0*DIFFphi_a[idx])*(DER(DIFFgamma12_a)*DIFFgamma23_a[idx] + DIFFgamma12_a[idx]*DER(DIFFgamma23_a) - DER(DIFFgamma13_a)*(1.0 + DIFFgamma22_a[idx]) - DIFFgamma13_a[idx]*DER(DIFFgamma22_a));
          pp.dgi[a][aIDX(2,3)] = -4.0*DER(DIFFphi_a)*pp.gi[aIDX(2,3)]
            + std::exp(-4.0*DIFFphi_a[idx])*(DER(DIFFgamma12_a)*DIFFgamma13_a[idx] + DIFFgamma12_a[idx]*DER(DIFFgamma13_a) - DER(DIFFgamma23_a)*(1.0 + DIFFgamma11_a[idx]) - DIFFgamma23_a[idx]*DER(DIFFgamma11_a));
        }

        corner_pp[i][j][k] = pp;
      }
  /* end i, j, k loops */

  real_t x_d[3];
  setX_d(p->X, x_d);
  return interpolatePrimitivesFromCorners(corner_pp, x_d);
}

/**
 * @brief      Implementation for evaluating a single RK step
 * @details    General implementation for computing one of the RK steps:
 *  y_p; y_a = y_c = yp; y_f = 0 (TODO)
 *  y_a = y_p
 *  
 *  y_c = y_p + dt/2 * f(y_a)    // RK1 step
 *  y_f += y_c
 *  y_a <-> y_c
 *  
 *  y_c = y_p + dt/2 * f(y_a)    // RK2 step
 *  y_f += 2 y_c
 *  y_a <-> y_c
 *  
 *  y_c = y_p + dt * f(y_a)      // RK3 step
 *  y_f += y_c
 *  y_a <-> y_c
 *  
 *  y_c = y_p + dt/2 * f(y_a)    // RK4 step
 *  y_f += y_c
 *  y_a <-> y_c
 *  
 *  (y_f = 5y_p + 1/2 K1 + K2 + K3 + 1/2 K4)
 *  y_f = y_f / 3 - 2/3 y_p      // finalize
 *  y_f <-> y_p
 *  
 *  Using Baumgarte & Shapiro, 5.223-5.225
 *
 *  For a single RK step, we can just do: y_c = y_p + h * f(y_a)
 *  and y_f += RK_sum_coeff * y_c
 *
 * @param      pr            Register of particles to integrate
 * @param[in]  h             step size for particular RK step
 * @param[in]  RK_sum_coeff  coefficient when adding RK step
 * @param      bssn_fields   Map from bssn class to fields.
 */
void Particles::RKStep(ParticleRegister<real_t> * pr, real_t h, real_t RK_sum_coeff,
  map_t & bssn_fields)
{
  Particle<real_t> & p_p = pr->p_p;
  Particle<real_t> & p_a = pr->p_a;
  Particle<real_t> & p_c = pr->p_c;
  Particle<real_t> & p_f = pr->p_f;

  ParticleMetricPrimitives<real_t> pp_a = getInterpolatedPrimitives(& p_a, bssn_fields);

// TODO: check!
  real_t W = std::sqrt( 1.0 + 
      pp_a.gi[aIDX(1,1)]*p_a.U[0]*p_a.U[0] + pp_a.gi[aIDX(2,2)]*p_a.U[1]*p_a.U[1] + pp_a.gi[aIDX(3,3)]*p_a.U[2]*p_a.U[2]
      + 2.0*( pp_a.gi[aIDX(1,2)]*p_a.U[0]*p_a.U[1] + pp_a.gi[aIDX(1,3)]*p_a.U[0]*p_a.U[2] + pp_a.gi[aIDX(2,3)]*p_a.U[1]*p_a.U[2] )
    );
  real_t U0 = W/pp_a.alpha;

  for(int i=1; i<=3; i++)
  {
    p_c.X[iIDX(i)] = p_p.X[iIDX(i)] + h*( (pp_a.gi[aIDX(i,1)]*p_a.U[iIDX(1)] + pp_a.gi[aIDX(i,2)]*p_a.U[iIDX(2)] + pp_a.gi[aIDX(i,3)]*p_a.U[iIDX(3)]) / U0 - pp_a.beta[iIDX(i)]);
    p_c.U[iIDX(i)] = p_p.U[iIDX(i)] + h*(
      -1.0*W*pp_a.dalpha[iIDX(i)] + p_a.U[iIDX(1)]*pp_a.dbeta[iIDX(i)][iIDX(1)] + p_a.U[iIDX(2)]*pp_a.dbeta[iIDX(i)][iIDX(2)] + p_a.U[iIDX(3)]*pp_a.dbeta[iIDX(i)][iIDX(3)]
      -1.0/2.0/U0*(
        pp_a.dgi[iIDX(i)][aIDX(1,1)]*p_a.U[0]*p_a.U[0] + pp_a.dgi[iIDX(i)][aIDX(2,2)]*p_a.U[1]*p_a.U[1] + pp_a.dgi[iIDX(i)][aIDX(3,3)]*p_a.U[2]*p_a.U[2]
        + 2.0*( pp_a.dgi[iIDX(i)][aIDX(1,2)]*p_a.U[0]*p_a.U[1] + pp_a.dgi[iIDX(i)][aIDX(1,3)]*p_a.U[0]*p_a.U[2] + pp_a.dgi[iIDX(i)][aIDX(2,3)]*p_a.U[1]*p_a.U[2] )
      )
    );

    p_f.X[iIDX(i)] += RK_sum_coeff*p_c.X[iIDX(i)];
    p_f.U[iIDX(i)] += RK_sum_coeff*p_c.U[iIDX(i)];
  }
}

void Particles::RK1Step(map_t & bssn_fields)
{
  _timer["Particles::RKCalcs"].start();
  PARTICLES_PARALLEL_LOOP(pr)
  {
    RKStep(& (*pr), dt/2.0, 1.0, bssn_fields);
  }
  _timer["Particles::RKCalcs"].stop();
}

void Particles::RK2Step(map_t & bssn_fields)
{
  _timer["Particles::RKCalcs"].start();
  PARTICLES_PARALLEL_LOOP(pr)
  {
    RKStep(& (*pr), dt/2.0, 2.0, bssn_fields);
  }
  _timer["Particles::RKCalcs"].stop();
}

void Particles::RK3Step(map_t & bssn_fields)
{
  _timer["Particles::RKCalcs"].start();
  PARTICLES_PARALLEL_LOOP(pr)
  {
    RKStep(& (*pr), dt, 1.0, bssn_fields);
  }
  _timer["Particles::RKCalcs"].stop();
}

void Particles::RK4Step(map_t & bssn_fields)
{
  _timer["Particles::RKCalcs"].start();
  PARTICLES_PARALLEL_LOOP(pr)
  {
    RKStep(& (*pr), dt/2.0, 1.0, bssn_fields);
  }
  _timer["Particles::RKCalcs"].stop();
}

void Particles::stepInit(map_t & bssn_fields)
{
  _timer["Particles::RKCalcs"].start();
  PARTICLES_PARALLEL_LOOP(pr)
  {
    pr->p_c = pr->p_p;
    pr->p_a = pr->p_p;
    pr->p_f = {0};
  }
  _timer["Particles::RKCalcs"].stop();
}

void Particles::regSwap_c_a()
{
  _timer["Particles::RKCalcs"].start();
  PARTICLES_PARALLEL_LOOP(pr)
  {
    std::swap(pr->p_c, pr->p_a);
  }
  _timer["Particles::RKCalcs"].stop();
}

void Particles::stepTerm()
{
  _timer["Particles::RKCalcs"].start();
  PARTICLES_PARALLEL_LOOP(pr)
  {
    Particle<real_t> & p_p = pr->p_p;
    Particle<real_t> & p_f = pr->p_f;

    for(int i=0; i<3; i++)
    {
      p_p.X[i] = p_f.X[i]/3.0 - 2.0/3.0*p_p.X[i];
      p_p.U[i] = p_f.U[i]/3.0 - 2.0/3.0*p_p.U[i];
    }
  }
  _timer["Particles::RKCalcs"].stop();
}

real_t Particles::getKernelWeight(real_t r, real_t r_s)
{
  real_t weight = 0.0;
  real_t q = r/r_s;
  if(r > 2.0*r_s)
  {
    weight = 0;
  }
  else if(r > r_s)
  {
    weight = 1.0/PI/std::pow(r_s,3.0)*(1.0/4.0)*std::pow((2.0-q), 3.0);
  }
  else
  {
    weight = 1.0/PI/std::pow(r_s,3.0)*(1.0 - 3.0/2.0*pw2(q) + 3.0/4.0*std::pow(q, 3.0) );
  }

  return weight;
}

/**
 * Set bssn _a source registers
 * using data from particle _c register
 */
void Particles::addParticlesToBSSNSrc(BSSN * bssnSim)
{
  _timer["Particles::addToBSSNSrc"].start();

  // matter / source fields
  // will always be setting _a register from _a register
  arr_t & DIFFr_a = *bssnSim->fields["DIFFr_a"];
  arr_t & DIFFS_a = *bssnSim->fields["DIFFS_a"];
  arr_t & S1_a = *bssnSim->fields["S1_a"];
  arr_t & S2_a = *bssnSim->fields["S2_a"];
  arr_t & S3_a = *bssnSim->fields["S3_a"];
  arr_t & STF11_a = *bssnSim->fields["STF11_a"];
  arr_t & STF12_a = *bssnSim->fields["STF12_a"];
  arr_t & STF13_a = *bssnSim->fields["STF13_a"];
  arr_t & STF22_a = *bssnSim->fields["STF22_a"];
  arr_t & STF23_a = *bssnSim->fields["STF23_a"];
  arr_t & STF33_a = *bssnSim->fields["STF33_a"];

  // smoothing radius
  real_t r_s = std::stod(_config("smoothing_radius", "1.5")); // units of dx

  PARTICLES_PARALLEL_LOOP(pr)
  {
    Particle<real_t> & p_a = pr->p_a;
    ParticleMetricPrimitives<real_t> pp_a = getInterpolatedPrimitivesIncomplete(& p_a, bssnSim->fields);

    real_t W = std::sqrt( 1.0 + 
        pp_a.gi[aIDX(1,1)]*p_a.U[0]*p_a.U[0] + pp_a.gi[aIDX(2,2)]*p_a.U[1]*p_a.U[1] + pp_a.gi[aIDX(3,3)]*p_a.U[2]*p_a.U[2]
        + 2.0*( pp_a.gi[aIDX(1,2)]*p_a.U[0]*p_a.U[1] + pp_a.gi[aIDX(1,3)]*p_a.U[0]*p_a.U[2] + pp_a.gi[aIDX(2,3)]*p_a.U[1]*p_a.U[2] )
      );
    real_t MnA = p_a.M / W / dx/dx/dx / pp_a.rootdetg;

    real_t rho = MnA*W*W;
    real_t S = rho - MnA;
    real_t S_1 = MnA*W*p_a.U[0];
    real_t S_2 = MnA*W*p_a.U[1];
    real_t S_3 = MnA*W*p_a.U[2];

    // cubic interpolant with kernel of characteristic "softening" radius r_s and maximum width w_k
    // eg, Eq. 12.2: http://www.ita.uni-heidelberg.de/~dullemond/lectures/num_fluid_2011/Chapter_12.pdf
    real_t w_k = r_s*2.0;
    idx_t w_idx = (idx_t) (w_k + 1.0);
    
    // distribute mass to nearby gridpoints
    idx_t x_idx = (idx_t) (p_a.X[0]/dx);
    idx_t y_idx = (idx_t) (p_a.X[1]/dx);
    idx_t z_idx = (idx_t) (p_a.X[2]/dx);

    // ensure conservation of mass
    real_t total_weight = 0.0;
    for(idx_t x=x_idx-w_idx; x<=x_idx+w_idx+1; ++x)
      for(idx_t y=y_idx-w_idx; y<=y_idx+w_idx+1; ++y)
        for(idx_t z=z_idx-w_idx; z<=z_idx+w_idx+1; ++z)
        {
          real_t r = std::sqrt( pw2(x - p_a.X[0]/dx) + pw2(y - p_a.X[1]/dx) + pw2(z - p_a.X[2]/dx) );
          total_weight += getKernelWeight(r, r_s);
        }

#   pragma omp critical
    {
      for(idx_t x=x_idx-w_idx; x<=x_idx+w_idx+1; ++x)
        for(idx_t y=y_idx-w_idx; y<=y_idx+w_idx+1; ++y)
          for(idx_t z=z_idx-w_idx; z<=z_idx+w_idx+1; ++z)
      {
        idx_t idx = NP_INDEX( idx_t_mod(x,NX), idx_t_mod(y,NY), idx_t_mod(z,NZ) );

        real_t r = std::sqrt( pw2(x - p_a.X[0]/dx) + pw2(y - p_a.X[1]/dx) + pw2(z - p_a.X[2]/dx) );
        real_t weight = getKernelWeight(r, r_s) / total_weight;

        DIFFr_a[idx] += weight*rho;
        DIFFS_a[idx] += weight*S;

        S1_a[idx] += weight*S_1;
        S2_a[idx] += weight*S_2;
        S3_a[idx] += weight*S_3;

        // Not yet trace-free
        STF11_a[idx] += weight*(MnA*p_a.U[0]*p_a.U[0]);
        STF12_a[idx] += weight*(MnA*p_a.U[0]*p_a.U[1]);
        STF13_a[idx] += weight*(MnA*p_a.U[0]*p_a.U[2]);
        STF22_a[idx] += weight*(MnA*p_a.U[1]*p_a.U[1]);
        STF23_a[idx] += weight*(MnA*p_a.U[1]*p_a.U[2]);
        STF33_a[idx] += weight*(MnA*p_a.U[2]*p_a.U[2]);
      }
    } // end loop over nearby indexes 
  } // end particles loop

  // ensure STF is trace-free
  idx_t i, j, k;
# pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i, j, k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    BSSNData bd = {0};
    bssnSim->set_bd_values(i, j, k, &bd);
    real_t trS = exp(-4.0*bd.phi)*(
        STF11_a[idx]*bd.gammai11 + STF22_a[idx]*bd.gammai22 + STF33_a[idx]*bd.gammai33
        + 2.0*(STF12_a[idx]*bd.gammai12 + STF13_a[idx]*bd.gammai13 + STF23_a[idx]*bd.gammai23)
      );

    STF11_a[idx] -= (1.0/3.0)*exp(4.0*bd.phi)*bd.gamma11*trS;
    STF12_a[idx] -= (1.0/3.0)*exp(4.0*bd.phi)*bd.gamma12*trS;
    STF13_a[idx] -= (1.0/3.0)*exp(4.0*bd.phi)*bd.gamma13*trS;
    STF22_a[idx] -= (1.0/3.0)*exp(4.0*bd.phi)*bd.gamma22*trS;
    STF23_a[idx] -= (1.0/3.0)*exp(4.0*bd.phi)*bd.gamma23*trS;
    STF33_a[idx] -= (1.0/3.0)*exp(4.0*bd.phi)*bd.gamma33*trS;
  }


  _timer["Particles::addToBSSNSrc"].stop();
}

} /* namespace cosmo */
