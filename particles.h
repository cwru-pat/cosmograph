#ifndef COSMO_PARTICLES
#define COSMO_PARTICLES

#include "cosmo.h"

#define PARTICLES_INTERPOLATION(var) \
  linearInterpolation( \
    corner_pp_in[0][0][0].var, corner_pp_in[0][0][1].var, corner_pp_in[0][1][0].var, corner_pp_in[0][1][1].var, \
    corner_pp_in[1][0][0].var, corner_pp_in[1][0][1].var, corner_pp_in[1][1][0].var, corner_pp_in[1][1][1].var, \
    x_d \
  )

#define PARTICLES_ROUND(val) ((IT)( (val) + 0.5))

#define DER(field) derivative(x_idx, y_idx, z_idx, a+1, field)


namespace cosmo
{

/**
 * @brief Data structure for particles
 * @details Data structure storing particle position, velocity, and mass
 * 
 * @tparam RT simulation real type
 */
template<typename RT>
struct Particle {
  RT X[3];     // particle position
  RT U[3];     // Particle velocity (covariant / lowered)
  RT M;        // Particle Mass
};

/**
 * @brief Data structure containing 4 needed RK4 registers
 * @details Data structure used for an RK4 integration; contains
 * registers p_x for a Particle<RT> : _p, _a, _c, _f
 * 
 * @tparam RT simulation real type
 */
template<typename RT>
struct ParticleRegister {
  // 4 standard "registers" for each particle
  Particle<RT> p_p;
  Particle<RT> p_a;
  Particle<RT> p_c;
  Particle<RT> p_f;
};

/**
 * @brief Data structure for storing metric quantities ("primitives") at an
 * arbitrary point
 * @details Data structure for storing metric quantities: lapse, shift, inverse
 * 3-metric (full metric, not conformal!); and derivatives of the lapse, shift,
 * and the 3-metric.
 * 
 * @tparam RT simulation real type
 */
template<typename RT>
struct ParticleMetricPrimitives {
  // metric components
  RT rootdetg;      // rooted metric determinant
  RT alpha;         // lapse
  RT beta[3];       // shift
  RT dalpha[3];     // derivative of lapse
  RT dbeta[3][3];   // derivative of shift
  RT gi[6];         // inverse spatial metric
  RT dgi[3][6];     // Inverse metric derivative
};

/**
 * @brief Class for evolving non-interacting matter particles
 * @details Class stores a (private) vector of particles, and contains
 * routines for particle evolution, interpolation, and interactions with
 * GR (BSSN) fields.
 * 
 * @tparam RT simulation real type
 * @tparam IT simulation index type
 */
template <typename RT, typename IT>
class Particles
{
  // list of particle registers
  std::vector<ParticleRegister<RT>> particles;

public:
  
  Particles()
  {
    // anything to do?
  }

  ~Particles()
  {
    // anything to do?
  }

  /**
   * @brief Create particles
   * @details Add particles with randomized positions, zero velocity,
   * and unit mass to internal particles vector.
   * 
   * @param n_particles number of particles to create
   */
  void init(IT n_particles)
  {
    IT i = 0;

    std::mt19937 gen(7.0);
    std::uniform_real_distribution<RT> dist(0.0, 1.0);
    dist(gen);

    for(i=0; i<n_particles; i++)
    {
      Particle<RT> particle = {0};

      // Randomized position
      particle.X[0] = dist(gen)*N*dx;
      particle.X[1] = dist(gen)*N*dx;
      particle.X[2] = dist(gen)*N*dx;
      // Mass in untis TBD
      particle.M = 1.0;

      ParticleRegister<RT> particle_register;
      particle_register.p_p = particle;
      particle_register.p_a = particle;
      particle_register.p_c = particle;
      particle_register.p_f = particle;

      particles.push_back( particle_register );
    }
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
  void setX_d(RT X[3], RT x_d[3])
  {
    for(IT i=0; i<3; i++)
    {
      if(X[i] < 0.0)
      {
        x_d[i] = 1.0 - std::fmod(-X[i], dx)/dx;
      }
      else
      {
        x_d[i] = std::fmod(X[i], dx)/dx;
      }
    }
  }

  /**
   * @brief Returns the index nearest to a point
   * @details Returns the nearest index to a point in a particular dimension.
   * Assumes periodic boundary conditions.
   * 
   * @param X Point (3-vector) of interest.
   * @param dir dimension of interest (0=x, 1=y, 2=z)
   * @return index associated with position
   */
  IT getINDEX(RT X[3], IT dir)
  {
    if(X[dir] < 0)
      return N - (-1.0*PARTICLES_ROUND(X[dir]/dx) % N);

    return (PARTICLES_ROUND(X[dir]/dx) % N);
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
  ParticleMetricPrimitives<RT> interpolatePrimitivesFromCorners(
    ParticleMetricPrimitives<RT> corner_pp_in[2][2][2], RT x_d[3])
  {
    ParticleMetricPrimitives<RT> pp = {0};

    pp.rootdetg = PARTICLES_INTERPOLATION(rootdetg);
    pp.alpha = PARTICLES_INTERPOLATION(alpha);

    for(IT i=0; i<3; i++)
    {
      pp.beta[3] = PARTICLES_INTERPOLATION(beta[3]);
      pp.dalpha[3] = PARTICLES_INTERPOLATION(dalpha[3]);
      for(IT j=0; j<3; j++)
      {
        pp.dbeta[3][3] = PARTICLES_INTERPOLATION(dbeta[3][3]);
      }
    }

    for(IT i=0; i<6; i++)
    {
      pp.gi[6] = PARTICLES_INTERPOLATION(gi[6]);
      for(IT j=0; j<3; j++)
      {
        pp.dgi[3][6] = PARTICLES_INTERPOLATION(dgi[3][6]);
      }
    }
  }

  /**
   * Only set some metric components:
   * gamma_ij, rootdetg, and alpha
   * see also: getInterpolatedPrimitivesIncomplete
   */
  ParticleMetricPrimitives<RT> interpolatePrimitivesFromCornersIncomplete(
    ParticleMetricPrimitives<RT> corner_pp_in[2][2][2], RT x_d[3])
  {
    ParticleMetricPrimitives<RT> pp = {0};
    pp.rootdetg = PARTICLES_INTERPOLATION(rootdetg);
    pp.alpha = PARTICLES_INTERPOLATION(alpha);
    for(IT i=0; i<6; i++)
    {
      pp.gi[6] = PARTICLES_INTERPOLATION(gi[6]);
    }
  }

  RT linearInterpolation(
    RT C000, RT C001, RT C010, RT C011, /* values at "C"orners, C_{x,y,z} */
    RT C100, RT C101, RT C110, RT C111, /* "binary" order */
    RT x_d[3] /* normalized position within cube */
    )
  {
    RT C00 = C000*(1.0 - x_d[0]) + C100*x_d[0];
    RT C01 = C001*(1.0 - x_d[0]) + C101*x_d[0];
    RT C10 = C010*(1.0 - x_d[0]) + C110*x_d[0];
    RT C11 = C011*(1.0 - x_d[0]) + C111*x_d[0];

    RT C0 = C00*(1.0 - x_d[1]) + C10*x_d[1];
    RT C1 = C01*(1.0 - x_d[1]) + C11*x_d[1];

    // return "C"
    return C0*(1.0 - x_d[2]) + C1*x_d[2];
  }

  /**
   * Only set some metric components:
   * gamma_ij, rootdetg, and alpha
   * see also: interpolatePrimitivesFromCornersIncomplete
   * uses _c register.
   */
  ParticleMetricPrimitives<RT> getInterpolatedPrimitivesIncomplete(Particle<RT> * p,
    std::map <std::string, RT *> * bssn_fields)
  {
    RT * const DIFFalpha_a = bssn_fields["DIFFalpha_c"];
    RT * const DIFFphi_a = bssn_fields["DIFFphi_c"];
    RT * const DIFFgamma11_a = bssn_fields["DIFFgamma11_c"];
    RT * const DIFFgamma12_a = bssn_fields["DIFFgamma12_c"];
    RT * const DIFFgamma13_a = bssn_fields["DIFFgamma13_c"];
    RT * const DIFFgamma22_a = bssn_fields["DIFFgamma22_c"];
    RT * const DIFFgamma23_a = bssn_fields["DIFFgamma23_c"];
    RT * const DIFFgamma33_a = bssn_fields["DIFFgamma33_c"];

    // NGP interpolant
    IT x_idx = getINDEX(p.X, 0);
    IT y_idx = getINDEX(p.X, 1);
    IT z_idx = getINDEX(p.X, 2);

    ParticleMetricPrimitives<RT> corner_pp[2][2][2];
    for(IT i=0; i<2; i++)
      for(IT j=0; j<2; j++)
        for(IT k=0; k<2; k++)
        {
          ParticleMetricPrimitives<RT> pp;
          IT idx = INDEX(x_idx + i, y_idx + j, z_idx + k);

          pp.rootdetg = std::exp(6.0*DIFFphi_a[idx]);
          pp.alpha = DIFFalpha_a[idx];

          pp.gi[aIDX(1,1)] = std::exp(-4.0*DIFFphi_a[idx])*(1.0 + DIFFgamma22_a[idx] + DIFFgamma33_a[idx] - pw2(DIFFgamma23_a[idx]) + DIFFgamma22_a[idx]*DIFFgamma33_a[idx]);
          pp.gi[aIDX(2,2)] = std::exp(-4.0*DIFFphi_a[idx])*(1.0 + DIFFgamma11_a[idx] + DIFFgamma33_a[idx] - pw2(DIFFgamma13_a[idx]) + DIFFgamma11_a[idx]*DIFFgamma33_a[idx]);
          pp.gi[aIDX(3,3)] = std::exp(-4.0*DIFFphi_a[idx])*(1.0 + DIFFgamma11_a[idx] + DIFFgamma22_a[idx] - pw2(DIFFgamma12_a[idx]) + DIFFgamma11_a[idx]*DIFFgamma22_a[idx]);
          pp.gi[aIDX(1,2)] = std::exp(-4.0*DIFFphi_a[idx])*(DIFFgamma13_a[idx]*DIFFgamma23_a[idx] - DIFFgamma12_a[idx]*(1.0 + DIFFgamma33_a[idx]));
          pp.gi[aIDX(1,3)] = std::exp(-4.0*DIFFphi_a[idx])*(DIFFgamma12_a[idx]*DIFFgamma23_a[idx] - DIFFgamma13_a[idx]*(1.0 + DIFFgamma22_a[idx]));
          pp.gi[aIDX(2,3)] = std::exp(-4.0*DIFFphi_a[idx])*(DIFFgamma12_a[idx]*DIFFgamma13_a[idx] - DIFFgamma23_a[idx]*(1.0 + DIFFgamma11_a[idx]));
        }
    /* end i, j, k loops */

    RT x_d[3];
    setX_d(p.X, x_d);
    return interpolatePrimitivesFromCornersIncomplete(&corner_pp, &x_d);
  }

  ParticleMetricPrimitives<RT> getInterpolatedPrimitives(Particle<RT> * p,
    std::map <std::string, RT *> * bssn_fields)
  {
    RT * const DIFFalpha_a = bssn_fields["DIFFalpha_a"];
    
    RT * const beta1_a = bssn_fields["beta1_a"];
    RT * const beta2_a = bssn_fields["beta2_a"];
    RT * const beta3_a = bssn_fields["beta3_a"];

    RT * const DIFFphi_a = bssn_fields["DIFFphi_a"];

    RT * const DIFFgamma11_a = bssn_fields["DIFFgamma11_a"];
    RT * const DIFFgamma12_a = bssn_fields["DIFFgamma12_a"];
    RT * const DIFFgamma13_a = bssn_fields["DIFFgamma13_a"];
    RT * const DIFFgamma22_a = bssn_fields["DIFFgamma22_a"];
    RT * const DIFFgamma23_a = bssn_fields["DIFFgamma23_a"];
    RT * const DIFFgamma33_a = bssn_fields["DIFFgamma33_a"];

    // NGP interpolant
    IT x_idx = getINDEX(p.X, 0);
    IT y_idx = getINDEX(p.X, 1);
    IT z_idx = getINDEX(p.X, 2);

    ParticleMetricPrimitives<RT> corner_pp[2][2][2];
    for(IT i=0; i<2; i++)
      for(IT j=0; j<2; j++)
        for(IT k=0; k<2; k++)
        {
          ParticleMetricPrimitives<RT> pp;
          IT idx = INDEX(x_idx + i, y_idx + j, z_idx + k);

          pp.rootdetg = std::exp(6.0*DIFFphi_a[idx]);
          pp.alpha = DIFFalpha_a[idx];

          pp.beta[0] = beta1_a[idx];
          pp.beta[1] = beta2_a[idx];
          pp.beta[2] = beta3_a[idx];

          pp.gi[aIDX(1,1)] = std::exp(-4.0*DIFFphi_a[idx])*(1.0 + DIFFgamma22_a[idx] + DIFFgamma33_a[idx] - pw2(DIFFgamma23_a[idx]) + DIFFgamma22_a[idx]*DIFFgamma33_a[idx]);
          pp.gi[aIDX(2,2)] = std::exp(-4.0*DIFFphi_a[idx])*(1.0 + DIFFgamma11_a[idx] + DIFFgamma33_a[idx] - pw2(DIFFgamma13_a[idx]) + DIFFgamma11_a[idx]*DIFFgamma33_a[idx]);
          pp.gi[aIDX(3,3)] = std::exp(-4.0*DIFFphi_a[idx])*(1.0 + DIFFgamma11_a[idx] + DIFFgamma22_a[idx] - pw2(DIFFgamma12_a[idx]) + DIFFgamma11_a[idx]*DIFFgamma22_a[idx]);
          pp.gi[aIDX(1,2)] = std::exp(-4.0*DIFFphi_a[idx])*(DIFFgamma13_a[idx]*DIFFgamma23_a[idx] - DIFFgamma12_a[idx]*(1.0 + DIFFgamma33_a[idx]));
          pp.gi[aIDX(1,3)] = std::exp(-4.0*DIFFphi_a[idx])*(DIFFgamma12_a[idx]*DIFFgamma23_a[idx] - DIFFgamma13_a[idx]*(1.0 + DIFFgamma22_a[idx]));
          pp.gi[aIDX(2,3)] = std::exp(-4.0*DIFFphi_a[idx])*(DIFFgamma12_a[idx]*DIFFgamma13_a[idx] - DIFFgamma23_a[idx]*(1.0 + DIFFgamma11_a[idx]));

          for(IT a=0; a<3; a++)
          {
            pp.dalpha[a] = DER(DIFFalpha_a);

            pp.dbeta[a][0] = DER(beta1_a);
            pp.dbeta[a][1] = DER(beta2_a);
            pp.dbeta[a][2] = DER(beta3_a);

            pp.dgi[a][6];

            pp.dgi[a][aIDX(1,1)] = -4.0*DER(DIFFphi_a)*pp.gi[aIDX(1,1)]
              + std::exp(-4.0*DIFFphi_a[idx])*(DER(DIFFgamma22_a) + DER(DIFFgamma33_a) - 2.0*DER(DIFFgamma23_a) + DER(DIFFgamma22_a)*DIFFgamma33_a[idx] + DIFFgamma22_a[idx]*DER(DIFFgamma33_a));
            pp.dgi[a][aIDX(2,2)] = -4.0*DER(DIFFphi_a)*pp.gi[aIDX(2,2)]
              + std::exp(-4.0*DIFFphi_a[idx])*(DER(DIFFgamma11_a) + DER(DIFFgamma33_a) - 2.0*DER(DIFFgamma13_a) + DER(DIFFgamma11_a)*DIFFgamma33_a[idx] + DIFFgamma11_a[idx]*DER(DIFFgamma33_a));
            pp.dgi[a][aIDX(3,3)] = -4.0*DER(DIFFphi_a)*pp.gi[aIDX(3,3)]
              + std::exp(-4.0*DIFFphi_a[idx])*(DER(DIFFgamma11_a) + DER(DIFFgamma22_a) - 2.0*DER(DIFFgamma12_a) + DER(DIFFgamma11_a)*DIFFgamma22_a[idx] + DIFFgamma11_a[idx]*DER(DIFFgamma22_a));
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

    RT x_d[3];
    setX_d(p.X, x_d);
    return interpolatePrimitivesFromCorners(&corner_pp, &x_d);
  }

  /**
   * RK4 implementation:
   * y_p; y_a, y_c, y_f = 0 (TODO)
   * y_a = y_p
   * 
   * y_c = y_p + dt/2 * f(y_a)    // RK1 step
   * y_f += y_c
   * y_a <-> y_c
   * 
   * y_c = y_p + dt/2 * f(y_a)    // RK2 step
   * y_f += 2 y_c
   * y_a <-> y_c
   * 
   * y_c = y_p + dt * f(y_a)      // RK3 step
   * y_f += y_c
   * y_a <-> y_c
   * 
   * y_c = y_p + dt/2 * f(y_a)    // RK4 step
   * y_f += y_c
   * y_a <-> y_c
   * 
   * (y_f = 5y_p + 1/2 K1 + K2 + K3 + 1/2 K4)
   * y_f = y_f / 3 - 2/3 y_p      // finalize
   * y_f <-> y_p
   * 
   * Using Baumgarte & Shapiro, 5.223-5.225
   *
   * For a single RK step, we can just do: y_c = y_p + h * f(y_a)
   * and y_f += RK_sum_coeff * y_c
   * */
  void RKStep(ParticleRegister<RT> * pr, RT h, RT RK_sum_coeff,
    std::map <std::string, RT *> * bssn_fields)
  {
    Particle<RT> & p_p = pr->p_p;
    Particle<RT> & p_a = pr->p_a;
    Particle<RT> & p_c = pr->p_c;
    Particle<RT> & p_f = pr->p_f;

    ParticleMetricPrimitives<RT> pp_a = getInterpolatedPrimitives(& p_a, bssn_fields);

    RT W = std::sqrt( 1.0 + 
        pp_a.gi[aIDX(1,1)]*p_a->U[0]*p_a->U[0] + pp_a.gi[aIDX(2,2)]*p_a->U[1]*p_a->U[1] + pp_a.gi[aIDX(3,3)]*p_a->U[2]*p_a->U[2]
        + 2.0*( pp_a.gi[aIDX(1,2)]*p_a->U[0]*p_a->U[1] + pp_a.gi[aIDX(1,3)]*p_a->U[0]*p_a->U[2] + pp_a.gi[aIDX(2,3)]*p_a->U[1]*p_a->U[2] )
      );
    RT U0 = W/pp_a.alpha;

    for(int i=0; i<3; i++)
    {
      p_c->X[i] = p_p->X[i] + h*(pp_a.gi[aIDX(i,1)]*p_a->U[1] + pp_a.gi[aIDX(i,2)]*p_a->U[2] + pp_a.gi[aIDX(i,3)]*p_a->U[3] - pp_a.beta[i]);
      p_f->X[i] += RK_sum_coeff*p_c->X[i];

      p_c->U[i] = p_p->U[i] + h*(
        -1.0*W*pp_a.dalpha[i] + p_a->U[1]*pp_a.dbeta[i][1] + p_a->U[2]*pp_a.dbeta[i][2] + p_a->U[3]*pp_a.dbeta[i][3]
        -1.0/2.0/U0*(
          pp_a.dgi[i][aIDX(1,1)]*p_a->U[0]*p_a->U[0] + pp_a.dgi[i][aIDX(2,2)]*p_a->U[1]*p_a->U[1] + pp_a.dgi[i][aIDX(3,3)]*p_a->U[2]*p_a->U[2]
          + 2.0*( pp_a.dgi[i][aIDX(1,2)]*p_a->U[0]*p_a->U[1] + pp_a.dgi[i][aIDX(1,3)]*p_a->U[0]*p_a->U[2] + pp_a.dgi[i][aIDX(2,3)]*p_a->U[1]*p_a->U[2] )
        )
      );
      p_f->U[i] += RK_sum_coeff*p_c->U[i];
    }

// TODO: careful about pragma with this statement (racey):
    addParticleToBSSNSrc(& pr->p_c, bssn_fields);
    std::swap(pr->p_c, pr->p_a);
  }

  void RK1Step(std::map <std::string, RT *> * bssn_fields)
  {
    typename std::vector<ParticleRegister<RT>>::iterator pr;
// TODO... #pragma
    for(pr = particles.begin(); pr < particles.end(); ++pr)
    {
      RKStep(pr, dt/2.0, 1.0, bssn_fields);
    } 
  }

  void RK2Step(std::map <std::string, RT *> * bssn_fields)
  {
    typename std::vector<ParticleRegister<RT>>::iterator pr;
// TODO... #pragma 
    for(pr = particles.begin(); pr < particles.end(); ++pr)
    {
      RKStep(pr, dt/2.0, 2.0, bssn_fields);
    } 
  }

  void RK3Step(std::map <std::string, RT *> * bssn_fields)
  {
    typename std::vector<ParticleRegister<RT>>::iterator pr;
// TODO... #pragma 
    for(pr = particles.begin(); pr < particles.end(); ++pr)
    {
      RKStep(pr, dt, 1.0, bssn_fields);
    } 
  }

  void RK4Step(std::map <std::string, RT *> * bssn_fields)
  {
    typename std::vector<ParticleRegister<RT>>::iterator pr;
// TODO... #pragma 
    for(pr = particles.begin(); pr < particles.end(); ++pr)
    {
      RKStep(pr, dt/2.0, 1.0, bssn_fields);
    } 
  }


  /**
   * Set bssn _a register
   * using data from particle _c register
   */
  void addParticleToBSSNSrc(Particle<RT> * p_c,
    std::map <std::string, RT *> * bssn_fields)
  {
    // matter / source fields
    // will always be setting _a register from _a register
    RT * const DIFFr_a = bssn_fields["DIFFr_a"];
    RT * const S_a = bssn_fields["S_a"];
    RT * const S1_a = bssn_fields["S1_a"];
    RT * const S2_a = bssn_fields["S2_a"];
    RT * const S3_a = bssn_fields["S3_a"];
    RT * const S11_a = bssn_fields["S11_a"];
    RT * const S12_a = bssn_fields["S12_a"];
    RT * const S13_a = bssn_fields["S13_a"];
    RT * const S22_a = bssn_fields["S22_a"];
    RT * const S23_a = bssn_fields["S23_a"];
    RT * const S33_a = bssn_fields["S33_a"];

    // NGP interpolant
    // should do something better eventually
    IT x_idx = getINDEX(p_c.X, 0);
    IT y_idx = getINDEX(p_c.X, 1);
    IT z_idx = getINDEX(p_c.X, 2);
    IT idx = INDEX(x_idx, y_idx, z_idx);

    ParticleMetricPrimitives<RT> pp_c = getInterpolatedPrimitives(& p_c, bssn_fields);

    RT W = std::sqrt( 1.0 + 
        pp_c.gi[aIDX(1,1)]*p_c->U[0]*p_c->U[0] + pp_c.gi[aIDX(2,2)]*p_c->U[1]*p_c->U[1] + pp_c.gi[aIDX(3,3)]*p_c->U[2]*p_c->U[2]
        + 2.0*( pp_c.gi[aIDX(1,2)]*p_c->U[0]*p_c->U[1] + pp_c.gi[aIDX(1,3)]*p_c->U[0]*p_c->U[2] + pp_c.gi[aIDX(2,3)]*p_c->U[1]*p_c->U[2] )
      );

    RT nA = 1.0 / W / dx/dx/dx / std::sqrt(pp_c.rootdetg);

    // Eq . 5.226 in Baumgarte & Shapiro
    DIFFr_a[idx] += p_c.M*nA*W*W;
    S_a[idx] += p_c.M*nA*W*W - p_c.M*nA;
    S1_a[idx] += p_c.M*nA*W*p_c.U[0];
    S2_a[idx] += p_c.M*nA*W*p_c.U[1];
    S3_a[idx] += p_c.M*nA*W*p_c.U[2];
    S11_a[idx] += p_c.M*nA*p_c.U[0]*p_c.U[0];
    S12_a[idx] += p_c.M*nA*p_c.U[0]*p_c.U[1];
    S13_a[idx] += p_c.M*nA*p_c.U[0]*p_c.U[2];
    S22_a[idx] += p_c.M*nA*p_c.U[1]*p_c.U[1];
    S23_a[idx] += p_c.M*nA*p_c.U[1]*p_c.U[2];
    S33_a[idx] += p_c.M*nA*p_c.U[2]*p_c.U[2];

  }
};

}

#endif
