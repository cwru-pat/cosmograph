#ifndef COSMO_PARTICLES
#define COSMO_PARTICLES

#include "../../cosmo_types.h"
#include "../../cosmo_includes.h"
#include "particles_data.h"
#include "../bssn/bssn.h"

namespace cosmo
{

/**
 * @brief Class for evolving non-interacting matter particles
 * @details Class stores a (private) vector of particles, and contains
 * routines for particle evolution, interpolation, and interactions with
 * GR (BSSN) fields.
 */
class Particles
{
  // list of particle registers
  particle_vec * particles;

public:
  
  Particles();
  ~Particles();

  void init(idx_t n_particles);
  void addParticle(Particle<real_t> particle);
  
  particle_vec * getParticleVec();

  real_t getFractionalIndex(real_t x);
  idx_t getIndexBelow(real_t x);
  idx_t getNearestIndex(real_t x);
  void setX_d(real_t X[3], real_t x_d[3]);

  ParticleMetricPrimitives<real_t> interpolatePrimitivesFromCorners(
    ParticleMetricPrimitives<real_t> corner_pp_in[2][2][2], real_t x_d[3]);

  ParticleMetricPrimitives<real_t> interpolatePrimitivesFromCornersIncomplete(
    ParticleMetricPrimitives<real_t> corner_pp_in[2][2][2], real_t x_d[3]);

  real_t linearInterpolation(
    real_t C000, real_t C001, real_t C010, real_t C011, /* values at "C"orners, C_{x,y,z} */
    real_t C100, real_t C101, real_t C110, real_t C111, /* "binary" order */
    real_t x_d[3] /* normalized position within cube */
    );

  ParticleMetricPrimitives<real_t> getInterpolatedPrimitivesIncomplete(
    Particle<real_t> * p, map_t & bssn_fields);

  ParticleMetricPrimitives<real_t> getInterpolatedPrimitives(
    Particle<real_t> * p, map_t & bssn_fields);

  void RKStep(ParticleRegister<real_t> * pr, real_t h, real_t RK_sum_coeff,
    map_t & bssn_fields);

  void RK1Step(map_t & bssn_fields);
  void RK2Step(map_t & bssn_fields);
  void RK3Step(map_t & bssn_fields);
  void RK4Step(map_t & bssn_fields);

  void stepInit(map_t & bssn_fields);
  void regSwap_c_a();
  void stepTerm();

  void addParticlesToBSSNSrc(BSSN * bssnSim);
  void addParticleToBSSNSrc(Particle<real_t> * p_c, map_t & bssn_fields);
};

}

#endif
