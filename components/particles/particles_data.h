#ifndef COSMO_PARTICLES_DATA
#define COSMO_PARTICLES_DATA

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

typedef std::vector<ParticleRegister<real_t>> particle_vec;

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

}

#endif
