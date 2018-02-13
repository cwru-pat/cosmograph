/** @file particle_ic.h
 * @brief Functions to set initial conditions for the particle class.
 * Functions should be made callable via a config setting in the ParticleSim
 * class.
 */

#ifndef COSMO_PARTICLE_ICS
#define COSMO_PARTICLE_ICS

#include "../bssn/bssn.h"
#include "../../IO/IOData.h"
#include "../../utils/Fourier.h"
#include "particles.h"

namespace cosmo
{

void particle_ic_set_random(BSSN * bssnSim, Particles * particles, Fourier * fourier,
  IOData * iodata);
void particle_ic_set_sinusoid(BSSN * bssnSim, Particles * particles, IOData * iodata);
void particle_ic_set_sinusoid_to_compare(BSSN * bssnSim, Particles * particles, IOData * iodata);
void particle_ic_set_vectorpert(BSSN * bssnSim, Particles * particles, IOData * iodata);

}

#endif
