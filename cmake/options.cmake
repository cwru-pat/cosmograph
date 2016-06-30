
SET(COSMO_N FALSE CACHE STRING "Simulation resolution")
if(COSMO_N)
  add_definitions(-DN=${COSMO_N})
  message(STATUS "${Cyan}Setting N=${COSMO_N}.${ColorReset}")
endif()

SET(COSMO_NX FALSE CACHE STRING "Simulation resolution in X-direction (overrides N)")
if(COSMO_NX)
  add_definitions(-DNX=${COSMO_NX})
  message(STATUS "${Cyan}Setting NX=${COSMO_NX}.${ColorReset}")
endif()

SET(COSMO_NY FALSE CACHE STRING "Simulation resolution in Y-direction (overrides N)")
if(COSMO_NY)
  add_definitions(-DNY=${COSMO_NY})
  message(STATUS "${Cyan}Setting NY=${COSMO_NY}.${ColorReset}")
endif()

SET(COSMO_NZ FALSE CACHE STRING "Simulation resolution in Z-direction (overrides N)")
if(COSMO_NZ)
  add_definitions(-DNZ=${COSMO_NZ})
  message(STATUS "${Cyan}Setting NZ=${COSMO_NZ}.${ColorReset}")
endif()

SET(COSMO_H_LEN_FRAC FALSE CACHE STRING "Simulation box length L (so L = N*dx)")
if(COSMO_H_LEN_FRAC)
  add_definitions(-DH_LEN_FRAC=${COSMO_H_LEN_FRAC})
  message(STATUS "${Cyan}Setting H_LEN_FRAC=${COSMO_H_LEN_FRAC}.${ColorReset}")
endif()

SET(COSMO_USE_REFERENCE_FRW FALSE CACHE STRING
	"Use reference integrator? (required for lcdm sims; disallowed for scalar and particle sims)")
if(COSMO_USE_REFERENCE_FRW)
  add_definitions(-DUSE_REFERENCE_FRW=${COSMO_USE_REFERENCE_FRW})
  message(STATUS "${Cyan}Setting USE_REFERENCE_FRW=${COSMO_USE_REFERENCE_FRW}.${ColorReset}")
endif()

SET(COSMO_STENCIL_ORDER FALSE CACHE STRING "Set stencil order? Must be 2, 4, 6, 8.")
if(COSMO_STENCIL_ORDER)
	add_definitions(-DSTENCIL_ORDER=${COSMO_STENCIL_ORDER})
	message(STATUS "${Cyan}Setting STENCIL_ORDER=${COSMO_STENCIL_ORDER}.${ColorReset}")
endif()

SET(COSMO_USE_HARMONIC_ALPHA FALSE CACHE STRING "Use Harmonic gauge? (Disallowed for dust sims)")
if(COSMO_USE_HARMONIC_ALPHA)
	add_definitions(-DUSE_HARMONIC_ALPHA=${COSMO_USE_HARMONIC_ALPHA})
	message(STATUS "${Cyan}Setting USE_HARMONIC_ALPHA=${COSMO_USE_HARMONIC_ALPHA}.${ColorReset}")
endif()

SET(COSMO_USE_BSSN_SHIFT FALSE CACHE STRING "Use BSSN shift? Required for scalar and particle sims.")
if(COSMO_USE_BSSN_SHIFT)
	add_definitions(-DUSE_BSSN_SHIFT=${COSMO_USE_BSSN_SHIFT})
	message(STATUS "${Cyan}Setting USE_BSSN_SHIFT=${COSMO_USE_BSSN_SHIFT}.${ColorReset}")
endif()

SET(COSMO_NORMALIZE_GAMMAIJ_AIJ FALSE CACHE STRING "Normalize conformal metric and trace-free extrinsic curvature?")
if(COSMO_NORMALIZE_GAMMAIJ_AIJ)
	add_definitions(-DNORMALIZE_GAMMAIJ_AIJ=${COSMO_NORMALIZE_GAMMAIJ_AIJ})
	message(STATUS "${Cyan}Setting NORMALIZE_GAMMAIJ_AIJ=${COSMO_NORMALIZE_GAMMAIJ_AIJ}.${ColorReset}")
endif()

