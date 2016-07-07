# Simulation resolution
if(DEFINED COSMO_N)
  add_definitions(-DN=${COSMO_N})
  message(STATUS "${Cyan}Setting N=${COSMO_N}.${ColorReset}")
endif()

# Simulation resolution in X-direction (overrides N)
if(DEFINED COSMO_NX)
  add_definitions(-DNX=${COSMO_NX})
  message(STATUS "${Cyan}Setting NX=${COSMO_NX}.${ColorReset}")
endif()

# Simulation resolution in Y-direction (overrides N)
if(DEFINED COSMO_NY)
  add_definitions(-DNY=${COSMO_NY})
  message(STATUS "${Cyan}Setting NY=${COSMO_NY}.${ColorReset}")
endif()

# Simulation resolution in Z-direction (overrides N)
if(DEFINED COSMO_NZ)
  add_definitions(-DNZ=${COSMO_NZ})
  message(STATUS "${Cyan}Setting NZ=${COSMO_NZ}.${ColorReset}")
endif()

# Simulation box length L (so L = N*dx)
if(DEFINED COSMO_H_LEN_FRAC)
  add_definitions(-DH_LEN_FRAC=${COSMO_H_LEN_FRAC})
  message(STATUS "${Cyan}Setting H_LEN_FRAC=${COSMO_H_LEN_FRAC}.${ColorReset}")
endif()

# Use reference integrator? (required for lcdm sims; disallowed for scalar and particle sims)
if(DEFINED COSMO_USE_REFERENCE_FRW)
  add_definitions(-DUSE_REFERENCE_FRW=${COSMO_USE_REFERENCE_FRW})
  message(STATUS "${Cyan}Setting USE_REFERENCE_FRW=${COSMO_USE_REFERENCE_FRW}.${ColorReset}")
endif()

# Set stencil order? Must be 2, 4, 6, 8.
if(DEFINED COSMO_STENCIL_ORDER)
  add_definitions(-DSTENCIL_ORDER=${COSMO_STENCIL_ORDER})
  message(STATUS "${Cyan}Setting STENCIL_ORDER=${COSMO_STENCIL_ORDER}.${ColorReset}")
endif()

# Use Harmonic gauge? (Disallowed for dust sims)
if(DEFINED COSMO_USE_HARMONIC_ALPHA)
  add_definitions(-DUSE_HARMONIC_ALPHA=${COSMO_USE_HARMONIC_ALPHA})
  message(STATUS "${Cyan}Setting USE_HARMONIC_ALPHA=${COSMO_USE_HARMONIC_ALPHA}.${ColorReset}")
endif()

# Use BSSN shift? Required for scalar and particle sims.
if(DEFINED COSMO_USE_BSSN_SHIFT)
  add_definitions(-DUSE_BSSN_SHIFT=${COSMO_USE_BSSN_SHIFT})
  message(STATUS "${Cyan}Setting USE_BSSN_SHIFT=${COSMO_USE_BSSN_SHIFT}.${ColorReset}")
endif()

# Normalize conformal metric and trace-free extrinsic curvature?
if(DEFINED COSMO_NORMALIZE_GAMMAIJ_AIJ)
  add_definitions(-DNORMALIZE_GAMMAIJ_AIJ=${COSMO_NORMALIZE_GAMMAIJ_AIJ})
  message(STATUS "${Cyan}Setting NORMALIZE_GAMMAIJ_AIJ=${COSMO_NORMALIZE_GAMMAIJ_AIJ}.${ColorReset}")
endif()

# Exclude second-order terms close to zero?
if(DEFINED COSMO_EXCLUDE_SECOND_ORDER_SMALL)
  add_definitions(-DEXCLUDE_SECOND_ORDER_SMALL=${COSMO_EXCLUDE_SECOND_ORDER_SMALL})
  message(STATUS "${Cyan}Setting EXCLUDE_SECOND_ORDER_SMALL=${COSMO_EXCLUDE_SECOND_ORDER_SMALL}.${ColorReset}")
endif()

# Exclude second-order terms around FRW expansion?
if(DEFINED COSMO_EXCLUDE_SECOND_ORDER_FRW)
  add_definitions(-DEXCLUDE_SECOND_ORDER_FRW=${COSMO_EXCLUDE_SECOND_ORDER_FRW})
  message(STATUS "${Cyan}Setting EXCLUDE_SECOND_ORDER_FRW=${COSMO_EXCLUDE_SECOND_ORDER_FRW}.${ColorReset}")
endif()

# Remove these from cache
unset(COSMO_N CACHE)
unset(COSMO_NX CACHE)
unset(COSMO_NY CACHE)
unset(COSMO_NZ CACHE)
unset(COSMO_H_LEN_FRAC CACHE)
unset(COSMO_USE_REFERENCE_FRW CACHE)
unset(COSMO_STENCIL_ORDER CACHE)
unset(COSMO_USE_HARMONIC_ALPHA CACHE)
unset(COSMO_USE_BSSN_SHIFT CACHE)
unset(COSMO_NORMALIZE_GAMMAIJ_AIJ CACHE)
unset(COSMO_EXCLUDE_SECOND_ORDER_SMALL CACHE)
unset(COSMO_EXCLUDE_SECOND_ORDER_FRW CACHE)
