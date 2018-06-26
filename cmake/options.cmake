# Simulation resolution
if(DEFINED COSMO_N)
  add_definitions(-DCOSMO_N=${COSMO_N})
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

# Use Gamma Driver shift?
if(DEFINED COSMO_USE_GAMMA_DRIVER)
  add_definitions(-DUSE_GAMMA_DRIVER=${COSMO_USE_GAMMA_DRIVER})
  message(STATUS "${Cyan}Setting USE_GAMMA_DRIVER=${COSMO_USE_GAMMA_DRIVER}.${ColorReset}")
endif()

# Use BSSN shift? Required for scalar and particle sims.
if(DEFINED COSMO_USE_BSSN_SHIFT)
  add_definitions(-DUSE_BSSN_SHIFT=${COSMO_USE_BSSN_SHIFT})
  message(STATUS "${Cyan}Setting USE_BSSN_SHIFT=${COSMO_USE_BSSN_SHIFT}.${ColorReset}")
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

#set potential type
if(DEFINED COSMO_USE_COSMO_CONST_POTENTIAL)
  add_definitions(-DUSE_COSMO_CONST_POTENTIAL=${COSMO_USE_COSMO_CONST_POTENTIAL})
  message(STATUS "${Cyan}Setting USE_COSMO_CONST_POTENTIAL=${COSMO_USE_COSMO_CONST_POTENTIAL}.${ColorReset}")
endif()

if(DEFINED COSMO_COSMO_CONST)
  add_definitions(-DCOSMO_CONST=${COSMO_COSMO_CONST})
  message(STATUS "${Cyan}Setting COSMO_CONST=${COSMO_COSMO_CONST}.${ColorReset}")
endif()

# Z4c terms?
if(DEFINED COSMO_USE_Z4c_DAMPING)
  add_definitions(-DUSE_Z4c_DAMPING=${COSMO_USE_Z4c_DAMPING})
  message(STATUS "${Cyan}Setting USE_Z4c_DAMPING=${COSMO_USE_Z4c_DAMPING}.${ColorReset}")
endif()

# Generalized Newton gauge? see notes
if(DEFINED COSMO_USE_GENERALIZED_NEWTON)
  add_definitions(-DUSE_Z4c_DAMPING=${COSMO_USE_GENERALIZED_NEWTON})
  message(STATUS "${Cyan}Setting USE_GENERALIZED_NEWTON=${COSMO_USE_GENERALIZED_NEWTON}.${ColorReset}")
endif()


# Remove these from cache
unset(COSMO_N CACHE)
unset(COSMO_NX CACHE)
unset(COSMO_NY CACHE)
unset(COSMO_NZ CACHE)
unset(COSMO_H_LEN_FRAC CACHE)
unset(COSMO_USE_REFERENCE_FRW CACHE)
unset(COSMO_STENCIL_ORDER CACHE)
unset(COSMO_USE_GAMMA_DRIVER CACHE)
unset(COSMO_USE_BSSN_SHIFT CACHE)
unset(COSMO_EXCLUDE_SECOND_ORDER_SMALL CACHE)
unset(COSMO_EXCLUDE_SECOND_ORDER_FRW CACHE)
unset(COSMO_USE_COSMO_CONST_POTENTIAL CACHE)
unset(COSMO_COSMO_CONST CACHE)
unset(COSMO_USE_Z4c_DAMPING CACHE)
