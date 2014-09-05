#ifndef COSMO_HYDRO
#define COSMO_HYDRO

#include "cosmo.h"
#include "globals.h"

#include "hydro_macros.h"

namespace cosmo
{

class BSSN;

typedef struct {

  idx_t i, j, k, idx;



} HydroData;


/** Hydro class **/
class Hydro
{
  /* local values */
  HydroData paq;

  /* Reference to sim of GR Fields */
  BSSN *bssnSim;

  /* Fluid fields */
  HYDRO_APPLY_TO_FIELDS(RK4_ARRAY_CREATE)

public:
  std::map <std::string, real_t *> fields;

  Hydro(BSSN *bssnSimRef)
  {
    HYDRO_APPLY_TO_FIELDS(RK4_ARRAY_ALLOC)
    HYDRO_APPLY_TO_FIELDS(RK4_ARRAY_ADDMAP)

    bssnSim = bssnSimRef;
  }
  ~Hydro()
  {
    HYDRO_APPLY_TO_FIELDS(RK4_ARRAY_DELETE)
  }

  /* set current local field values */
  void set_local_vals()
  {

  }

  void extract_primitive_variables()
  {
    
  }

  /* set conserved values (constructed from other quantities) */
  void set_flux_src_vals()
  {
    idx_t idx;

    // flux values
    LOOP3(i, j, k)
    {
      // idx = INDEX(i, j, k);
      
      // F0_1[idx] = U0_1[idx]*v_1


    }


  }

  // calculate needed quantities
  void set_paq_values(idx_t i, idx_t j, idx_t k)
  {
    
  }

  void step()
  {

  }

  void init()
  {
    // initialize values

  }

};

}

#endif
