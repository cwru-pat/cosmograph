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

  // local copies of current field values
  HYDRO_APPLY_TO_FIELDS(DECLARE_REAL_T)

  // local copies of adjacent current field values for fast derivatives
  HYDRO_APPLY_TO_FIELDS(DECLARE_ADJACENT_REAL_T)


} HydroData;


/** Hydro class **/
class Hydro
{
  /* local values */
  HydroData paq;

  /* Reference to sim of GR Fields */
  BSSN *bssnSim;

  /* Fluid fields */
  HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_CREATE)

public:
  std::map <std::string, real_t *> fields;

  Hydro(BSSN *bssnSimRef)
  {
    HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_ALLOC)
    HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_ADDMAP)

    bssnSim = bssnSimRef;
  }
  ~Hydro()
  {
    HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_DELETE)
  }

  /* set current local field values */
void set_paq_values(idx_t i, idx_t j, idx_t k, HydroData *paq)
  {
    paq->i = i;
    paq->j = j;
    paq->k = k;
    paq->idx = INDEX(i,j,k);

    // draw data from cache
    //set_local_vals(paq);

    // pre-compute re-used quantities

  }

  void WENO_step()
  {
    set_flux_src_terms();
    // actual step
  }

  void set_flux_src_terms()
  {
//    LOOP
    {

    }
  }

  /* set conserved values (constructed from other quantities) */
  void set_bssn_src_vals(idx_t idx)
  {

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
