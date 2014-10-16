#ifndef COSMO_HYDRO
#define COSMO_HYDRO

#include "cosmo.h"
#include "globals.h"

#include "hydro_macros.h"

namespace cosmo
{

typedef struct {

  idx_t i, j, k, idx;

} HydroData;


/** Hydro class **/
class Hydro
{
  /* local values */
  HydroData paq;

  /* equation of state "w" value */
  real_t w_EOS;

  /* Fluid fields */
  HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_CREATE)
  HYDRO_APPLY_TO_FLUXES(FLUX_ARRAY_CREATE)
  HYDRO_APPLY_TO_SOURCES(GEN1_ARRAY_CREATE)

public:
  std::map <std::string, real_t *> fields;

  Hydro()
  {
    HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_ALLOC)
    HYDRO_APPLY_TO_FLUXES(FLUX_ARRAY_ALLOC)
    HYDRO_APPLY_TO_SOURCES(GEN1_ARRAY_ALLOC)

    HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_ADDMAP)
    HYDRO_APPLY_TO_FLUXES(FLUX_ARRAY_ADDMAP)
    HYDRO_APPLY_TO_SOURCES(GEN1_ARRAY_ADDMAP)

  }
  ~Hydro()
  {
    HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_DELETE)
    HYDRO_APPLY_TO_FLUXES(FLUX_ARRAY_DELETE)
    HYDRO_APPLY_TO_SOURCES(GEN1_ARRAY_DELETE)
  }

  // assumes 
  void WENOStepPt()
  {
    //setPrimitivesPt();
    //set_WENO_fluxes();
    //set_sources();

    // next - perform actual evolution (finite volume method)

  }

  void setPrimitivesPt(cosmo::BSSNData *paq)
  {
    //pq = (*bssnSim).paq;

    // printf("%g poop\n", (*paq).phi);

    //real_t W; /* lorentz factor */
    //real_t g; /* metric determinant */

    //g = exp(4*paq->phi);

    /* W = (gamma^ij S_j S_j / D^2 / (1+w^2) + 1)^1/2 */
    // W = sqrt(
    //     1 + (
    //        paq->gammai11*US1[idx]*US1[idx] + paq->gammai22*US2[idx]*US2[idx] + paq->gammai33*US3[idx]*US3[idx]
    //         + 2*paq->gammai12*US1[idx]*US2[idx] + 2*paq->gammai13*US1[idx]*US3[idx] + 2*paq->gammai23*US2[idx]*US3[idx]
    //       ) / UD[idx] / (1 + w_EOS*w_EOS)
    //     );

    /* \rho = D / */


  }

  void set_WENO_fluxes()
  {
    LOOP3(i, j, k)
    {

    }
  }

  void set_sources()
  {
    LOOP3(i, j, k)
    {

    }
  }

  void init()
  {
    // initialize values

  }

};

}

#endif
