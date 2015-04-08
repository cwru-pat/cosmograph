#ifndef COSMO_HYDRO
#define COSMO_HYDRO

#include "cosmo.h"
#include "globals.h"

#include "hydro_data.h"
#include "hydro_macros.h"

namespace cosmo
{

/** Hydro class **/
class Hydro
{
  /* Fluid fields */
  HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_CREATE)
  HYDRO_APPLY_TO_FLUXES(FLUX_ARRAY_CREATE)
  HYDRO_APPLY_TO_SOURCES(GEN1_ARRAY_CREATE)
  /* primitives */
  HYDRO_APPLY_TO_PRIMITIVES(GEN1_ARRAY_CREATE)
  HYDRO_APPLY_TO_FLUXES_INT(FLUX_ARRAY_CREATE)

public:
  std::map <std::string, real_t *> fields;

  Hydro();
  ~Hydro();

  void setQuantitiesCell(BSSNData *paq, HydroData *hdp);
  void setFluxesCell(BSSNData *paq, HydroData *hdp);

  void setOneFluxInt(idx_t i, idx_t j, idx_t k, int d, real_t *U_ARR,
      real_t *F_ARR, real_t *F_ARR_INT);
  void setAllFluxInt(idx_t i, idx_t j, idx_t k);

  void evolveFluid(idx_t i, idx_t j, idx_t k);

  void addBSSNSrc(std::map <std::string, real_t *> & bssn_fields);

  void init();

  void stepTerm();

};

}

#endif
