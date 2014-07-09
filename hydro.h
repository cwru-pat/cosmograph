#ifndef COSMO_HYDRO
#define COSMO_HYDRO

#include "cosmo.h"

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

public:
  std::map <std::string, real_t *> fields;

  BSSN();
  ~BSSN();

  real_t der(real_t field_adj[3][3][3], int d)
  {
   // return 0;
    
    switch (d) {
      case 1:
        return field_adj[2][1][1] - field_adj[0][1][1];
        break;
      case 2:
        return field_adj[1][2][1] - field_adj[1][0][1];
        break;
      case 3:
        return field_adj[1][1][2] - field_adj[1][1][0];
        break;
    }

    /* XXX */
    return 0;
  }

  /* set current local field values */
  void set_local_vals(PointData *paq)
  {

  }

  // calculate needed quantities (need the inverse metric set everywhere first)
  void set_paq_values(idx_t i, idx_t j, idx_t k, PointData *paq)
  {

  }

  void step();
  void init();

};

}

#endif
