#ifndef COSMO_DUST
#define COSMO_DUST

#include "../../cosmo_types.h"
#include "../bssn/bssn.h"


namespace cosmo
{

typedef struct {

  // conformal field values (D = tilde(D) from notes)
  // or D = sqrt(gamma) D From Rezzolla
  real_t D, S1, S2, S3;

  // Derived stuff
  real_t W;
  real_t v1, v2, v3; // contravariant velocities

} DustData;


/**
 * @brief Class implementing functionality for a dust fluid that relies on a
 * BSSN instance.
 */
class Dust
{

public:
  register_t D;
  register_t S1;
  register_t S2;
  register_t S3;

  arr_t aDv1, aDv2, aDv3;
  
  arr_t aS1v1, aS1v2, aS1v3;
  arr_t aS2v1, aS2v2, aS2v3;
  arr_t aS3v1, aS3v2, aS3v3;

  arr_t S1src, S2src, S3src;

  arr_t detg, g11, g12, g13, g22, g23, g33;
  arr_t W;

  Dust();
  ~Dust();

  void setDt(real_t dt);

  void stepInit(BSSN *bssn);
  void K1Finalize();
  void K2Finalize();
  void K3Finalize();
  void K4Finalize();
  void populateDerivedFields(BSSN *bssn);
  void RKEvolve(BSSN *bssn);

  real_t dt_D(idx_t i, idx_t j, idx_t k);
  real_t dt_S1(idx_t i, idx_t j, idx_t k);
  real_t dt_S2(idx_t i, idx_t j, idx_t k);
  real_t dt_S3(idx_t i, idx_t j, idx_t k);

  DustData getDustData(BSSNData *bd);

  void addBSSNSrc(BSSN *bssn);
};

} // namespace cosmo

#endif
