#ifndef COSMO_SCALAR
#define COSMO_SCALAR

#include "../../cosmo_types.h"
#include "../bssn/bssn.h"


namespace cosmo
{

typedef struct {

  // field values
  real_t phi, Pi, psi1, psi2, psi3;

  // derivatives of fields
  real_t d1phi, d2phi, d3phi;
  real_t d1Pi, d2Pi, d3Pi;
  real_t d1psi1, d2psi1, d3psi1;
  real_t d1psi2, d2psi2, d3psi2;
  real_t d1psi3, d2psi3, d3psi3;

} ScalarData;

/**
 * @brief Class implementing functionality for a scalar field that relies on a
 * BSSN instance.
 */
class Scalar
{
  // real_t (*_dV)(real_t);
  // real_t (*_V)(real_t);

public:
  register_t phi;
  register_t Pi;
  register_t psi1;
  register_t psi2;
  register_t psi3;

  Scalar();
  ~Scalar();

  void stepInit();
  void K1Finalize();
  void K2Finalize();
  void K3Finalize();
  void K4Finalize();
  void RKEvolve(BSSNData *bd);
  void RKEvolvePt(BSSNData *bd);

  ScalarData getScalarData(BSSNData *bd);

  real_t dt_phi(BSSNData *bd, ScalarData *sd);
  real_t dt_Pi(BSSNData *bd, ScalarData *sd);
  real_t dt_psi1(BSSNData *bd, ScalarData *sd);
  real_t dt_psi2(BSSNData *bd, ScalarData *sd);
  real_t dt_psi3(BSSNData *bd, ScalarData *sd);

  // Determine potential using an anonymous function?
  // template<typename F>
  // set_dV()
  // {
  // }

  real_t dV(real_t phi_in);
  real_t V(real_t phi_in);
  void addBSSNSource(BSSN * bssn);

  real_t scalarConstraint(idx_t i, idx_t j, idx_t k, idx_t dir);

};

} // namespace cosmo

#endif
