#ifndef COSMO_SCALAR
#define COSMO_SCALAR

#include "../cosmo_includes.h"
#include "../cosmo_types.h"
#include "../globals.h"

#include "../bssn/bssn.h"
#include "../utils/math.h"
#include "../utils/Array.h"
#include "../utils/RK4Register.h"

namespace cosmo
{

typedef struct {

  real_t phi, Pi, psi1, psi2, psi3;

  real_t d1phi, d2phi, d3phi;
  real_t d1Pi, d2Pi, d3Pi;
  real_t d1psi1, d2psi1, d3psi1;
  real_t d1psi2, d2psi2, d3psi2;
  real_t d1psi3, d2psi3, d3psi3;

} ScalarData;

/** Scalar class **/
class Scalar
{
public:
  RK4Register<idx_t, real_t> phi;
  RK4Register<idx_t, real_t> Pi;
  RK4Register<idx_t, real_t> psi1;
  RK4Register<idx_t, real_t> psi2;
  RK4Register<idx_t, real_t> psi3;

  Scalar() :
    phi(), Pi(), psi1(), psi2(), psi3()
  {
    if(!USE_BSSN_SHIFT)
    {
      std::cout << "BSSN shift must be enabled.";
      throw -1;
    }
    if(USE_REFERENCE_FRW)
    {
      std::cout << "Reference FRW not supported.";
      throw -1;
    }

    phi.init(NX, NY, NZ, dt);
    Pi.init(NX, NY, NZ, dt);
    psi1.init(NX, NY, NZ, dt);
    psi2.init(NX, NY, NZ, dt);
    psi3.init(NX, NY, NZ, dt);
  }

  ~Scalar()
  {
    phi.~RK4Register();
    Pi.~RK4Register();
    psi1.~RK4Register();
    psi2.~RK4Register();
    psi3.~RK4Register();
  }

  void stepInit()
  {
    phi.stepInit();
    Pi.stepInit();
    psi1.stepInit();
    psi2.stepInit();
    psi3.stepInit();
  }

  void RK1Finalize()
  {
    phi.RK1Finalize();
    Pi.RK1Finalize();
    psi1.RK1Finalize();
    psi2.RK1Finalize();
    psi3.RK1Finalize();
  }

  void RK2Finalize()
  {
    phi.RK2Finalize();
    Pi.RK2Finalize();
    psi1.RK2Finalize();
    psi2.RK2Finalize();
    psi3.RK2Finalize();
  }

  void RK3Finalize()
  {
    phi.RK3Finalize();
    Pi.RK3Finalize();
    psi1.RK3Finalize();
    psi2.RK3Finalize();
    psi3.RK3Finalize();
  }

  void RK4Finalize()
  {
    phi.RK4Finalize();
    Pi.RK4Finalize();
    psi1.RK4Finalize();
    psi2.RK4Finalize();
    psi3.RK4Finalize();
  }

  ScalarData getScalarData(idx_t i, idx_t j, idx_t k)
  {
    ScalarData sd = {0};
    idx_t idx = NP_INDEX(i,j,k);

    sd.phi = phi._array_a[idx];
    sd.Pi = Pi._array_a[idx];
    sd.psi1 = psi1._array_a[idx];
    sd.psi2 = psi2._array_a[idx];
    sd.psi3 = psi3._array_a[idx];

    sd.d1phi = derivative(i, j, k, 1, phi._array_a);
    sd.d2phi = derivative(i, j, k, 2, phi._array_a);
    sd.d3phi = derivative(i, j, k, 3, phi._array_a);

    sd.d1Pi = derivative(i, j, k, 1, Pi._array_a);
    sd.d2Pi = derivative(i, j, k, 2, Pi._array_a);
    sd.d3Pi = derivative(i, j, k, 3, Pi._array_a);

    sd.d1psi1 = derivative(i, j, k, 1, psi1._array_a);
    sd.d2psi1 = derivative(i, j, k, 2, psi1._array_a);
    sd.d3psi1 = derivative(i, j, k, 3, psi1._array_a);

    sd.d1psi2 = derivative(i, j, k, 1, psi2._array_a);
    sd.d2psi2 = derivative(i, j, k, 2, psi2._array_a);
    sd.d3psi2 = derivative(i, j, k, 3, psi2._array_a);

    sd.d1psi3 = derivative(i, j, k, 1, psi3._array_a);
    sd.d2psi3 = derivative(i, j, k, 2, psi3._array_a);
    sd.d3psi3 = derivative(i, j, k, 3, psi3._array_a);

    return sd;
  }

  void RKEvolvePt(BSSNData *paq)
  {
    ScalarData sd = getScalarData(paq->i, paq->j, paq->k);

    phi._array_c = ev_phi(paq, &sd);
    Pi._array_c = ev_Pi(paq, &sd);
    psi1._array_c = ev_psi1(paq, &sd);
    psi2._array_c = ev_psi2(paq, &sd);
    psi3._array_c = ev_psi3(paq, &sd);
  }

  real_t ev_phi(BSSNData *paq, ScalarData *sd)
  {
    return 0;
  }

  real_t ev_Pi(BSSNData *paq, ScalarData *sd)
  {
    return 0;
  }

  real_t ev_psi1(BSSNData *paq, ScalarData *sd)
  {
    return 0;
  }

  real_t ev_psi2(BSSNData *paq, ScalarData *sd)
  {
    return 0;
  }

  real_t ev_psi3(BSSNData *paq, ScalarData *sd)
  {
    return 0;
  }

  void addBSSNSource()
  {
    // todo
  }

};

} // namespace cosmo

#endif
