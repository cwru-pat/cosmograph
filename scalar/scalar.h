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

  // field values
  real_t phi, Pi, psi1, psi2, psi3;

  // derivatives of fields
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

  ScalarData getScalarData(BSSNData *bd)
  {
    ScalarData sd = {0};
    idx_t i = bd->i, j = bd->j, k = bd->k;
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

  void RKEvolvePt(BSSNData *bd)
  {
    ScalarData sd = getScalarData(bd);

    phi._array_c = dt_phi(bd, &sd);
    Pi._array_c = dt_Pi(bd, &sd);
    psi1._array_c = dt_psi1(bd, &sd);
    psi2._array_c = dt_psi2(bd, &sd);
    psi3._array_c = dt_psi3(bd, &sd);
  }

  real_t dt_phi(BSSNData *bd, ScalarData *sd)
  {
    return (
      bd->beta1*sd->d1phi + bd->beta2*sd->d2phi + bd->beta3*sd->d3phi
      - bd->alpha*sd->Pi
    );
  }

  real_t dt_Pi(BSSNData *bd, ScalarData *sd)
  {
    return (
      bd->beta1*sd->d1Pi + bd->beta2*sd->d2Pi + bd->beta3*sd->d3Pi
      -exp(4.0*bd->phi)*(
        bd->gammai11*(bd->alpha*sd->d1psi1 + sd->psi1*bd->d1a) + bd->gammai12*(bd->alpha*sd->d1psi2 + sd->psi1*bd->d2a) + bd->gammai13*(bd->alpha*sd->d1psi3 + sd->psi1*bd->d3a)
        + bd->gammai21*(bd->alpha*sd->d2psi1 + sd->psi2*bd->d1a) + bd->gammai22*(bd->alpha*sd->d2psi2 + sd->psi2*bd->d2a) + bd->gammai23*(bd->alpha*sd->d2psi3 + sd->psi2*bd->d3a)
        + bd->gammai31*(bd->alpha*sd->d3psi1 + sd->psi3*bd->d1a) + bd->gammai32*(bd->alpha*sd->d3psi2 + sd->psi3*bd->d2a) + bd->gammai33*(bd->alpha*sd->d3psi3 + sd->psi3*bd->d3a)
      ) + bd->alpha*( (
          bd->Gamma1 - 2.0*(bd->gammai11*bd->d1phi + bd->gammai12*bd->d2phi + bd->gammai13*bd->d3phi)
        )*sd->psi1 + (
          bd->Gamma2 - 2.0*(bd->gammai21*bd->d1phi + bd->gammai22*bd->d2phi + bd->gammai23*bd->d3phi)
        )*sd->psi2 + (
          bd->Gamma3 - 2.0*(bd->gammai31*bd->d1phi + bd->gammai32*bd->d2phi + bd->gammai33*bd->d3phi)
        )*sd->psi3
        + bd->K*sd->Pi
        // TODO: + V'(sd->phi)
      )
    );
  }

  real_t dt_psi1(BSSNData *bd, ScalarData *sd)
  {
    return (
      bd->beta1*sd->d1psi1 + bd->beta2*sd->d2psi1 + bd->beta3*sd->d3psi1
      + sd->psi1*bd->d1beta1 + sd->psi2*bd->d1beta2 + sd->psi3*bd->d1beta3
      - bd->alpha*sd->d1Pi
      - sd->Pi*bd->d1a
    );
  }

  real_t dt_psi2(BSSNData *bd, ScalarData *sd)
  {
    return (
      bd->beta1*sd->d1psi2 + bd->beta2*sd->d2psi2 + bd->beta3*sd->d3psi2
      + sd->psi1*bd->d2beta1 + sd->psi2*bd->d2beta2 + sd->psi3*bd->d2beta3
      - bd->alpha*sd->d2Pi
      - sd->Pi*bd->d2a
    );
  }

  real_t dt_psi3(BSSNData *bd, ScalarData *sd)
  {
    return (
      bd->beta1*sd->d1psi3 + bd->beta2*sd->d2psi3 + bd->beta3*sd->d3psi3
      + sd->psi1*bd->d3beta1 + sd->psi2*bd->d3beta2 + sd->psi3*bd->d3beta3
      - bd->alpha*sd->d3Pi
      - sd->Pi*bd->d3a
    );
  }

  void addBSSNSource(BSSN * bssnSim)
  {
    arr_t & DIFFr_a = *bssnSim->fields["DIFFr_a"];
    arr_t & DIFFS_a = *bssnSim->fields["DIFFS_a"];
    arr_t & S1_a = *bssnSim->fields["S1_a"];
    arr_t & S2_a = *bssnSim->fields["S2_a"];
    arr_t & S3_a = *bssnSim->fields["S3_a"];
    arr_t & S11_a = *bssnSim->fields["STF11_a"];
    arr_t & S12_a = *bssnSim->fields["STF12_a"];
    arr_t & S13_a = *bssnSim->fields["STF13_a"];
    arr_t & S22_a = *bssnSim->fields["STF22_a"];
    arr_t & S23_a = *bssnSim->fields["STF23_a"];
    arr_t & S33_a = *bssnSim->fields["STF33_a"];

    idx_t i, j, k;

    #pragma omp parallel for default(shared) private(i, j, k)
    LOOP3(i, j, k)
    {
      idx_t idx = INDEX(i,j,k);

      BSSNData bd = {0};
      // TODO: remove redundant computations here?
      bssnSim->set_paq_values(i, j, k, &bd);
      ScalarData sd = getScalarData(&bd);

      // additional quantities
      // n^mu d_mu phi
      real_t nmudmuphi = (
        dt_phi(&bd, &sd)
        - bd.beta1*sd.d1phi - bd.beta2*sd.d2phi - bd.beta3*sd.d3phi
      )/bd.alpha;
      // gammai^ij d_j phi d_i phi
      real_t diphidiphi = (
        bd.gammai11*sd.d1phi*sd.d1phi + bd.gammai22*sd.d2phi*sd.d2phi + bd.gammai33*sd.d3phi*sd.d3phi
        + 2.0*(bd.gammai12*sd.d1phi*sd.d2phi + bd.gammai13*sd.d1phi*sd.d3phi + bd.gammai23*sd.d2phi*sd.d3phi)
      );
      // g^munu d_mu phi d_nu phi
      real_t dmuphidmuphi = (
        0.0 // TODO
      );

      DIFFr_a[idx] = 0.5*nmudmuphi*nmudmuphi
        + 0.5*diphidiphi*diphidiphi; // TODO: + V(phi)

      S1_a[idx] = -sd.d1phi*nmudmuphi;
      S2_a[idx] = -sd.d2phi*nmudmuphi;
      S3_a[idx] = -sd.d3phi*nmudmuphi;

      S11_a[idx] = sd.d1phi*sd.d1phi - bd.gamma11*(-0.5*dmuphidmuphi); // TODO: V(phi)
      S12_a[idx] = sd.d1phi*sd.d2phi - bd.gamma12*(-0.5*dmuphidmuphi); // TODO: V(phi)
      S13_a[idx] = sd.d1phi*sd.d3phi - bd.gamma13*(-0.5*dmuphidmuphi); // TODO: V(phi)
      S22_a[idx] = sd.d2phi*sd.d2phi - bd.gamma22*(-0.5*dmuphidmuphi); // TODO: V(phi)
      S23_a[idx] = sd.d2phi*sd.d3phi - bd.gamma23*(-0.5*dmuphidmuphi); // TODO: V(phi)
      S33_a[idx] = sd.d3phi*sd.d3phi - bd.gamma33*(-0.5*dmuphidmuphi); // TODO: V(phi)
    }

  }

};

} // namespace cosmo

#endif
