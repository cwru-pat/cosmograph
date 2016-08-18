#include "scalar.h"
#include "../../utils/math.h"
#include "../../cosmo_includes.h"
#include "../../cosmo_globals.h"

namespace cosmo
{

/**
 * @brief Constructor: initialize fields needed for scalar evolution,
 * set timestep according to global `dt`.
 */
Scalar::Scalar():
  phi(), Pi(), psi1(), psi2(), psi3()
{
  std::cout << "Creating scalar class with dt=" << dt << "\n";

  phi.init(NX, NY, NZ, dt);
  Pi.init(NX, NY, NZ, dt);
  psi1.init(NX, NY, NZ, dt);
  psi2.init(NX, NY, NZ, dt);
  psi3.init(NX, NY, NZ, dt);
}

Scalar::~Scalar()
{
  phi.~RK4Register();
  Pi.~RK4Register();
  psi1.~RK4Register();
  psi2.~RK4Register();
  psi3.~RK4Register();
}

/**
 * @brief Call RK4Register::stepInit for fields.
 */
void Scalar::stepInit()
{
  phi.stepInit();
  Pi.stepInit();
  psi1.stepInit();
  psi2.stepInit();
  psi3.stepInit();
}

/**
 * @brief Call RK4Register::K1Finalize for fields.
 */
void Scalar::K1Finalize()
{
  phi.K1Finalize();
  Pi.K1Finalize();
  psi1.K1Finalize();
  psi2.K1Finalize();
  psi3.K1Finalize();
}

void Scalar::K2Finalize()
{
  phi.K2Finalize();
  Pi.K2Finalize();
  psi1.K2Finalize();
  psi2.K2Finalize();
  psi3.K2Finalize();
}

void Scalar::K3Finalize()
{
  phi.K3Finalize();
  Pi.K3Finalize();
  psi1.K3Finalize();
  psi2.K3Finalize();
  psi3.K3Finalize();
}

void Scalar::K4Finalize()
{
  phi.K4Finalize();
  Pi.K4Finalize();
  psi1.K4Finalize();
  psi2.K4Finalize();
  psi3.K4Finalize();
}

void Scalar::RKEvolve(BSSNData *bd)
{
  idx_t i, j, k;
  
  #pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    RKEvolvePt(bd);
  }
}

void Scalar::RKEvolvePt(BSSNData *bd)
{
  ScalarData sd = getScalarData(bd);
  idx_t idx = bd->idx;

  phi._array_c[idx] = dt_phi(bd, &sd);
  Pi._array_c[idx] = dt_Pi(bd, &sd);
  psi1._array_c[idx] = dt_psi1(bd, &sd);
  psi2._array_c[idx] = dt_psi2(bd, &sd);
  psi3._array_c[idx] = dt_psi3(bd, &sd);
}

ScalarData Scalar::getScalarData(BSSNData *bd)
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

real_t Scalar::dt_phi(BSSNData *bd, ScalarData *sd)
{
  return (
    #if(USE_BSSN_SHIFT)
      bd->beta1*sd->d1phi + bd->beta2*sd->d2phi + bd->beta3*sd->d3phi
    #endif
    - bd->alpha*sd->Pi
  );
}

real_t Scalar::dt_Pi(BSSNData *bd, ScalarData *sd)
{
  return (
    #if(USE_BSSN_SHIFT)
      bd->beta1*sd->d1Pi + bd->beta2*sd->d2Pi + bd->beta3*sd->d3Pi
    #endif
    -exp(-4.0*bd->phi)*(
      bd->gammai11*(bd->alpha*sd->d1psi1 + sd->psi1*bd->d1a) + bd->gammai12*(bd->alpha*sd->d1psi2 + sd->psi1*bd->d2a) + bd->gammai13*(bd->alpha*sd->d1psi3 + sd->psi1*bd->d3a)
      + bd->gammai21*(bd->alpha*sd->d2psi1 + sd->psi2*bd->d1a) + bd->gammai22*(bd->alpha*sd->d2psi2 + sd->psi2*bd->d2a) + bd->gammai23*(bd->alpha*sd->d2psi3 + sd->psi2*bd->d3a)
      + bd->gammai31*(bd->alpha*sd->d3psi1 + sd->psi3*bd->d1a) + bd->gammai32*(bd->alpha*sd->d3psi2 + sd->psi3*bd->d2a) + bd->gammai33*(bd->alpha*sd->d3psi3 + sd->psi3*bd->d3a)
    ) + bd->alpha*( (
        bd->Gamma1 - 2.0*(bd->gammai11*bd->d1phi + bd->gammai12*bd->d2phi + bd->gammai13*bd->d3phi)
      )*sd->psi1*exp(-4.0*bd->phi) + (
        bd->Gamma2 - 2.0*(bd->gammai21*bd->d1phi + bd->gammai22*bd->d2phi + bd->gammai23*bd->d3phi)
      )*sd->psi2*exp(-4.0*bd->phi) + (
        bd->Gamma3 - 2.0*(bd->gammai31*bd->d1phi + bd->gammai32*bd->d2phi + bd->gammai33*bd->d3phi)
      )*sd->psi3* exp(-4.0*bd->phi)
      + bd->K*sd->Pi
      + dV(sd->phi)
    )
  );
}

real_t Scalar::dt_psi1(BSSNData *bd, ScalarData *sd)
{
  return (
    #if(USE_BSSN_SHIFT)
      bd->beta1*sd->d1psi1 + bd->beta2*sd->d2psi1 + bd->beta3*sd->d3psi1
      + sd->psi1*bd->d1beta1 + sd->psi2*bd->d1beta2 + sd->psi3*bd->d1beta3
    #endif
    - bd->alpha*sd->d1Pi
    - sd->Pi*bd->d1a
  );
}

real_t Scalar::dt_psi2(BSSNData *bd, ScalarData *sd)
{
  return (
    #if(USE_BSSN_SHIFT)
      bd->beta1*sd->d1psi2 + bd->beta2*sd->d2psi2 + bd->beta3*sd->d3psi2
      + sd->psi1*bd->d2beta1 + sd->psi2*bd->d2beta2 + sd->psi3*bd->d2beta3
    #endif
    - bd->alpha*sd->d2Pi
    - sd->Pi*bd->d2a
  );
}

real_t Scalar::dt_psi3(BSSNData *bd, ScalarData *sd)
{
  return (
    #if(USE_BSSN_SHIFT)
      bd->beta1*sd->d1psi3 + bd->beta2*sd->d2psi3 + bd->beta3*sd->d3psi3
      + sd->psi1*bd->d3beta1 + sd->psi2*bd->d3beta2 + sd->psi3*bd->d3beta3
    #endif
    - bd->alpha*sd->d3Pi
    - sd->Pi*bd->d3a
  );
}

void Scalar::addBSSNSource(BSSN * bssn)
{
  arr_t & DIFFr_a = *bssn->fields["DIFFr_a"];
  arr_t & DIFFS_a = *bssn->fields["DIFFS_a"];
  arr_t & S1_a = *bssn->fields["S1_a"];
  arr_t & S2_a = *bssn->fields["S2_a"];
  arr_t & S3_a = *bssn->fields["S3_a"];
  arr_t & STF11_a = *bssn->fields["STF11_a"];
  arr_t & STF12_a = *bssn->fields["STF12_a"];
  arr_t & STF13_a = *bssn->fields["STF13_a"];
  arr_t & STF22_a = *bssn->fields["STF22_a"];
  arr_t & STF23_a = *bssn->fields["STF23_a"];
  arr_t & STF33_a = *bssn->fields["STF33_a"];

  idx_t i, j, k;

  #pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i, j, k)
  {
    idx_t idx = INDEX(i,j,k);

    BSSNData bd = {0};
    // TODO: remove redundant computations here?
    bssn->set_bd_values(i, j, k, &bd);
    ScalarData sd = getScalarData(&bd);

    // n^mu d_mu phi
    real_t nmudmuphi = (
      dt_phi(&bd, &sd)
      #if(USE_BSSN_SHIFT)
        - bd.beta1*sd.d1phi - bd.beta2*sd.d2phi - bd.beta3*sd.d3phi
      #endif
    )/bd.alpha;
    
    // gammai^ij d_j phi d_i phi
    real_t diphidiphi = (
      bd.gammai11*sd.d1phi*sd.d1phi + bd.gammai22*sd.d2phi*sd.d2phi + bd.gammai33*sd.d3phi*sd.d3phi
      + 2.0*(bd.gammai12*sd.d1phi*sd.d2phi + bd.gammai13*sd.d1phi*sd.d3phi + bd.gammai23*sd.d2phi*sd.d3phi)
    );

    DIFFr_a[idx] += 0.5*nmudmuphi*nmudmuphi
      + 0.5*exp(-4.0*bd.phi)*diphidiphi + V(sd.phi);

    DIFFS_a[idx] += 3.0/2.0*nmudmuphi*nmudmuphi
      - 0.5*exp(-4.0*bd.phi)*diphidiphi - 3.0*V(sd.phi);

    S1_a[idx] += -nmudmuphi*sd.d1phi;
    S2_a[idx] += -nmudmuphi*sd.d2phi;
    S3_a[idx] += -nmudmuphi*sd.d3phi;

    STF11_a[idx] += sd.d1phi*sd.d1phi - bd.gamma11/3.0*diphidiphi;
    STF12_a[idx] += sd.d1phi*sd.d2phi - bd.gamma12/3.0*diphidiphi;
    STF13_a[idx] += sd.d1phi*sd.d3phi - bd.gamma13/3.0*diphidiphi;
    STF22_a[idx] += sd.d2phi*sd.d2phi - bd.gamma22/3.0*diphidiphi;
    STF23_a[idx] += sd.d2phi*sd.d3phi - bd.gamma23/3.0*diphidiphi;
    STF33_a[idx] += sd.d3phi*sd.d3phi - bd.gamma33/3.0*diphidiphi;
  }

  return;
}

real_t Scalar::dV(real_t phi_in)
{
  return 0;
}

real_t Scalar::V(real_t phi_in)
{
  // TODO: discuss how to set this better
#if USE_COSMO_CONST_POTENTIAL
  return COSMO_CONST;
#endif
  return 1.0;
}

real_t Scalar::scalarConstraint(idx_t i, idx_t j, idx_t k, idx_t dir)
{
  switch(dir)
  {
    case 1:
      return derivative(i, j, k, dir, phi._array_a) - psi1._array_a[INDEX(i,j,k)];
    case 2:
      return derivative(i, j, k, dir, phi._array_a) - psi2._array_a[INDEX(i,j,k)];
    case 3:
      return derivative(i, j, k, dir, phi._array_a) - psi3._array_a[INDEX(i,j,k)];
  }

  throw -1;
  return 0;
}

} // namespace cosmo
