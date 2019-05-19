#include "dust.h"
#include "../../utils/math.h"
#include "../../cosmo_includes.h"
#include "../../cosmo_globals.h"

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
 * @brief Constructor: initialize fields needed for dust evolution,
 * set timestep according to global `dt`.
 * Doesn't work with a shift! (For now?)
 */
Dust::Dust():
  D(), S1(), S2(), S3(), S4(),
  aDv1(), aDv2(), aDv3(),
  aS1v1(), aS1v2(), aS1v3(),
  aS2v1(), aS2v2(), aS2v3(),
  aS3v1(), aS3v2(), aS3v3(),
  S1src(), S2src(), S3src(),
  detg(), g11(), g12(), g13(), g22(), g23(), g33(),
  W()
{
  std::cout << "Creating dust class with dt=" << dt << "\n";

  D.init(NX, NY, NZ, dt);
  S1.init(NX, NY, NZ, dt);
  S2.init(NX, NY, NZ, dt);
  S3.init(NX, NY, NZ, dt);

  aDv1.init(NX, NY, NZ); aDv2.init(NX, NY, NZ); aDv3.init(NX, NY, NZ);

  aS1v1.init(NX, NY, NZ); aS1v2.init(NX, NY, NZ); aS1v3.init(NX, NY, NZ);
  aS2v1.init(NX, NY, NZ); aS2v2.init(NX, NY, NZ); aS2v3.init(NX, NY, NZ);
  aS3v1.init(NX, NY, NZ); aS3v2.init(NX, NY, NZ); aS3v3.init(NX, NY, NZ);

  S1src.init(NX, NY, NZ); S2src.init(NX, NY, NZ); S3src.init(NX, NY, NZ);
  detg.init(NX, NY, NZ);
  g11.init(NX, NY, NZ); g12.init(NX, NY, NZ); g13.init(NX, NY, NZ);
  g22.init(NX, NY, NZ); g23.init(NX, NY, NZ); g33.init(NX, NY, NZ);

  W.init(NX, NY, NZ);
}

Dust::~Dust()
{
  D.~RK4Register();
  S1.~RK4Register();
  S2.~RK4Register();
  S3.~RK4Register();
}

void Dust::setDt(real_t dt)
{
  D.setDt(dt);
  S1.setDt(dt);
  S2.setDt(dt);
  S3.setDt(dt);
}

/**
 * @brief Call RK4Register::stepInit for fields.
 */
void Dust::stepInit(BSSN *bssn)
{
  D.stepInit();
  S1.stepInit();
  S2.stepInit();
  S3.stepInit();
  populateDerivedFields(bssn);
}

/**
 * @brief Call RK4Register::K1Finalize for fields.
 */
void Dust::K1Finalize()
{
  D.K1Finalize();
  S1.K1Finalize();
  S2.K1Finalize();
  S3.K1Finalize();
}

void Dust::K2Finalize()
{
  D.K2Finalize();
  S1.K2Finalize();
  S2.K2Finalize();
  S3.K2Finalize();
}

void Dust::K3Finalize()
{
  D.K3Finalize();
  S1.K3Finalize();
  S2.K3Finalize();
  S3.K3Finalize();
}

void Dust::K4Finalize()
{
  D.K4Finalize();
  S1.K4Finalize();
  S2.K4Finalize();
  S3.K4Finalize();
}

void Dust::populateDerivedFields(BSSN *bssn)
{
  idx_t i, j, k;
  #pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);

    BSSNData bd = {0};
    bssn->set_bd_values(i, j, k, &bd); // TODO: get only relevant variables here
    DustData dd = getDustData(&bd);

    aDv1[idx] = bd.alpha*D[idx]*dd.v1;
    aDv2[idx] = bd.alpha*D[idx]*dd.v2;
    aDv3[idx] = bd.alpha*D[idx]*dd.v3;

    aS1v1[idx] = bd.alpha*S1[idx]*dd.v1;
    aS1v2[idx] = bd.alpha*S1[idx]*dd.v2;
    aS1v3[idx] = bd.alpha*S1[idx]*dd.v3;
    aS2v1[idx] = bd.alpha*S2[idx]*dd.v1;
    aS2v2[idx] = bd.alpha*S2[idx]*dd.v2;
    aS2v3[idx] = bd.alpha*S2[idx]*dd.v3;
    aS3v1[idx] = bd.alpha*S3[idx]*dd.v1;
    aS3v2[idx] = bd.alpha*S3[idx]*dd.v2;
    aS3v3[idx] = bd.alpha*S3[idx]*dd.v3;

    S1src[idx] = 1.0/2.0 * bd.alpha * dd.W * D[idx] * std::exp(4.0*bd->phi) *(
        dd.v1*dd.v1*(4.0*bd.gamma11*bd.d1phi + bd.d1g11)  + dd.v2*dd.v2*(4.0*bd.gamma22*bd.d1phi + bd.d1g22) + dd.v3*dd.v3*(4.0*bd.gamma33*bd.d1phi + bd.d1g33)
        + 2.0*( dd.v1*dd.v2*(4.0*bd.gamma12*bd.d1phi + bd.d1g12)  + dd.v1*dd.v3*(4.0*bd.gamma13*bd.d1phi + bd.d1g13) + dd.v2*dd.v3*(4.0*bd.gamma23*bd.d1phi + bd.d1g23) )
      ) - dd.W*D[idx]*bd.d1a;
    S2src[idx] = 1.0/2.0 * bd.alpha * dd.W * D[idx] * std::exp(4.0*bd->phi) *(
        dd.v1*dd.v1*(4.0*bd.gamma11*bd.d2phi + bd.d2g11)  + dd.v2*dd.v2*(4.0*bd.gamma22*bd.d2phi + bd.d2g22) + dd.v3*dd.v3*(4.0*bd.gamma33*bd.d2phi + bd.d2g33)
        + 2.0*( dd.v1*dd.v2*(4.0*bd.gamma12*bd.d2phi + bd.d2g12)  + dd.v1*dd.v3*(4.0*bd.gamma13*bd.d2phi + bd.d2g13) + dd.v2*dd.v3*(4.0*bd.gamma23*bd.d2phi + bd.d2g23) )
      ) - dd.W*D[idx]*bd.d2a;
    S3src[idx] = 1.0/2.0 * bd.alpha * dd.W * D[idx] * std::exp(4.0*bd->phi) *(
        dd.v1*dd.v1*(4.0*bd.gamma11*bd.d3phi + bd.d3g11)  + dd.v2*dd.v2*(4.0*bd.gamma22*bd.d3phi + bd.d3g22) + dd.v3*dd.v3*(4.0*bd.gamma33*bd.d3phi + bd.d3g33)
        + 2.0*( dd.v1*dd.v2*(4.0*bd.gamma12*bd.d3phi + bd.d3g12)  + dd.v1*dd.v3*(4.0*bd.gamma13*bd.d3phi + bd.d3g13) + dd.v2*dd.v3*(4.0*bd.gamma23*bd.d3phi + bd.d3g23) )
      ) - dd.W*D[idx]*bd.d3a;

    real_t e4phi = std::exp(4.0*bd.phi);
    detg[idx] = std::exp(6.0*bd.phi);
    g11[idx] = e4phi*bd.gamma11;
    g12[idx] = e4phi*bd.gamma12;
    g13[idx] = e4phi*bd.gamma13;
    g22[idx] = e4phi*bd.gamma22;
    g23[idx] = e4phi*bd.gamma23;
    g33[idx] = e4phi*bd.gamma33;
    W[idx] = dd.W;
  }
}

void Dust::RKEvolve(BSSN *bssn)
{
  idx_t i, j, k;
  populateDerivedFields(bssn);

  #pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    D._array_c[idx] = dt_D(i,j,k);
    S1._array_c[idx] = dt_S1(i,j,k);
    S2._array_c[idx] = dt_S2(i,j,k);
    S3._array_c[idx] = dt_S3(i,j,k);
  }
}


DustData Dust::getDustData(BSSNData *bd)
{
  DustData dd = {0};
  idx_t i = bd->i, j = bd->j, k = bd->k;
  idx_t idx = NP_INDEX(i,j,k);

  dd.D = D._array_a[idx];
  dd.S1 = S1._array_a[idx];
  dd.S2 = S2._array_a[idx];
  dd.S3 = S3._array_a[idx];

  real_t em4phi = std::exp(-4.0*bd->phi);
  real_t e6phi = std::exp(6.0*bd->phi);

  real_t giIJSISJ = em4phi*( bd->gammai11*dd.S1*dd.S1 + bd->gammai22*dd.S2*dd.S2 + bd->gammai33*dd.S3*dd.S3
    + 2.0*(bd->gammai12*dd.S1*dd.S2 + bd->gammai13*dd.S1*dd.S3 + bd->gammai23*dd.S2*dd.S3) );
  dd.W = std::sqrt(1.0 + giIJSISJ/dd.D/dd.D);

  dd.v1 = em4phi / e6phi * ( bd->gammai11*dd.S1 + bd->gammai12*dd.S2 + bd->gammai13*dd.S3 );
  dd.v2 = em4phi / e6phi * ( bd->gammai12*dd.S1 + bd->gammai22*dd.S2 + bd->gammai23*dd.S3 );
  dd.v3 = em4phi / e6phi * ( bd->gammai13*dd.S1 + bd->gammai23*dd.S2 + bd->gammai33*dd.S3 );

  return dd;
}


real_t Dust::dt_D(idx_t i, idx_t j, idx_t k)
{
  return derivative(i,j,k,1,aDv1)+derivative(i,j,k,2,aDv2)+derivative(i,j,k,3,aDv3);
}

real_t Dust::dt_S1(idx_t i, idx_t j, idx_t k)
{
  return derivative(i,j,k,1,aS1v1)+derivative(i,j,k,2,aS1v2)+derivative(i,j,k,3,aS1v3)
  + S1src[NP_INDEX(i,j,k)];
}

real_t Dust::dt_S2(idx_t i, idx_t j, idx_t k)
{
  return derivative(i,j,k,1,aS2v1)+derivative(i,j,k,2,aS2v2)+derivative(i,j,k,3,aS2v3)
  + S2src[NP_INDEX(i,j,k)];
}

real_t Dust::dt_S3(idx_t i, idx_t j, idx_t k)
{
  return derivative(i,j,k,1,aS3v1)+derivative(i,j,k,2,aS3v2)+derivative(i,j,k,3,aS3v3)
  + S3src[NP_INDEX(i,j,k)];
}

void Dust::addBSSNSource(BSSN * bssn)
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

    DIFFr_a[idx] += W[idx]*D[idx]/detg[idx];

    DIFFS_a[idx] += D[idx]/detg[idx] * (W[idx]*W[idx]-1.0)/W[idx];

    S1_a[idx] += S1[idx]/detg[idx];
    S2_a[idx] += S2[idx]/detg[idx];
    S3_a[idx] += S3[idx]/detg[idx];

    STF11_a[idx] += S1[idx]*S1[idx]/W[idx]/D[idx]/detg[idx] - 1.0/3.0*g11[idx]*DIFFS_a[idx];
    STF12_a[idx] += S1[idx]*S2[idx]/W[idx]/D[idx]/detg[idx] - 1.0/3.0*g12[idx]*DIFFS_a[idx];
    STF13_a[idx] += S1[idx]*S3[idx]/W[idx]/D[idx]/detg[idx] - 1.0/3.0*g13[idx]*DIFFS_a[idx];
    STF22_a[idx] += S2[idx]*S2[idx]/W[idx]/D[idx]/detg[idx] - 1.0/3.0*g22[idx]*DIFFS_a[idx];
    STF23_a[idx] += S2[idx]*S3[idx]/W[idx]/D[idx]/detg[idx] - 1.0/3.0*g23[idx]*DIFFS_a[idx];
    STF33_a[idx] += S3[idx]*S3[idx]/W[idx]/D[idx]/detg[idx] - 1.0/3.0*g33[idx]*DIFFS_a[idx];
  }

  return;
}


} // namespace cosmo
