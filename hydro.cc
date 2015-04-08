#include "hydro.h"

namespace cosmo
{

Hydro::Hydro()
{
  // w=0 fluid *only* for now.
  // only works in synchronous gauge too.

  HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_ALLOC)
  HYDRO_APPLY_TO_FLUXES(FLUX_ARRAY_ALLOC)
  HYDRO_APPLY_TO_FLUXES_INT(FLUX_ARRAY_ALLOC)

  HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_ADDMAP)
  HYDRO_APPLY_TO_FLUXES(FLUX_ARRAY_ADDMAP)
  HYDRO_APPLY_TO_FLUXES_INT(FLUX_ARRAY_ADDMAP)
}

Hydro::~Hydro()
{
  HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_DELETE)
  HYDRO_APPLY_TO_FLUXES(FLUX_ARRAY_DELETE)
  HYDRO_APPLY_TO_FLUXES_INT(FLUX_ARRAY_DELETE)
}

void Hydro::setQuantitiesCell(BSSNData *paq, HydroData *hdp)
{
  setFluxesCell(paq, hdp);
}

void Hydro::setFluxesCell(BSSNData *paq, HydroData *hdp)
{
  idx_t idx = paq->idx;
  idx_t f_idx_1 = F_NP_INDEX(paq->i, paq->j, paq->k, 1);
  idx_t f_idx_2 = F_NP_INDEX(paq->i, paq->j, paq->k, 2);
  idx_t f_idx_3 = F_NP_INDEX(paq->i, paq->j, paq->k, 3);

  real_t gijsisj = exp(4.0*paq->phi)*( // gammas from BSSN sim are conformal gammas
        paq->gamma11*US1_a[idx]*US1_a[idx] + paq->gamma22*US2_a[idx]*US2_a[idx] + paq->gamma33*US3_a[idx]*US3_a[idx] 
       + 2.0*(paq->gamma12*US1_a[idx]*US2_a[idx] + paq->gamma13*US1_a[idx]*US3_a[idx] + paq->gamma23*US2_a[idx]*US3_a[idx] )
      );

  // common factor
  real_t sqrt_D2_S2 = sqrt(
      pw2(UD_a[idx])
      + gijsisj 
    );

  real_t US1cont = exp(4.0*paq->phi)*(paq->gammai11*US1_a[idx] + paq->gammai12*US2_a[idx] + paq->gammai13*US3_a[idx]);
  real_t US2cont = exp(4.0*paq->phi)*(paq->gammai12*US1_a[idx] + paq->gammai22*US2_a[idx] + paq->gammai23*US3_a[idx]);
  real_t US3cont = exp(4.0*paq->phi)*(paq->gammai13*US1_a[idx] + paq->gammai23*US2_a[idx] + paq->gammai33*US3_a[idx]);

  // Calculate fluxes in each cell
  // D
    FD_a[f_idx_1] = UD_a[idx]*US1cont/sqrt_D2_S2;
    FD_a[f_idx_2] = UD_a[idx]*US2cont/sqrt_D2_S2;
    FD_a[f_idx_3] = UD_a[idx]*US3cont/sqrt_D2_S2;
  // S1
    FS1_a[f_idx_1] = US1_a[idx]*US1cont/sqrt_D2_S2;
    FS1_a[f_idx_2] = US1_a[idx]*US2cont/sqrt_D2_S2;
    FS1_a[f_idx_3] = US1_a[idx]*US3cont/sqrt_D2_S2;
  // S2
    FS2_a[f_idx_1] = US2_a[idx]*US1cont/sqrt_D2_S2;
    FS2_a[f_idx_2] = US2_a[idx]*US2cont/sqrt_D2_S2;
    FS2_a[f_idx_3] = US2_a[idx]*US3cont/sqrt_D2_S2;
  // S3
    FS3_a[f_idx_1] = US3_a[idx]*US1cont/sqrt_D2_S2;
    FS3_a[f_idx_2] = US3_a[idx]*US2cont/sqrt_D2_S2;
    FS3_a[f_idx_3] = US3_a[idx]*US3cont/sqrt_D2_S2;

if(idx == 32*32*15)
{
  std::cout << "\nD flux: " << FD_a[f_idx_1] << "\n"
            << "S flux: " << FS1_a[f_idx_1] << "\n"
            << "sqrt_D2_S2: " << sqrt_D2_S2 << "\n"
            << "gijsisj: " << gijsisj << "\n";
}

}

void Hydro::setOneFluxInt(idx_t i, idx_t j, idx_t k, int d, real_t *U_ARR,
    real_t *F_ARR, real_t *F_ARR_INT)
{
  // Calculate fluxes at cell interfaces using a WENO scheme
  // http://lsec.cc.ac.cn/lcfd/DEWENO/weno.pdf
  // => interfaces stored are at "+1/2" to i,j,k cell

  real_t g1, g2, g0;
  real_t p1, p2, p0;
  real_t b1, b2, b0, ebtot2;
  real_t w1, w2, w0, wtot;

  // construct indexes depending on direction of flux
  idx_t f_idx_m2=0, f_idx_m1=0, f_idx=0, f_idx_p1=0, f_idx_p2=0, f_idx_p3=0;
  idx_t u_idx = INDEX(i,j,k), u_idx_p1=0;
  switch (d)
  {
    case 1:
      f_idx_m2 = F_INDEX(i-2,j,k,1);
      f_idx_m1 = F_INDEX(i-1,j,k,1);
      f_idx    = F_INDEX(i  ,j,k,1);
      f_idx_p1 = F_INDEX(i+1,j,k,1);
      f_idx_p2 = F_INDEX(i+2,j,k,1);
      f_idx_p3 = F_INDEX(i+3,j,k,1);
      u_idx_p1 = INDEX(i+1,j,k);
      break;
    case 2:
      f_idx_m2 = F_INDEX(i,j-2,k,2);
      f_idx_m1 = F_INDEX(i,j-1,k,2);
      f_idx    = F_INDEX(i,j  ,k,2);
      f_idx_p1 = F_INDEX(i,j+1,k,2);
      f_idx_p2 = F_INDEX(i,j+2,k,2);
      f_idx_p3 = F_INDEX(i,j+3,k,2);
      u_idx_p1 = INDEX(i,j+1,k);
      break;
    case 3:
      f_idx_m2 = F_INDEX(i,j,k-2,3);
      f_idx_m1 = F_INDEX(i,j,k-1,3);
      f_idx    = F_INDEX(i,j,k  ,3);
      f_idx_p1 = F_INDEX(i,j,k+1,3);
      f_idx_p2 = F_INDEX(i,j,k+2,3);
      f_idx_p3 = F_INDEX(i,j,k+3,3);
      u_idx_p1 = INDEX(i,j,k+1);
      break;
  }

  // 1) 
  // determine direction of prop. across interface (Roe speed)
  // => determines direction to take stencils
  real_t a = SIGN( (F_ARR[f_idx] - F_ARR[f_idx_p1])
             / (U_ARR[u_idx] - U_ARR[u_idx_p1]) );

  // 2) 
  // stencils are determined by the sign of `a' (mirror images of each other)
  if(a > 0.0)
  {
    p0 = F_ARR[f_idx_m2]/3.0 - 7.0/6.0*F_ARR[f_idx_m1] + 11.0/6.0*F_ARR[f_idx];
    p1 = -F_ARR[f_idx_m1]/6.0 + 5.0/6.0*F_ARR[f_idx] + F_ARR[f_idx_p1]/3.0;
    p2 = F_ARR[f_idx]/3.0 + 5.0/6.0*F_ARR[f_idx_p1] - F_ARR[f_idx_p2]/6.0;
  }
  else
  {
    p0 = F_ARR[f_idx_p3]/3.0 - 7.0/6.0*F_ARR[f_idx_p2] + 11.0/6.0*F_ARR[f_idx_p1];
    p1 = -F_ARR[f_idx_p2]/6.0 + 5.0/6.0*F_ARR[f_idx_p1] + F_ARR[f_idx]/3.0;
    p2 = F_ARR[f_idx_p1]/3.0 + 5.0/6.0*F_ARR[f_idx] - F_ARR[f_idx_m1]/6.0;
  }

  // 3) calculate "smoothness" indicators
  if(a > 0.0)
  {
    b0 = 13.0/12.0 * pw2(F_ARR[f_idx_m2] - 2.0*F_ARR[f_idx_m1] + F_ARR[f_idx])
         + 1.0/4.0 * pw2(3.0*F_ARR[f_idx_m2] - 4.0*F_ARR[f_idx_m1] + F_ARR[f_idx]);
    b1 = 13.0/12.0 * pw2(F_ARR[f_idx_m1] - 2.0*F_ARR[f_idx] + F_ARR[f_idx_p1])
         + 1.0/4.0 * pw2(3.0*F_ARR[f_idx_m1] - F_ARR[f_idx_p1]);
    b2 = 13.0/12.0 * pw2(F_ARR[f_idx] - 2.0*F_ARR[f_idx_p1] + F_ARR[f_idx_p2])
         + 1.0/4.0 * pw2(F_ARR[f_idx] - 4.0*F_ARR[f_idx_p1] + F_ARR[f_idx_p2]);
  }
  else
  {
    b0 = 13.0/12.0 * pw2(F_ARR[f_idx_p3] - 2.0*F_ARR[f_idx_p2] + F_ARR[f_idx_p1])
         + 1.0/4.0 * pw2(3.0*F_ARR[f_idx_p3] - 4.0*F_ARR[f_idx_p2] + F_ARR[f_idx_p1]);
    b1 = 13.0/12.0 * pw2(F_ARR[f_idx_p2] - 2.0*F_ARR[f_idx_p1] + F_ARR[f_idx])
         + 1.0/4.0 * pw2(3.0*F_ARR[f_idx_p2] - F_ARR[f_idx]);
    b2 = 13.0/12.0 * pw2(F_ARR[f_idx_p1] - 2.0*F_ARR[f_idx] + F_ARR[f_idx_m1])
         + 1.0/4.0 * pw2(F_ARR[f_idx_p1] - 4.0*F_ARR[f_idx] + F_ARR[f_idx_m1]);
  }

  // initial 5th-order accurate weights
  g0 = 0.1;
  g1 = 0.6;
  g2 = 0.3;

  // nonlinear weights:
  ebtot2 = pw2(EPS + b0) + pw2(EPS + b1) + pw2(EPS + b2);
  w0 = g0 / ebtot2;
  w1 = g1 / ebtot2;
  w2 = g2 / ebtot2;
  wtot = w0 + w1 + w2;
  w0 /= wtot;
  w1 /= wtot;
  w2 /= wtot;

  // Final flux
  real_t flux;
  flux = w0*p0 + w1*p1 + w2*p2;

  F_ARR_INT[f_idx] = flux;
}
  
void Hydro::setAllFluxInt(idx_t i, idx_t j, idx_t k)
{
  // "D" flux
  setOneFluxInt(i, j, k, 1, UD_a, FD_a, FD_int_a);
  setOneFluxInt(i, j, k, 2, UD_a, FD_a, FD_int_a);
  setOneFluxInt(i, j, k, 3, UD_a, FD_a, FD_int_a);

  // "S1" flux
  setOneFluxInt(i, j, k, 1, US1_a, FS1_a, FS1_int_a);
  setOneFluxInt(i, j, k, 2, US1_a, FS1_a, FS1_int_a);
  setOneFluxInt(i, j, k, 3, US1_a, FS1_a, FS1_int_a);

  // "S2" flux
  setOneFluxInt(i, j, k, 1, US2_a, FS2_a, FS2_int_a);
  setOneFluxInt(i, j, k, 2, US2_a, FS2_a, FS2_int_a);
  setOneFluxInt(i, j, k, 3, US2_a, FS2_a, FS2_int_a);

  // "S3" flux
  setOneFluxInt(i, j, k, 1, US3_a, FS3_a, FS3_int_a);
  setOneFluxInt(i, j, k, 2, US3_a, FS3_a, FS3_int_a);
  setOneFluxInt(i, j, k, 3, US3_a, FS3_a, FS3_int_a);
}

void Hydro::evolveFluid(idx_t i, idx_t j, idx_t k)
{
  idx_t idx = INDEX(i,j,k);

  /* semi-log scheme for density variable (enforces positivity) */
  /* May need to check for UD_a being too small or zero? */
  /* http://www.wccm-eccm-ecfd2014.org/admin/files/filePaper/p2390.pdf */
  if(UD_a[idx] < 1e-10)
  {
    UD_f[idx] = UD_a[idx];
  }
  else
  {
    UD_f[idx] = UD_a[idx]*exp(-dt/dx/dx/UD_a[idx]*(
        FD_int_a[F_NP_INDEX(i,j,k,1)] + FD_int_a[F_NP_INDEX(i,j,k,2)] + FD_int_a[F_NP_INDEX(i,j,k,3)]
        - FD_int_a[F_INDEX(i-1,j,k,1)] - FD_int_a[F_INDEX(i,j-1,k,2)] - FD_int_a[F_INDEX(i,j,k-1,3)]
      ));
  }

  US1_f[idx] = US1_a[idx] - dt/dx/dx*(
      FS1_int_a[F_NP_INDEX(i,j,k,1)] + FS1_int_a[F_NP_INDEX(i,j,k,2)] + FS1_int_a[F_NP_INDEX(i,j,k,3)]
      - FS1_int_a[F_INDEX(i-1,j,k,1)] - FS1_int_a[F_INDEX(i,j-1,k,2)] - FS1_int_a[F_INDEX(i,j,k-1,3)]
    );

  US2_f[idx] = US2_a[idx] - dt/dx/dx*(
      FS2_int_a[F_NP_INDEX(i,j,k,1)] + FS2_int_a[F_NP_INDEX(i,j,k,2)] + FS2_int_a[F_NP_INDEX(i,j,k,3)]
      - FS2_int_a[F_INDEX(i-1,j,k,1)] - FS2_int_a[F_INDEX(i,j-1,k,2)] - FS2_int_a[F_INDEX(i,j,k-1,3)]
    );

  US3_f[idx] = US3_a[idx] - dt/dx/dx*(
      FS3_int_a[F_NP_INDEX(i,j,k,1)] + FS3_int_a[F_NP_INDEX(i,j,k,2)] + FS3_int_a[F_NP_INDEX(i,j,k,3)]
      - FS3_int_a[F_INDEX(i-1,j,k,1)] - FS3_int_a[F_INDEX(i,j-1,k,2)] - FS3_int_a[F_INDEX(i,j,k-1,3)]
    );

}

void Hydro::addBSSNSrc(std::map <std::string, real_t *> & bssn_fields)
{

 // use const_iterator to walk through elements of pairs
 // for ( std::map< std::string, real_t *>::const_iterator iter = bssn_fields.begin();
 //    iter != bssn_fields.end(); ++iter )
 //    std::cout << iter->first << '\n';

  idx_t i, j, k;
  #pragma omp parallel for default(shared) private(i, j, k)
  LOOP3(i,j,k)
  {
    idx_t idx = INDEX(i,j,k);

    real_t g11 = bssn_fields["gamma11_a"][idx];
    real_t g12 = bssn_fields["gamma12_a"][idx];
    real_t g13 = bssn_fields["gamma13_a"][idx];
    real_t g22 = bssn_fields["gamma22_a"][idx];
    real_t g23 = bssn_fields["gamma23_a"][idx];
    real_t g33 = bssn_fields["gamma33_a"][idx];

    real_t gi11 = g22*g33 - g23*g23;
    real_t gi12 = g13*g23 - g12*g33;
    real_t gi13 = g12*g23 - g13*g22;
    real_t gi22 = g11*g33 - g13*g13;
    real_t gi23 = g12*g13 - g23*g11;
    real_t gi33 = g11*g22 - g12*g12;

    real_t a = exp(4.0*bssn_fields["phi_a"][idx]); // ~ scale factor a
    real_t g = exp(12.0*bssn_fields["phi_a"][idx]);
    real_t rg = exp(6.0*bssn_fields["phi_a"][idx]);

    real_t W_rel, r, u1, u2, u3;
    W_rel = 1.0;
    r = 0.0;
    u1 = 0.0;
    u2 = 0.0;
    u3 = 0.0;

    if(UD_a[idx] > 0.0)
    {

      W_rel = sqrt(
          1.0 + a * a * (
             gi11*US1_a[idx]*US1_a[idx] + gi22*US2_a[idx]*US2_a[idx] + gi33*US3_a[idx]*US3_a[idx]
              + 2.0*(gi12*US1_a[idx]*US2_a[idx] + gi13*US1_a[idx]*US3_a[idx] + gi23*US2_a[idx]*US3_a[idx])
            ) / UD_a[idx] / UD_a[idx]
          );

      r =  UD_a[idx] / rg / W_rel;
      u1 = US1_a[idx] / W_rel / rg / r;
      u2 = US2_a[idx] / W_rel / rg / r;
      u3 = US3_a[idx] / W_rel / rg / r;
    }

    bssn_fields["r_a"][idx] += r*W_rel*W_rel;

    bssn_fields["S1_a"][idx] += r*W_rel*u1;
    bssn_fields["S2_a"][idx] += r*W_rel*u2;
    bssn_fields["S3_a"][idx] += r*W_rel*u3;

    real_t S = r*( W_rel*W_rel - 1.0 );
    bssn_fields["S_a"][idx] += S;

    bssn_fields["STF11_a"][idx] += r*(u1*u1 - a*g11*(W_rel*W_rel - 1.0));
    bssn_fields["STF12_a"][idx] += r*(u1*u2 - a*g12*(W_rel*W_rel - 1.0));
    bssn_fields["STF13_a"][idx] += r*(u1*u3 - a*g13*(W_rel*W_rel - 1.0));
    bssn_fields["STF22_a"][idx] += r*(u2*u2 - a*g22*(W_rel*W_rel - 1.0));
    bssn_fields["STF23_a"][idx] += r*(u2*u3 - a*g23*(W_rel*W_rel - 1.0));
    bssn_fields["STF33_a"][idx] += r*(u3*u3 - a*g33*(W_rel*W_rel - 1.0));
  }

}

void Hydro::init()
{
  // initialize values
  idx_t i, j, k;
  LOOP3(i,j,k)
  {
    idx_t idx = INDEX(i,j,k);

    // empty universe
    UD_a[idx]  = UD_f[idx]  = 0.0;
    US1_a[idx] = US1_f[idx] = 0.0;
    US2_a[idx] = US2_f[idx] = 0.0;
    US3_a[idx] = US3_f[idx] = 0.0;
  }

}

void Hydro::stepTerm()
{
  std::swap(UD_a, UD_f);
  std::swap(fields["UD_a"], fields["UD_f"]);
  std::swap(US1_a, US1_f);
  std::swap(fields["US1_a"], fields["US1_f"]);
  std::swap(US2_a, US2_f);
  std::swap(fields["US2_a"], fields["US2_f"]);
  std::swap(US3_a, US3_f);
  std::swap(fields["US3_a"], fields["US3_f"]);
}


} /* namespace cosmo */
