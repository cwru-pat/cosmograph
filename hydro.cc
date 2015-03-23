#include "hydro.h"

namespace cosmo
{

Hydro::Hydro(real_t w)
{
  w_EOS = w;

  HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_ALLOC)
  HYDRO_APPLY_TO_FLUXES(FLUX_ARRAY_ALLOC)
  HYDRO_APPLY_TO_SOURCES(GEN1_ARRAY_ALLOC)
  HYDRO_APPLY_TO_PRIMITIVES(GEN1_ARRAY_ALLOC)
  HYDRO_APPLY_TO_FLUXES_INT(FLUX_ARRAY_ALLOC)

  HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_ADDMAP)
  HYDRO_APPLY_TO_FLUXES(FLUX_ARRAY_ADDMAP)
  HYDRO_APPLY_TO_SOURCES(GEN1_ARRAY_ADDMAP)
  HYDRO_APPLY_TO_PRIMITIVES(GEN1_ARRAY_ADDMAP)
  HYDRO_APPLY_TO_FLUXES_INT(FLUX_ARRAY_ADDMAP)

}

Hydro::~Hydro()
{
  HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_DELETE)
  HYDRO_APPLY_TO_FLUXES(FLUX_ARRAY_DELETE)
  HYDRO_APPLY_TO_SOURCES(GEN1_ARRAY_DELETE)
  HYDRO_APPLY_TO_PRIMITIVES(GEN1_ARRAY_DELETE)
  HYDRO_APPLY_TO_FLUXES_INT(FLUX_ARRAY_DELETE)
}

void Hydro::setQuantitiesCell(BSSNData *paq, HydroData *hdp)
{
  setPrimitivesCell(paq, hdp);
  setFluxesCell(paq, hdp);
  setSourcesCell(paq, hdp);
}

void Hydro::setPrimitivesCell(BSSNData *paq, HydroData *hdp)
{
  idx_t idx = paq->idx;

  /* root of metric determinant */
  hdp->rg = exp(6.0*paq->phi);
  hdp->g  = exp(12.0*paq->phi);
  hdp->a  = exp(2.0*paq->phi);

  if(UD_a[idx] <= 0.0)
  {
    hdp->W = 1.0;
    r_a[idx] = 0.0;
    v1_a[idx] = 0.0;
    v2_a[idx] = 0.0;
    v3_a[idx] = 0.0;
  }
  else
  {
    /* lorentz factor */
    /* W = (gamma^ij S_j S_j / D^2 / (1+w^2) + 1)^1/2 */
    hdp->W = sqrt(
        1.0 + hdp->a * hdp->a * (
           paq->gammai11*US1_a[idx]*US1_a[idx] + paq->gammai22*US2_a[idx]*US2_a[idx] + paq->gammai33*US3_a[idx]*US3_a[idx]
            + 2.0*paq->gammai12*US1_a[idx]*US2_a[idx] + 2.0*paq->gammai13*US1_a[idx]*US3_a[idx] + 2.0*paq->gammai23*US2_a[idx]*US3_a[idx]
          ) / UD_a[idx] / UD_a[idx] / pw2(1.0 + w_EOS)
        );

    /* \rho = D / gamma^1/2 / W */
    r_a[idx] = UD_a[idx] / hdp->rg / hdp->W;

    /* fluid 4-velocities (covariant) */
    real_t u1, u2, u3;
    u1 = US1_a[idx] / hdp->W / hdp->rg / r_a[idx] / (1.0 + w_EOS);
    u2 = US2_a[idx] / hdp->W / hdp->rg / r_a[idx] / (1.0 + w_EOS);
    u3 = US3_a[idx] / hdp->W / hdp->rg / r_a[idx] / (1.0 + w_EOS);

    /* velocities (contravariant) */
    v1_a[idx] = paq->alpha / hdp->W * hdp->g * ( paq->gammai11*u1 + paq->gammai12*u2 + paq->gammai13*u3 ) - paq->beta1;
    v2_a[idx] = paq->alpha / hdp->W * hdp->g * ( paq->gammai12*u1 + paq->gammai22*u2 + paq->gammai23*u3 ) - paq->beta2;
    v3_a[idx] = paq->alpha / hdp->W * hdp->g * ( paq->gammai13*u1 + paq->gammai23*u2 + paq->gammai33*u3 ) - paq->beta3;
  }

if(idx==640 || idx==INDEX(14,13,10)|| idx==INDEX(14,13,9))
{
  std::cout << "\n W  = " << hdp->W
            << "\n UD = " << UD_a[idx]
            << "\n";
}

}

void Hydro::setFluxesCell(BSSNData *paq, HydroData *hdp)
{
  idx_t idx = paq->idx;
  idx_t f_idx_1 = F_NP_INDEX(paq->i, paq->j, paq->k, 1);
  idx_t f_idx_2 = F_NP_INDEX(paq->i, paq->j, paq->k, 2);
  idx_t f_idx_3 = F_NP_INDEX(paq->i, paq->j, paq->k, 3);

  // Calculate fluxes in each cell
  // D
    FD_a[f_idx_1] = UD_a[idx]*v1_a[idx];
    FD_a[f_idx_2] = UD_a[idx]*v2_a[idx];
    FD_a[f_idx_3] = UD_a[idx]*v3_a[idx];
  // S1
    FS1_a[f_idx_1] = US1_a[idx]*v1_a[idx]/paq->alpha
                      + w_EOS * r_a[idx] * hdp->rg * paq->alpha;
    FS1_a[f_idx_2] = US1_a[idx]*v2_a[idx]/paq->alpha;
    FS1_a[f_idx_3] = US1_a[idx]*v3_a[idx]/paq->alpha;
  // S2
    FS2_a[f_idx_1] = US2_a[idx]*v1_a[idx]/paq->alpha;
    FS2_a[f_idx_2] = US2_a[idx]*v2_a[idx]/paq->alpha
                      + w_EOS * r_a[idx] * hdp->rg * paq->alpha;
    FS2_a[f_idx_3] = US2_a[idx]*v3_a[idx]/paq->alpha;
  // S3
    FS3_a[f_idx_1] = US3_a[idx]*v1_a[idx]/paq->alpha;
    FS3_a[f_idx_2] = US3_a[idx]*v2_a[idx]/paq->alpha;
    FS3_a[f_idx_3] = US3_a[idx]*v3_a[idx]/paq->alpha
                      + w_EOS * r_a[idx] * hdp->rg * paq->alpha;
}

void Hydro::setSourcesCell(BSSNData *paq, HydroData *hdp)
{
  idx_t idx = paq->idx;
  real_t pref = 0.5*paq->alpha*hdp->rg*r_a[idx];
  real_t uu_fac = (w_EOS+1.0)*pw2(hdp->W)/pw2(paq->alpha);

  // 0.5*Alpha*gamma^1/2*T*g,j
  SS1_a[idx] = pref*TABGAB_J(1);
  SS2_a[idx] = pref*TABGAB_J(2);
  SS3_a[idx] = pref*TABGAB_J(3);
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
  // => determines direction to take stencils across
  real_t a = SIGN( (F_ARR[f_idx] - F_ARR[f_idx_p1])
             / (U_ARR[u_idx] - U_ARR[u_idx_p1]) );

  // 2) 
  // Two stencils are always the same; then one upwind stencil
  // is determined depending on the sign of `a'

  p1 = -F_ARR[f_idx_m1]/6.0 + 5.0/6.0*F_ARR[f_idx] + F_ARR[f_idx_p1]/3.0;
  p2 = F_ARR[f_idx]/3.0 + 5.0/6.0*F_ARR[f_idx_p1] + -F_ARR[f_idx_p2]/3.0;

  if(a > 0.0)
  {
    p0 = F_ARR[f_idx_m2]/3.0 - 7.0/6.0*F_ARR[f_idx_m1] + 11.0/6.0*F_ARR[f_idx];
  }
  else
  {
    p0 = F_ARR[f_idx_p3]/3.0 - 7.0/6.0*F_ARR[f_idx_p2] + 11.0/6.0*F_ARR[f_idx_p1];
  }

  // 3) calculate "smoothness" indicators
  b1 = 13.0/12.0 * pw2(F_ARR[f_idx_m1] - 2.0*F_ARR[f_idx] + F_ARR[f_idx_p1])
       + 1.0/4.0 * pw2(3.0*F_ARR[f_idx_m1] - F_ARR[f_idx_p1]);
  b2 = 13.0/12.0 * pw2(F_ARR[f_idx] - 2.0*F_ARR[f_idx_p1] + F_ARR[f_idx_p2])
       + 1.0/4.0 * pw2(F_ARR[f_idx] - 4.0*F_ARR[f_idx_p1] - F_ARR[f_idx_p2]);

  if(a > 0.0)
  {
    b0 = 13.0/12.0 * pw2(F_ARR[f_idx_m2] - 2.0*F_ARR[f_idx_m1] + F_ARR[f_idx])
       + 1.0/4.0 * pw2(3.0*F_ARR[f_idx_m2] - 4.0*F_ARR[f_idx_m1] + F_ARR[f_idx]);
  }
  else
  {
     b0 = 13.0/12.0 * pw2(F_ARR[f_idx_p3] - 2.0*F_ARR[f_idx_p2] + F_ARR[f_idx_p1])
       + 1.0/4.0 * pw2(3.0*F_ARR[f_idx_p3] - 4.0*F_ARR[f_idx_p2] + F_ARR[f_idx_p1]);
  }

  // initial 5th-order accurate weights
  if(a > 0.0)
  {
    g1 = 0.6;
    g2 = 0.3;
  }
  else
  {
    g1 = 0.3;
    g2 = 0.6;
  }
  g0 = 0.1;

  // nonlinear weights:
  ebtot2 = pw2(EPS + b1) + pw2(EPS + b2) + pw2(EPS + b0);
  w1 = g1 / ebtot2;
  w2 = g2 / ebtot2;
  w0 = g0 / ebtot2;
  wtot = w1 + w2 + w0;
  w1 /= wtot;
  w2 /= wtot;
  w0 /= wtot;

  // Final flux
  real_t flux;
  flux = w1*p1 + w1*p2 + w0*p0;

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
  UD_f[idx] = exp(log(UD_a[idx]) + dt/UD_a[idx]*(
      FD_int_a[F_NP_INDEX(i,j,k,1)] + FD_int_a[F_NP_INDEX(i,j,k,2)] + FD_int_a[F_NP_INDEX(i,j,k,3)]
      - FD_int_a[F_INDEX(i-1,j,k,1)] - FD_int_a[F_INDEX(i,j-1,k,2)] - FD_int_a[F_INDEX(i,j,k-1,3)]
    ));

  US1_f[idx] = US1_a[idx] + dt*(
      FS1_int_a[F_NP_INDEX(i,j,k,1)] + FS1_int_a[F_NP_INDEX(i,j,k,2)] + FS1_int_a[F_NP_INDEX(i,j,k,3)]
      - FS1_int_a[F_INDEX(i-1,j,k,1)] - FS1_int_a[F_INDEX(i,j-1,k,2)] - FS1_int_a[F_INDEX(i,j,k-1,3)]
    );

  US2_f[idx] = US2_a[idx] + dt*(
      FS2_int_a[F_NP_INDEX(i,j,k,1)] + FS2_int_a[F_NP_INDEX(i,j,k,2)] + FS2_int_a[F_NP_INDEX(i,j,k,3)]
      - FS2_int_a[F_INDEX(i-1,j,k,1)] - FS2_int_a[F_INDEX(i,j-1,k,2)] - FS2_int_a[F_INDEX(i,j,k-1,3)]
    );

  US3_f[idx] = US3_a[idx] + dt*(
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

    real_t a = exp(4.0*bssn_fields["phi_a"][idx]);
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
            ) / UD_a[idx] / UD_a[idx] / pw2(1.0 + w_EOS)
          );

      r =  UD_a[idx] / rg / W_rel;
      u1 = US1_a[idx] / W_rel / rg / r / (1.0 + w_EOS);
      u2 = US2_a[idx] / W_rel / rg / r / (1.0 + w_EOS);
      u3 = US3_a[idx] / W_rel / rg / r / (1.0 + w_EOS);
    }

    bssn_fields["r_a"][idx] += r*((1.0 + w_EOS)*W_rel*W_rel - w_EOS);

    bssn_fields["S1_a"][idx] += r*(1.0 + w_EOS)*W_rel*u1;
    bssn_fields["S2_a"][idx] += r*(1.0 + w_EOS)*W_rel*u2;
    bssn_fields["S3_a"][idx] += r*(1.0 + w_EOS)*W_rel*u3;

    real_t S = r*( 3.0*w_EOS + (1.0 + w_EOS)*(W_rel*W_rel-1.0) );
    bssn_fields["S_a"][idx] += S;

    bssn_fields["STF11_a"][idx] += r*(1.0 + w_EOS)*(u1*u1 - a*g11*(W_rel*W_rel - 1.0));
    bssn_fields["STF12_a"][idx] += r*(1.0 + w_EOS)*(u1*u2 - a*g12*(W_rel*W_rel - 1.0));
    bssn_fields["STF13_a"][idx] += r*(1.0 + w_EOS)*(u1*u3 - a*g13*(W_rel*W_rel - 1.0));
    bssn_fields["STF22_a"][idx] += r*(1.0 + w_EOS)*(u2*u2 - a*g22*(W_rel*W_rel - 1.0));
    bssn_fields["STF23_a"][idx] += r*(1.0 + w_EOS)*(u2*u3 - a*g23*(W_rel*W_rel - 1.0));
    bssn_fields["STF33_a"][idx] += r*(1.0 + w_EOS)*(u3*u3 - a*g33*(W_rel*W_rel - 1.0));

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
