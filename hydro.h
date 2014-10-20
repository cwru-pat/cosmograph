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
  /* equation of state "w" value */
  real_t w_EOS;

  /* Fluid fields */
  HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_CREATE)
  HYDRO_APPLY_TO_FLUXES(FLUX_ARRAY_CREATE)
  HYDRO_APPLY_TO_SOURCES(GEN1_ARRAY_CREATE)
  /* primitives */
  HYDRO_APPLY_TO_PRIMITIVES(GEN1_ARRAY_CREATE)
  HYDRO_APPLY_TO_FLUXES_INT(FLUX_ARRAY_CREATE)

public:
  std::map <std::string, real_t *> fields;

  Hydro()
  {
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
  ~Hydro()
  {
    HYDRO_APPLY_TO_FIELDS(GEN2_ARRAY_DELETE)
    HYDRO_APPLY_TO_FLUXES(FLUX_ARRAY_DELETE)
    HYDRO_APPLY_TO_SOURCES(GEN1_ARRAY_DELETE)
    HYDRO_APPLY_TO_PRIMITIVES(GEN1_ARRAY_DELETE)
    HYDRO_APPLY_TO_FLUXES_INT(FLUX_ARRAY_DELETE)
  }

  void setPrimitivesPt(BSSNData *paq, HydroData *hdp)
  {
    idx_t idx = paq->idx;

    /* root of metric determinant */
    hdp->rg = exp(2.0*paq->phi);

    /* lorentz factor */
    /* W = (gamma^ij S_j S_j / D^2 / (1+w^2) + 1)^1/2 */
    hdp->W = sqrt(
        1.0 + (
           paq->gammai11*US1_a[idx]*US1_a[idx] + paq->gammai22*US2_a[idx]*US2_a[idx] + paq->gammai33*US3_a[idx]*US3_a[idx]
            + 2.0*paq->gammai12*US1_a[idx]*US2_a[idx] + 2.0*paq->gammai13*US1_a[idx]*US3_a[idx] + 2.0*paq->gammai23*US2_a[idx]*US3_a[idx]
          ) / UD_a[idx] / (1.0 + w_EOS*w_EOS)
        );

    /* \rho = D / gamma^1/2 / W */
    r_a[idx] = UD_a[idx] / hdp->rg / hdp->W;

    /* fluid 4-velocities (covariant) */
    real_t u1, u2, u3;
    u1 = US1_a[idx] / hdp->W / hdp->rg / r_a[idx] / (1.0 + w_EOS);
    u2 = US2_a[idx] / hdp->W / hdp->rg / r_a[idx] / (1.0 + w_EOS);
    u3 = US3_a[idx] / hdp->W / hdp->rg / r_a[idx] / (1.0 + w_EOS);
    /* velocities (contravariant) */
    v1_a[idx] = paq->alpha / hdp->W * ( paq->gammai11*u1 + paq->gammai12*u2 + paq->gammai13*u3 ) - paq->beta1;
    v2_a[idx] = paq->alpha / hdp->W * ( paq->gammai12*u1 + paq->gammai22*u2 + paq->gammai23*u3 ) - paq->beta2;
    v3_a[idx] = paq->alpha / hdp->W * ( paq->gammai13*u1 + paq->gammai23*u2 + paq->gammai33*u3 ) - paq->beta3;

  }

  void setFluxesPt(BSSNData *paq, HydroData *hdp)
  {
    idx_t idx = paq->idx;
    idx_t f_idx_1 = F_INDEX(paq->i, paq->j, paq->k, 0);
    idx_t f_idx_2 = F_INDEX(paq->i, paq->j, paq->k, 1);
    idx_t f_idx_3 = F_INDEX(paq->i, paq->j, paq->k, 2);

    // Calculate fluxes in each cell
    // D
      FD_a[f_idx_1] = FD_a[idx]*v1_a[idx];
      FD_a[f_idx_2] = FD_a[idx]*v2_a[idx];
      FD_a[f_idx_3] = FD_a[idx]*v3_a[idx];
    // S1
      FS1_a[f_idx_1] = FS1_a[idx]*v1_a[idx]/paq->alpha
                        + hdp->W * r_a[idx] * hdp->rg * paq->alpha;
      FS1_a[f_idx_2] = FS1_a[idx]*v2_a[idx]/paq->alpha;
      FS1_a[f_idx_3] = FS1_a[idx]*v3_a[idx]/paq->alpha;
    // S2
      FS2_a[f_idx_1] = FS2_a[idx]*v1_a[idx]/paq->alpha;
      FS2_a[f_idx_2] = FS2_a[idx]*v2_a[idx]/paq->alpha
                        + hdp->W * r_a[idx] * hdp->rg * paq->alpha;
      FS2_a[f_idx_3] = FS2_a[idx]*v3_a[idx]/paq->alpha;
    // S3
      FS3_a[f_idx_1] = FS3_a[idx]*v1_a[idx]/paq->alpha;
      FS3_a[f_idx_2] = FS3_a[idx]*v2_a[idx]/paq->alpha;
      FS3_a[f_idx_3] = FS3_a[idx]*v3_a[idx]/paq->alpha
                        + hdp->W * r_a[idx] * hdp->rg * paq->alpha;
  }

  void setSourcesPt(BSSNData *paq, HydroData *hdp)
  {
    idx_t idx = paq->idx;
    real_t pref = 0.5*paq->alpha*hdp->rg;

    // 0.5*Alpha*gamma^1/2*T*g,j
    SS1[idx] = 0;
    SS2[idx] = 0;
    SS3[idx] = 0;
  }

  void setOneFluxInt(idx_t i, idx_t j, idx_t k, int d, real_t *U_ARR,
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
    
  void setAllFluxInt(idx_t i, idx_t j, idx_t k)
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

  void init()
  {
    // initialize values

  }

};

}

#endif
