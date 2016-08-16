#include "static_ic.h"
#include "../../cosmo_includes.h"
#include "../../cosmo_types.h"
#include "../../cosmo_globals.h"
#include "../../ICs/ICs.h"
#include "../../utils/math.h"

namespace cosmo
{

/**
 * @brief Gaussian random field ICs
 */
void dust_ic_set_random(BSSN * bssn, Static * dust, Fourier * fourier,
  IOData * iodata)
{
  idx_t i, j, k;

  arr_t & DIFFr_a = *bssn->fields["DIFFr_a"];
  arr_t & DIFFphi_p = *bssn->fields["DIFFphi_p"];
  arr_t & DIFFphi_a = *bssn->fields["DIFFphi_a"];
  arr_t & DIFFphi_f = *bssn->fields["DIFFphi_f"];

  arr_t & DIFFD_a = *dust->fields["DIFFD_a"];

  auto & frw = bssn->frw;

  ICsData icd = cosmo_get_ICsData();
  iodata->log( "Generating ICs with peak at k = " + stringify(icd.peak_k) );
  iodata->log( "Generating ICs with peak amp. = " + stringify(icd.peak_amplitude) );

  // the conformal factor in front of metric is the solution to
  // d^2 exp(\phi) = -2*pi exp(5\phi) * \delta_rho
  // generate gaussian random field 1 + xi = exp(phi) (store xi in phi_p temporarily):
  set_gaussian_random_field(DIFFphi_p, fourier, &icd);

  // delta_rho = -lap(phi)/(1+xi)^5/2pi
# pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k) {
    DIFFr_a[NP_INDEX(i,j,k)] = -0.5/PI/(
      pow(1.0 + DIFFphi_p[NP_INDEX(i,j,k)], 5.0)
    )*(
      double_derivative(i, j, k, 1, 1, DIFFphi_p)
      + double_derivative(i, j, k, 2, 2, DIFFphi_p)
      + double_derivative(i, j, k, 3, 3, DIFFphi_p)
    );
  }

  // phi = ln(xi)
# pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k) {
    idx_t idx = NP_INDEX(i,j,k);
    DIFFphi_a[idx] = log1p(DIFFphi_p[idx]);
    DIFFphi_f[idx] = log1p(DIFFphi_p[idx]);
    DIFFphi_p[idx] = log1p(DIFFphi_p[idx]);
  }

  // Make sure min density value > 0
  // Set conserved density variable field
  real_t min = icd.rho_K_matter;
  real_t max = min;
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    real_t rho_FRW = icd.rho_K_matter;
    real_t DIFFr = DIFFr_a[idx];
    real_t rho = rho_FRW + DIFFr;
    // phi_FRW = 0
    real_t DIFFphi = DIFFphi_a[idx];
    // phi = DIFFphi
    // DIFFK = 0

    DIFFD_a[idx] =
      rho_FRW*expm1(6.0*DIFFphi) + exp(6.0*DIFFphi)*DIFFr;

    if(rho < min)
    {
      min = rho;
    }
    if(rho > max)
    {
      max = rho;
    }
    if(rho != rho)
    {
      iodata->log("Error: NaN energy density.");
      throw -1;
    }
  }

  iodata->log( "Minimum fluid density: " + stringify(min) );
  iodata->log( "Maximum fluid density: " + stringify(max) );
  iodata->log( "Average fluctuation density: " + stringify(average(DIFFD_a)) );
  iodata->log( "Std.dev fluctuation density: " + stringify(standard_deviation(DIFFD_a)) );
  if(min < 0.0)
  {
    iodata->log("Error: negative density in some regions.");
    throw -1;
  }

# if USE_REFERENCE_FRW
  // Set values in reference FRW integrator
  real_t rho_FRW = icd.rho_K_matter;
  real_t K_frw = -sqrt(24.0*PI*rho_FRW);

  frw->set_phi(0.0);
  frw->set_K(K_frw);
  frw->addFluid(rho_FRW, 0.0 /* w=0 */);
# else
  arr_t & DIFFK_p = *bssn->fields["DIFFK_p"];
  arr_t & DIFFK_a = *bssn->fields["DIFFK_a"];
  // add in FRW pieces to ICs
  // phi is unchanged
  // rho (D) and K get contribs
  // w=0 fluid only
# pragma omp parallel for default(shared) private(i,j,k)
  LOOP3(i,j,k)
  {
    idx_t idx = NP_INDEX(i,j,k);
    real_t rho_FRW = icd.rho_K_matter;
    real_t D_FRW = rho_FRW; // on initial slice

    DIFFr_a[idx] += rho_FRW;

    DIFFK_a[idx] = -sqrt(24.0*PI*rho_FRW);
    DIFFK_p[idx] = -sqrt(24.0*PI*rho_FRW);

    DIFFD_a[idx] += D_FRW;
  }
# endif
}

}
