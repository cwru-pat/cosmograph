#include "ICs.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

namespace cosmo
{

void set_BH_ICs(std::map <std::string, real_t *> & bssn_fields)
{
  real_t M = 5.0;
  const idx_t spline_res = 20;
  const idx_t spline_size = spline_res*N;
  idx_t i, j, k;

  // first need to get r_s in terms of r
  real_t r_vals[spline_size], rs_vals[spline_size];

  for(i = 0; i<spline_size; ++i)
  {
    rs_vals[i] = 3.0/2.0*M*(1.0 + sqrt(3.0)*i/spline_res + 0.0000001);
    r_vals[i] = (2.0*rs_vals[i] + M + sqrt(4.0*rs_vals[i]*rs_vals[i] + 4.0*M*rs_vals[i] + 3.0*M*M))/4.0
           * pow(
              (4.0+3.0*sqrt(2.0))*(2.0*rs_vals[i] - 3.0*M)
              /(8.0*rs_vals[i] + 6.0*M + 3.0*sqrt(8.0*rs_vals[i]*rs_vals[i] + 8.0*M*rs_vals[i] + 6.0*M*M))
            , 1.0/sqrt(2.0));
  }

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, spline_size);
  gsl_spline_init(spline, r_vals, rs_vals, spline_size);

  std::cout << "First spline vals are: r[0]=" << r_vals[0] << ", rs[0]=" << rs_vals[0] << "; "
            << "r[1]=" << r_vals[1] << ", rs[1]=" << rs_vals[1] << "\n"
            << "Spline value is: " << gsl_spline_eval(spline, 10.0, acc) << "\n";

  std::cout << "psi slice is: ";
  idx_t i0 = N/2, j0 = N/2, k0 = N/2;
  LOOP3(i, j, k)
  {
    // set psi, \beta, and \alpha values
    // (everything else is 0)
    real_t r = sqrt((real_t) ((i-i0)*(i-i0) + (j-j0)*(j-j0) + (k-k0)*(k-k0)))+0.0001;
    real_t rs = gsl_spline_eval(spline, r, acc);

    real_t psi = sqrt(
        4.0*rs/(2.0*rs + M + sqrt(4.0*rs*rs + 4.0*rs*M + 3.0*M*M))
      ) * pow(
        (8.0*rs + 6.0*M + 3.0*sqrt(8.0*rs*rs + 8.0*M*rs + 6.0*M*M))/(4.0+3.0*sqrt(2.0))/(2.0*rs - 3.0*M)
      , 1.0/2.0/sqrt(2));
    if(i==N/2&&j==N/2){
      std::cout << psi << ", ";
    }
    bssn_fields["phi_p"][NP_INDEX(i,j,k)] = log(psi);

    real_t alpha = sqrt(1.0-2.0*M/rs + 27.0/16.0*pow(M/rs,4));
    bssn_fields["alpha_p"][NP_INDEX(i,j,k)] = alpha;

    real_t beta_r = 3.0*sqrt(3)*M*M/4.0*r/pow(rs,3);
    bssn_fields["beta1_p"][NP_INDEX(i,j,k)] = (i-i0)/r*beta_r;
    bssn_fields["beta2_p"][NP_INDEX(i,j,k)] = (j-j0)/r*beta_r;
    bssn_fields["beta3_p"][NP_INDEX(i,j,k)] = (k-k0)/r*beta_r;
  }
  std::cout << "\n";

}

} // namespace cosmo
