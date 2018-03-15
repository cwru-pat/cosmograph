/** @file bardeen.h
 * @brief compute Bardeen potentials for metric.
 * @ Currently works in synchronous gauge only; ie. with dust field.
 */

#ifndef COSMO_BARDEEN
#define COSMO_BARDEEN

#include "../../cosmo_includes.h"
#include "../../cosmo_macros.h"
#include "../../cosmo_types.h"
#include "../../utils/Fourier.h"
#include "bssn.h"

namespace cosmo
{

/**
 * Compute bardeen potentials
 * ** Restricted to synchronous gauge evolution
 * ** Reference metric not supported
 */
class Bardeen
{
  BSSN * bssn;
  Fourier * fourier;

  bool use_matter_scale_factor;

public:  
  // perturbed metric & time derivatives
  arr_t h11, h12, h13, h22, h23, h33;
  arr_t dt_h11, dt_h12, dt_h13, dt_h22, dt_h23, dt_h33;
  arr_t d2t_h11, d2t_h12, d2t_h13, d2t_h22, d2t_h23, d2t_h33;
  arr_t h01, h02, h03;
  arr_t dt_h01, dt_h02, dt_h03;

  // storage for intermediate metric calc's
  arr_t dt_g11, dt_g12, dt_g13,
        dt_g22, dt_g23, dt_g33;
  arr_t d2t_g11, d2t_g12, d2t_g13,
        d2t_g22, d2t_g23, d2t_g33;
  arr_t dt_beta1, dt_beta2, dt_beta3;
  arr_t dt_phi, d2t_phi;

  // metric scalar potentials
  arr_t A, dt_A, d2t_A;
  arr_t B, dt_B, d2t_B;
  arr_t F, dt_F, // higher-order time derivs. of these
        E;       // aren't needed (for now).

  // Vector potentials
  arr_t G1, G2, G3;
  arr_t C1, C2, C3;

  // Tensor
  arr_t D11, D12, D13, D22, D23, D33;

  // Bardeen potentials
  arr_t Phi, Psi;

  // constraint violation
  arr_t lin_viol, lin_viol_der_mag, lin_viol_der;

  real_t * viols;

  Bardeen(BSSN * bssn_in, Fourier * fourier_in)
  {
    bssn = bssn_in;
    fourier = fourier_in;

    use_matter_scale_factor = true;

    h11.init(NX, NY, NZ); h12.init(NX, NY, NZ); h13.init(NX, NY, NZ);
    h22.init(NX, NY, NZ); h23.init(NX, NY, NZ); h33.init(NX, NY, NZ);
    dt_h11.init(NX, NY, NZ); dt_h12.init(NX, NY, NZ); dt_h13.init(NX, NY, NZ);
    dt_h22.init(NX, NY, NZ); dt_h23.init(NX, NY, NZ); dt_h33.init(NX, NY, NZ);
    d2t_h11.init(NX, NY, NZ); d2t_h12.init(NX, NY, NZ); d2t_h13.init(NX, NY, NZ);
    d2t_h22.init(NX, NY, NZ); d2t_h23.init(NX, NY, NZ); d2t_h33.init(NX, NY, NZ);

    h01.init(NX, NY, NZ); h02.init(NX, NY, NZ); h03.init(NX, NY, NZ);
    dt_h01.init(NX, NY, NZ); dt_h02.init(NX, NY, NZ); dt_h03.init(NX, NY, NZ);

    dt_g11.init(NX, NY, NZ); dt_g12.init(NX, NY, NZ); dt_g13.init(NX, NY, NZ);
    dt_g22.init(NX, NY, NZ); dt_g23.init(NX, NY, NZ); dt_g33.init(NX, NY, NZ);
    d2t_g11.init(NX, NY, NZ); d2t_g12.init(NX, NY, NZ); d2t_g13.init(NX, NY, NZ);
    d2t_g22.init(NX, NY, NZ); d2t_g23.init(NX, NY, NZ); d2t_g33.init(NX, NY, NZ);
    dt_beta1.init(NX, NY, NZ); dt_beta2.init(NX, NY, NZ); dt_beta3.init(NX, NY, NZ);
    dt_phi.init(NX, NY, NZ); d2t_phi.init(NX, NY, NZ);

    A.init(NX, NY, NZ); dt_A.init(NX, NY, NZ); d2t_A.init(NX, NY, NZ);
    B.init(NX, NY, NZ); dt_B.init(NX, NY, NZ); d2t_B.init(NX, NY, NZ);
    F.init(NX, NY, NZ); dt_F.init(NX, NY, NZ);
    E.init(NX, NY, NZ);

    G1.init(NX, NY, NZ); G2.init(NX, NY, NZ); G3.init(NX, NY, NZ);
    C1.init(NX, NY, NZ); C2.init(NX, NY, NZ); C3.init(NX, NY, NZ);

    D11.init(NX, NY, NZ);
    D12.init(NX, NY, NZ);
    D13.init(NX, NY, NZ);
    D22.init(NX, NY, NZ);
    D23.init(NX, NY, NZ);
    D33.init(NX, NY, NZ);

    Phi.init(NX, NY, NZ); Psi.init(NX, NY, NZ);

    lin_viol.init(NX, NY, NZ);
    lin_viol_der_mag.init(NX, NY, NZ);
    lin_viol_der.init(NX, NY, NZ);

    viols = new real_t[7];
    for(int i=0; i<7; ++i)
      viols[0] = 0;

    // add Bardeen potentials to BSSN fields map
    bssn->fields["Bardeen_Phi"] = & Phi;
    bssn->fields["Bardeen_Psi"] = & Psi;
    bssn->fields["Bardeen_A"] = & A;
    bssn->fields["Bardeen_dt_A"] = & dt_A;
    bssn->fields["Bardeen_d2t_A"] = & d2t_A;
    bssn->fields["Bardeen_B"] = & B;
    bssn->fields["Bardeen_dt_B"] = & dt_B;
    bssn->fields["Bardeen_d2t_B"] = & d2t_B;

    bssn->fields["Bardeen_E"] = & E;
    bssn->fields["Bardeen_F"] = & F;

    bssn->fields["Bardeen_G1"] = & G1;
    bssn->fields["Bardeen_G2"] = & G2;
    bssn->fields["Bardeen_G3"] = & G3;
    
    bssn->fields["Bardeen_C1"] = & C1;
    bssn->fields["Bardeen_C2"] = & C2;
    bssn->fields["Bardeen_C3"] = & C3;

    bssn->fields["Bardeen_D11"] = & D11;
    bssn->fields["Bardeen_D12"] = & D12;
    bssn->fields["Bardeen_D13"] = & D13;
    bssn->fields["Bardeen_D22"] = & D22;
    bssn->fields["Bardeen_D23"] = & D23;
    bssn->fields["Bardeen_D33"] = & D33;
  }
  
  ~Bardeen()
  {
    // anything to do?
  }

  void setUseMatterScaleFactor(bool use)
  {
    use_matter_scale_factor = use;
  }

  void setPotentials(real_t elapsed_sim_time);

  void getSVTViolations(real_t * viols_copyto);

}; // Bardeen class

}

#endif
