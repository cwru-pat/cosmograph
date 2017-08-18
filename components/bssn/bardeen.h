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

public:  
  // perturbed metric & time derivatives
  arr_t h11, h12, h13, h22, h23, h33;
  arr_t dt_h11, dt_h12, dt_h13, dt_h22, dt_h23, dt_h33;
  arr_t d2t_h11, d2t_h12, d2t_h13, d2t_h22, d2t_h23, d2t_h33;

  // metric scalar potentials
  arr_t A, dt_A, d2t_A;
  arr_t B, dt_B, d2t_B;
  // Bardeen potentials
  arr_t Phi, Psi;

  arr_t tmp;

  Bardeen(BSSN * bssn_in, Fourier * fourier_in)
  {
    bssn = bssn_in;
    fourier = fourier_in;

    h11.init(NX, NY, NZ); h12.init(NX, NY, NZ); h13.init(NX, NY, NZ);
    h22.init(NX, NY, NZ); h23.init(NX, NY, NZ); h33.init(NX, NY, NZ);
    dt_h11.init(NX, NY, NZ); dt_h12.init(NX, NY, NZ); dt_h13.init(NX, NY, NZ);
    dt_h22.init(NX, NY, NZ); dt_h23.init(NX, NY, NZ); dt_h33.init(NX, NY, NZ);
    d2t_h11.init(NX, NY, NZ); d2t_h12.init(NX, NY, NZ); d2t_h13.init(NX, NY, NZ);
    d2t_h22.init(NX, NY, NZ); d2t_h23.init(NX, NY, NZ); d2t_h33.init(NX, NY, NZ);

    A.init(NX, NY, NZ); dt_A.init(NX, NY, NZ); d2t_A.init(NX, NY, NZ);
    B.init(NX, NY, NZ); dt_B.init(NX, NY, NZ); d2t_B.init(NX, NY, NZ);
    Phi.init(NX, NY, NZ); Psi.init(NX, NY, NZ);

    tmp.init(NX, NY, NZ);

    // add Bardeen potentials to BSSN fields map
    bssn->fields["Bardeen_Phi"] = & Phi;
    bssn->fields["Bardeen_Psi"] = & Psi;
    bssn->fields["Bardeen_A"] = & A;
    bssn->fields["Bardeen_dt_A"] = & dt_A;
    bssn->fields["Bardeen_d2t_A"] = & d2t_A;
    bssn->fields["Bardeen_B"] = & B;
    bssn->fields["Bardeen_dt_B"] = & dt_B;
    bssn->fields["Bardeen_d2t_B"] = & d2t_B;
  }
  
  ~Bardeen()
  {
    // anything to do?
  }

  void setPotentials();

}; // Bardeen class

}

#endif
