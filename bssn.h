#ifndef COSMO_BSSN
#define COSMO_BSSN

#include BSSN_macros.h

/** BSSN class **/
class BSSN
{
  /* create arrays for storing fields */
    BSSN_APPLY_TO_FIELDS(RK4_ARRAY_CREATE);

  /* Local variables: */

      // generic var for misc. expressions
      real_t trace, expression;

      // trace-free ricci tensor components
      real_t ricciTF11, ricciTF12, ricciTF13, ricciTF22, ricciTF23, ricciTF33;

      // derivatives of \alpha
        // covariant double-derivatives 
        real_t D1D1a, D1D2a, D1D3a, D2D2a, D2D3a, D3D3a;
        // normal derivatives of
        real_t d1a, d2a, d3a;

      // derivatives of phi
        // covariant double-derivatives 
        real_t D1D1phi, D1D2phi, D1D3phi, D2D2phi, D2D3phi, D3D3phi;
        // normal derivatives of
        real_t d1phi, d2phi, d3phi;

      // Christoffel symbols
      real_t G111, G112, G113, G122, G123, G133,
             G211, G212, G213, G222, G223, G233,
             G311, G312, G313, G322, G323, G333;

      // derivatives of the metric, d_i g_jk
      real_t d1g11, d1g12, d1g13, d1g22, d1g23, d1g33,
             d2g11, d2g12, d2g13, d2g22, d2g23, d2g33,
             d3g11, d3g12, d3g13, d3g22, d3g23, d3g33;

      // derivatives of the inverse metric d_i g^jk
      real_t d1gi11, d1gi12, d1gi13, d1gi22, d1gi23, d1gi33,
             d2gi11, d2gi12, d2gi13, d2gi22, d2gi23, d2gi33,
             d3gi11, d3gi12, d3gi13, d3gi22, d3gi23, d3gi33;

      // second derivatives of the metric d_i d_j g_kl
      real_t d1d1g11, d1d1g12, d1d1g13, d1d1g22, d1d1g23, d1d1g33,
             d1d2g11, d1d2g12, d1d2g13, d1d2g22, d1d2g23, d1d2g33,
             d1d3g11, d1d3g12, d1d3g13, d1d3g22, d1d3g23, d1d3g33,
             d2d2g11, d2d2g12, d2d2g13, d2d2g22, d2d2g23, d2d2g33,
             d2d3g11, d2d3g12, d2d3g13, d2d3g22, d2d3g23, d2d3g33,
             d3d3g11, d3d3g12, d3d3g13, d3d3g22, d3d3g23, d3d3g33;

      // local copies of current field values
      BSSN_APPLY_TO_FIELDS(DECLARE_REAL_T);

    /* End local variable declarations */

public:
  BSSN()
  {
    BSSN_APPLY_TO_FIELDS(RK4_ARRAY_ALLOC);
    BSSN_APPLY_TO_FIELDS(RK4_ARRAY_ADDMAP);
  }

  ~BSSN()
  {
    BSSN_APPLY_TO_FIELDS(RK4_ARRAY_DELETE);
  }
  
  void step() {
    LOOP3(i, j, k)
    {
      /* evolve stuff... */
      
    }

    BSSN_APPLY_TO_FIELDS(RK4_ARRAY_CYCLE);
  }

  /* set current local field values */
  inline void calculate_gammai(idx_t idx)
  {
    BSSN_APPLY_TO_FIELDS(SET_LOCAL_VALUES);
  }

  /* Calculate inverse metric components */
  inline void calculate_gammai(idx_t idx)
  {
    // gamma is unitary, so this is actually simple:
    gammai11 = gamma22*gamma33 - gamma23*gamma23;
    gammai12 = gamma13*gamma23 - gamma12*gamma33;
    gammai13 = gamma12*gamma23 - gamma13*gamma22;
    gammai22 = gamma11*gamma33 - gamma13*gamma13;
    gammai23 = gamma12*gamma13 - gamma23*gamma11;
    gammai33 = gamma11*gamma22 - gamma12*gamma12;
  }

  /* Calculate metric derivatives */
  inline void calculate_dgamma(idx_t idx)
  {
    BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_DGAMMA);

    d1d1g11 = dder(gamma11, 1, 1, idx); 
    d1d1g12 = dder(gamma12, 1, 1, idx); 
    d1d1g13 = dder(gamma13, 1, 1, idx); 
    d1d1g22 = dder(gamma22, 1, 1, idx); 
    d1d1g23 = dder(gamma23, 1, 1, idx); 
    d1d1g33 = dder(gamma33, 1, 1, idx); 
             
    d1d2g11 = dder(gamma11, 1, 2, idx); 
    d1d2g12 = dder(gamma12, 1, 2, idx); 
    d1d2g13 = dder(gamma13, 1, 2, idx); 
    d1d2g22 = dder(gamma22, 1, 2, idx); 
    d1d2g23 = dder(gamma23, 1, 2, idx); 
    d1d2g33 = dder(gamma33, 1, 2, idx); 
             
    d1d3g11 = dder(gamma11, 1, 3, idx); 
    d1d3g12 = dder(gamma12, 1, 3, idx); 
    d1d3g13 = dder(gamma13, 1, 3, idx); 
    d1d3g22 = dder(gamma22, 1, 3, idx); 
    d1d3g23 = dder(gamma23, 1, 3, idx); 
    d1d3g33 = dder(gamma33, 1, 3, idx); 
             
    d2d2g11 = dder(gamma11, 2, 2, idx); 
    d2d2g12 = dder(gamma12, 2, 2, idx); 
    d2d2g13 = dder(gamma13, 2, 2, idx); 
    d2d2g22 = dder(gamma22, 2, 2, idx); 
    d2d2g23 = dder(gamma23, 2, 2, idx); 
    d2d2g33 = dder(gamma33, 2, 2, idx); 
             
    d2d3g11 = dder(gamma11, 2, 3, idx); 
    d2d3g12 = dder(gamma12, 2, 3, idx); 
    d2d3g13 = dder(gamma13, 2, 3, idx); 
    d2d3g22 = dder(gamma22, 2, 3, idx); 
    d2d3g23 = dder(gamma23, 2, 3, idx); 
    d2d3g33 = dder(gamma33, 2, 3, idx); 
             
    d3d3g11 = dder(gamma11, 3, 3, idx); 
    d3d3g12 = dder(gamma12, 3, 3, idx); 
    d3d3g13 = dder(gamma13, 3, 3, idx); 
    d3d3g22 = dder(gamma22, 3, 3, idx); 
    d3d3g23 = dder(gamma23, 3, 3, idx); 
    d3d3g33 = dder(gamma33, 3, 3, idx);

  }

  /* Calculate metric derivatives */
  inline void calculate_dgammai(idx_t idx)
  {
    BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_DGAMMAI);
  }

  inline void calculate_christoffels(idx_t idx)
  {
    // christoffel symbols: \Gamma^i_{jk} = Gijk
    BSSN_APPLY_TO_IJK_PERMS(BSSN_CALCULATE_CHRISTOFFEL);
  }

  /* Calculate trace-free ricci tensor components */
  void calculateRicciTF(idx_t idx)
  {
    ricciTF11 = (
        - 0.5*(
          gammai11*d1d1g11 + gammai22*d2d2g11 + gammai33*d3d3g11
          + 2.0*(gammai12*d1d2g11 + gammai13*d1d3g11 + gammai23*d2d3g11)
        )
        + 0.5*(
          gamma11*der(Gamma1, 1, idx) + gamma21*der(Gamma2, 1, idx) + gamma31*der(Gamma3, 1, idx)
          + gamma11*der(Gamma1, 1, idx) + gamma21*der(Gamma2, 1, idx) + gamma31*der(Gamma3, 1, idx)
        )
        - 0.5*(
          d1g11*d1gi11 + d2g12*d1gi22 + d3g13*d1gi33 + 2.0*(d1g12*d1gi12 + d1g13*d1gi13 + d2g13*d1gi23)
          + d1g11*d1gi11 + d2g12*d1gi22 + d3g13*d1gi33 + 2.0*(d1g12*d1gi12 + d1g13*d1gi13 + d2g13*d1gi23)
          - Gamma1*d1g11 - Gamma2*d2g11 - Gamma3*d3g11
        )
        - (
          G111*G111 + G212*G212 + G313*G313 + 2.0*(G112*G211 + G113*G311 + G213*G312)
        )
      );
    ricciTF12 = (
        - 0.5*(
          gammai11*d1d1g12 + gammai22*d2d2g12 + gammai33*d3d3g12
          + 2.0*(gammai12*d1d2g12 + gammai13*d1d3g12 + gammai23*d2d3g12)
        )
        + 0.5*(
          gamma11*der(Gamma1, 2, idx) + gamma21*der(Gamma2, 2, idx) + gamma31*der(Gamma3, 2, idx)
          + gamma12*der(Gamma1, 1, idx) + gamma22*der(Gamma2, 1, idx) + gamma32*der(Gamma3, 1, idx)
        )
        - 0.5*(
          d1g11*d2gi11 + d2g12*d2gi22 + d3g13*d2gi33 + 2.0*(d1g12*d2gi12 + d1g13*d2gi13 + d2g13*d2gi23)
          + d1g21*d1gi11 + d2g22*d1gi22 + d3g23*d1gi33 + 2.0*(d1g22*d1gi12 + d1g23*d1gi13 + d2g23*d1gi23)
          - Gamma1*d1g12 - Gamma2*d2g12 - Gamma3*d3g12
        )
        - (
          G111*G121 + G212*G222 + G313*G323 + 2.0*(G112*G221 + G113*G321 + G213*G322)
        )
      );
    ricciTF13 = (
        - 0.5*(
          gammai11*d1d1g13 + gammai22*d2d2g13 + gammai33*d3d3g13
          + 2.0*(gammai12*d1d2g13 + gammai13*d1d3g13 + gammai23*d2d3g13)
        )
        + 0.5*(
          gamma11*der(Gamma1, 3, idx) + gamma21*der(Gamma2, 3, idx) + gamma31*der(Gamma3, 3, idx)
          + gamma13*der(Gamma1, 1, idx) + gamma23*der(Gamma2, 1, idx) + gamma33*der(Gamma3, 1, idx)
        )
        - 0.5*(
          d1g11*d3gi11 + d2g12*d3gi22 + d3g13*d3gi33 + 2.0*(d1g12*d3gi12 + d1g13*d3gi13 + d2g13*d3gi23)
          + d1g31*d1gi11 + d2g32*d1gi22 + d3g33*d1gi33 + 2.0*(d1g32*d1gi12 + d1g33*d1gi13 + d2g33*d1gi23)
          - Gamma1*d1g13 - Gamma2*d2g13 - Gamma3*d3g13
        )
        - (
          G111*G131 + G212*G232 + G313*G333 + 2.0*(G112*G231 + G113*G331 + G213*G332)
        )
      );
    ricciTF22 = (
        - 0.5*(
          gammai11*d1d1g22 + gammai22*d2d2g22 + gammai33*d3d3g22
          + 2.0*(gammai12*d1d2g22 + gammai13*d1d3g22 + gammai23*d2d3g22)
        )
        + 0.5*(
          gamma12*der(Gamma1, 2, idx) + gamma22*der(Gamma2, 2, idx) + gamma32*der(Gamma3, 2, idx)
          + gamma12*der(Gamma1, 2, idx) + gamma22*der(Gamma2, 2, idx) + gamma32*der(Gamma3, 2, idx)
        )
        - 0.5*(
          d1g21*d2gi11 + d2g22*d2gi22 + d3g23*d2gi33 + 2.0*(d1g22*d2gi12 + d1g23*d2gi13 + d2g23*d2gi23)
          + d1g21*d2gi11 + d2g22*d2gi22 + d3g23*d2gi33 + 2.0*(d1g22*d2gi12 + d1g23*d2gi13 + d2g23*d2gi23)
          - Gamma1*d1g22 - Gamma2*d2g22 - Gamma3*d3g22
        )
        - (
          G121*G121 + G222*G222 + G323*G323 + 2.0*(G122*G221 + G123*G321 + G223*G322)
        )
      );
    ricciTF23 = (
        - 0.5*(
          gammai11*d1d1g23 + gammai22*d2d2g23 + gammai33*d3d3g23
          + 2.0*(gammai12*d1d2g23 + gammai13*d1d3g23 + gammai23*d2d3g23)
        )
        + 0.5*(
          gamma12*der(Gamma1, 3, idx) + gamma22*der(Gamma2, 3, idx) + gamma32*der(Gamma3, 3, idx)
          + gamma13*der(Gamma1, 2, idx) + gamma23*der(Gamma2, 2, idx) + gamma33*der(Gamma3, 2, idx)
        )
        - 0.5*(
          d1g21*d3gi11 + d2g22*d3gi22 + d3g23*d3gi33 + 2.0*(d1g22*d3gi12 + d1g23*d3gi13 + d2g23*d3gi23)
          + d1g31*d2gi11 + d2g32*d2gi22 + d3g33*d2gi33 + 2.0*(d1g32*d2gi12 + d1g33*d2gi13 + d2g33*d2gi23)
          - Gamma1*d1g23 - Gamma2*d2g23 - Gamma3*d3g23
        )
        - (
          G121*G131 + G222*G232 + G323*G333 + 2.0*(G122*G231 + G123*G331 + G223*G332)
        )
      );
    ricciTF33 = (
        - 0.5*(
          gammai11*d1d1g33 + gammai22*d2d2g33 + gammai33*d3d3g33
          + 2.0*(gammai12*d1d2g33 + gammai13*d1d3g33 + gammai23*d2d3g33)
        )
        + 0.5*(
          gamma13*der(Gamma1, 3, idx) + gamma23*der(Gamma2, 3, idx) + gamma33*der(Gamma3, 3, idx)
          + gamma13*der(Gamma1, 3, idx) + gamma23*der(Gamma2, 3, idx) + gamma33*der(Gamma3, 3, idx)
        )
        - 0.5*(
          d1g31*d3gi11 + d2g32*d3gi22 + d3g33*d3gi33 + 2.0*(d1g32*d3gi12 + d1g33*d3gi13 + d2g33*d3gi23)
          + d1g31*d3gi11 + d2g32*d3gi22 + d3g33*d3gi33 + 2.0*(d1g32*d3gi12 + d1g33*d3gi13 + d2g33*d3gi23)
          - Gamma1*d1g33 - Gamma2*d2g33 - Gamma3*d3g33
        )
        - (
          G131*G131 + G232*G232 + G333*G333 + 2.0*(G132*G231 + G133*G331 + G233*G332)
        )
      );


    real_t expression = (
      gammai11*(D1D1phi - 2.0*d1phi*d1phi) + gammai22*(D2D2phi - 2.0*d2phi*d2phi) + gammai33*(D3D3phi - 2.0*d3phi*d3phi)
      + 2.0*(gammai12*(D1D2phi - 2.0*d1phi*d2phi) + gammai13*(D1D3phi - 2.0*d1phi*d3phi) + gammai23*(D2D3phi - 2.0*d2phi*d3phi))
    );

    /* phi-piece */
    ricciTF11 += -2.0*( D1D1phi -2.0*d1phi*d1phi + gamma11*(expression) );
    ricciTF12 += -2.0*( D1D2phi -2.0*d1phi*d2phi + gamma12*(expression) );
    ricciTF13 += -2.0*( D1D3phi -2.0*d1phi*d3phi + gamma13*(expression) );
    ricciTF22 += -2.0*( D2D2phi -2.0*d2phi*d2phi + gamma22*(expression) );
    ricciTF23 += -2.0*( D2D3phi -2.0*d2phi*d3phi + gamma23*(expression) );
    ricciTF33 += -2.0*( D3D3phi -2.0*d3phi*d3phi + gamma33*(expression) );

    /* remove trace... */ 
    trace = gammai11*ricciTF11 + gammai22*ricciTF22 + gammai33*ricciTF33
        + 2.0*(gammai12*ricciTF12 + gammai13*ricciTF13 + gammai23*ricciTF23);
    
    ricciTF11 -= (1/3.0)*gamma11*trace;
    ricciTF12 -= (1/3.0)*gamma12*trace;
    ricciTF13 -= (1/3.0)*gamma13*trace;
    ricciTF22 -= (1/3.0)*gamma22*trace;
    ricciTF23 -= (1/3.0)*gamma23*trace;
    ricciTF33 -= (1/3.0)*gamma33*trace;

  }

  void calculateDDphi(idx_t idx)
  {
    // normal derivatives of alpha
    d1phi = der(phi, 1, idx);
    d2phi = der(phi, 2, idx);
    d3phi = der(phi, 3, idx);

    // double covariant derivatives, using normal metric
    D1D1phi = dder(phi, 1, 1, idx) - (G111*d1phi + G211*d2phi + G311*d3phi);
    D1D2phi = dder(phi, 1, 2, idx) - (G112*d1phi + G212*d2phi + G312*d3phi);
    D1D3phi = dder(phi, 1, 3, idx) - (G113*d1phi + G213*d2phi + G313*d3phi);
    D2D2phi = dder(phi, 2, 2, idx) - (G122*d1phi + G222*d2phi + G322*d3phi);
    D2D3phi = dder(phi, 2, 3, idx) - (G123*d1phi + G223*d2phi + G323*d3phi);
    D3D3phi = dder(phi, 3, 3, idx) - (G133*d1phi + G233*d2phi + G333*d3phi);
  }

  void calculateDDalpha(idx_t idx)
  {
    // normal derivatives of alpha
    d1a = der(alpha, 1, idx);
    d2a = der(alpha, 2, idx);
    d3a = der(alpha, 3, idx);

    // double covariant derivatives - use non-unitary metric - extra pieces that depend on phi!
    // these are needed for the BSSN_CALCULATE_DIDJALPHA macro
    real_t gamma1ldlphi = gammai11*d1phi + gammai12*d2phi + gammai13*d3phi;
    real_t gamma2ldlphi = gammai21*d1phi + gammai22*d2phi + gammai23*d3phi;
    real_t gamma3ldlphi = gammai31*d1phi + gammai32*d2phi + gammai33*d3phi;
    BSSN_APPLY_TO_IJ_PERMS(BSSN_CALCULATE_DIDJALPHA);

  }


};

#endif
