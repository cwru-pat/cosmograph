#ifndef COSMO_BSSN
#define COSMO_BSSN

/** BSSN class **/
class BSSN
{
  /* create arrays for storing fields */
    BSSN_APPLY_TO_FIELDS(RK4_ARRAY_CREATE);

  /* Local variables: */

      // generic var for tracing
      real_t trace;

      // trace-free ricci tensor components
      real_t ricciTF11;
      real_t ricciTF12;
      real_t ricciTF13;
      real_t ricciTF22;
      real_t ricciTF23;
      real_t ricciTF33;

      // derivatives of \alpha
        // covariant double-derivatives 
        real_t D1D1a, D1D2a, D1D3a, D2D2a, D2D3a, D3D3a;
        // normal derivatives of
        real_t d1a, d2a, d3a;

      // derivatives of the metric, d_i g_jk
      real_t d1g11, d1g12, d1g13, d1g22, d1g23, d1g33,
             d2g11, d2g12, d2g13, d2g22, d2g23, d2g33,
             d3g11, d3g12, d3g13, d3g22, d3g23, d3g33;

      // derivatives of the inverse metric d_i g^jk
      real_t d1gi11, d1gi12, d1gi13, d1gi22, d1gi23, d1gi33,
             d2gi11, d2gi12, d2gi13, d2gi22, d2gi23, d2gi33,
             d3gi11, d3gi12, d3gi13, d3gi22, d3gi23, d3gi33;

    /* End Local variables */

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

  /* Calculate inverse metric components */
  inline void calculate_gammai(real_t idx)
  {
    gammai11[idx] = gamma22[idx]*gamma33[idx] - gamma23[idx]*gamma23[idx];
    gammai12[idx] = gamma13[idx]*gamma23[idx] - gamma12[idx]*gamma33[idx];
    gammai13[idx] = gamma12[idx]*gamma23[idx] - gamma13[idx]*gamma22[idx];
    gammai22[idx] = gamma11[idx]*gamma33[idx] - gamma13[idx]*gamma13[idx];
    gammai23[idx] = gamma12[idx]*gamma13[idx] - gamma23[idx]*gamma11[idx];
    gammai33[idx] = gamma11[idx]*gamma22[idx] - gamma12[idx]*gamma12[idx];
  }

  /* Calculate metric derivatives */
  inline void calculate_dgamma(real_t idx)
  {
    d1g11 = der(d1g11, 1, idx);
    d1g12 = der(d1g12, 1, idx);
    d1g13 = der(d1g13, 1, idx);
    d1g22 = der(d1g22, 1, idx);
    d1g23 = der(d1g23, 1, idx);
    d1g33 = der(d1g33, 1, idx);
    d2g11 = der(d2g11, 2, idx);
    d2g12 = der(d2g12, 2, idx);
    d2g13 = der(d2g13, 2, idx);
    d2g22 = der(d2g22, 2, idx);
    d2g23 = der(d2g23, 2, idx);
    d2g33 = der(d2g33, 2, idx);
    d3g11 = der(d3g11, 3, idx);
    d3g12 = der(d3g12, 3, idx);
    d3g13 = der(d3g13, 3, idx);
    d3g22 = der(d3g22, 3, idx);
    d3g23 = der(d3g23, 3, idx);
    d3g33 = der(d3g33, 3, idx);
  }

  /* Calculate metric derivatives */
  inline void calculate_dgammai(real_t idx)
  {
    d1gi11 = der(d1gi11, 1, idx);
    d1gi12 = der(d1gi12, 1, idx);
    d1gi13 = der(d1gi13, 1, idx);
    d1gi22 = der(d1gi22, 1, idx);
    d1gi23 = der(d1gi23, 1, idx);
    d1gi33 = der(d1gi33, 1, idx);
    d2gi11 = der(d2gi11, 2, idx);
    d2gi12 = der(d2gi12, 2, idx);
    d2gi13 = der(d2gi13, 2, idx);
    d2gi22 = der(d2gi22, 2, idx);
    d2gi23 = der(d2gi23, 2, idx);
    d2gi33 = der(d2gi33, 2, idx);
    d3gi11 = der(d3gi11, 3, idx);
    d3gi12 = der(d3gi12, 3, idx);
    d3gi13 = der(d3gi13, 3, idx);
    d3gi22 = der(d3gi22, 3, idx);
    d3gi23 = der(d3gi23, 3, idx);
    d3gi33 = der(d3gi33, 3, idx);
  }

  inline void calculate_christoffels(real_t idx)
  {
    // christoffel symbols: \Gamma^i_{jk} = Gijk  
    // Perhaps optimize this some after making sure it works
    G111 = 0.5*(
        gammai11 * (d1g11 + d1g11 - d1g11) +
        gammai12 * (d1g12 + d1g12 - d2g11) +
        gammai13 * (d1g13 + d1g13 - d3g11)
      );
    G112 = 0.5*(
        gammai11 * (d2g11 + d1g21 - d1g12) +
        gammai12 * (d2g12 + d1g22 - d2g12) +
        gammai13 * (d2g13 + d1g23 - d3g12)
      );
    G113 = 0.5*(
        gammai11 * (d3g11 + d1g31 - d1g13) +
        gammai12 * (d3g12 + d1g32 - d2g13) +
        gammai13 * (d3g13 + d1g33 - d3g13)
      );
    G122 = 0.5*(
        gammai11 * (d2g21 + d2g21 - d1g22) +
        gammai12 * (d2g22 + d2g22 - d2g22) +
        gammai13 * (d2g23 + d2g23 - d3g22)
      );
    G123 = 0.5*(
        gammai11 * (d3g21 + d2g31 - d1g23) +
        gammai12 * (d3g22 + d2g32 - d2g23) +
        gammai13 * (d3g23 + d2g33 - d3g23)
      );
    G133 = 0.5*(
        gammai11 * (d3g31 + d3g31 - d1g33) +
        gammai12 * (d3g32 + d3g32 - d2g33) +
        gammai13 * (d3g33 + d3g33 - d3g33)
      );
    G211 = 0.5*(
        gammai21 * (d1g11 + d1g11 - d1g11) +
        gammai22 * (d1g12 + d1g12 - d2g11) +
        gammai23 * (d1g13 + d1g13 - d3g11)
      );
    G212 = 0.5*(
        gammai21 * (d2g11 + d1g21 - d1g12) +
        gammai22 * (d2g12 + d1g22 - d2g12) +
        gammai23 * (d2g13 + d1g23 - d3g12)
      );
    G213 = 0.5*(
        gammai21 * (d3g11 + d1g31 - d1g13) +
        gammai22 * (d3g12 + d1g32 - d2g13) +
        gammai23 * (d3g13 + d1g33 - d3g13)
      );
    G222 = 0.5*(
        gammai21 * (d2g21 + d2g21 - d1g22) +
        gammai22 * (d2g22 + d2g22 - d2g22) +
        gammai23 * (d2g23 + d2g23 - d3g22)
      );
    G223 = 0.5*(
        gammai21 * (d3g21 + d2g31 - d1g23) +
        gammai22 * (d3g22 + d2g32 - d2g23) +
        gammai23 * (d3g23 + d2g33 - d3g23)
      );
    G233 = 0.5*(
        gammai21 * (d3g31 + d3g31 - d1g33) +
        gammai22 * (d3g32 + d3g32 - d2g33) +
        gammai23 * (d3g33 + d3g33 - d3g33)
      );
    G311 = 0.5*(
        gammai31 * (d1g11 + d1g11 - d1g11) +
        gammai32 * (d1g12 + d1g12 - d2g11) +
        gammai33 * (d1g13 + d1g13 - d3g11)
      );
    G312 = 0.5*(
        gammai31 * (d2g11 + d1g21 - d1g12) +
        gammai32 * (d2g12 + d1g22 - d2g12) +
        gammai33 * (d2g13 + d1g23 - d3g12)
      );
    G313 = 0.5*(
        gammai31 * (d3g11 + d1g31 - d1g13) +
        gammai32 * (d3g12 + d1g32 - d2g13) +
        gammai33 * (d3g13 + d1g33 - d3g13)
      );
    G322 = 0.5*(
        gammai31 * (d2g21 + d2g21 - d1g22) +
        gammai32 * (d2g22 + d2g22 - d2g22) +
        gammai33 * (d2g23 + d2g23 - d3g22)
      );
    G323 = 0.5*(
        gammai31 * (d3g21 + d2g31 - d1g23) +
        gammai32 * (d3g22 + d2g32 - d2g23) +
        gammai33 * (d3g23 + d2g33 - d3g23)
      );
    G333 = 0.5*(
        gammai31 * (d3g31 + d3g31 - d1g33) +
        gammai32 * (d3g32 + d3g32 - d2g33) +
        gammai33 * (d3g33 + d3g33 - d3g33)
      );
  }

  /* Calculate trace-free ricci tensor components */
  void calculateRicciTF(real_t idx)
  {

  }

  void calculateDDalpha(real_t idx)
  {

  }

  /*
   * Time-evolution functions for all the variables
   */

  // gamma (unit determinant metric)
  real_t dt_gamma11(real_t idx) { return BSSN_GAMMA_EQUATION(1, 1, idx); }
  real_t dt_gamma12(real_t idx) { return BSSN_GAMMA_EQUATION(1, 2, idx); }
  real_t dt_gamma13(real_t idx) { return BSSN_GAMMA_EQUATION(1, 3, idx); }
  real_t dt_gamma22(real_t idx) { return BSSN_GAMMA_EQUATION(2, 2, idx); }
  real_t dt_gamma23(real_t idx) { return BSSN_GAMMA_EQUATION(2, 3, idx); }
  real_t dt_gamma33(real_t idx) { return BSSN_GAMMA_EQUATION(3, 3, idx); }

  real_t dt_A11(real_t idx) { return BSSN_A_EQUATION(1, 1, idx); }
  real_t dt_A12(real_t idx) { return BSSN_A_EQUATION(1, 2, idx); }
  real_t dt_A13(real_t idx) { return BSSN_A_EQUATION(1, 3, idx); }
  real_t dt_A22(real_t idx) { return BSSN_A_EQUATION(2, 2, idx); }
  real_t dt_A23(real_t idx) { return BSSN_A_EQUATION(2, 3, idx); }
  real_t dt_A33(real_t idx) { return BSSN_A_EQUATION(3, 3, idx); }




};

#endif
