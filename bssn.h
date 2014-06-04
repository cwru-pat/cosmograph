#ifndef COSMO_BSSN
#define COSMO_BSSN

/** BSSN class **/
class BSSN
{
  /* arrays for storing fields */

  // unitary metric, gamma_ij
  real_t gammaxx;
  real_t gammaxy;
  real_t gammayy;
  real_t gammayz;
  real_t gammazz;
  real_t gammaxz;

  // conformal factor
  real_t W;

  // A
  real_t Axx;
  real_t Axy;
  real_t Ayy;
  real_t Ayz;
  real_t Azz;
  real_t Axz;

  // ext. curv.
  real_t K;

  // christoffel derivative
  real_t Gammax;
  real_t Gammay;
  real_t Gammaz;

public:
  person() : age(5) { }
  void print() const;
};


public:
  BSSN() {};

  

};

#endif
