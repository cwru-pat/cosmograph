#ifndef COSMO_BSSN
#define COSMO_BSSN

/** BSSN class **/
class BSSN
{
  /* create arrays for storing fields */
  BSSN_APPLY_TO_FIELDS(RK4_ARRAY_CREATE);

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

};

#endif
