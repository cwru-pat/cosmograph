// g++ -g --std=c++11 -lrt scalar.cc ../bssn/bssn.cc ../scalar/scalar.cc ../utils/Timer.cc ../utils/ConfigParser.cc

#include "../cosmo_includes.h"
#include "../cosmo_types.h"
#include "../cosmo_globals.h"

#include "../scalar/scalar.h"
#include "../bssn/bssn.h"

using namespace std;
using namespace cosmo;

/* global definitions */
TimerManager _timer;
ConfigParser _config;
#ifndef dt
  real_t dt;
#endif
#ifndef dx
  real_t dx;
#endif

int main()
{
  idx_t i=0, j=0, k=0;
  dx = 1.0;
  dt = 0.1;

  std::cout << "Creating BSSN and Scalar classes..." << std::endl;

  // initialize GR sim; Scalar needs to reference this
  // but - don't evolve, remaining Minkowski/flat
  BSSN bssnSim;
  BSSNData bd;
  bssnSim.set_bd_values(0,0,0,&bd);

  // make a scalar field class
  Scalar scalarSim;

  // initialize with \phi = x + y + z
  std::cout << "Initializing scalar field..." << std::endl;
  #pragma omp parallel for default(shared) private(i, j, k, b_data)
  LOOP3(i,j,k)
  {
    scalarSim.phi._array_p[NP_INDEX(i,j,k)] = i + j + k;
    scalarSim.psi1._array_p[NP_INDEX(i,j,k)] = 1;
    scalarSim.psi2._array_p[NP_INDEX(i,j,k)] = 1;
    scalarSim.psi3._array_p[NP_INDEX(i,j,k)] = 1;
    scalarSim.Pi._array_p[NP_INDEX(i,j,k)] = 1;
  }
  scalarSim.stepInit(); // copy over to _a register 

  std::cout << "Scalar field constraint residual: "
    << scalarSim.scalarConstraint(NX/2, NY/2, NZ/2, 1) << std::endl;

  bssnSim.clearSrc();
  scalarSim.addBSSNSource(&bssnSim);

  std::cout << "BSSN 'DIFFr_a' Source: " << (*bssnSim.fields["DIFFr_a"])[NP_INDEX(NX/2, NY/2, NZ/2)] << std::endl;
  std::cout << "BSSN 'DIFFS_a' Source: " << (*bssnSim.fields["DIFFS_a"])[NP_INDEX(NX/2, NY/2, NZ/2)] << std::endl;
  std::cout << "BSSN 'S1_a' Source: " << (*bssnSim.fields["S1_a"])[NP_INDEX(NX/2, NY/2, NZ/2)] << std::endl;
  std::cout << "BSSN 'S2_a' Source: " << (*bssnSim.fields["S2_a"])[NP_INDEX(NX/2, NY/2, NZ/2)] << std::endl;
  std::cout << "BSSN 'S3_a' Source: " << (*bssnSim.fields["S3_a"])[NP_INDEX(NX/2, NY/2, NZ/2)] << std::endl;
  std::cout << "BSSN 'STF11_a' Source: " << (*bssnSim.fields["STF11_a"])[NP_INDEX(NX/2, NY/2, NZ/2)] << std::endl;
  std::cout << "BSSN 'STF12_a' Source: " << (*bssnSim.fields["STF12_a"])[NP_INDEX(NX/2, NY/2, NZ/2)] << std::endl;
  std::cout << "BSSN 'STF13_a' Source: " << (*bssnSim.fields["STF13_a"])[NP_INDEX(NX/2, NY/2, NZ/2)] << std::endl;
  std::cout << "BSSN 'STF22_a' Source: " << (*bssnSim.fields["STF22_a"])[NP_INDEX(NX/2, NY/2, NZ/2)] << std::endl;
  std::cout << "BSSN 'STF23_a' Source: " << (*bssnSim.fields["STF23_a"])[NP_INDEX(NX/2, NY/2, NZ/2)] << std::endl;
  std::cout << "BSSN 'STF33_a' Source: " << (*bssnSim.fields["STF33_a"])[NP_INDEX(NX/2, NY/2, NZ/2)] << std::endl;


  std::cout << "Running RK step (evolving scalar, not evolving metric)..." << std::endl;
  // run a test simulation step (dynamic field only)
  scalarSim.stepInit();

  // First RK step
  scalarSim.RKEvolve(&bd);
  scalarSim.K1Finalize();

  // Second RK step
  scalarSim.RKEvolve(&bd);
  scalarSim.K2Finalize();

  // Third RK step
  scalarSim.RKEvolve(&bd);
  scalarSim.K3Finalize();

  // Fourth RK step
  scalarSim.RKEvolve(&bd);
  scalarSim.K4Finalize();

  std::cout << "Done!" << std::endl;
  exit(EXIT_SUCCESS);
}
