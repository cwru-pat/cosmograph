
#include "cosmo.h"
#include "globals.h"
#include <cmath>
#include <cfloat>

using namespace std;
using namespace cosmo;

/* global definitions */
TimerManager _timer;
ConfigParser _config;

int main(int argc, char **argv)
{
  _timer["MAIN"].start();
  idx_t i, j, k, s, steps;

  // read in config file
  if(argc != 2)
  {
    cout << "Error: please supply exactly one config filename as an argument.\n";
    return EXIT_FAILURE;
  }
  else
  {
    _config.parse(argv[1]);
  }

  // IO init - will use this for logging.
  IOData iodata;
  io_init(&iodata, _config["output_dir"]);
  // save a copy of config.txt
  io_config_backup(&iodata, argv[1]);

  steps = stoi(_config["steps"]);
  omp_set_num_threads(stoi(_config["omp_num_threads"]));

  // Create simulation
  _timer["init"].start();
    LOG(iodata.log, "Creating initial conditions...\n");

    // Fluid fields
    // Static matter (w=0)
    Static staticSim;
    staticSim.init();
    // DE
    Lambda lambdaSim;

    // GR Fields
    BSSN bssnSim;
    bssnSim.init();

    // generic reusable fourier class for N^3 arrays
    Fourier fourier;
    fourier.Initialize(N, staticSim.fields["D_a"] /* just any N^3 array for planning */);

    // "conformal" initial conditions:
    LOG(iodata.log, "Using conformal initial conditions...\n");
    set_conformal_ICs(bssnSim.fields, staticSim.fields, &fourier, &iodata);

    // Trial FRW class
    staticSim.addBSSNSrc(bssnSim.fields);
    real_t frw_rho = average(bssnSim.fields["r_a"]);
    FRW<real_t> frw (0.0, -sqrt(24.0*PI*frw_rho) /* K */);
    frw.addFluid(frw_rho /* rho */, 0.0 /* 'w' */);

  _timer["init"].stop();

  // evolve simulation
  LOG(iodata.log, "Running simulation...\n");
  _timer["loop"].start();
  for(s=0; s < steps; ++s) {

    _timer["Reference FRW"].start();
    // LOG(iodata.log,    "\n"
    //                 << frw.get_phi()
    //                 << "\n"
    //                 << bssnSim.fields["phi_a"][0]
    //                 << "\n"
    //               );
    int subiters = 100;
    for(int p=0; p<subiters; ++p)
    {
      frw.step(dt/((real_t) subiters));
    }
    _timer["Reference FRW"].stop();

    // Init arrays and calculate source term for next step
      // _p is copied to _a here, which hydro uses
      bssnSim.stepInit();
      // clear existing data
      bssnSim.clearSrc();
      // add hydro source to bssn sim
      staticSim.addBSSNSrc(bssnSim.fields);
      lambdaSim.addBSSNSrc(bssnSim.fields);


#pragma omp parallel for default(shared) private(i, j, k)
LOOP3(i,j,k)
{
  BSSNData b_paq = {0}; // data structure associated with bssn sim
  bssnSim.set_paq_values(i,j,k,&b_paq);
}


    // output simulation information
    // these generally output any data in the _a registers.
    _timer["output"].start();
      io_show_progress(s, steps);
      io_data_dump(bssnSim.fields, staticSim.fields, &iodata, s, &fourier);
    _timer["output"].stop();

    // Run RK steps explicitly here (ties together BSSN + Hydro stuff).
    // See bssn class or hydro class for more comments.
    _timer["RK_steps"].start();


// DEBUG
// if(0) { //N==64 || (N==128 && (s+1)%2) || (N==256 && (s+1)%4)) {

// BSSNData tmp_paq = {0};
// bssnSim.clearSrc();
// staticSim.addBSSNSrc(bssnSim.fields);
// lambdaSim.addBSSNSrc(bssnSim.fields);
// bssnSim.set_paq_values(0, 0, 0, &tmp_paq);

// real_t term1 = 0.5*(
//   tmp_paq.gammai11*tmp_paq.d1d1g11 + tmp_paq.gammai22*tmp_paq.d2d2g11 + tmp_paq.gammai33*tmp_paq.d3d3g11
//   + 2.0*(tmp_paq.gammai12*tmp_paq.d1d2g11 + tmp_paq.gammai13*tmp_paq.d1d3g11 + tmp_paq.gammai23*tmp_paq.d2d3g11)
// );
// real_t term2 = 0.5*(
//   tmp_paq.gamma11*derivative(tmp_paq.i, tmp_paq.j, tmp_paq.k, 1, bssnSim.fields["Gamma1_a"])
//     + tmp_paq.gamma21*derivative(tmp_paq.i, tmp_paq.j, tmp_paq.k, 1, bssnSim.fields["Gamma2_a"])
//     + tmp_paq.gamma31*derivative(tmp_paq.i, tmp_paq.j, tmp_paq.k, 1, bssnSim.fields["Gamma3_a"]) +
//   tmp_paq.gamma11*derivative(tmp_paq.i, tmp_paq.j, tmp_paq.k, 1, bssnSim.fields["Gamma1_a"])
//     + tmp_paq.gamma21*derivative(tmp_paq.i, tmp_paq.j, tmp_paq.k, 1, bssnSim.fields["Gamma2_a"])
//     + tmp_paq.gamma31*derivative(tmp_paq.i, tmp_paq.j, tmp_paq.k, 1, bssnSim.fields["Gamma3_a"])
// );
// real_t term3 = 0.5*(
//   tmp_paq.d1g11*tmp_paq.d1gi11 + tmp_paq.d1g12*tmp_paq.d2gi11 + tmp_paq.d1g13*tmp_paq.d3gi11
//     + tmp_paq.d2g11*tmp_paq.d1gi12 + tmp_paq.d2g12*tmp_paq.d2gi12 + tmp_paq.d2g13*tmp_paq.d3gi12
//     + tmp_paq.d3g11*tmp_paq.d1gi13 + tmp_paq.d3g12*tmp_paq.d2gi13 + tmp_paq.d3g13*tmp_paq.d3gi13
//   + tmp_paq.d1g11*tmp_paq.d1gi11 + tmp_paq.d1g12*tmp_paq.d2gi11 + tmp_paq.d1g13*tmp_paq.d3gi11
//     + tmp_paq.d2g11*tmp_paq.d1gi12 + tmp_paq.d2g12*tmp_paq.d2gi12 + tmp_paq.d2g13*tmp_paq.d3gi12
//     + tmp_paq.d3g11*tmp_paq.d1gi13 + tmp_paq.d3g12*tmp_paq.d2gi13 + tmp_paq.d3g13*tmp_paq.d3gi13
//   - tmp_paq.Gamma1*tmp_paq.d1g11 - tmp_paq.Gamma2*tmp_paq.d2g11 - tmp_paq.Gamma3*tmp_paq.d3g11
// );
// real_t term4 = 0.5*(
//     tmp_paq.G111*tmp_paq.G111 + tmp_paq.G112*tmp_paq.G211 + tmp_paq.G113*tmp_paq.G311
//   + tmp_paq.G211*tmp_paq.G112 + tmp_paq.G212*tmp_paq.G212 + tmp_paq.G213*tmp_paq.G312
//   + tmp_paq.G311*tmp_paq.G113 + tmp_paq.G312*tmp_paq.G213 + tmp_paq.G313*tmp_paq.G313
// );
// real_t expr = (
//     tmp_paq.gammai11*(tmp_paq.D1D1phi + 2.0*tmp_paq.d1phi*tmp_paq.d1phi)
//     + tmp_paq.gammai22*(tmp_paq.D2D2phi + 2.0*tmp_paq.d2phi*tmp_paq.d2phi)
//     + tmp_paq.gammai33*(tmp_paq.D3D3phi + 2.0*tmp_paq.d3phi*tmp_paq.d3phi)
//     + 2.0*(
//       tmp_paq.gammai12*(tmp_paq.D1D2phi + 2.0*tmp_paq.d1phi*tmp_paq.d2phi)
//       + tmp_paq.gammai13*(tmp_paq.D1D3phi + 2.0*tmp_paq.d1phi*tmp_paq.d3phi)
//       + tmp_paq.gammai23*(tmp_paq.D2D3phi + 2.0*tmp_paq.d2phi*tmp_paq.d3phi)
//     )
// );

// real_t sumterm1 = (
//       derivative(tmp_paq.i, tmp_paq.j, tmp_paq.k, 1, bssnSim.fields["Gamma1_a"])
//       + derivative(tmp_paq.i, tmp_paq.j, tmp_paq.k, 2, bssnSim.fields["Gamma2_a"])
//       + derivative(tmp_paq.i, tmp_paq.j, tmp_paq.k, 3, bssnSim.fields["Gamma3_a"])
//     );
// real_t sumterm2 = (
//   tmp_paq.gammai11*tmp_paq.d1d1g11 + tmp_paq.gammai22*tmp_paq.d2d2g11 + tmp_paq.gammai33*tmp_paq.d3d3g11
//       + tmp_paq.gammai11*tmp_paq.d1d1g22 + tmp_paq.gammai22*tmp_paq.d2d2g22 + tmp_paq.gammai33*tmp_paq.d3d3g22
//       + tmp_paq.gammai11*tmp_paq.d1d1g33 + tmp_paq.gammai22*tmp_paq.d2d2g33 + tmp_paq.gammai33*tmp_paq.d3d3g33
//   );

// LOG(iodata.log, "\n Some values: { ");
// // LOG(iodata.log, " " << tmp_paq.K);
// // LOG(iodata.log, ", " << tmp_paq.K * tmp_paq.K);
// // LOG(iodata.log, ", " << tmp_paq.phi);
// // LOG(iodata.log, ", " << tmp_paq.d1phi);
// // LOG(iodata.log, ", " << tmp_paq.D1D1phi);
// // LOG(iodata.log, ", " << tmp_paq.rho);
// // LOG(iodata.log, ", " << tmp_paq.AijAij);
// // LOG(iodata.log, ", " << tmp_paq.A11);
// // LOG(iodata.log, ", " << tmp_paq.gamma11);
// // LOG(iodata.log, ", " << tmp_paq.gammai11);
// // LOG(iodata.log, ", " << tmp_paq.G111);
// // LOG(iodata.log, ", " << tmp_paq.GL111);
// // LOG(iodata.log, ", " << tmp_paq.d1g11);
// // LOG(iodata.log, ", " << tmp_paq.d1d1g11);

// LOG(iodata.log, ", " << tmp_paq.A11);
// LOG(iodata.log, ", " << tmp_paq.A12);
// LOG(iodata.log, ", " << tmp_paq.A13);

// LOG(iodata.log, ", " << derivative(tmp_paq.i, tmp_paq.j, tmp_paq.k, 1, bssnSim.fields["K_a"]));
// LOG(iodata.log, ", " << derivative(tmp_paq.i, tmp_paq.j, tmp_paq.k, 2, bssnSim.fields["K_a"]));
// LOG(iodata.log, ", " << derivative(tmp_paq.i, tmp_paq.j, tmp_paq.k, 3, bssnSim.fields["K_a"]));

// LOG(iodata.log, ", " << tmp_paq.d1phi);
// LOG(iodata.log, ", " << tmp_paq.d2phi);
// LOG(iodata.log, ", " << tmp_paq.d3phi);

// LOG(iodata.log, ", " << tmp_paq.D1D1phi);
// LOG(iodata.log, ", " << tmp_paq.D2D2phi);
// LOG(iodata.log, ", " << tmp_paq.D3D3phi);

// LOG(iodata.log, ", " << tmp_paq.Gamma1);
// LOG(iodata.log, ", " << tmp_paq.Gamma2);
// LOG(iodata.log, ", " << tmp_paq.Gamma3);

// LOG(iodata.log, ", " << tmp_paq.Uricci11);
// LOG(iodata.log, ", " << tmp_paq.ricciTF11);
// LOG(iodata.log, ", " << tmp_paq.gammai11);
// LOG(iodata.log, ", " << tmp_paq.Uricci22);
// LOG(iodata.log, ", " << tmp_paq.ricciTF22);
// LOG(iodata.log, ", " << tmp_paq.gammai22);
// LOG(iodata.log, ", " << tmp_paq.Uricci33);
// LOG(iodata.log, ", " << tmp_paq.ricciTF33);
// LOG(iodata.log, ", " << tmp_paq.gammai33);
// LOG(iodata.log, ", " << tmp_paq.Uricci12);
// LOG(iodata.log, ", " << tmp_paq.ricciTF12);
// LOG(iodata.log, ", " << tmp_paq.gammai12);
// LOG(iodata.log, ", " << tmp_paq.Uricci23);
// LOG(iodata.log, ", " << tmp_paq.ricciTF23);
// LOG(iodata.log, ", " << tmp_paq.gammai23);
// LOG(iodata.log, ", " << tmp_paq.Uricci13);
// LOG(iodata.log, ", " << tmp_paq.ricciTF13);
// LOG(iodata.log, ", " << tmp_paq.gammai13);

// LOG(iodata.log, ", " << tmp_paq.ricci);
// LOG(iodata.log, ", " << tmp_paq.unitRicci);

// LOG(iodata.log, ", " << 
//     tmp_paq.ricciTF11*tmp_paq.gammai11 + tmp_paq.ricciTF22*tmp_paq.gammai22 + tmp_paq.ricciTF33*tmp_paq.gammai33
//     + 2.0*(tmp_paq.ricciTF12*tmp_paq.gammai12 + tmp_paq.ricciTF13*tmp_paq.gammai13 + tmp_paq.ricciTF23*tmp_paq.gammai23)
//   );

// LOG(iodata.log, ", " << tmp_paq.trace);
// LOG(iodata.log, ", " << term1);
// LOG(iodata.log, ", " << term2);
// LOG(iodata.log, ", " << term3);
// LOG(iodata.log, ", " << term4);
// LOG(iodata.log, ", " << expr);
// LOG(iodata.log, ", " << tmp_paq.D1D1phi - 2.0*tmp_paq.d1phi*tmp_paq.d1phi);
// LOG(iodata.log, ", " << sumterm1);
// LOG(iodata.log, ", " << sumterm2);
// LOG(iodata.log, "}, \n");

// }

// DEBUG
// if( (s>100 && s<150) || (s>250 && s<300) || (s>400 && s<450) || (s>550 && s<600) ) {
//   LOOP3(i,j,k) {
//     BSSNData b_paq = {0}; // data structure associated with bssn sim
//     bssnSim.set_paq_values(i,j,k,&b_paq);
//     staticSim.fields["D_a"][b_paq.idx] = exp(6.0*b_paq.phi)/16.0/PI*(
//       b_paq.ricci + 2.0/3.0*pw2(b_paq.K) - b_paq.AijAij
//     );
//   }
// }

    #pragma omp parallel for default(shared) private(i, j, k)
    LOOP3(i, j, k) {
      bssnSim.set_detgamma(i,j,k);
      bssnSim.set_gammai_values(i, j, k);
    }

    // First RK step, Set Hydro Vars, & calc. constraint
    #pragma omp parallel for default(shared) private(i, j, k)
    LOOP3(i, j, k) {
      BSSNData b_paq = {0}; // data structure associated with bssn sim
      bssnSim.set_paq_values(i,j,k,&b_paq);
    }
    #pragma omp parallel for default(shared) private(i, j, k)
    LOOP3(i, j, k)
    {
      BSSNData b_paq = {0}; // data structure associated with bssn sim
      bssnSim.K1CalcPt(i, j, k, &b_paq);
    }

    // reset source using new metric
    bssnSim.clearSrc();
    // add hydro source to bssn sim
    staticSim.addBSSNSrc(bssnSim.fields);
    bssnSim.regSwap_c_a();


    // Subsequent BSSN steps
      // Second RK step
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k) {
        bssnSim.set_detgamma(i,j,k);
        bssnSim.set_gammai_values(i, j, k);
      }
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k) {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.set_paq_values(i,j,k,&b_paq);
      }
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k)
      {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.K2CalcPt(i, j, k, &b_paq);
      }

      // reset source using new metric
      bssnSim.clearSrc();
      // add hydro source to bssn sim
      staticSim.addBSSNSrc(bssnSim.fields);
      lambdaSim.addBSSNSrc(bssnSim.fields);

      bssnSim.regSwap_c_a();


      // Third RK step
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k) {
        bssnSim.set_detgamma(i,j,k);
        bssnSim.set_gammai_values(i, j, k);
      }
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k) {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.set_paq_values(i,j,k,&b_paq);
      }
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k)
      {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.K3CalcPt(i, j, k, &b_paq);
      }

      // reset source using new metric
      bssnSim.clearSrc();
      // add hydro source to bssn sim
      staticSim.addBSSNSrc(bssnSim.fields);
      lambdaSim.addBSSNSrc(bssnSim.fields); 

      bssnSim.regSwap_c_a();


      // Fourth RK step
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k) {
        bssnSim.set_detgamma(i,j,k);
        bssnSim.set_gammai_values(i, j, k);
      }
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k) {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.set_paq_values(i,j,k,&b_paq);
      }
      #pragma omp parallel for default(shared) private(i, j, k)
      LOOP3(i, j, k)
      {
        BSSNData b_paq = {0}; // data structure associated with bssn sim
        bssnSim.K4CalcPt(i, j, k, &b_paq);
      }

    // Wrap up
      // bssn _f <-> _p
      bssnSim.stepTerm();
      // "current" data is in the _p array.

    _timer["RK_steps"].stop();


    _timer["meta_output_interval"].start();

      if(s%iodata.meta_output_interval == 0)
      {
        idx_t isNaN = 0;
        bssnSim.stepInit();
        bssnSim.clearSrc();
        staticSim.addBSSNSrc(bssnSim.fields);

        real_t mean_hamiltonian_constraint = bssnSim.hamiltonianConstraintMagMean();
        real_t stdev_hamiltonian_constraint = bssnSim.hamiltonianConstraintMagStDev(mean_hamiltonian_constraint);
        real_t max_hamiltonian_constraint = bssnSim.hamiltonianConstraintMagMax();
        real_t mean_momentum_constraint = bssnSim.momentumConstraintMagMean();
        real_t stdev_momentum_constraint = bssnSim.momentumConstraintMagStDev(mean_momentum_constraint);
        real_t max_momentum_constraint = bssnSim.momentumConstraintMagMax();
        io_dump_data(mean_hamiltonian_constraint, &iodata, "H_violations");
        io_dump_data(stdev_hamiltonian_constraint, &iodata, "H_violations");
        io_dump_data(max_hamiltonian_constraint, &iodata, "H_violations");
        io_dump_data(mean_momentum_constraint, &iodata, "M_violations");
        io_dump_data(stdev_momentum_constraint, &iodata, "M_violations");
        io_dump_data(max_momentum_constraint, &iodata, "M_violations");

        // LOG(iodata.log,
        //   "\nFractional Mean / St. Dev Momentum constraint violation is: " <<
        //   mean_momentum_constraint << " / " << stdev_momentum_constraint <<
        //   "\nFractional Mean / St. Dev Hamiltonian constraint violation is: " <<
        //   mean_hamiltonian_constraint << " / " << stdev_hamiltonian_constraint <<
        //   "\n"
        // );

        if(numNaNs(bssnSim.fields["phi_a"]) > 0)
        {
          LOG(iodata.log, "\nNAN detected!\n");
          _timer["meta_output_interval"].stop();
          break;
        }
      }
    _timer["meta_output_interval"].stop();

  }
  _timer["loop"].stop();

  _timer["output"].start();
  LOG(iodata.log, "\nAverage conformal factor reached " << average(bssnSim.fields["phi_p"]) << "\n");
  LOG(iodata.log, "Ending simulation.\n");
  _timer["output"].stop();

  _timer["MAIN"].stop();

  LOG(iodata.log, endl << _timer << endl);

  return EXIT_SUCCESS;
}
