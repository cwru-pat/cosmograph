/** @file sheet.h
 * @brief Class implementing a "sheet" / interpolated phase-space
 * method for solving the Vlasov equation.
 *  Vocabulary:
 *   - Phase-space: Reference to 6-d (x, y, z, vx, vy, vz) space in which
 *     matter lives.
 *   - Sheet: 3-dimensional slice through 6-d phase-space on which all matter
 *      content lives.
 *   - Sheet coordinates: coordinates parameterizing the 3-d slice (s1, s2, s3)
 *   - Metric coordinates: x, y, z
 *   - Position fields: \vec{x}(\vec{s})
 *
 *  Two coordinate systems:
 *  - Sheet coordinates
 *    position and velocity fields live on these coordinates
 *  - Metric coordinates
 *    Density and metric
 */
#ifndef COSMO_SHEET_H
#define COSMO_SHEET_H

#include "../../utils/RK4Register.h"
#include "../../utils/Array.h"
#include "../../utils/Timer.h"
#include "../../cosmo_types.h"
#include "../bssn/bssn.h"
#include "../../utils/TriCubicInterpolator.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>
#include "zlib.h"
#include <limits>
#include <sstream>
#include <iomanip>
#include <fstream>

namespace cosmo
{


/**
 * Class used to run a sheet sim.
 */
class Sheet
{
public:
  // Simulation information
  idx_t nx, ny, nz;
  idx_t ns1, ns2, ns3; ///< Phase-space sheet resolution
  real_t lx, ly, lz; ///< Metric grid physical dimensions
  real_t dx, dy, dz; ///< Metric grid physical spacing
  
  register_t Dx, Dy, Dz; ///< Metric-space density
  register_t vx, vy, vz; ///< Phase-space velocity fields
  
  arr_t tmp; ///< Array for misc. tmp storage (such as deconvolving)

  idx_t step;

  // internal types
  enum carrierCountScheme { per_dx = 0, per_ds = 1};
  carrierCountScheme carrier_count_scheme;

  enum depositScheme { CIC = 0, PCS = 1, CINT = 2 };
  depositScheme deposit;

  idx_t carriers_per_dx,
        carriers_per_dy,
        carriers_per_dz;

  Sheet();
  ~Sheet();

  /**
   * Functions to convert s-indices to non-displaced coordinates
   */
  real_t _S1IDXtoX0(idx_t s1) { return (real_t)s1*lx/ns1; }
  real_t _S2IDXtoY0(idx_t s2) { return (real_t)s2*ly/ns2; }
  real_t _S3IDXtoZ0(idx_t s3) { return (real_t)s3*lz/ns3; }
  real_t _S1IDXtoX0(real_t s1) {  return (real_t)s1*lx/(real_t)ns1; }
  real_t _S2IDXtoY0(real_t s2) {  return (real_t)s2*ly/(real_t)ns2; }
  real_t _S3IDXtoZ0(real_t s3) {  return (real_t)s3*lz/(real_t)ns3; }

  void _MassDeposit(real_t weight, real_t x_idx, real_t y_idx, real_t z_idx,
                    arr_t &rho);

  void _CICDeposit(real_t weight, real_t x_idx, real_t y_idx, real_t z_idx,
                   arr_t &rho);

  void _PCSDeposit(real_t weight, real_t x_idx, real_t y_idx, real_t z_idx,
                   arr_t &rho);

  void _CINTDeposit(real_t weight, real_t x_idx, real_t y_idx, real_t z_idx,
                   arr_t &rho);

  void _deconvolve(arr_t &field);

  /**
   * Compute conribution to rho(x) from data in a phase-space
   * sheet voxel and add to rho(x) grid. Do so via 1501.01959 mass
   * deposition scheme.
   * TODO: improve; consider higher-order or analytic versions of this
   */
  void _pushSheetMassToRho(idx_t s1, idx_t s2, idx_t s3);

  /**
   * Set metric potentials (or really just derivatives thereof)
   */
  void _setMetricPotentials();

  /**
   * Intermediate RK4 calculations for phase space fields
   */
  void _RK4Calc();

  void _stepInit();

  void _K1Finalize();

  void _K2Finalize();

  void _K3Finalize();

  void _K4Finalize();

  real_t _getXRangeInSVoxel(register_t & DX, idx_t s1_idx, idx_t s2_idx,
                            idx_t s3_idx, real_t X0_lower, real_t X0_upper);

  void addBSSNSource(BSSN *bssn, real_t tot_mass);

  void RKStep(BSSN *bssn);

  void stepInit();

  void K1Finalize();
  void K2Finalize();
  void K3Finalize();
  void K4Finalize();
  
  arr_t d1alpha_a, d2alpha_a, d3alpha_a;

  arr_t d1gammai11_a, d1gammai22_a, d1gammai33_a,
    d1gammai12_a, d1gammai13_a, d1gammai23_a;

  arr_t d2gammai11_a, d2gammai22_a, d2gammai33_a,
    d2gammai12_a, d2gammai13_a, d2gammai23_a;

  arr_t d3gammai11_a, d3gammai22_a, d3gammai33_a,
    d3gammai12_a, d3gammai13_a, d3gammai23_a;

  arr_t d1beta1_a, d1beta2_a, d1beta3_a;
  arr_t d2beta1_a, d2beta2_a, d2beta3_a;
  arr_t d3beta1_a, d3beta2_a, d3beta3_a;

};

} //namespace cosmo
#endif // include guard
