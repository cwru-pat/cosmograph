/**
 * @file BSSNGaugeHandler.h
 * @brief Functions to determine gauge evolution for the BSSN class.
 * Functions are determined via a config setting in the CosmoSim class.
 */

#ifndef COSMO_BSSN_GAUGE_FNS
#define COSMO_BSSN_GAUGE_FNS

#include "bssn_data.h"
#include <string>
#include <iostream>
#include "../../utils/ConfigParser.h"

namespace cosmo
{

class BSSN; // forward declaration of BSSN class

class BSSNGaugeHandler
{
private:
  typedef real_t (BSSNGaugeHandler::*bssn_gauge_func_t)(BSSNData *bd); ///< internal function pointer type

  // Reference to BSSN instance
  BSSN * bssn;

  // Maps to available functions
  std::map<std::string, bssn_gauge_func_t> lapse_gauge_map;
  std::map<std::string, std::map<std::string, bssn_gauge_func_t>> shift_gauge_map;

  // pointers to functions being used
  bssn_gauge_func_t lapse_fn; ///< Lapse evolution function
  bssn_gauge_func_t shift_fn1; ///< Shift evolution function
  bssn_gauge_func_t shift_fn2; ///< Shift evolution function
  bssn_gauge_func_t shift_fn3; ///< Shift evolution function

  // Generic, not evolving gauge
  real_t Static(BSSNData *bd);

  // Harmonic gauge lapse
  real_t HarmonicLapse(BSSNData *bd);

  // 1+log gauge slicing
  real_t gd_c; ///< Tunable gauge parameter
  real_t OnePlusLogLapse(BSSNData *bd);

  // Untested/experimental lapses
  real_t AnharmonicLapse(BSSNData *bd);
  real_t ConformalSyncLapse(BSSNData *bd);

  // Gamma driver shift function
  real_t GammaDriverShift1(BSSNData *bd);
  real_t GammaDriverShift2(BSSNData *bd);
  real_t GammaDriverShift3(BSSNData *bd);

  // Damped wave gauge
  real_t dw_mu_l; ///< damped wave "mu_l" parameter
  real_t dw_mu_s; ///< damped wave "mu_s" parameter
  real_t dw_p; ///< damped wave "p" parameter
  real_t DampedWaveLapse(BSSNData *bd);
  real_t DampedWaveShift1(BSSNData *bd);
  real_t DampedWaveShift2(BSSNData *bd);
  real_t DampedWaveShift3(BSSNData *bd);

  // AwA Gauge Wave test lapse
  real_t gauge_wave_dir; ///< wave direction of prop. (\in {1,2,3})
  real_t AwAGaugeWaveLapse(BSSNData *bd);

  // AwA Shifted Gauge Wave test gauge
  // also uses gauge_wave_dir
  real_t AwAShiftedWaveLapse(BSSNData *bd);
  real_t AwAShiftedWaveShift1(BSSNData *bd);
  real_t AwAShiftedWaveShift2(BSSNData *bd);
  real_t AwAShiftedWaveShift3(BSSNData *bd);

  real_t k_driver_coeff;
  real_t TestKDriverLapse(BSSNData *bd);
  real_t TestAijDriverLapse(BSSNData *bd);
  real_t AijDriverShift1(BSSNData *bd);
  real_t AijDriverShift2(BSSNData *bd);
  real_t AijDriverShift3(BSSNData *bd);

  // Map of strings to functions
  void _initGaugeMaps()
  {
    // Lapse functions
    lapse_gauge_map["Static"] = &BSSNGaugeHandler::Static;
    lapse_gauge_map["Harmonic"] = &BSSNGaugeHandler::HarmonicLapse;
    lapse_gauge_map["Anharmonic"] = &BSSNGaugeHandler::AnharmonicLapse;
    lapse_gauge_map["OnePlusLog"] = &BSSNGaugeHandler::OnePlusLogLapse;
    lapse_gauge_map["DampedWave"] = &BSSNGaugeHandler::DampedWaveLapse;
    lapse_gauge_map["ConformalSync"] = &BSSNGaugeHandler::ConformalSyncLapse;
    lapse_gauge_map["AwAGaugeWave"] = &BSSNGaugeHandler::AwAGaugeWaveLapse;
    lapse_gauge_map["AwAShiftedWave"] = &BSSNGaugeHandler::AwAShiftedWaveLapse;

    lapse_gauge_map["TestKDriverLapse"] = &BSSNGaugeHandler::TestKDriverLapse;
    lapse_gauge_map["TestAijDriverLapse"] = &BSSNGaugeHandler::TestAijDriverLapse;

    // Shift functions
    // Static gauge
    shift_gauge_map["Static"]["1"] = &BSSNGaugeHandler::Static;
    shift_gauge_map["Static"]["2"] = &BSSNGaugeHandler::Static;
    shift_gauge_map["Static"]["3"] = &BSSNGaugeHandler::Static;
    // gamma driver
    shift_gauge_map["GammaDriver"]["1"] = &BSSNGaugeHandler::GammaDriverShift1;
    shift_gauge_map["GammaDriver"]["2"] = &BSSNGaugeHandler::GammaDriverShift2;
    shift_gauge_map["GammaDriver"]["3"] = &BSSNGaugeHandler::GammaDriverShift3;
    // Damped wave
    shift_gauge_map["DampedWave"]["1"] = &BSSNGaugeHandler::DampedWaveShift1;
    shift_gauge_map["DampedWave"]["2"] = &BSSNGaugeHandler::DampedWaveShift2;
    shift_gauge_map["DampedWave"]["3"] = &BSSNGaugeHandler::DampedWaveShift3;
    // AwA shifted wave test
    shift_gauge_map["AwAShiftedWave"]["1"] = &BSSNGaugeHandler::AwAShiftedWaveShift1;
    shift_gauge_map["AwAShiftedWave"]["2"] = &BSSNGaugeHandler::AwAShiftedWaveShift2;
    shift_gauge_map["AwAShiftedWave"]["3"] = &BSSNGaugeHandler::AwAShiftedWaveShift3;


    shift_gauge_map["AijDriverShift"]["1"] = &BSSNGaugeHandler::AijDriverShift1;
    shift_gauge_map["AijDriverShift"]["2"] = &BSSNGaugeHandler::AijDriverShift2;
    shift_gauge_map["AijDriverShift"]["3"] = &BSSNGaugeHandler::AijDriverShift3;
  }

  void _initDefaultParameters(ConfigParser *config)
  {
    gauge_wave_dir = std::stoi((*config)("gauge_wave_dir", "1"));
    dw_mu_l = std::stod((*config)("dw_mu_l", "0.0"));
    dw_mu_s = std::stod((*config)("dw_mu_s", "0.0"));
    dw_p = std::stod((*config)("dw_p", "0.0"));
    gd_c = std::stod((*config)("gd_c", "0.0"));

    k_driver_coeff = std::stod((*config)("k_driver_coeff", "0.04"));
  }

public:

  /**
   * @brief Initialize with static, non-evolving gauge
   */
  BSSNGaugeHandler()
  {
    ConfigParser emptyConfig;
    _initGaugeMaps();
    _initDefaultParameters(&emptyConfig);
    setLapseFn("Static");
    setShiftFn("Static");
  }

  /**
   * @brief Initialize with gauge determined by config file (default to a "static", non-evolving gauge)
   */
  BSSNGaugeHandler(ConfigParser *config, BSSN *bssnSim)
  {
    _initGaugeMaps();
    _initDefaultParameters(config);
    setLapseFn((*config)("lapse", "Static"));
    setShiftFn((*config)("shift", "Static"));

    bssn = bssnSim;
  }

  /**
   * @brief Set the lapse function
   */
  void setLapseFn(std::string name)
  {
    if ( lapse_gauge_map.find(name) == lapse_gauge_map.end() )
    {
      std::cout << "Error: Lapse gauge not found: `" << name << "`!\n";
      throw -1;
    }

    std::cout << "Using lapse: `" << name << "`.\n";
    lapse_fn = lapse_gauge_map[name];
  }

  /**
   * @brief Set the shift function
   */
  void setShiftFn(std::string name)
  {
    // Shift needs to be enabled for non-trivial evolution
    if(name != "Static" && !USE_BSSN_SHIFT)
    {
      std::cerr << "Code must be compiled with shift enabled to use non-Static shift!";
      throw -1;
    }

    // Gamma driver needs an extra field compiled in
    if(name == "gammadriver" && !USE_GAMMA_DRIVER)
    {
      std::cerr << "Code must be compiled with gamma driver enabled to use the gamma driver gauge.";
      throw -1;
    }

    
    if ( shift_gauge_map.find(name) == shift_gauge_map.end() )
    {
      std::cout << "Error: Shift gauge not found: `" << name << "`!\n";
      throw -1;
    }

    shift_fn1 = shift_gauge_map[name]["1"];
    shift_fn2 = shift_gauge_map[name]["2"];
    shift_fn3 = shift_gauge_map[name]["3"];
  }

  /**
   * @brief Lapse evolution function for BSSN class to call
   */
  real_t ev_lapse(BSSNData *bd)
  {
    return (*this.*lapse_fn)(bd);
  }

  /**
   * @brief Shift in x-dir evolution function for BSSN class to call
   */
  real_t ev_shift1(BSSNData *bd)
  {
    return (*this.*shift_fn1)(bd);
  }

  /**
   * @brief Shift in y-dir evolution function for BSSN class to call
   */
  real_t ev_shift2(BSSNData *bd)
  {
    return (*this.*shift_fn2)(bd);
  }

  /**
   * @brief Shift in z-dir evolution function for BSSN class to call
   */
  real_t ev_shift3(BSSNData *bd)
  {
    return (*this.*shift_fn3)(bd);
  }

};

}

#endif
