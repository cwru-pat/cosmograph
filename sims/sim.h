#ifndef COSMO_SIM_H
#define COSMO_SIM_H

#include "../cosmo_includes.h"
#include "../cosmo_types.h"

#include "../utils/Fourier.h"
#include "../utils/FRW.h"

#include "../IO/io.h"
#include "../components/bssn/bssn.h"
#include "../components/bssn/bardeen.h"

namespace cosmo
{

class CosmoSim
{
protected:
  idx_t step;
  idx_t num_steps;
  bool dt_flip;
  idx_t dt_flip_step;
  real_t t; ///< Time @ current step

  std::string simulation_type;
  IOData * iodata;
  Fourier * fourier;
  
  BSSN * bssnSim;

  Bardeen * bardeen;
  bool use_bardeen;

  int verbosity;

# if USE_COSMOTRACE
  bool ray_integrate;
  idx_t ray_flip_step;
  std::vector<RayTrace<real_t, idx_t> *> rays;
  bool simple_raytrace;
# endif

public:
  CosmoSim();
  virtual ~CosmoSim() {};

  // These functions will be called in main();
  // Each derived class should implement them.
  virtual void init() = 0;
  virtual void runStep() = 0;
  virtual void setICs() = 0;

  void simInit();
  void run();
  void runCommonStepTasks();

# if USE_COSMOTRACE
  void runRayTraceStep();
  void outputRayTraceStep();
# endif

  void prepBSSNOutput();
  void outputStateInformation();

  idx_t simNumNaNs();
  void setVerbosity(int verbosity_in);

};

} /* namespace cosmo */

#endif
