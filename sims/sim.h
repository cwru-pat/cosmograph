#ifndef COSMO_SIM_H
#define COSMO_SIM_H

#include "../cosmo_includes.h"
#include "../cosmo_types.h"
#include "../cosmo_globals.h"

#include "../utils/Fourier.h"
#include "../utils/FRW.h"

#include "../ICs/ICs.h"
#include "../IO/io.h"

#include "../bssn/bssn.h"


namespace cosmo
{

class CosmoSim
{
protected:
  idx_t step;
  idx_t num_steps;
  bool dt_flip;
  idx_t dt_flip_step;

  std::string simulation_type;
  IOData * iodata;
  Fourier * fourier;
  
  BSSN * bssnSim;

  int verbosity;

  bool ray_integrate;
  idx_t ray_flip_step;
  std::vector<RayTrace<real_t, idx_t> *> rays;

public:
  CosmoSim();
  ~CosmoSim();

  // functions will be called in main()
  virtual void init() = 0;
  virtual void runStep() = 0;

  void simInit();
  void run();
  void runCommonStepTasks();

  void runRayTraceStep();
  void outputRayTraceStep();

  void prepBSSNOutput();
  void outputStateInformation();

  idx_t simNumNaNs();
  void setVerbosity(int verbosity_in);

};

} /* namespace cosmo */

#endif
