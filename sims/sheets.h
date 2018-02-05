#ifndef COSMO_SHEETS_SIM_H
#define COSMO_SHEETS_SIM_H

#include "sim.h"
#include "../components/phase_space_sheet/sheets.h"


namespace cosmo
{

/**
 * derived class based on CosmoSim class (sim.h)
 */
class SheetSim : public CosmoSim
{
protected:
  Sheet * sheetSim;
  real_t tot_mass;
public:
  SheetSim();
  ~SheetSim(){}

  void init();
  void setICs();
  void initSheetStep();
  void outputSheetStep();
  void runSheetStep();
  void runStep();
};

} /* namespace cosmo */

#endif
