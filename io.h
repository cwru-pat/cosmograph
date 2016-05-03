#ifndef COSMO_UTILS_IO_H
#define COSMO_UTILS_IO_H

#include "cosmo_includes.h"
#include "cosmo_types.h"
#include "globals.h"

#include "utils/Fourier.h"
#include "utils/reference_frw.h"
#include "utils/math.h"

#include "cosmotrace/raytrace.h"

#include "io_data.h"
#include "bssn.h"

#define LOG(log_fout, message) \
  std::cout << message << std::flush; \
  log_fout << message; \
  log_fout.flush();

namespace cosmo
{

void io_init(IOData *iodata, std::string output_dir);
void io_config_backup(IOData *iodata, std::string config_file);
void io_show_progress(idx_t s, idx_t maxs);

void io_bssn_fields_snapshot(IOData *iodata, idx_t step,
  std::map <std::string, real_t *> & bssn_fields);
void io_bssn_fields_powerdump(IOData *iodata, idx_t step,
  std::map <std::string, real_t *> & bssn_fields, Fourier *fourier);
void io_bssn_constraint_violation(IOData *iodata, idx_t step, BSSN * bssnSim);
void io_bssn_dump_statistics(IOData *iodata, idx_t step,
  std::map <std::string, real_t *> & bssn_fields, FRW<real_t> *frw);

void io_raytrace_dump(IOData *iodata, idx_t step,
  std::vector<RayTrace<real_t, idx_t> *> const * rays);

void io_dump_2dslice(IOData *iodata, real_t *field, std::string filename);
void io_dump_3dslice(IOData *iodata, real_t *field, std::string filename);
void io_dump_strip(IOData *iodata, real_t *field, std::string file,
  int axis, idx_t n1, idx_t n2);
void io_dump_value(IOData *iodata, real_t value, std::string filename,
  std::string delimiter);
void io_dump_2d_array(IOData *iodata, real_t * array, idx_t n_x, idx_t n_y,
  std::string filename, std::string dataset_name);

}

#endif
