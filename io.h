#ifndef COSMO_UTILS_IO_H
#define COSMO_UTILS_IO_H

#include "cosmo.h"
#include "io_data.h"

#define LOG(log_fout, message) \
  std::cout << message << std::flush; \
  log_fout << message; \
  log_fout.flush();

namespace cosmo
{

void io_init(IOData *iodata, std::string output_dir);

void io_config_backup(IOData *iodata, std::string config_file);

void io_show_progress(idx_t s, idx_t maxs);

void io_data_dump(std::map <std::string, real_t *> & bssn_fields,
                  std::map <std::string, real_t *> & hydro_fields,
                  IOData *iodata, idx_t step, Fourier *fourier, FRW<real_t> *frw);

void io_dump_strip(real_t *field, int axis, idx_t n1, idx_t n2, IOData *iodata);

void io_dump_statistics(std::map <std::string, real_t *> & bssn_fields,
                        std::map <std::string, real_t *> & hydro_fields,
                        IOData *iodata, FRW<real_t> *frw);

void io_dump_data(real_t value, IOData *iodata, std::string filename);

void io_dump_2dslice(real_t *field, std::string filename, IOData *iodata);
void io_dump_3dslice(real_t *field, std::string filename, IOData *iodata);

}

#endif
