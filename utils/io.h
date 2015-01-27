#ifndef COSMO_UTILS_IO_H
#define COSMO_UTILS_IO_H

#include "../cosmo.h"

namespace cosmo
{

typedef struct {
  std::string output_dir;
} IOData;

void io_init(IOData *iodata);

void io_dump_strip(real_t *field, int axis, idx_t n1, idx_t n2, IOData *iodata);

void io_dump_quantities(std::map <std::string, real_t *> & bssn_fields,
                        std::map <std::string, real_t *> & hydro_fields,
                        std::string filename, IOData *iodata);

void io_dump_2dslice(real_t *field, std::string filename, IOData *iodata);
void io_dump_3dslice(real_t *field, std::string filename, IOData *iodata);

}

#endif
