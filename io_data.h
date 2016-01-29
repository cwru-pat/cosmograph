#ifndef COSMO_IO_DATA
#define COSMO_IO_DATA

namespace cosmo
{

typedef struct {
  std::string output_dir;
  std::string dump_file;
  std::ofstream log;
  idx_t slice_output_interval;
  idx_t grid_output_interval;
  idx_t meta_output_interval;
  idx_t spec_output_interval;
} IOData;

} /* namespace cosmo */

#endif
