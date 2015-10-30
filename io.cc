
#include "io.h"
#include <hdf5.h>
#include <sys/stat.h>
#include <sstream>

#define DETAILS(field) \
  real_t avg_##field = conformal_average(bssn_fields[#field "_a"], bssn_fields["phi_a"]); \
  sprintf(data, "%g\t", (double) avg_##field); \
  gzwrite(datafile, data, strlen(data)); \
  real_t std_##field = conformal_standard_deviation(bssn_fields[#field "_a"], bssn_fields["phi_a"], avg_##field); \
  sprintf(data, "%g\t", (double) std_##field); \
  gzwrite(datafile, data, strlen(data));


namespace cosmo
{

void io_init(IOData *iodata, std::string output_dir)
{
  iodata->output_dir = output_dir;

  /* ensure output_dir ends with '/', unless empty string is specified. */
  size_t len_dir_name = iodata->output_dir.length();
  if((iodata->output_dir)[len_dir_name - 1] != '/' && len_dir_name != 0)
  {
    iodata->output_dir += '/';
  }

  /* create data_dir */
  if(len_dir_name != 0)
    mkdir(iodata->output_dir.c_str(), 0755);

  iodata->slice_output_interval = stoi(_config["slice_output_interval"]);
  iodata->grid_output_interval = stoi(_config["grid_output_interval"]);
  iodata->meta_output_interval = stoi(_config["meta_output_interval"]);
  iodata->spec_output_interval = stoi(_config["spec_output_interval"]);
  iodata->dump_file = _config["dump_file"];

  iodata->log.open(iodata->output_dir + "log.txt");
  LOG(iodata->log, "Log file open.\n");
  LOG(iodata->log, "Running with NX = " << NX
                           << ", NY = " << NY
                           << ", NZ = " << NZ
                           << "\n");
  LOG(iodata->log, "Running with dt = " << dt << "\n");
  LOG(iodata->log, "Running with dx = " << dx << "\n");
}

void io_config_backup(IOData *iodata, std::string config_file)
{
  std::ifstream source(config_file, std::ios::binary);
  std::ofstream dest(iodata->output_dir + "config.txt", std::ios::binary);
  dest << source.rdbuf();
  source.close();
  dest.close();
}

void io_data_dump(std::map <std::string, real_t *> & bssn_fields,
                  std::map <std::string, real_t *> & static_field,
                  IOData *iodata, idx_t step, Fourier *fourier)
{
  if(step % iodata->slice_output_interval == 0)
  {
    io_dump_2dslice(bssn_fields["K_p"], "K_slice." + std::to_string(step), iodata);
    io_dump_2dslice(bssn_fields["phi_p"], "phi_slice." + std::to_string(step), iodata);
    io_dump_2dslice(static_field["D_a"], "UD_slice."  + std::to_string(step), iodata);
    io_dump_2dslice(bssn_fields["gamma11_p"], "gamma11." + std::to_string(step), iodata);
  }
  if(step % iodata->grid_output_interval == 0)
  {
    io_dump_3dslice(bssn_fields["gamma11_p"], "gamma11." + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["gamma12_p"], "gamma12." + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["gamma13_p"], "gamma13." + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["gamma22_p"], "gamma22." + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["gamma23_p"], "gamma23." + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["gamma33_p"], "gamma33." + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["phi_p"],     "phi."     + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["A11_p"],     "A11."     + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["A12_p"],     "A12."     + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["A13_p"],     "A13."     + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["A22_p"],     "A22."     + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["A23_p"],     "A23."     + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["A33_p"],     "A33."     + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["K_p"],       "K."       + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["ricci_a"],   "ricci."   + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["AijAij_a"],  "AijAij."  + std::to_string(step), iodata);
    io_dump_3dslice(static_field["D_a"],      "D."       + std::to_string(step), iodata);

    io_dump_3dslice(bssn_fields["r_a"],       "r_a."     + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["K_a"],       "K_a."     + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["phi_a"],     "phi_a."   + std::to_string(step), iodata);

    io_dump_3dslice(bssn_fields["KDx_a"],     "KDx_a."   + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["KDy_a"],     "KDy_a."   + std::to_string(step), iodata);
    io_dump_3dslice(bssn_fields["KDz_a"],     "KDz_a."   + std::to_string(step), iodata);
  }
  if(step % iodata->spec_output_interval == 0)
  {
    fourier->powerDump(bssn_fields["phi_p"], iodata);
    fourier->powerDump(bssn_fields["r_a"], iodata);
  }
  if(step % iodata->meta_output_interval == 0)
  {
    // some statistical values
    io_dump_statistics(bssn_fields, static_field, iodata);
  }
}

void io_show_progress(idx_t s, idx_t maxs) // terminal output only
{
  if(s > maxs)
  {
    std::cout << "Simulation has run more steps than allowed.";
    throw -1;
    return;
  }

  idx_t s_digits = (int) log10 ((double) s+1.0) + 1;
  if(s==0)
  {
    s_digits = 1;
  }
  idx_t maxs_digits = (int) log10 ((double) maxs) + 1;

  std::cout << " Running step " << s+1;
  for(int i=s_digits; i<=maxs_digits; ++i)
  {
    std::cout << " ";
  }
  std::cout << " / " << maxs;

  idx_t ndots = 20;
  std::cout << " [";
  for(int i=1; i<=ndots; ++i)
  {
    if(i <= ndots*(s+1)/maxs)
    {
      std::cout << "=";
    }
    else
    {
      std::cout << " ";
    }
  }
  std::cout << "]\r" << std::flush;
  return;
}

void io_dump_strip(real_t *field, int axis, idx_t n1, idx_t n2, IOData *iodata)
{
  std::string filename = iodata->output_dir + "strip.dat.gz";
  char data[20];

  gzFile datafile = gzopen(filename.c_str(), "ab");
  if(datafile == Z_NULL) {
    LOG(iodata->log, "Error opening file: " << filename.c_str() << "\n");
    return;
  }

  switch (axis)
  {
    case 1:
      for(idx_t i=0; i<NX; i++)
      {
        sprintf(data, "%g\t", (double) field[INDEX(i,n1,n2)]);
        gzwrite(datafile, data, strlen(data));
      }
      break;
    case 2:
      for(idx_t j=0; j<NY; j++)
      {
        sprintf(data, "%g\t", (double) field[INDEX(n1,j,n2)]);
        gzwrite(datafile, data, strlen(data));
      }
      break;
    case 3:
      for(idx_t k=0; k<NZ; k++)
      {
        sprintf(data, "%g\t", (double) field[INDEX(n1,n2,k)]);
        gzwrite(datafile, data, strlen(data));
      }
      break;
  }

  gzwrite(datafile, "\n", strlen("\n"));

  gzclose(datafile);

  return;
}


void io_dump_2dslice(real_t *field, std::string filename, IOData *iodata)
{
  // dump the first NY*NZ points (a 2-d slice on a boundary)
  std::string dump_filename = iodata->output_dir + filename + ".2d_grid.h5.gz";

  hid_t       file, space, dset, dcpl;  /* Handles */
  herr_t      status;
  htri_t      avail;
  H5Z_filter_t  filter_type;
  hsize_t     dims[2] = {NY, NZ},
              maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED},
              chunk[2] = {6, 6};

  file = H5Fcreate (dump_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  space = H5Screate_simple (2, dims, maxdims);
  dcpl = H5Pcreate (H5P_DATASET_CREATE);
  status = H5Pset_deflate (dcpl, 9);
  status = H5Pset_chunk (dcpl, 2, chunk);
  dset = H5Dcreate2 (file, "Dataset1", H5T_IEEE_F64LE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);

  status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, field);

  status = H5Pclose (dcpl);
  status = H5Dclose (dset);
  status = H5Sclose (space);
  status = H5Fclose (file);

  status = status; // suppress "unused" warning
  return;
}


/* 
 * Write full 3D slice to a file.
 */
void io_dump_3dslice(real_t *field, std::string filename, IOData *iodata)
{
  // dump all NX*NY*NZ points
  std::string dump_filename = iodata->output_dir + filename + ".3d_grid.h5.gz";

  hid_t       file, space, dset, dcpl;  /* Handles */
  herr_t      status;
  htri_t      avail;
  H5Z_filter_t  filter_type;
  hsize_t     dims[3] = {NX, NY, NZ},
              maxdims[3] = {H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED},
              chunk[3] = {6, 6, 6};

  file = H5Fcreate (dump_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  space = H5Screate_simple (3, dims, maxdims);
  dcpl = H5Pcreate (H5P_DATASET_CREATE);
  status = H5Pset_deflate (dcpl, 9);
  status = H5Pset_chunk (dcpl, 3, chunk);
  dset = H5Dcreate2 (file, "DS1", H5T_IEEE_F64LE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);

  status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, field);

  status = H5Pclose (dcpl);
  status = H5Dclose (dset);
  status = H5Sclose (space);
  status = H5Fclose (file);

  status = status; // suppress "unused" warning
  return;
}


void io_dump_statistics(std::map <std::string, real_t *> & bssn_fields,
                        std::map <std::string, real_t *> & static_field,
                        IOData *iodata)
{
  std::string filename = iodata->dump_file;

  // output misc. info about simulation here.
  char data[20];
  std::string dump_filename = iodata->output_dir + filename + ".dat.gz";

  gzFile datafile = gzopen(dump_filename.c_str(), "ab");
  if(datafile == Z_NULL) {
    LOG(iodata->log, "Error opening file: " << dump_filename << "\n");
    return;
  }

  // phi output
  real_t avg_phi = conformal_average(bssn_fields["phi_a"], bssn_fields["phi_a"]);
  sprintf(data, "%g\t", (double) avg_phi);
  gzwrite(datafile, data, strlen(data));

  real_t std_phi = conformal_standard_deviation(bssn_fields["phi_a"], bssn_fields["phi_a"], avg_phi);
  sprintf(data, "%g\t", (double) std_phi);
  gzwrite(datafile, data, strlen(data));

  // K output
  real_t avg_K = conformal_average(bssn_fields["K_a"], bssn_fields["phi_a"]);
  sprintf(data, "%g\t", (double) avg_K);
  gzwrite(datafile, data, strlen(data));

  real_t std_K = conformal_standard_deviation(bssn_fields["K_a"], bssn_fields["phi_a"], avg_K);
  sprintf(data, "%g\t", (double) std_K);
  gzwrite(datafile, data, strlen(data));

  // rho output
  real_t avg_rho = conformal_average(bssn_fields["r_a"], bssn_fields["phi_a"]);
  sprintf(data, "%g\t", (double) avg_rho);
  gzwrite(datafile, data, strlen(data));

  real_t std_rho = conformal_standard_deviation(bssn_fields["r_a"], bssn_fields["phi_a"], avg_rho);
  sprintf(data, "%g\t", (double) std_rho);
  gzwrite(datafile, data, strlen(data));

  // ricci output
  real_t avg_ricci = conformal_average(bssn_fields["ricci_a"], bssn_fields["phi_a"]);
  sprintf(data, "%g\t", (double) avg_ricci);
  gzwrite(datafile, data, strlen(data));

  real_t std_ricci = conformal_standard_deviation(bssn_fields["ricci_a"], bssn_fields["phi_a"], avg_ricci);
  sprintf(data, "%g\t", (double) std_ricci);
  gzwrite(datafile, data, strlen(data));

  // average volume
  sprintf(data, "%g\t", (double) volume_average(bssn_fields["phi_a"]));
  gzwrite(datafile, data, strlen(data));

/*
  DETAILS(A11)
  DETAILS(A12)
  DETAILS(A13)
  DETAILS(A22)
  DETAILS(A23)
  DETAILS(A33)
  DETAILS(gamma11)
  DETAILS(gamma12)
  DETAILS(gamma13)
  DETAILS(gamma22)
  DETAILS(gamma23)
  DETAILS(gamma33)
  DETAILS(Gamma1)
  DETAILS(Gamma2)
  DETAILS(Gamma3)
*/

  gzwrite(datafile, "\n", strlen("\n"));
  gzclose(datafile);
}

void io_dump_data(real_t value, IOData *iodata, std::string filename)
{
  // output misc. info about simulation here.
  char data[20];
  std::string dump_filename = iodata->output_dir + filename + ".dat.gz";

  gzFile datafile = gzopen(dump_filename.c_str(), "ab");
  if(datafile == Z_NULL) {
    LOG(iodata->log, "Error opening file: " << dump_filename << "\n");
    return;
  }

  sprintf(data, "%g", (double) value);
  gzwrite(datafile, data, strlen(data));

  gzwrite(datafile, "\n", strlen("\n"));
  gzclose(datafile);
}


} /* namespace cosmo */
