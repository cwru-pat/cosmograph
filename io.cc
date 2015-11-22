
#include "io.h"
#include <hdf5.h>
#include <sys/stat.h>
#include <sstream>

#define DETAILS(field) \
  real_t avg_##field = conformal_average(bssn_fields[#field "_a"], bssn_fields["DIFFphi_a"], phi_FRW); \
  sprintf(data, "%.15g\t", (double) avg_##field); \
  gzwrite(datafile, data, strlen(data)); \
  real_t std_##field = conformal_standard_deviation(bssn_fields[#field "_a"], bssn_fields["DIFFphi_a"], phi_FRW, avg_##field); \
  sprintf(data, "%.15g\t", (double) std_##field); \
  gzwrite(datafile, data, strlen(data));


#define STRINGIFY_STRINGIFIER(function) #function
#define STRINGIFY_EVALUATOR(function) STRINGIFY_STRINGIFIER(function)
#define STRINGIFY(function) STRINGIFY_EVALUATOR(function)


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

  /* Create data_dir if needed */
  if(len_dir_name != 0)
    mkdir(iodata->output_dir.c_str(), 0755);

  /* check for conflicting data in directory */
  if(std::ifstream(iodata->output_dir + "log.txt"))
  {
    LOG(iodata->log, "Data files in output directory seem to already exist!");
    int s=1;
    while(std::ifstream(
        (iodata->output_dir).substr(0,iodata->output_dir.length()-1) + "." + std::to_string(s) + "/log.txt"
      ))
      s += 1;

    iodata->output_dir = (iodata->output_dir).substr(0,iodata->output_dir.length()-1) + "." + std::to_string(s);
    LOG(iodata->log, " Using '" << iodata->output_dir << "' instead.\n");
    iodata->output_dir += '/';
    mkdir(iodata->output_dir.c_str(), 0755);
  }

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

  LOG(iodata->log, "Other parameters: \n");
  LOG(iodata->log, "  H_LEN_FRAC = " << H_LEN_FRAC << "\n");
  LOG(iodata->log, "  USE_REFERENCE_FRW = " << USE_REFERENCE_FRW << "\n");
  LOG(iodata->log, "  NORMALIZE_GAMMAIJ_AIJ = " << NORMALIZE_GAMMAIJ_AIJ << "\n");
  LOG(iodata->log, "  KO_ETA = " << KO_ETA << "\n");
  LOG(iodata->log, "  BS_H_DAMPING_AMPLITUDE = " << BS_H_DAMPING_AMPLITUDE << "\n");
  LOG(iodata->log, "  JM_K_DAMPING_AMPLITUDE = " << JM_K_DAMPING_AMPLITUDE << "\n");
  LOG(iodata->log, "  USE_Z4c_DAMPING = " << USE_Z4c_DAMPING << "\n");
  LOG(iodata->log, "  Z4c_K1_DAMPING_AMPLITUDE = " << Z4c_K1_DAMPING_AMPLITUDE << "\n");
  LOG(iodata->log, "  Z4c_K2_DAMPING_AMPLITUDE = " << Z4c_K2_DAMPING_AMPLITUDE << "\n");
  LOG(iodata->log, "  STENCIL_ORDER = " << STRINGIFY(STENCIL_ORDER_FUNCTION(Odx)) << "\n");
  LOG(iodata->log, "  USE_HARMONIC_ALPHA = " << USE_HARMONIC_ALPHA << "\n");
  LOG(iodata->log, "  USE_BSSN_SHIFT = " << USE_BSSN_SHIFT << "\n");
}

void io_config_backup(IOData *iodata, std::string config_file)
{
  std::ifstream source(config_file, std::ios::binary);
  std::ofstream dest(iodata->output_dir + "config.txt", std::ios::binary);
  dest << source.rdbuf();
  source.close();
  dest.close();
}

void io_data_dump(std::map <std::string, periodicArray<idx_t, real_t> *> & bssn_fields,
                  std::map <std::string, periodicArray<idx_t, real_t> *> & static_field,
                  IOData *iodata, idx_t step, Fourier *fourier, FRW<real_t> *frw)
{
  if(step % iodata->slice_output_interval == 0)
  {
    // io_dump_2dslice(bssn_fields["DIFFK_a"], "DIFFK_slice." + std::to_string(step), iodata);
    // io_dump_2dslice(bssn_fields["DIFFphi_a"], "DIFFphi_slice." + std::to_string(step), iodata);
    // io_dump_2dslice(static_field["DIFFD_a"], "DIFFUD_slice."  + std::to_string(step), iodata);
    // io_dump_2dslice(bssn_fields["DIFFgamma11_a"], "DIFFgamma11." + std::to_string(step), iodata);
  }
  if(step % iodata->grid_output_interval == 0)
  {
    // io_dump_3dslice(bssn_fields["DIFFgamma11_a"], "DIFFgamma11." + std::to_string(step), iodata);
    // io_dump_3dslice(bssn_fields["DIFFgamma12_a"], "DIFFgamma12." + std::to_string(step), iodata);
    // io_dump_3dslice(bssn_fields["DIFFgamma13_a"], "DIFFgamma13." + std::to_string(step), iodata);
    // io_dump_3dslice(bssn_fields["DIFFgamma22_a"], "DIFFgamma22." + std::to_string(step), iodata);
    // io_dump_3dslice(bssn_fields["DIFFgamma23_a"], "DIFFgamma23." + std::to_string(step), iodata);
    // io_dump_3dslice(bssn_fields["DIFFgamma33_a"], "DIFFgamma33." + std::to_string(step), iodata);
    // io_dump_3dslice(bssn_fields["DIFFphi_a"],     "DIFFphi."     + std::to_string(step), iodata);
    // io_dump_3dslice(bssn_fields["DIFFK_a"],       "DIFFK."       + std::to_string(step), iodata);
    // io_dump_3dslice(bssn_fields["ricci_a"],   "ricci."   + std::to_string(step), iodata);
    // io_dump_3dslice(bssn_fields["AijAij_a"],  "AijAij."  + std::to_string(step), iodata);
    // io_dump_3dslice(static_field["DIFFD_a"],      "DIFFUD."       + std::to_string(step), iodata);

    // io_dump_3dslice(bssn_fields["KDx_a"],     "KDx_a."   + std::to_string(step), iodata);
    // io_dump_3dslice(bssn_fields["KDy_a"],     "KDy_a."   + std::to_string(step), iodata);
    // io_dump_3dslice(bssn_fields["KDz_a"],     "KDz_a."   + std::to_string(step), iodata);
  }
  if(step % iodata->spec_output_interval == 0)
  {
    // fourier->powerDump(bssn_fields["DIFFphi_a"], iodata);
    // fourier->powerDump(bssn_fields["DIFFr_a"], iodata);
  }
  if(step % iodata->meta_output_interval == 0)
  {
    // some statistical values
    io_dump_statistics(bssn_fields, static_field, iodata, frw);
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

void io_dump_strip(periodicArray<idx_t, real_t> *field, int axis, idx_t n1, idx_t n2, IOData *iodata)
{
  std::string filename = iodata->output_dir + "strip.dat.gz";
  char data[35];

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
        sprintf(data, "%.15g\t", (double) (*field)[INDEX(i,n1,n2)]);
        gzwrite(datafile, data, strlen(data));
      }
      break;
    case 2:
      for(idx_t j=0; j<NY; j++)
      {
        sprintf(data, "%.15g\t", (double) (*field)[INDEX(n1,j,n2)]);
        gzwrite(datafile, data, strlen(data));
      }
      break;
    case 3:
      for(idx_t k=0; k<NZ; k++)
      {
        sprintf(data, "%.15g\t", (double) (*field)[INDEX(n1,n2,k)]);
        gzwrite(datafile, data, strlen(data));
      }
      break;
  }

  gzwrite(datafile, "\n", strlen("\n"));

  gzclose(datafile);

  return;
}


void io_dump_2dslice(periodicArray<idx_t, real_t> *field, std::string filename, IOData *iodata)
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

  status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, field->_array);

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
void io_dump_3dslice(periodicArray<idx_t, real_t> *field, std::string filename, IOData *iodata)
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

  status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, field->_array);

  status = H5Pclose (dcpl);
  status = H5Dclose (dset);
  status = H5Sclose (space);
  status = H5Fclose (file);

  status = status; // suppress "unused" warning
  return;
}


void io_dump_statistics(std::map <std::string, periodicArray<idx_t, real_t> *> & bssn_fields,
                        std::map <std::string, periodicArray<idx_t, real_t> *> & static_field,
                        IOData *iodata, FRW<real_t> *frw)
{
  std::string filename = iodata->dump_file;
  char data[35];

  // dump FRW quantities
  std::string dump_filename = iodata->output_dir + filename + ".frwdat.gz";
  gzFile datafile = gzopen(dump_filename.c_str(), "ab");
  if(datafile == Z_NULL) {
    LOG(iodata->log, "Error opening file: " << dump_filename << "\n");
    return;
  }

  real_t phi_FRW = frw->get_phi();
  sprintf(data, "%.15g\t", (double) phi_FRW);
  gzwrite(datafile, data, strlen(data));
  real_t K_FRW = frw->get_K();
  sprintf(data, "%.15g\t", (double) K_FRW);
  gzwrite(datafile, data, strlen(data));
  real_t rho_FRW = frw->get_rho();
  sprintf(data, "%.15g\t", (double) rho_FRW);
  gzwrite(datafile, data, strlen(data));
  real_t S_FRW = frw->get_S();
  sprintf(data, "%.15g\t", (double) S_FRW);
  gzwrite(datafile, data, strlen(data));

  gzwrite(datafile, "\n", strlen("\n"));
  gzclose(datafile);


  // output misc. info about simulation here.
  dump_filename = iodata->output_dir + filename + ".dat.gz";
  datafile = gzopen(dump_filename.c_str(), "ab");
  if(datafile == Z_NULL) {
    LOG(iodata->log, "Error opening file: " << dump_filename << "\n");
    return;
  }

  // phi output
  DETAILS(DIFFphi)
  // K output
  DETAILS(DIFFK)
  // rho output
  DETAILS(DIFFr)
  // ricci output
  DETAILS(ricci)
  // average volume
  sprintf(data, "%.15g\t", (double) volume_average(bssn_fields["DIFFphi_a"], phi_FRW));
  gzwrite(datafile, data, strlen(data));

  gzwrite(datafile, "\n", strlen("\n"));
  gzclose(datafile);
}

void io_dump_data(real_t value, IOData *iodata, std::string filename)
{
  // output misc. info about simulation here.
  char data[35];
  std::string dump_filename = iodata->output_dir + filename + ".dat.gz";

  gzFile datafile = gzopen(dump_filename.c_str(), "ab");
  if(datafile == Z_NULL) {
    LOG(iodata->log, "Error opening file: " << dump_filename << "\n");
    return;
  }

  sprintf(data, "%.15g", (double) value);
  gzwrite(datafile, data, strlen(data));

  gzwrite(datafile, "\n", strlen("\n"));
  gzclose(datafile);
}


} /* namespace cosmo */
