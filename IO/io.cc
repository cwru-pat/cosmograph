  
#include "io.h"
#include <hdf5.h>
#include <zlib.h>
#include <sys/stat.h>
#include <sstream>

#define DETAILS(field) \
  real_t avg_##field = conformal_average(*bssn_fields[#field "_a"], *bssn_fields["DIFFphi_a"], phi_FRW); \
  sprintf(data, "%.15g\t", (double) avg_##field); \
  gzwrite(datafile, data, strlen(data)); \
  real_t std_##field = conformal_standard_deviation(*bssn_fields[#field "_a"], *bssn_fields["DIFFphi_a"], phi_FRW, avg_##field); \
  sprintf(data, "%.15g\t", (double) std_##field); \
  gzwrite(datafile, data, strlen(data));


#define STRINGIFY_STRINGIFIER(function) #function
#define STRINGIFY_EVALUATOR(function) STRINGIFY_STRINGIFIER(function)
#define STRINGIFY(function) (STRINGIFY_EVALUATOR(function))


namespace cosmo
{

void log_defines(IOData *iodata)
{
  iodata->log( "Running with NX = " + stringify(NX)
                        + ", NY = " + stringify(NY)
                        + ", NZ = " + stringify(NZ) );
  iodata->log( "Running with dt = " + stringify(dt) );
  iodata->log( "Running with dx = " + stringify(dx) );

  iodata->log( "Other parameters:");
  iodata->log( "  H_LEN_FRAC = " + stringify(H_LEN_FRAC) );
  iodata->log( "  USE_REFERENCE_FRW = " + stringify(USE_REFERENCE_FRW) );
  iodata->log( "  NORMALIZE_GAMMAIJ_AIJ = " + stringify(NORMALIZE_GAMMAIJ_AIJ) );
  iodata->log( "  KO_ETA = " + stringify(KO_ETA) );
  iodata->log( "  BS_H_DAMPING_AMPLITUDE = " + stringify(BS_H_DAMPING_AMPLITUDE) );
  iodata->log( "  JM_K_DAMPING_AMPLITUDE = " + stringify(JM_K_DAMPING_AMPLITUDE) );
  iodata->log( "  USE_Z4c_DAMPING = " + stringify(USE_Z4c_DAMPING) );
  iodata->log( "  Z4c_K1_DAMPING_AMPLITUDE = " + stringify(Z4c_K1_DAMPING_AMPLITUDE) );
  iodata->log( "  Z4c_K2_DAMPING_AMPLITUDE = " + stringify(Z4c_K2_DAMPING_AMPLITUDE) );
  std::string stencil_order = STRINGIFY(STENCIL_ORDER_FUNCTION(Odx));
  iodata->log( "  STENCIL_ORDER = " + stencil_order);
  iodata->log( "  USE_HARMONIC_ALPHA = " + stringify(USE_HARMONIC_ALPHA) );
  iodata->log( "  USE_BSSN_SHIFT = " + stringify(USE_BSSN_SHIFT) );
}

/**
 * @brief      Print out a progress bar in the terminal (not log file)
 *
 * @param[in]  s     step number
 * @param[in]  maxs  maximum number of steps
 */
void io_show_progress(idx_t s, idx_t maxs)
{
  if(s > maxs)
  {
    std::cout << "Simulation has run more steps than allowed.";
    throw -1;
    return;
  }

  idx_t s_digits = (int) log10 ((double) s) + 1;
  if(s==0)
  {
    s_digits = 1;
  }
  idx_t maxs_digits = (int) log10 ((double) maxs);

  std::cout << " Running step " << s;
  for(int i=s_digits; i<=maxs_digits; ++i)
  {
    std::cout << " ";
  }
  std::cout << " / " << maxs;

  idx_t ndots = 20;
  std::cout << " [";
  for(int i=1; i<=ndots; ++i)
  {
    if(i <= ndots*s/maxs)
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

/**
 * @brief      Output 3d snapshot of BSSN fields in _a register
 *
 * @param      bssn_fields   map to bssn fields
 * @param[in]  step          step number (part of file names)
 * @param[in]  dim           # dimensions to output (1, 2, or 3)
 */
void io_bssn_fields_snapshot(IOData *iodata, idx_t step,
  map_t & bssn_fields)
{
  std::string step_str = std::to_string(step);
  if( step % std::stoi(_config["IO_3D_grid_interval"]) == 0 )
  {
    io_dump_3dslice(iodata, *bssn_fields["DIFFphi_a"], "3D_DIFFphi." + step_str);
    io_dump_3dslice(iodata, *bssn_fields["DIFFK_a"],   "3D_DIFFK."   + step_str);
    io_dump_3dslice(iodata, *bssn_fields["ricci_a"],   "3D_ricci."   + step_str);
    io_dump_3dslice(iodata, *bssn_fields["DIFFr_a"],   "3D_DIFFr."   + step_str);
  }
  
  if( step % std::stoi(_config["IO_2D_grid_interval"]) == 0 )
  {
    io_dump_2dslice(iodata, *bssn_fields["DIFFphi_a"], "2D_DIFFphi." + step_str);
    io_dump_2dslice(iodata, *bssn_fields["DIFFK_a"],   "2D_DIFFK."   + step_str);
    io_dump_2dslice(iodata, *bssn_fields["ricci_a"],   "2D_ricci."   + step_str);
    io_dump_2dslice(iodata, *bssn_fields["DIFFr_a"],   "2D_DIFFr."   + step_str);
  }
  
  if( step % std::stoi(_config["IO_1D_grid_interval"]) == 0 )
  {
    io_dump_strip(iodata, *bssn_fields["DIFFgamma11_a"], "1D_DIFFgamma11", 1, 0, 0);
    io_dump_strip(iodata, *bssn_fields["DIFFphi_a"],     "1D_DIFFphi", 1, 0, 0);
    io_dump_strip(iodata, *bssn_fields["DIFFK_a"],       "1D_DIFFK",   1, 0, 0);
    io_dump_strip(iodata, *bssn_fields["ricci_a"],       "1D_ricci",   1, 0, 0);
    io_dump_strip(iodata, *bssn_fields["DIFFr_a"],       "1D_DIFFr",   1, 0, 0);
  }
}

/**
 * @brief      Output power spectrum of some bssn fields
 *
 * @param      bssn_fields  map to bssn fields
 * @param      fourier      initialized Fourier instance
 */
void io_bssn_fields_powerdump(IOData *iodata, idx_t step,
  map_t & bssn_fields, Fourier *fourier)
{
  if( step % std::stoi(_config["IO_powerspec_interval"]) == 0 )
  {
    fourier->powerDump(bssn_fields["DIFFphi_a"]->_array, iodata);
    fourier->powerDump(bssn_fields["DIFFr_a"]->_array, iodata);
  }
}

/**
 * @brief      Output constraint violation amplitudes
 *
 * @param      iodata   { parameter_description }
 * @param      bssnSim  { parameter_description }
 */
void io_bssn_constraint_violation(IOData *iodata, idx_t step, BSSN * bssnSim)
{
  if( step % std::stoi(_config["IO_constraint_interval"]) == 0 )
  {
    real_t H_calcs[7] = {0}, M_calcs[7] = {0};

    // Constraint Violation Calculations
    bssnSim->setHamiltonianConstraintCalcs(H_calcs, false);
    io_dump_value(iodata, H_calcs[4], "H_violations", "\t"); // mean(H/[H])
    io_dump_value(iodata, H_calcs[5], "H_violations", "\t"); // stdev(H/[H])
    io_dump_value(iodata, H_calcs[6], "H_violations", "\t"); // max(H/[H])
    io_dump_value(iodata, H_calcs[2], "H_violations", "\n"); // max(H)

    bssnSim->setMomentumConstraintCalcs(M_calcs);
    io_dump_value(iodata, M_calcs[4], "M_violations", "\t"); // mean(M/[M])
    io_dump_value(iodata, M_calcs[5], "M_violations", "\t"); // stdev(M/[M])
    io_dump_value(iodata, M_calcs[6], "M_violations", "\t"); // max(M/[M])
    io_dump_value(iodata, M_calcs[2], "M_violations", "\n"); // max(M)
  }
}

/**
 * @brief      Output statistical information about simulation:
 *   averaged fields (FRW-analogues), standard deviations, FRW values, etc
 * 
 * @param      bssn_fields  map to bssn sim fields
 * @param      iodata       initialized IOData struct
 * @param      frw          frw instance from bssn sim
 */
void io_bssn_dump_statistics(IOData *iodata, idx_t step,
  map_t & bssn_fields, FRW<real_t> *frw)
{
  /* no output if not @ correct interval */
  if( step % std::stoi(_config["IO_bssnstats_interval"]) != 0 )
    return;

  std::string filename = _config["dump_file"];
  char data[35];

  // dump FRW quantities
  std::string dump_filename = iodata->dir() + filename + ".frwdat.gz";
  gzFile datafile = gzopen(dump_filename.c_str(), "ab");
  if(datafile == Z_NULL) {
    iodata->log("Error opening file: " + dump_filename);
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
  dump_filename = iodata->dir() + filename + ".dat.gz";
  datafile = gzopen(dump_filename.c_str(), "ab");
  if(datafile == Z_NULL) {
    iodata->log("Error opening file: " + dump_filename);
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
  sprintf(data, "%.15g\t", (double) volume_average(*bssn_fields["DIFFphi_a"], phi_FRW));
  gzwrite(datafile, data, strlen(data));

  gzwrite(datafile, "\n", strlen("\n"));
  gzclose(datafile);
}

#if USE_COSMOTRACE
/**
 * @brief      Write ray information to file
 *
 * @param      iodata  Initialized IOData
 * @param      rays    Vector of rays
 */
void io_raytrace_dump(IOData *iodata, idx_t step,
  std::vector<RayTrace<real_t, idx_t> *> const * rays)
{
  /* no output if not @ correct interval */
  if( step % std::stoi(_config["IO_raytrace_interval"]) != 0 )
    return;
  
  idx_t num_values = 9;
  idx_t num_rays = rays->size();

  real_t total_E = 0.0;
  real_t total_Phi = 0.0;
  real_t total_ell = 0.0;
  real_t total_ellrho = 0.0;
  real_t total_rho = 0.0;

  RaytraceData<real_t> tmp_rd = {0};
  real_t * ray_dump_values = new real_t[num_rays * num_values];
  
  for(idx_t n=0; n<num_rays; n++)
  {
    tmp_rd = (*rays)[n]->getRaytraceData();
    ray_dump_values[n*num_values + 0] = tmp_rd.E;
    ray_dump_values[n*num_values + 1] = tmp_rd.x[0];
    ray_dump_values[n*num_values + 2] = tmp_rd.x[1];
    ray_dump_values[n*num_values + 3] = tmp_rd.x[2];
    ray_dump_values[n*num_values + 4] = tmp_rd.Phi;
    ray_dump_values[n*num_values + 5] = tmp_rd.ell;
    ray_dump_values[n*num_values + 6] = tmp_rd.sig_Re;
    ray_dump_values[n*num_values + 7] = tmp_rd.sig_Im;
    ray_dump_values[n*num_values + 8] = tmp_rd.rho;

    total_E += tmp_rd.E;
    total_Phi += tmp_rd.Phi;
    total_ell += tmp_rd.ell;
    total_ellrho += tmp_rd.ell*tmp_rd.rho;
    total_rho += tmp_rd.E*tmp_rd.rho;
  }

  // Write data from individual rays to file
  std::string dataset_name = "step_" + std::to_string(step);
  std::string file_name = "raytracedata";
  io_dump_2d_array(iodata, ray_dump_values, num_rays, num_values,
    file_name, dataset_name);

  // write average, weighted average (biased) values to file
  if(total_rho == 0.0)
    total_rho = 1.0;
  io_dump_value(iodata, total_E / (real_t) num_rays, "raytracedata_avg", "\t");
  io_dump_value(iodata, total_Phi / (real_t) num_rays, "raytracedata_avg", "\t");
  io_dump_value(iodata, total_ell / (real_t) num_rays, "raytracedata_avg", "\t");
  io_dump_value(iodata, total_rho / (real_t) num_rays, "raytracedata_avg", "\t");
  io_dump_value(iodata, total_ellrho / total_rho, "raytracedata_avg", "\n");

  delete[] ray_dump_values;

  return;
}
#endif

/**
 * @brief      Output 3d snapshot of BSSN fields in _a register
 *
 * @param      bssn_fields   map to bssn fields
 * @param[in]  step          step number (part of file names)
 * @param[in]  dim           # dimensions to output (1, 2, or 3)
 */
void io_scalar_snapshot(IOData *iodata, idx_t step, Scalar * scalar)
{
  std::string step_str = std::to_string(step);
  if( step % std::stoi(_config["IO_3D_grid_interval"]) == 0 )
  {
    io_dump_3dslice(iodata, scalar->phi._array_a, "3D_scalar_phi." + step_str);
    io_dump_3dslice(iodata, scalar->Pi._array_a, "3D_scalar_Pi." + step_str);
  }
  
  if( step % std::stoi(_config["IO_2D_grid_interval"]) == 0 )
  {
    io_dump_2dslice(iodata, scalar->phi._array_a, "2D_scalar_phi." + step_str);
    io_dump_2dslice(iodata, scalar->Pi._array_a, "2D_scalar_Pi." + step_str);
  }
  
  if( step % std::stoi(_config["IO_1D_grid_interval"]) == 0 )
  {
    io_dump_strip(iodata, scalar->phi._array_a, "1D_scalar_phi", 1, NY/2, NZ/2);
    io_dump_strip(iodata, scalar->Pi._array_a, "1D_scalar_Pi", 1, NY/2, NZ/2);
  }
}

/**
 * @brief      Write full 3D slice to a file.
 *
 * @param      iodata    initialized IOData
 * @param      field     Field to write
 * @param[in]  filename  filename to write to (minus suffix)
 */
void io_dump_3dslice(IOData *iodata, arr_t & field, std::string filename)
{
  // dump all NX*NY*NZ points
  std::string dump_filename = iodata->dir() + filename + ".3d_grid.h5.gz";

  hid_t       file, space, dset, dcpl;  /* Handles */
  herr_t      status;
  htri_t      avail;
  H5Z_filter_t  filter_type;
  hsize_t     dims[3] = {(hsize_t) NX, (hsize_t) NY, (hsize_t) NZ},
              maxdims[3] = {H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED},
              chunk[3] = {6, 6, 6};

  file = H5Fcreate (dump_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  space = H5Screate_simple (3, dims, maxdims);
  dcpl = H5Pcreate (H5P_DATASET_CREATE);
  status = H5Pset_deflate (dcpl, 9);
  status = H5Pset_chunk (dcpl, 3, chunk);
  dset = H5Dcreate2 (file, "DS1", H5T_IEEE_F64LE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);

  status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, field._array);

  status = H5Pclose (dcpl);
  status = H5Dclose (dset);
  status = H5Sclose (space);
  status = H5Fclose (file);

  status = status; // suppress "unused" warning
  return;
}

/**
 * @brief      Write a 2D slice of a field to a file.
 *
 * @param      iodata    initialized IOData
 * @param      field     Field to write
 * @param[in]  filename  filename to write to (minus suffix)
 */
void io_dump_2dslice(IOData *iodata, arr_t & field, std::string filename)
{
  // dump the first NY*NZ points (a 2-d slice on a boundary)
  std::string dump_filename = iodata->dir() + filename + ".2d_grid.h5.gz";

  hid_t       file, space, dset, dcpl;  /* Handles */
  herr_t      status;
  htri_t      avail;
  H5Z_filter_t  filter_type;
  hsize_t     dims[2] = {(hsize_t) NY, (hsize_t) NZ},
              maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED},
              chunk[2] = {6, 6};

  file = H5Fcreate (dump_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  space = H5Screate_simple (2, dims, maxdims);
  dcpl = H5Pcreate (H5P_DATASET_CREATE);
  status = H5Pset_deflate (dcpl, 9);
  status = H5Pset_chunk (dcpl, 2, chunk);
  dset = H5Dcreate2 (file, "Dataset1", H5T_IEEE_F64LE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);

  status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, field._array);

  status = H5Pclose (dcpl);
  status = H5Dclose (dset);
  status = H5Sclose (space);
  status = H5Fclose (file);

  status = status; // suppress "unused" warning
  return;
}

/**
 * @brief      Output a 1-d slice of a simulation
 *
 * @param      field     field to output
 * @param[in]  file      filename to output strip in (minus suffix)
 * @param[in]  axis      axis along which to output data (1, 2, or 3)
 * @param[in]  n1        grid point number along first direction orthogonal to axis
 * @param[in]  n2        grid point number along second direction orthogonal to axis
 * @param      iodata    initialized IOData struct
 */
void io_dump_strip(IOData *iodata, arr_t & field, std::string file,
  int axis, idx_t n1, idx_t n2)
{
  std::string filename = iodata->dir() + file + ".strip.dat.gz";
  char data[35];

  gzFile datafile = gzopen(filename.c_str(), "ab");
  if(datafile == Z_NULL) {
    iodata->log("Error opening file: " + filename);
    return;
  }

  switch (axis)
  {
    case 1:
      for(idx_t i=0; i<NX; i++)
      {
        sprintf(data, "%.15g\t", (double) field[INDEX(i,n1,n2)]);
        gzwrite(datafile, data, strlen(data));
      }
      break;
    case 2:
      for(idx_t j=0; j<NY; j++)
      {
        sprintf(data, "%.15g\t", (double) field[INDEX(n1,j,n2)]);
        gzwrite(datafile, data, strlen(data));
      }
      break;
    case 3:
      for(idx_t k=0; k<NZ; k++)
      {
        sprintf(data, "%.15g\t", (double) field[INDEX(n1,n2,k)]);
        gzwrite(datafile, data, strlen(data));
      }
      break;
  }

  gzwrite(datafile, "\n", strlen("\n"));

  gzclose(datafile);

  return;
}

/**
 * @brief      Write a single (real) value (plus delimiter) to a file
 *
 * @param[in]  value      numerical value to write
 * @param      iodata     initialized IOData
 * @param[in]  filename   filename to output to (minus suffix)
 * @param[in]  delimiter  delimited (tab, newline, etc)
 */
void io_dump_value(IOData *iodata, real_t value, std::string filename,
  std::string delimiter)
{
  // output misc. info about simulation here.
  char data[35];
  std::string dump_filename = iodata->dir() + filename + ".dat.gz";

  gzFile datafile = gzopen(dump_filename.c_str(), "ab");
  if(datafile == Z_NULL) {
    iodata->log("Error opening file: " + dump_filename);
    return;
  }

  sprintf(data, "%.15g", (double) value);
  gzwrite(datafile, data, strlen(data));

  gzwrite(datafile, delimiter.c_str(), strlen(delimiter.c_str()));
  gzclose(datafile);
}


void io_dump_2d_array(IOData *iodata, real_t * array, idx_t n_x, idx_t n_y,
  std::string filename, std::string dataset_name)
{
  std::string dump_filename = iodata->dir() + filename + ".values.h5.gz";

  hid_t       file, space, dset, dcpl;  /* Handles */
  herr_t      status;
  htri_t      avail;
  H5Z_filter_t  filter_type;
  hsize_t     dims[2] = {(hsize_t) n_x, (hsize_t) n_y},
              maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED},
              chunk[2] = {6, 6};

  if(std::ifstream(dump_filename))
  {
    file = H5Fopen(dump_filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  }
  else
  {
    file = H5Fcreate(dump_filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
  }

  space = H5Screate_simple (2, dims, maxdims);
  dcpl = H5Pcreate (H5P_DATASET_CREATE);
  status = H5Pset_deflate (dcpl, 9);
  status = H5Pset_chunk (dcpl, 2, chunk);
  dset = H5Dcreate2 (file, dataset_name.c_str(), H5T_IEEE_F64LE, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);

  status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);

  status = H5Pclose (dcpl);
  status = H5Dclose (dset);
  status = H5Sclose (space);
  status = H5Fclose (file);

  status = status; // suppress "unused" warning
  return;
}

} /* namespace cosmo */
