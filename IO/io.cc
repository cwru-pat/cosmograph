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
  iodata->log( "  USE_Z4c_DAMPING = " + stringify(USE_Z4c_DAMPING) );
  iodata->log( "  Z4c_K1_DAMPING_AMPLITUDE = " + stringify(Z4c_K1_DAMPING_AMPLITUDE) );
  iodata->log( "  Z4c_K2_DAMPING_AMPLITUDE = " + stringify(Z4c_K2_DAMPING_AMPLITUDE) );
  std::string stencil_order = STRINGIFY(STENCIL_ORDER_FUNCTION(Odx));
  iodata->log( "  STENCIL_ORDER = " + stencil_order);
  iodata->log( "  USE_GAMMA_DRIVER = " + stringify(USE_GAMMA_DRIVER) );
  iodata->log( "  USE_COSMO_CONST_POTENTIAL = " + stringify(USE_COSMO_CONST_POTENTIAL) );
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
  bool output_step = false;
  bool output_this_step = false;

  output_step = ( std::stoi(_config("IO_3D_grid_interval", "0")) > 0 );
  output_this_step = ( 0 == step % std::stoi(_config("IO_3D_grid_interval", "1")) );
  if( output_step && output_this_step )
  {
    for ( const auto &field_reg : bssn_fields ) {
      // look for names of the form: "IO_3D_field_r"
      if( std::stoi(_config( "IO_3D_" + field_reg.first , "0")) )
      {
        io_dump_3dslice(iodata, *bssn_fields[field_reg.first],
          "3D_" + (field_reg.first) + step_str);
      }
    }
  }
  
  output_step = ( std::stoi(_config("IO_2D_grid_interval", "0")) > 0 );
  output_this_step = ( 0 == step % std::stoi(_config("IO_2D_grid_interval", "1")) );
  if( output_step && output_this_step )
  {
    for ( const auto &field_reg : bssn_fields ) {
      // look for names of the form: "IO_2D_field_r"
      if( std::stoi(_config( "IO_2D_" + field_reg.first , "0")) )
      {
        io_dump_2dslice(iodata, *bssn_fields[field_reg.first],
          "2D_" + (field_reg.first) + step_str);
      }
    }
  }
  
  output_step = ( std::stoi(_config("IO_1D_grid_interval", "0")) > 0 );
  output_this_step = ( 0 == step % std::stoi(_config("IO_1D_grid_interval", "1")) );
  if( output_step && output_this_step )
  {
    for ( const auto &field_reg : bssn_fields ) {
      // look for names of the form: "IO_2D_field_r"
      if( std::stoi(_config( "IO_1D_" + field_reg.first , "0")) )
      {
        io_dump_strip(iodata, *bssn_fields[field_reg.first],
          "1D_" + (field_reg.first),
          std::stoi(_config( "IO_1D_" + field_reg.first + "_axis" , "1")),
          std::stoi(_config( "IO_1D_" + field_reg.first + "_xoffset" , "0")),
          std::stoi(_config( "IO_1D_" + field_reg.first + "_yoffset" , "0"))
        );
      }
    }
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
  bool output_step = ( std::stoi(_config("IO_powerspec_interval", "0")) > 0 );
  bool output_this_step = (0 == step % std::stoi(_config("IO_powerspec_interval", "1")));
  if( output_step && output_this_step )
  {
    fourier->powerDump(bssn_fields["DIFFphi_a"]->_array, iodata);
    fourier->powerDump(bssn_fields["DIFFr_a"]->_array, iodata);
  }
}

/**
 * @brief      Output constraint violation amplitudes
 */
void io_bssn_constraint_violation(IOData *iodata, idx_t step, BSSN * bssnSim)
{
  bool output_step = ( std::stoi(_config("IO_constraint_interval", "0")) > 0 );
  bool output_this_step = (0 == step % std::stoi(_config("IO_constraint_interval", "1")));
  if( output_step && output_this_step )
  {
    real_t H_calcs[7] = {0}, M_calcs[7] = {0}, G_calcs[7] = {0},
           A_calcs[7] = {0}, S_calcs[7] = {0};

    // Constraint Violation Calculations
    bssnSim->setConstraintCalcs(H_calcs, M_calcs, G_calcs,
                                A_calcs, S_calcs);

    io_dump_value(iodata, H_calcs[4], "H_violations", "\t"); // mean(H/[H])
    io_dump_value(iodata, H_calcs[5], "H_violations", "\t"); // stdev(H/[H])
    io_dump_value(iodata, H_calcs[6], "H_violations", "\t"); // max(H/[H])
    io_dump_value(iodata, H_calcs[2], "H_violations", "\n"); // max(H)

    io_dump_value(iodata, M_calcs[4], "M_violations", "\t"); // mean(M/[M])
    io_dump_value(iodata, M_calcs[5], "M_violations", "\t"); // stdev(M/[M])
    io_dump_value(iodata, M_calcs[6], "M_violations", "\t"); // max(M/[M])
    io_dump_value(iodata, M_calcs[2], "M_violations", "\n"); // max(M)

    io_dump_value(iodata, G_calcs[4], "G_violations", "\t"); // mean(G/[G])
    io_dump_value(iodata, G_calcs[5], "G_violations", "\t"); // stdev(G/[G])
    io_dump_value(iodata, G_calcs[6], "G_violations", "\t"); // max(G/[G])
    io_dump_value(iodata, G_calcs[2], "G_violations", "\n"); // max(G)

    io_dump_value(iodata, A_calcs[4], "A_violations", "\t"); // mean(A/[A])
    io_dump_value(iodata, A_calcs[5], "A_violations", "\t"); // stdev(A/[A])
    io_dump_value(iodata, A_calcs[6], "A_violations", "\t"); // max(A/[A])
    io_dump_value(iodata, A_calcs[2], "A_violations", "\n"); // max(A)

    io_dump_value(iodata, S_calcs[4], "S_violations", "\t"); // mean(S/[S])
    io_dump_value(iodata, S_calcs[5], "S_violations", "\t"); // stdev(S/[S])
    io_dump_value(iodata, S_calcs[6], "S_violations", "\t"); // max(S/[S])
    io_dump_value(iodata, S_calcs[2], "S_violations", "\n"); // max(S)
  }

  bool output_g11m1 = ( std::stoi(_config("IO_constraint_g11m1", "0")) > 0 );
  if( output_step && output_this_step && output_g11m1 )
  {
    arr_t & Dg11 = *bssnSim->fields["DIFFphi_a"];
    idx_t i, j, k;
    real_t mean = 0, stdev = 0, max = 0;
#   pragma omp parallel for default(shared) private(i, j, k) reduction(+:mean)
    LOOP3(i, j, k)
    {
      idx_t idx = NP_INDEX(i,j,k);
      mean += Dg11[idx];
#     pragma omp critical
      {
        max = max > std::fabs(Dg11[idx]) ? max : std::fabs(Dg11[idx]);
      }
    }
    mean /= POINTS;
#   pragma omp parallel for default(shared) private(i, j, k) reduction(+:stdev)
    LOOP3(i, j, k)
    {
      idx_t idx = NP_INDEX(i,j,k);
      stdev += std::pow(mean - Dg11[idx], 2.0);
    }
    stdev = std::sqrt(stdev/(POINTS-1));
    io_dump_value(iodata, stdev, "g11_violations", "\t");
    io_dump_value(iodata, mean, "g11_violations", "\t");
    io_dump_value(iodata, max, "g11_violations", "\n");
  }

}

/**
 * @brief Print constraint violation (useful for debugging)
 */
void io_print_constraint_violation(IOData *iodata, BSSN * bssnSim)
{
  real_t H_calcs[7] = {0}, M_calcs[7] = {0}, G_calcs[7] = {0},
         A_calcs[7] = {0}, S_calcs[7] = {0};
  bssnSim->setConstraintCalcs(H_calcs, M_calcs, G_calcs,
                              A_calcs, S_calcs);

  iodata->log( "\nConstraint Violation: " );

  // Constraint Violation Calculations
  iodata->log( "Max |H/[H]|: " + stringify(H_calcs[6]) );
  iodata->log( "Max |H|: " + stringify(H_calcs[2]) );

  iodata->log( "Max |M/[M]|: " + stringify(M_calcs[6]) );
  iodata->log( "Max |M|: " + stringify(M_calcs[2]) );
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
  bool output_step = ( std::stoi(_config("IO_bssnstats_interval", "0")) > 0 );
  bool output_this_step = (0 == step % std::stoi(_config("IO_bssnstats_interval", "1")));
  if( !output_step || !output_this_step )
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
  // average expantion output
#if USE_BSSN_SHIFT
  DETAILS(expN)
#endif
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
  bool output_step = false;
  bool output_this_step = false;

  output_step = ( std::stoi(_config("IO_3D_grid_interval", "0")) > 0 );
  output_this_step = ( 0 == step % std::stoi(_config("IO_3D_grid_interval", "1")) );
  if( output_step && output_this_step )
  {
    io_dump_3dslice(iodata, scalar->phi._array_a, "3D_scalar_phi." + step_str);
    io_dump_3dslice(iodata, scalar->Pi._array_a, "3D_scalar_Pi." + step_str);
  }
  
  output_step = ( std::stoi(_config("IO_2D_grid_interval", "0")) > 0 );
  output_this_step = ( 0 == step % std::stoi(_config("IO_2D_grid_interval", "1")) );
  if( output_step && output_this_step )
  {
    io_dump_2dslice(iodata, scalar->phi._array_a, "2D_scalar_phi." + step_str);
    io_dump_2dslice(iodata, scalar->Pi._array_a, "2D_scalar_Pi." + step_str);
  }
  
  output_step = ( std::stoi(_config("IO_1D_grid_interval", "0")) > 0 );
  output_this_step = ( 0 == step % std::stoi(_config("IO_1D_grid_interval", "1")) );
  if( output_step && output_this_step )
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

void io_print_strip(IOData *iodata, arr_t & field,
  int axis, idx_t n1, idx_t n2)
{
  std::cout << "\nField values:\n{";
  switch (axis)
  {
    case 1:
      for(idx_t i=0; i<NX; i++)
      {
        std::cout << std::setprecision(19) << field[INDEX(i,n1,n2)] << ", ";
      }
      break;
    case 2:
      for(idx_t j=0; j<NY; j++)
      {
        std::cout << std::setprecision(19) << field[INDEX(n1,j,n2)] << ", ";
      }
      break;
    case 3:
      for(idx_t k=0; k<NZ; k++)
      {
        std::cout << std::setprecision(19) << field[INDEX(n1,n2,k)] << ", ";
      }
      break;
  }

  std::cout << "}\n";

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

void io_print_particles(IOData *iodata, idx_t step, Particles *particles)
{
  bool output_step = ( std::stoi(_config("IO_particles", "0")) > 0 );
  bool output_this_step = (0 == step % std::stoi(_config("IO_particles", "1")));
  if( output_step && output_this_step )
  {
    // output misc. info about simulation here.
    char data[35];
    std::string dump_filename = iodata->dir() + "particles.dat.gz";

    gzFile datafile = gzopen(dump_filename.c_str(), "ab");
    if(datafile == Z_NULL) {
      iodata->log("Error opening file: " + dump_filename);
      return;
    }

// modified for 1-d output
    particle_vec * p_vec = particles->getParticleVec();
    for(particle_vec::iterator it = p_vec->begin(); it != p_vec->end(); ++it) {
      if(std::abs(it->p_a.X[1]) < 1e-10 && std::abs(it->p_a.X[2]) < 1e-10)
      {
        sprintf(data, "%.15g\t", (double) it->p_a.X[0]);
        gzwrite(datafile, data, strlen(data));
      }
      // sprintf(data, "%.15g\t", (double) it->p_a.X[1]);
      // gzwrite(datafile, data, strlen(data));
      // sprintf(data, "%.15g\t", (double) it->p_a.X[2]);
      // gzwrite(datafile, data, strlen(data));
    }

    sprintf(data, "\n");
    gzwrite(datafile, data, strlen(data));
    gzclose(datafile);
  }
}

} /* namespace cosmo */
