#include "io.h"
#include <hdf5.h>
#include <math.h>
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

#if USE_LONG_DOUBLES
# define H5T_TO_USE H5T_NATIVE_LDOUBLE
# define H5T_ALLOC H5T_NATIVE_LDOUBLE
#else
# define H5T_TO_USE H5T_NATIVE_DOUBLE
# define H5T_ALLOC H5T_IEEE_F64LE
#endif

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
      // look for names of the form: "IO_1D_field_r"
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
  if(!output_step) return;
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
  if(!output_step) return;

  bool output_this_step = (0 == step % std::stoi(_config("IO_constraint_interval", "1")));

  // whether dump 1D constraint everywhere
  //  bool dump_1d_hamiltonian_constraint = ( std::stoi(_config("IO_1D_hamiltonian_constraint", "0")) > 0 );
  bool output_constraint_snapshot = (0 == step % std::stoi(_config("IO_constraint_snapshot_interval", "999999999")));

  if(output_step && output_constraint_snapshot)
  {
      real_t *H_values, *M_values;
      H_values = new real_t[NX];
      M_values = new real_t[NX];
      
      bssnSim->set1DConstraintOutput(
        H_values, M_values, 1, 0, 0);


      std::string hfile = "1D_hamiltonian_constraints";
      std::string hfilename = iodata->dir() + hfile + ".strip.dat.gz";
      char hdata[35];

      gzFile hdatafile = gzopen(hfilename.c_str(), "ab");
      if(hdatafile == Z_NULL) {
        iodata->log("Error opening file: " + hfilename);
        return;
      }

      std::string mfile = "1D_momentum_constraints";
      std::string mfilename = iodata->dir() + mfile + ".strip.dat.gz";
      char mdata[35];

      gzFile mdatafile = gzopen(mfilename.c_str(), "ab");
      if(mdatafile == Z_NULL) {
        iodata->log("Error opening file: " + mfilename);
        return;
      }

      for(idx_t i=0; i<NX; i++)
      {
        sprintf(hdata, "%.15g\t", (double) H_values[i]);
        gzwrite(hdatafile, hdata, strlen(hdata));

        sprintf(mdata, "%.15g\t", (double) M_values[i]);
        gzwrite(mdatafile, mdata, strlen(mdata));

      }
      gzwrite(hdatafile, "\n", strlen("\n"));
      gzwrite(mdatafile, "\n", strlen("\n"));

      gzclose(hdatafile);
      gzclose(mdatafile);

  }
  
  if( output_step && output_this_step )
  {
    real_t H_calcs[8] = {0}, M_calcs[8] = {0}, G_calcs[7] = {0},
           A_calcs[7] = {0}, S_calcs[7] = {0};

    // Constraint Violation Calculations
    bssnSim->setConstraintCalcs(H_calcs, M_calcs, G_calcs,
                                A_calcs, S_calcs);

    io_dump_value(iodata, H_calcs[4], "H_violations", "\t"); // mean(H/[H])
    io_dump_value(iodata, H_calcs[5], "H_violations", "\t"); // stdev(H/[H])
    io_dump_value(iodata, H_calcs[6], "H_violations", "\t"); // max(H/[H])
    io_dump_value(iodata, H_calcs[2], "H_violations", "\t"); // max(H)
    io_dump_value(iodata, H_calcs[7], "H_violations", "\n"); // H_L2_norm

    io_dump_value(iodata, M_calcs[4], "M_violations", "\t"); // mean(M/[M])
    io_dump_value(iodata, M_calcs[5], "M_violations", "\t"); // stdev(M/[M])
    io_dump_value(iodata, M_calcs[6], "M_violations", "\t"); // max(M/[M])
    io_dump_value(iodata, M_calcs[2], "M_violations", "\t"); // max(M)
    io_dump_value(iodata, M_calcs[7], "M_violations", "\n"); // M_L2_norm

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
  real_t H_calcs[8] = {0}, M_calcs[8] = {0}, G_calcs[7] = {0},
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
  if(!output_step) return;

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
  
# pragma omp parallel for private(tmp_rd)
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

#   pragma omp atomic
    total_E += tmp_rd.E;
#   pragma omp atomic
    total_Phi += tmp_rd.Phi;
#   pragma omp atomic
    total_ell += tmp_rd.ell;
#   pragma omp atomic
    total_ellrho += tmp_rd.ell*tmp_rd.rho;
#   pragma omp atomic
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

void io_sheets_snapshot(IOData *iodata, idx_t step, Sheet * sheets)
{
  std::string step_str = std::to_string(step);
  bool output_step = false;
  bool output_this_step = false;

  output_step = ( std::stoi(_config("IO_3D_grid_interval", "0")) > 0 );
  output_this_step = ( 0 == step % std::stoi(_config("IO_3D_grid_interval", "1")) );
  bool output_Dx = (std::stoi(_config("IO_sheets_displacement_x", "0")) > 0);
  bool output_Dy = (std::stoi(_config("IO_sheets_displacement_y", "0")) > 0);
  bool output_Dz = (std::stoi(_config("IO_sheets_displacement_z", "0")) > 0);

  bool output_vx = (std::stoi(_config("IO_sheets_velocity_x", "0")) > 0);
  bool output_vy = (std::stoi(_config("IO_sheets_velocity_y", "0")) > 0);
  bool output_vz = (std::stoi(_config("IO_sheets_velocity_z", "0")) > 0);

  
  if( output_step && output_this_step )
  {
    if(output_Dx)
      io_dump_3dslice(iodata, sheets->Dx._array_a, "3D_sheets_Dx." + step_str);
    if(output_Dy)
      io_dump_3dslice(iodata, sheets->Dy._array_a, "3D_sheets_Dy." + step_str);
    if(output_Dz)
      io_dump_3dslice(iodata, sheets->Dz._array_a, "3D_sheets_Dz." + step_str);

    if(output_vx)
      io_dump_3dslice(iodata, sheets->vx._array_a, "3D_sheets_vx." + step_str);
    if(output_vy)
      io_dump_3dslice(iodata, sheets->vy._array_a, "3D_sheets_vy." + step_str);
    if(output_vz)
      io_dump_3dslice(iodata, sheets->vz._array_a, "3D_sheets_vz." + step_str);
  }
  
  output_step = ( std::stoi(_config("IO_2D_grid_interval", "0")) > 0 );
  output_this_step = ( 0 == step % std::stoi(_config("IO_2D_grid_interval", "1")) );
  if( output_step && output_this_step )
  {
    if(output_Dx)
      io_dump_2dslice(iodata, sheets->Dx._array_a, "2D_sheets_Dx." + step_str);
    if(output_Dy)
      io_dump_2dslice(iodata, sheets->Dy._array_a, "2D_sheets_Dy." + step_str);
    if(output_Dz)
      io_dump_2dslice(iodata, sheets->Dz._array_a, "2D_sheets_Dz." + step_str);

    if(output_vx)
      io_dump_2dslice(iodata, sheets->vx._array_a, "2D_sheets_vx." + step_str);
    if(output_vy)
      io_dump_2dslice(iodata, sheets->vy._array_a, "2D_sheets_vy." + step_str);
    if(output_vz)
      io_dump_2dslice(iodata, sheets->vz._array_a, "2D_sheets_vz." + step_str);

  }
  
  output_step = ( std::stoi(_config("IO_1D_grid_interval", "0")) > 0 );
  output_this_step = ( 0 == step % std::stoi(_config("IO_1D_grid_interval", "1")) );

  if( output_step && output_this_step )
  {
    int _axis = std::stoi(_config("axis", "1"));
    int _xoffset = std::stoi(_config("xoffset", "0"));
    int _yoffset = std::stoi(_config("yoffset", "0"));


    if(output_Dx)
      io_dump_strip(iodata, sheets->Dx._array_a,
                    "1D_sheets_Dx", _axis, _xoffset, _yoffset);
    if(output_Dy)
      io_dump_strip(iodata, sheets->Dy._array_a,
                    "1D_sheets_Dy", _axis, _xoffset, _yoffset);
    if(output_Dz)
      io_dump_strip(iodata, sheets->Dz._array_a,
                    "1D_sheets_Dz", _axis, _xoffset, _yoffset);

    if(output_vx)
      io_dump_strip(iodata, sheets->vx._array_a,
                    "1D_sheets_vx", _axis, _xoffset, _yoffset);
    if(output_vy)
      io_dump_strip(iodata, sheets->vy._array_a,
                    "1D_sheets_vy", _axis, _xoffset, _yoffset);
    if(output_vz)
      io_dump_strip(iodata, sheets->vz._array_a,
                    "1D_sheets_vz", _axis, _xoffset, _yoffset);

  }

  
}

/**
 * @brief      Read full 3D slice from a file.
 *
 * @param      iodata    initialized IOData
 * @param      field     Field to real
 * @param[in]  filename  filename to read (minus suffix)
 */
bool io_read_3dslice(IOData *iodata, arr_t & field, std::string filename)
{
  std::string dump_filename = filename + ".3d_grid.h5.gz";
  hid_t file_id = H5Fopen(dump_filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  if(file_id < 0) return false;

  hid_t temp;

  temp = H5Dopen(file_id, "DS1", H5P_DEFAULT);

  if( temp < 0 ) {
    std::cout << "Dataset DS1 does not exist." << std::endl;
    return false;
  }
  H5Dclose(temp);

  hid_t dset_id = H5Dopen(file_id, "DS1", H5P_DEFAULT);
  std::cout << "Dataset DS1 exists!" << std::endl << std::flush;

  herr_t status;

  // hsize_t dims[3] = {(hsize_t) field.nx, (hsize_t) field.ny, (hsize_t) field.nz},
  //   maxdims[3] = {H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED},
  //   chunk[3] = {6, 6, 6};

  status = H5Dread(dset_id, H5T_TO_USE, H5S_ALL, H5S_ALL, H5P_DEFAULT, field._array);
  status = status; // suppress unused variable warning
  
  return true;
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
  hsize_t     dims[3] = {(hsize_t) field.nx, (hsize_t) field.ny, (hsize_t) field.nz},
              maxdims[3] = {H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED},
              chunk[3] = {6, 6, 6};

  file = H5Fcreate (dump_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  space = H5Screate_simple (3, dims, maxdims);
  dcpl = H5Pcreate (H5P_DATASET_CREATE);
  status = H5Pset_deflate (dcpl, 9);
  status = H5Pset_chunk (dcpl, 3, chunk);
  dset = H5Dcreate2 (file, "DS1", H5T_ALLOC, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);

  status = H5Dwrite (dset, H5T_TO_USE, H5S_ALL, H5S_ALL, H5P_DEFAULT, field._array);

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
  hsize_t     dims[2] = {(hsize_t) field.ny, (hsize_t) field.nz},
              maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED},
              chunk[2] = {6, 6};

  file = H5Fcreate (dump_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  space = H5Screate_simple (2, dims, maxdims);
  dcpl = H5Pcreate (H5P_DATASET_CREATE);
  status = H5Pset_deflate (dcpl, 9);
  status = H5Pset_chunk (dcpl, 2, chunk);
  dset = H5Dcreate2 (file, "Dataset1", H5T_ALLOC, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);

  status = H5Dwrite (dset, H5T_TO_USE, H5S_ALL, H5S_ALL, H5P_DEFAULT, field._array);

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
      for(idx_t i=0; i<field.nx; i++)
      {
        sprintf(data, "%.15g\t", (double) field(i,n1,n2));
        gzwrite(datafile, data, strlen(data));
      }
      break;
    case 2:
      for(idx_t j=0; j<field.ny; j++)
      {
        sprintf(data, "%.15g\t", (double) field(n1,j,n2));
        gzwrite(datafile, data, strlen(data));
      }
      break;
    case 3:
      for(idx_t k=0; k<field.nz; k++)
      {
        sprintf(data, "%.15g\t", (double) field(n1,n2,k));
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
      for(idx_t i=0; i<field.nx; i++)
      {
        std::cout << std::setprecision(19) << field[INDEX(i,n1,n2)] << ", ";
      }
      break;
    case 2:
      for(idx_t j=0; j<field.ny; j++)
      {
        std::cout << std::setprecision(19) << field[INDEX(n1,j,n2)] << ", ";
      }
      break;
    case 3:
      for(idx_t k=0; k<field.nz; k++)
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
  dset = H5Dcreate2 (file, dataset_name.c_str(), H5T_ALLOC, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);

  status = H5Dwrite (dset, H5T_TO_USE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);

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
  if(!output_step) return;

  bool output_this_step = (0 == step % std::stoi(_config("IO_particles", "1")));
  bool output_phase_diagram = (0 == step % std::stoi(_config("IO_particles_diagram", "1")));

  bool output_x = (std::stoi(_config("IO_particles_x", "0")) > 0);
  bool output_y = (std::stoi(_config("IO_particles_y", "0")) > 0);
  bool output_z = (std::stoi(_config("IO_particles_z", "0")) > 0);

  bool output_vx = (std::stoi(_config("IO_particles_vx", "0")) > 0);
  bool output_vy = (std::stoi(_config("IO_particles_vy", "0")) > 0);
  bool output_vz = (std::stoi(_config("IO_particles_vz", "0")) > 0);

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

    particle_vec * p_vec = particles->getParticleVec();
    for(particle_vec::iterator it = p_vec->begin(); it != p_vec->end(); ++it) {
      sprintf(data, "%.15g\t", (double) it->p_a.X[0]);
      gzwrite(datafile, data, strlen(data));
      sprintf(data, "%.15g\t", (double) it->p_a.X[1]);
      gzwrite(datafile, data, strlen(data));
      sprintf(data, "%.15g\t", (double) it->p_a.X[2]);
      gzwrite(datafile, data, strlen(data));
    }

    sprintf(data, "\n");
    gzwrite(datafile, data, strlen(data));
    gzclose(datafile);
  }

  if( output_step  && output_phase_diagram)
  {
    const real_t TOL = 0.01;
    std::vector<double> x_cache, y_cache, z_cache, vx_cache, vy_cache, vz_cache;
    particle_vec * p_vec = particles->getParticleVec();
    for(particle_vec::iterator it = p_vec->begin(); it != p_vec->end(); ++it)
    {
      if(output_x && std::fabs(pw2(it->p_a.X[1]) + pw2(it->p_a.X[2])) < TOL)
        x_cache.push_back(it->p_a.X[0]);
      if(output_y && std::fabs(pw2(it->p_a.X[0]) + pw2(it->p_a.X[2])) < TOL)
        y_cache.push_back(it->p_a.X[1]);
      if(output_z && std::fabs(pw2(it->p_a.X[1]) + pw2(it->p_a.X[0])) < TOL)
        z_cache.push_back(it->p_a.X[2]);

      if(output_vx && std::fabs(pw2(it->p_a.X[1]) + pw2(it->p_a.X[2])) < TOL)
        vx_cache.push_back(it->p_a.U[0]);
      if(output_vy && std::fabs(pw2(it->p_a.X[0]) + pw2(it->p_a.X[2])) < TOL)
        vy_cache.push_back(it->p_a.U[1]);
      if(output_vz && std::fabs(pw2(it->p_a.X[1]) + pw2(it->p_a.X[0])) < TOL)
        vz_cache.push_back(it->p_a.U[2]);
    }

    if(output_x)
    {
      std::string filename = iodata->dir() + "particle_x.dat.gz";
      char data[35];

      gzFile datafile = gzopen(filename.c_str(), "ab");
      if(datafile == Z_NULL) {
        iodata->log("Error opening file: " + filename);
        return;
      }

      for(uint i=0; i<x_cache.size(); i++)
      {
        sprintf(data, "%.15g\t", x_cache[i]);
        gzwrite(datafile, data, strlen(data));
      }
      gzwrite(datafile, "\n", strlen("\n"));

      gzclose(datafile);
    }

    if(output_y)
    {
      std::string filename = iodata->dir() + "particle_y.dat.gz";
      char data[35];

      gzFile datafile = gzopen(filename.c_str(), "ab");
      if(datafile == Z_NULL) {
        iodata->log("Error opening file: " + filename);
        return;
      }

      for(uint i=0; i<y_cache.size(); i++)
      {
        sprintf(data, "%.15g\t", y_cache[i]);
        gzwrite(datafile, data, strlen(data));
      }
      gzwrite(datafile, "\n", strlen("\n"));

      gzclose(datafile);
    }

    if(output_z)
    {
      std::string filename = iodata->dir() + "particle_z.dat.gz";
      char data[35];

      gzFile datafile = gzopen(filename.c_str(), "ab");
      if(datafile == Z_NULL) {
        iodata->log("Error opening file: " + filename);
        return;
      }

      for(uint i=0; i<z_cache.size(); i++)
      {
        sprintf(data, "%.15g\t", z_cache[i]);
        gzwrite(datafile, data, strlen(data));
      }
      gzwrite(datafile, "\n", strlen("\n"));

      gzclose(datafile);

    }

    if(output_vx)
    {
      std::string filename = iodata->dir() + "particle_vx.dat.gz";
      char data[35];

      gzFile datafile = gzopen(filename.c_str(), "ab");
      if(datafile == Z_NULL) {
        iodata->log("Error opening file: " + filename);
        return;
      }

      for(uint i=0; i<vx_cache.size(); i++)
      {
        sprintf(data, "%.15g\t", vx_cache[i]);
        gzwrite(datafile, data, strlen(data));
      }
      gzwrite(datafile, "\n", strlen("\n"));

      gzclose(datafile);

    }

    if(output_vy)
    {
      std::string filename = iodata->dir() + "particle_vy.dat.gz";
      char data[35];

      gzFile datafile = gzopen(filename.c_str(), "ab");
      if(datafile == Z_NULL) {
        iodata->log("Error opening file: " + filename);
        return;
      }

      for(uint i=0; i<vy_cache.size(); i++)
      {
        sprintf(data, "%.15g\t", vy_cache[i]);
        gzwrite(datafile, data, strlen(data));
      }
      gzwrite(datafile, "\n", strlen("\n"));

      gzclose(datafile);

    }

    if(output_vz)
    {
      std::string filename = iodata->dir() + "particle_vz.dat.gz";
      char data[35];

      gzFile datafile = gzopen(filename.c_str(), "ab");
      if(datafile == Z_NULL) {
        iodata->log("Error opening file: " + filename);
        return;
      }

      for(uint i=0; i<vz_cache.size(); i++)
      {
        sprintf(data, "%.15g\t", vz_cache[i]);
        gzwrite(datafile, data, strlen(data));
      }
      gzwrite(datafile, "\n", strlen("\n"));

      gzclose(datafile);

    } 
  }
}

#if USE_COSMOTRACE
void io_raytrace_bardeen_dump(IOData *iodata, idx_t step,
  std::vector<RayTrace<real_t, idx_t> *> const * rays, Bardeen * bardeen,
  real_t t)
{
  bardeen->setPotentials(t);

  idx_t num_rays = rays->size();
  
  real_t * Phis, * Psis, * dt_Bs;
  Phis = new real_t [num_rays];
  Psis = new real_t [num_rays];
  dt_Bs = new real_t [num_rays];

# pragma omp parallel for
  for(idx_t n=0; n<num_rays; n++)
  {
    RaytraceData<real_t> tmp_rd = {0};
    tmp_rd = (*rays)[n]->getRaytraceData();
    Phis[n] = interp(tmp_rd.x[0]/dx, tmp_rd.x[1]/dx, tmp_rd.x[2]/dx,
      NX, NY, NZ, bardeen->Phi);
    Psis[n] = interp(tmp_rd.x[0]/dx, tmp_rd.x[1]/dx, tmp_rd.x[2]/dx,
      NX, NY, NZ, bardeen->Psi);
    dt_Bs[n] = interp(tmp_rd.x[0]/dx, tmp_rd.x[1]/dx, tmp_rd.x[2]/dx,
      NX, NY, NZ, bardeen->dt_B);
  }

  // Write data from individual rays to file
  std::string dataset_name = "step_" + std::to_string(step);
  std::string file_name = "bardeen_phi_data";
  io_dump_2d_array(iodata, Phis, num_rays, 1,
    file_name, dataset_name);
  file_name = "bardeen_psi_data";
  io_dump_2d_array(iodata, Psis, num_rays, 1,
    file_name, dataset_name);
  file_name = "bardeen_dt_B_data";
  io_dump_2d_array(iodata, dt_Bs, num_rays, 1,
    file_name, dataset_name);

  delete[] Phis;
  delete[] Psis;
  delete[] dt_Bs;

  return;
}
#endif


void io_svt_violation(IOData *iodata, idx_t step, Bardeen * bardeen, real_t t)
{
  // potentials should be set per sim call to prepBSSNOutput
  bool output_step = ( std::stoi(_config("SVT_constraint_interval", "0")) > 0 );
  if(!output_step) return;

  bool output_this_step = (0 == step % std::stoi(_config("SVT_constraint_interval", "1")));
  if( output_step && output_this_step )
  {
    real_t SVT_calcs[NUM_BARDEEN_VIOLS] = {0};
    bardeen->getSVTViolations(SVT_calcs);

    for(int i=0; i<NUM_BARDEEN_VIOLS-1; ++i)
      io_dump_value(iodata, SVT_calcs[i], "SVT_violations", "\t");
    io_dump_value(iodata, SVT_calcs[NUM_BARDEEN_VIOLS-1], "SVT_violations", "\n");
  }
}


void io_raysheet_dump(IOData *iodata, idx_t step,
  Sheet * raySheet, BSSN *bssnSim, Lambda * lambda)
{
  bool output_step = ( std::stoi(_config("IO_raysheet_interval", "0")) > 0 );
  if(!output_step) return;

  bool output_this_step = (0 == step % std::stoi(_config("IO_raysheet_interval", "1")));
  bool output_minimalwrite = !!std::stoi(_config("IO_raysheet_minimalwrite", "1"));
  if( output_step && output_this_step )
  {
    // output misc. info about simulation here.
    char data[35];
    std::string dump_filename = iodata->dir() + "raysheet.dat.gz";

    gzFile datafile = gzopen(dump_filename.c_str(), "ab");
    if(datafile == Z_NULL) {
      iodata->log("Error opening file: " + dump_filename);
      return;
    }

    for(idx_t r=0; r<raySheet->ns1; ++r)
    {
      std::vector<real_t> sheet_data = raySheet->getRayDataAtS(r,bssnSim,lambda);
      idx_t size = (idx_t) sheet_data.size();
      idx_t i;
      for(i=output_minimalwrite?6:0; i<size; ++i)
      {
        sprintf(data, "%.15g\t", (double) sheet_data[i]);
        gzwrite(datafile, data, strlen(data));        
      }
    }

    sprintf(data, "\n");
    gzwrite(datafile, data, strlen(data));
    gzclose(datafile);
  }

  return;
}


} /* namespace cosmo */
