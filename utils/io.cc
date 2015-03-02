
#include "io.h"
#include <hdf5.h>
#include <sys/stat.h>

namespace cosmo
{

void io_show_progress(idx_t s, idx_t maxs)
{
  if(s > maxs)
  {
    std::cout << "Simulation has ran more steps than allowed.";
    throw -1;
    return;
  }

  idx_t s_digits = (int) log10 ((double) s) + 1;
  if(s==0)
  {
    s_digits = 1;
  }
  idx_t maxs_digits = (int) log10 ((double) maxs) + 1;

  std::cout << " Running step " << s;
  for(int i=s_digits; i<maxs_digits; ++i)
  {
    std::cout << " ";
  }
  std::cout << " / " << maxs;

  idx_t ndots = 20;
  std::cout << "[";
  for(int i=0; i<ndots; ++i)
  {
    if(i < (ndots*s)/maxs)
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

void io_init(IOData *iodata)
{
  /* ensure data_dir ends with '/', unless empty string is specified. */
  size_t len_dir_name = iodata->output_dir.length();
  if((iodata->output_dir)[len_dir_name - 1] != '/' && len_dir_name != 0)
  {
    iodata->output_dir += '/';
  }
  /* create data_dir */
  if(len_dir_name != 0)
    mkdir(iodata->output_dir.c_str(), 0755);
}

void io_dump_strip(real_t *field, int axis, idx_t n1, idx_t n2, IOData *iodata)
{
  std::string filename = iodata->output_dir + "strip.dat.gz";
  char data[20];

  gzFile datafile = gzopen(filename.c_str(), "ab");
  if(datafile == Z_NULL) {
    printf("Error opening file: %s\n", filename.c_str());
    return;
  }

  switch (axis)
  {
    case 1:
      for(idx_t i=0; i<N; i++)
      {
        sprintf(data, "%g\t", field[INDEX(i,n1,n2)]);
        gzwrite(datafile, data, strlen(data));
      }
      break;
    case 2:
      for(idx_t j=0; j<N; j++)
      {
        sprintf(data, "%g\t", field[INDEX(n1,j,n2)]);
        gzwrite(datafile, data, strlen(data));
      }
      break;
    case 3:
      for(idx_t k=0; k<N; k++)
      {
        sprintf(data, "%g\t", field[INDEX(n1,n2,k)]);
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
  // dump the first N*N points (a 2-d slice on a boundary)
  std::string dump_filename = iodata->output_dir + filename + ".2d_grid.h5.gz";

  hid_t       file, space, dset, dcpl;  /* Handles */
  herr_t      status;
  htri_t      avail;
  H5Z_filter_t  filter_type;
  hsize_t     dims[2] = {N, N},
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
  // dump all N*N*N points
  std::string dump_filename = iodata->output_dir + filename + ".3d_grid.h5.gz";

  hid_t       file, space, dset, dcpl;  /* Handles */
  herr_t      status;
  htri_t      avail;
  H5Z_filter_t  filter_type;
  hsize_t     dims[3] = {N, N, N},
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


void io_dump_quantities(std::map <std::string, real_t *> & bssn_fields,
                      std::map <std::string, real_t *> & hydro_fields,
                      std::string filename, IOData *iodata)
{
  // output misc. info about simulation here.
  char data[20];
  std::string dump_filename = iodata->output_dir + filename + ".dat.gz";

  gzFile datafile = gzopen(dump_filename.c_str(), "ab");
  if(datafile == Z_NULL) {
    std::cout << "Error opening file: " << dump_filename << "\n";
    return;
  }

  // average phi
  sprintf(data, "%g\t", average(bssn_fields["phi_a"]));
  gzwrite(datafile, data, strlen(data));

  // average K
  sprintf(data, "%g\t", average(bssn_fields["K_a"]));
  gzwrite(datafile, data, strlen(data));

  // average UD
  sprintf(data, "%g\t", average(hydro_fields["UD_a"]));
  gzwrite(datafile, data, strlen(data));

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
    std::cout << "Error opening file: " << dump_filename << "\n";
    return;
  }

  sprintf(data, "%g", value);
  gzwrite(datafile, data, strlen(data));

  gzwrite(datafile, "\n", strlen("\n"));
  gzclose(datafile);
}


} /* namespace cosmo */
