
#include "io.h"
#include <hdf5.h>

namespace cosmo
{

void io_dump_strip(real_t *field, int axis, idx_t n1, idx_t n2)
{
  char filename[] = "strip.dat.gz";
  char data[20];

  gzFile datafile = gzopen(filename, "ab");
  if(datafile == Z_NULL) {
    printf("Error opening file: %s\n", filename);
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


void io_dump_2dslice(real_t *field, std::string filename)
{
  // dump the first N*N points (a 2-d slice on a boundary)
  std::string dump_filename = filename + ".grid.h5.gz";

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

  return;
}


void io_dump_quantities(std::map <std::string, real_t *> & bssn_fields,
                      std::map <std::string, real_t *> & hydro_fields,
                      std::string filename)
{
  // output misc. info about simulation here.
  char data[20];
  std::string dump_filename = filename + ".dat.gz";

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

  // FRW Piece of Hamiltonian constraint
  idx_t i, j, k;
  real_t sum = 0.0;
  LOOP3(i, j, k)
  {
    // 4/3*A^2 - K^2 + 24 pi e^(5phi)*\rho
    sum += 4.0/3.0*(
        0.0 // need A_ij^2
      ) - pw2(bssn_fields["K_a"][NP_INDEX(i,j,k)])
      + 24.0*PI*bssn_fields["r_a"][NP_INDEX(i,j,k)]
    ;
  }
  // std::cout << "FRW Const. Resid. :" << sum << "\n";
  // std::cout << "Aij contrib. Resid. :" << average(bssn_fields["A11_a"])*POINTS << "\n";

  gzwrite(datafile, "\n", strlen("\n"));
  gzclose(datafile);
}


}
