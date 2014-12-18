
#include "io.h"

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

  gzwrite(datafile, "\n", strlen("\n"));
  gzclose(datafile);
}


}
