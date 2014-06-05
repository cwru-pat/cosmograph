
#include "io.h"

namespace cosmo
{

void dump_strip(real_t *field, int axis, idx_t n1, idx_t n2)
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

}
