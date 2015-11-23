#ifndef PERIODARRAY_H
#define PERIODARRAY_H

namespace cosmo
{

template<typename IT, typename RT>
class periodicArray
{
  public:
    IT nx, ny, nz;
    IT pts;

    RT* _array;

    periodicArray()
    {
      // call init() to initialize data
    }

    ~periodicArray()
    {
      delete [] _array;
    }

    void init(IT nx_in, IT ny_in, IT nz_in)
    {
      nx = nx_in;
      ny = ny_in;
      nz = nz_in;

      pts = nx*ny*nz;

      _array = new RT[pts];
    }

    RT& operator()(IT i, IT j, IT k)
    {
      // indexing only works down to negative nx, ny, nz
      IT idx = ( ((i+nx)%(nx))*(ny)*(nz) + ((j+ny)%(ny))*(nz) + (k+nz)%(nz) );
      return _array[idx];
    }

    RT& operator[](IT idx)
    {
      return _array[idx];
    }
};

} // namespace cosmo

#endif
