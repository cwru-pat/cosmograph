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

    std::string name;

    RT* _array;

    periodicArray()
    {
      // call init() to initialize data
    }

    ~periodicArray()
    {
      deleteDev();
      delete [] _array;
    }

    void createDev()
    {
      #pragma acc enter data copyin( this[0:1], _array[0:pts] )
      #pragma acc update device(nx, ny, nz, pts, name)
    }

    void deleteDev()
    {
      #pragma acc exit data delete( _array[0:pts], this[0:1] )
    }

    void updateDev()
    {
      #pragma acc update device(_array[0:pts])
    }

    void updateHost()
    {
      #pragma acc update self(_array[0:pts])
    }


    void setName(std::string name_in)
    {
      name = name_in;
    }

    void init(IT nx_in, IT ny_in, IT nz_in)
    {
      nx = nx_in;
      ny = ny_in;
      nz = nz_in;

      pts = nx*ny*nz;

      _array = new RT[pts];
      IT i;
      for(i=0; i<pts; ++i)
      {
        _array[i] = 0.0;
      }

      createDev();
    }

    IT idx(IT i_in, IT j_in, IT k_in)
    {
      IT i=i_in, j=j_in, k=k_in;

      // indexing only works down to negative nx, ny, nz
      if(i_in < 0 || i_in >= nx) i = (i_in+100*nx)%nx;
      if(j_in < 0 || j_in >= ny) j = (j_in+100*ny)%ny;
      if(k_in < 0 || k_in >= nz) k = (k_in+100*nz)%nz;
      return ( i*ny*nz + j*nz + k );
    }

    RT& operator()(IT i_in, IT j_in, IT k_in)
    {
      IT i=i_in, j=j_in, k=k_in;

      // indexing only works down to negative nx, ny, nz
      if(i_in < 0 || i_in >= nx) i = (i_in+100*nx)%nx;
      if(j_in < 0 || j_in >= ny) j = (j_in+100*ny)%ny;
      if(k_in < 0 || k_in >= nz) k = (k_in+100*nz)%nz;
      IT idx = ( i*ny*nz + j*nz + k );

      return _array[idx];
    }

    periodicArray& operator=(const periodicArray& other) {
      // check for self-assignment
      if(&other == this)
        return *this;

      this->nx = other.nx;
      this->nz = other.ny;
      this->ny = other.nz;
      this->pts = other.pts;

      this->name = other.name;

      this->_array = other._array;

      return *this;
    }

    RT& operator[](IT idx)
    {
      return _array[idx];
    }
};

} // namespace cosmo

#endif
