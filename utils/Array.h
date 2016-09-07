#ifndef COSMO_UTILS_ARRAY_H
#define COSMO_UTILS_ARRAY_H

#include <string>
#include <utility>

namespace cosmo
{

template<typename IT, typename RT>
class CosmoArray
{
  public:
    IT nx, ny, nz;
    IT pts;

    std::string name;

    RT* _array;

    CosmoArray() {}

    CosmoArray(IT n_in)
    {
      init(n_in, n_in, n_in);
    }

    CosmoArray(IT nx_in, IT ny_in, IT nz_in)
    {
      init(nx_in, ny_in, nz_in);
    }

    ~CosmoArray()
    {
      delete [] _array;
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

      int status = posix_memalign((void **) &_array, 64, pts * sizeof(RT));
      if(status != 0)
      {
        throw -1;
      }

      #pragma omp parallel for
      for(IT i=0; i<pts; ++i)
      {
        _array[i] = 0.0;
      }
    }

    RT sum()
    {
      RT res = 0.0;
      #pragma omp parallel for reduction(+:res)
      for(IT i=0; i<pts; ++i)
      {
        res += _array[i];
      }
      return res;

    }
    
    RT avg()
    {
      return sum() / (RT)pts;
    }

    RT min()
    {
      RT min_res = 1e100;
      #pragma omp parallel for 
      for(IT i = 0; i<pts; i++)
      {
        #pragma omp critical
        if(_array[i] < min_res)
          min_res = _array[i];
      }
      return min_res;
    }

    RT max()
    {
      RT max_res = -1e100;
      #pragma omp parallel for
      for(IT i = 0; i<pts; i++)
      {
        #pragma omp critical
        {
          if(_array[i] > max_res)
            max_res = _array[i];
        }
      }

      return max_res;

    }
    IT idx(IT i_in, IT j_in, IT k_in)
    {
      IT i=i_in, j=j_in, k=k_in;

      // indexing only works down to negative 100*(nx, ny, nz)?
      // Using this is slow.
      if(i_in < 0 || i_in >= nx) i = (i_in+100*nx)%nx;
      if(j_in < 0 || j_in >= ny) j = (j_in+100*ny)%ny;
      if(k_in < 0 || k_in >= nz) k = (k_in+100*nz)%nz;
      return ( i*ny*nz + j*nz + k );
    }

    RT& operator()(IT i, IT j, IT k)
    {
      IT x = idx(i, j, k);
      return _array[x];
    }

    CosmoArray& operator=(const CosmoArray& other) {
      // check for self-assignment
      if(&other == this)
        return *this;

      this->nx = other.nx;
      this->nz = other.ny;
      this->ny = other.nz;
      this->pts = other.pts;

      this->name = other.name;

      #pragma omp parallel for
      for(IT i=0; i<pts; ++i)
      {
        this->_array[i] = other._array[i];
      }

      return *this;
    }

    RT& operator[](IT idx)
    {
      return _array[idx];
    }
};

template <class T>
void cosmoArraySwap(T & arr1, T & arr2)
{
  std::swap(arr1.nx, arr2.nx);
  std::swap(arr1.ny, arr2.ny);
  std::swap(arr1.nz, arr2.nz);
  std::swap(arr1.pts, arr2.pts);
  std::swap(arr1.name, arr2.name);
  std::swap(arr1._array, arr2._array);
}

}

#endif
