#ifndef COSMO_UTILS_ARRAY_H
#define COSMO_UTILS_ARRAY_H

#include <string>
#include <utility>
#include <iostream>

#include "TriCubicInterpolator.h"

namespace cosmo
{

template<typename IT, typename RT>
class CosmoArray
{
  public:
    IT nx, ny, nz;
    IT pts = 0;

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
      if(pts > 0)
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

      _array = new RT[pts];

#pragma omp parallel for
      for(IT i=0; i<pts; ++i)
      {
        _array[i] = 0.0;
      }
    }
    
    IT _IT_mod(IT n, IT d) const
    {
      IT mod = n % d;
      if(mod < 0)
        mod += d;
      return mod;
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
      // Using this is slow. Use a macro instead.
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

    // Weighted averaging / trilinear interpolation via
    // https://en.wikipedia.org/wiki/Trilinear_interpolation#Method
    RT getInterpolatedValue(RT i_in, RT j_in, RT k_in)
    {
      IT il = i_in < 0 ? (IT) i_in - 1 : (IT) i_in; // Index "left" of i
      RT id = i_in - il; // fractional difference
      IT jl = j_in < 0 ? (IT) j_in - 1 : (IT) j_in; // same as ^ but j
      RT jd = j_in - jl;
      IT kl = k_in < 0 ? (IT) k_in - 1 : (IT) k_in; // same as ^ but k
      RT kd = k_in - kl;

      RT c00 = _array[idx(il, jl, kl)]*(1-id) + _array[idx(il+1, jl, kl)]*id;
      RT c01 = _array[idx(il, jl, kl+1)]*(1-id) + _array[idx(il+1, jl, kl+1)]*id;
      RT c10 = _array[idx(il, jl+1, kl)]*(1-id) + _array[idx(il+1, jl+1, kl)]*id;
      RT c11 = _array[idx(il, jl+1, kl+1)]*(1-id) + _array[idx(il+1, jl+1, kl+1)]*id;
      RT c0 = c00*(1-jd) + c10*jd;
      RT c1 = c01*(1-jd) + c11*jd;

      return c0*(1-kd) + c1*kd;
    }

    // Catmull-Rom cubic spline
    RT CINT(RT u, RT p0, RT p1, RT p2, RT p3)
    {
      return 0.5*(
            (u*u*(2.0 - u) - u)*p0
          + (u*u*(3.0*u - 5.0) + 2)*p1
          + (u*u*(4.0 - 3.0*u) + u)*p2
          + u*u*(u - 1.0)*p3
        );
    }

    RT getTriCubicInterpolatedValue(RT i_in, RT j_in, RT k_in)
    {
      // TODO: need higher order interpolation methods
      // For now, weighted average (linear interpolation)
      IT il = i_in < 0 ? (IT) i_in - 1 : (IT) i_in; // Index "left" of i
      RT id = i_in - il; // fractional difference
      IT jl = j_in < 0 ? (IT) j_in - 1 : (IT) j_in; // same as ^ but j
      RT jd = j_in - jl;
      IT kl = k_in < 0 ? (IT) k_in - 1 : (IT) k_in; // same as ^ but k
      RT kd = k_in - kl;

      // interpolated value at (i*, j*, k_in)
      RT * F_i_j_kd = new RT[16];
      for(IT i=0; i<4; ++i)
        for(IT j=0; j<4; ++j)
          F_i_j_kd[i*4+j] = CINT(kd,
            _array[idx(il+i-1, jl+j-1, kl-1)], _array[idx(il+i-1, jl+j-1, kl+0)],
            _array[idx(il+i-1, jl+j-1, kl+1)], _array[idx(il+i-1, jl+j-1, kl+2)]);

      // interpolated value at (i*, j_in, k_in)
      RT * F_i_jd_kd = new RT[4];
      for(IT i=0; i<4; ++i)
        F_i_jd_kd[i] = CINT(jd, F_i_j_kd[i*4+0], F_i_j_kd[i*4+1], F_i_j_kd[i*4+2], F_i_j_kd[i*4+3]);

      // interpolated value at (i_in, j_in, k_in)
      RT Fijk = CINT(id, F_i_jd_kd[0], F_i_jd_kd[1], F_i_jd_kd[2], F_i_jd_kd[3]);

      delete [] F_i_j_kd;
      delete [] F_i_jd_kd;

      return Fijk;
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
