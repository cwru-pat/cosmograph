#ifndef COSMO_UTILS_RK4REGISTER_H
#define COSMO_UTILS_RK4REGISTER_H

#include <string>
#include <utility>
#include "Array.h"

namespace cosmo
{

/**
 * @brief RK4 Class for integration
 * @details See the docs/RK4_integration.pptx file.
 * 
 * @tparam IT Index type
 * @tparam RT Real type
 */
template<typename IT, typename RT>
class RK4Register
{
  private:
    std::string name;
    IT points;
    RT sim_dt;

  public:
    CosmoArray<IT, RT> _array_p; ///< "_p" register: contains data from _p_revious step
    CosmoArray<IT, RT> _array_a; ///< "_a" register: containes _a_ctive data needed for _c_omputations
    CosmoArray<IT, RT> _array_c; ///< "_c" register: contains _c_omputed values
    CosmoArray<IT, RT> _array_f; ///< "_f" register: containes final value of RK4 step

    /**
     * @brief Constructor calls CosmoArray constructors; call RK4Register::init
     * to initialize class
     */
    RK4Register() :
      _array_p(), _array_a(), _array_c(), _array_f()
    {
      // call init()
    }

 RK4Register(IT nx_in, IT ny_in, IT nz_in, RT sim_dt_in):
      _array_p(nx_in, ny_in, nz_in),
        _array_a(nx_in, ny_in, nz_in),
        _array_c(nx_in, ny_in, nz_in),
        _array_f(nx_in, ny_in, nz_in)
        {
          setDt(sim_dt_in);

          points = nx_in*ny_in*nz_in;

          // call init()
        }

      
    RK4Register(IT nx_in, IT ny_in, IT nz_in, RT lx_in, RT ly_in, RT lz_in, RT sim_dt_in):
      _array_p(nx_in, ny_in, nz_in, lx_in, ly_in, lz_in),
        _array_a(nx_in, ny_in, nz_in, lx_in, ly_in, lz_in),
        _array_c(nx_in, ny_in, nz_in, lx_in, ly_in, lz_in),
        _array_f(nx_in, ny_in, nz_in, lx_in, ly_in, lz_in)
      {
        setDt(sim_dt_in);

        points = nx_in*ny_in*nz_in;
      

        // call init()
      }

      
    /**
     * @brief Initialize class variables; call CosmoArray::init for array members
     * @details Set "dt" for class instance; grid dimensions
     * 
     * @param nx_in num. grid points in x-direction
     * @param ny_in num. grid points in y-direction
     * @param nz_in num. grid points in z-direction
     * @param sim_dt_in initial timestep
     */
    void init(IT nx_in, IT ny_in, IT nz_in, RT sim_dt_in)
    {
      setDt(sim_dt_in);

      points = nx_in*ny_in*nz_in;
      
      _array_p.init(nx_in, ny_in, nz_in);
      _array_a.init(nx_in, ny_in, nz_in);
      _array_c.init(nx_in, ny_in, nz_in);
      _array_f.init(nx_in, ny_in, nz_in);
    }

    /**
     * @brief Set "dt" for RK4Register instance
     */
    void setDt(RT sim_dt_in)
    {
      sim_dt = sim_dt_in;
    }

    ~RK4Register() {}

    /**
     * @brief Set "name" property of instance, CosmoArray member instances
     * 
     * @param name_in name
     */
    void setName(std::string name_in)
    {
      name = name_in;
      _array_p.setName(name_in + "_p");
      _array_a.setName(name_in + "_a");
      _array_c.setName(name_in + "_c");
      _array_f.setName(name_in + "_f");
    }

    /**
     * @brief Swap data in and names of _a and _c registers
     */
    void swap_a_c()
    {
      std::swap(_array_a.name, _array_c.name);
      std::swap(_array_a._array, _array_c._array);
    }

    /**
     * @brief Swap data in and names of _p and _f registers
     */
    void swap_p_f()
    {
      std::swap(_array_p.name, _array_f.name);
      std::swap(_array_p._array, _array_f._array);
    }

    void stepInit()
    {
      IT i;
      #pragma omp parallel for default(shared) private(i)
      for(i=0; i<points; ++i)
      {
        _array_a[i] = _array_p[i];
        _array_f[i] = 0;
      }
    }

    void K1Finalize()
    {
      #pragma omp parallel for
      for(IT i=0; i<points; ++i)
      {
        _array_f[i] += sim_dt*_array_c[i]/6.0;
        _array_c[i] = _array_p[i] + sim_dt*_array_c[i]/2.0;
      }

      swap_a_c();
    }

    void K2Finalize()
    {
      #pragma omp parallel for
      for(IT i=0; i<points; ++i)
      {
        _array_f[i] += sim_dt*_array_c[i]/3.0;
        _array_c[i] = _array_p[i] + sim_dt*_array_c[i]/2.0;
      }
      
      swap_a_c();
    }

    void K3Finalize()
    {
      #pragma omp parallel for
      for(IT i=0; i<points; ++i)
      {
        _array_f[i] += sim_dt*_array_c[i]/3.0;
        _array_c[i] = _array_p[i] + sim_dt*_array_c[i];
      }
      
      swap_a_c();
    }

    void K4Finalize()
    {
      #pragma omp parallel for
      for(IT i=0; i<points; ++i)
      {
        _array_f[i] += sim_dt*_array_c[i]/6.0 + _array_p[i];
        _array_p[i] = _array_f[i];
      }

      swap_a_c();
    }

    RT& _p(const IT & i, const IT & j, const IT & k) { return _array_p(i, j, k); }
    RT& _a(const IT & i, const IT & j, const IT & k) { return _array_a(i, j, k); }
    RT& _c(const IT & i, const IT & j, const IT & k) { return _array_c(i, j, k); }
    RT& _f(const IT & i, const IT & j, const IT & k) { return _array_f(i, j, k); }

    RT& operator()(const IT & i, const IT & j, const IT & k)
    {
      return _array_a(i, j, k);
      return _array_a[_array_a.idx(i, j, k)];
    }

    RT& operator[](IT idx)
    {
      return _array_a[idx];
    }

    
};

} // end namespace

#endif
