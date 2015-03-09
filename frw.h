#ifndef COSMO_FRW
#define COSMO_FRW

#define PI_L 3.141592653589793238462643383279502884L
#include <list>

namespace cosmo
{

/** FRW class **/
template<typename RT>
class FRW
{
  // metric variables in FRW
  RT phi;
  RT K;
  // densities and corresponding EOS "w"s
  std::list< std::pair<RT,RT> > fluids;

public:
  FRW(RT phi_in, RT K_in)
  {
    phi = phi_in;
    K = K_in;
  }

  void addFluid(RT rho_in, RT w_in)
  {
    fluids.emplace_back(rho_in, w_in);
  }
  void step(RT h);
};

template <typename RT>
void FRW<RT>::step(RT h)
{
  RT source = 0.0;
  for(auto& x: fluids)
  {
    source += x.first*(1.0+x.second);
  }

  RT new_phi = -K/6.0;
  RT new_K = 12*PI_L*source;

  phi += h*new_phi;
  K += h*new_K;
}

} // namespace cosmo

#endif
