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
  RT get_phi()
  {
    return phi;
  }
  RT get_K()
  {
    return K;
  }
};

template <typename RT>
void FRW<RT>::step(RT h)
{
  // phi evolution term
  RT new_phi = -K/6.0;

  // K evolution term
  // Also evolve the fluid here.
  RT K_source = 0.0;
  for(auto& x: fluids)
  {
    // source is 4*pi*(3*rho + S)
    // = 4*pi*(3*rho + 3*P)
    // = 12*pi*rho*(1 + w)
    K_source += 12*PI_L*x.first*(1.0 + x.second);

    // K is now sourced; evolve the fluid here so we
    // don't have to loop a bunch later.
    // d/dt rho = -3*H*rho*(1+w)
    RT new_rho = K*x.first*(1.0 + x.second);
    x.first += h*new_rho;
  }
  RT new_K = K_source;

  // Eulerian integration
  phi += h*new_phi;
  K += h*new_K;
}

} // namespace cosmo

#endif
