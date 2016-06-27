#ifndef COSMO_FRW
#define COSMO_FRW

#define PI_L 3.141592653589793238462643383279502884L
#include <vector>

namespace cosmo
{

/** FRW class **/
template<typename RT>
class FRW
{
  // metric variables in FRW
  RT phi;
  RT K;
  RT alpha;
  // densities and corresponding EOS "w"s
  std::vector< std::pair<RT,RT> > fluids;
  int num_fluids;

  // RK variables
  RT phi_K1, K_K1, alpha_K1,
     phi_K2, K_K2, alpha_K2,
     phi_K3, K_K3, alpha_K3,
     phi_K4, K_K4, alpha_K4;
  std::vector<RT> fluids_K1,
                  fluids_K2,
                  fluids_K3,
                  fluids_K4;

  // variables to return
  RT phi_get, K_get, alpha_get, rho_get, S_get;

public:
  FRW(RT phi_in, RT K_in)
  {
    phi_get = phi = phi_in;
    K_get = K = K_in;
    alpha_get = alpha = 1.0;
    num_fluids = 0.0;
    rho_get = 0.0;
    S_get = 0.0;
  }

  // "set" variables
  void set_phi(RT phi_in)
  {
    phi_get = phi = phi_in;
  }
  void set_K(RT K_in)
  {
    K_get = K = K_in;
  }
  void set_alpha(RT alpha_in)
  {
    alpha_get = alpha = alpha_in;
  }
  void addFluid(RT rho_in, RT w_in)
  {
    fluids.emplace_back(rho_in, w_in);
    num_fluids += 1;
    fluids_K1.resize(num_fluids);
    fluids_K2.resize(num_fluids);
    fluids_K3.resize(num_fluids);
    fluids_K4.resize(num_fluids);
    rho_get = 0.0;
    S_get = 0.0;
    for(int n=0; n<num_fluids; ++n)
    {
      RT rho = fluids[n].first;
      rho_get += rho;
      S_get += 3.0*rho*fluids[n].second;
    }
  }
  
  // get variables
  RT get_phi() { return phi_get; }
  RT get_K() { return K_get; }
  RT get_alpha() { return alpha_get; }
  RT get_rho() { return rho_get; }
  RT get_S() { return S_get; }

  // RK calculations
  void P1_step(RT h);
  void P2_step(RT h);
  void P3_step(RT h);
  void RK_total_step(RT h);

};

template <typename RT>
void FRW<RT>::P1_step(RT h)
{
  // Phi evolution
  phi_K1 = -1.0*K/6.0;

  // K evolution
  RT K_source = 0.0;
  for(int n=0; n<num_fluids; ++n)
  {
    std::pair<RT,RT> x = fluids[n];
    K_source += 4.0*PI_L*x.first*(1.0 + 3.0*x.second);
  }
  K_K1 = 1.0/3.0*K*K + K_source;

  // matter fields
  for(int n=0; n<num_fluids; ++n)
  {
    std::pair<RT,RT> x = fluids[n];
    // d/dt rho = -3*H*rho*(1+w)
    RT rho = x.first;
    fluids_K1[n] = K*rho*(1.0 + x.second);
  }

  // BSSNSim will expect these returned at this point
  phi_get = phi + (h/2.0)*phi_K1;
  K_get = K + (h/2.0)*K_K1;
  // matter fields
  rho_get = 0.0;
  S_get = 0.0;
  for(int n=0; n<num_fluids; ++n)
  {
    RT rho = fluids[n].first + (h/2.0)*fluids_K1[n];
    rho_get += rho;
    S_get += 3.0*rho*fluids[n].second;
  }
}

template <typename RT>
void FRW<RT>::P2_step(RT h)
{
  // Phi evolution
  RT K_interm = K + (h/2.0)*K_K1;
  phi_K2 = -1.0*K_interm/6.0;

  // K evolution
  RT K_source = 0.0;
  for(int n=0; n<num_fluids; ++n)
  {
    std::pair<RT,RT> x = fluids[n];
    RT rho_interm = fluids[n].first + (h/2.0)*fluids_K1[n];
    K_source += 4.0*PI_L*rho_interm*(1.0 + 3.0*x.second);
  }
  K_K2 = 1.0/3.0*K_interm*K_interm + K_source;

  // matter fields
  for(int n=0; n<num_fluids; ++n)
  {
    std::pair<RT,RT> x = fluids[n];
    // d/dt rho = -3*H*rho*(1+w)
    RT rho_interm = fluids[n].first + (h/2.0)*fluids_K1[n];
    fluids_K2[n] = K_interm*rho_interm*(1.0 + x.second);
  }

  // BSSNSim will expect these returned at this point
  phi_get = phi + (h/2.0)*phi_K2;
  K_get = K + (h/2.0)*K_K2;
  // matter fields
  rho_get = 0.0;
  S_get = 0.0;
  for(int n=0; n<num_fluids; ++n)
  {
    RT rho = fluids[n].first + (h/2.0)*fluids_K2[n];
    rho_get += rho;
    S_get += 3.0*rho*fluids[n].second;
  }
}

template <typename RT>
void FRW<RT>::P3_step(RT h)
{
  // Phi evolution
  RT K_interm = K + (h/2.0)*K_K2;
  phi_K3 = -1.0*K_interm/6.0;

  // K evolution
  RT K_source = 0.0;
  for(int n=0; n<num_fluids; ++n)
  {
    std::pair<RT,RT> x = fluids[n];
    RT rho_interm = fluids[n].first + (h/2.0)*fluids_K2[n];
    K_source += 4.0*PI_L*rho_interm*(1.0 + 3.0*x.second);
  }
  K_K3 = 1.0/3.0*K_interm*K_interm + K_source;

  // matter fields
  for(int n=0; n<num_fluids; ++n)
  {
    std::pair<RT,RT> x = fluids[n];
    // d/dt rho = -3*H*rho*(1+w)
    RT rho_interm = fluids[n].first + (h/2.0)*fluids_K2[n];
    fluids_K3[n] = K_interm*rho_interm*(1.0 + x.second);
  }

  // BSSNSim will expect these returned at this point
  phi_get = phi + h*phi_K3;
  K_get = K + h*K_K3;
  // matter fields
  rho_get = 0.0;
  S_get = 0.0;
  for(int n=0; n<num_fluids; ++n)
  {
    RT rho = fluids[n].first + h*fluids_K3[n];
    rho_get += rho;
    S_get += 3.0*rho*fluids[n].second;
  }
}

template <typename RT>
void FRW<RT>::RK_total_step(RT h)
{
  // Phi evolution
  RT K_interm = K + h*K_K3;
  phi_K4 = -1.0*K_interm/6.0;

  // K evolution
  RT K_source = 0.0;
  for(int n=0; n<num_fluids; ++n)
  {
    std::pair<RT,RT> x = fluids[n];
    RT rho_interm = fluids[n].first + h*fluids_K3[n];
    K_source += 4.0*PI_L*rho_interm*(1.0 + 3.0*x.second);
  }
  K_K4 = 1.0/3.0*K_interm*K_interm + K_source;

  // matter fields
  for(int n=0; n<num_fluids; ++n)
  {
    std::pair<RT,RT> x = fluids[n];
    // d/dt rho = -3*H*rho*(1+w)
    RT rho_interm = fluids[n].first + h*fluids_K3[n];
    fluids_K4[n] = K_interm*rho_interm*(1.0 + x.second);
  }

  // RK update step
  phi += h/6.0*(phi_K1 + 2.0*phi_K2 + 2.0*phi_K3 + phi_K4);
  K += h/6.0*(K_K1 + 2.0*K_K2 + 2.0*K_K3 + K_K4);
  for(int n=0; n<num_fluids; ++n)
    fluids[n].first += h/6.0*(fluids_K1[n] + 2.0*fluids_K2[n] + 2.0*fluids_K3[n] + fluids_K4[n]);

  // BSSNSim will expect these returned at this point
  phi_get = phi;
  K_get = K;
  rho_get = 0.0;
  S_get = 0.0;
  for(int n=0; n<num_fluids; ++n)
  {
    RT rho = fluids[n].first;
    rho_get += rho;
    S_get += 3.0*rho*fluids[n].second;
  }
}

} // namespace cosmo

#endif
