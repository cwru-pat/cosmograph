// g++ --std=c++11 reference_frw.cc && ./a.out

#include "../utils/FRW.h"
#include <iostream>
#include <string>

#define PW2(x) ((x)*(x))

typedef float real_t;

int main()
{
  real_t dt = 0.01; // seems close to optimal for floats
  std::cout << "Running test (dt = " << dt << ")..." << std::endl;

  real_t H0 = 2.89440502;

  cosmo::FRW<real_t> frw (0.0, -3.0*H0);
  // some sort of universe "content"
  frw.addFluid(0.25, 0.0);
  frw.addFluid(0.25, 1.0/3.0);
  frw.addFluid(0.25, -1.0);
  frw.addFluid(0.25, -0.5);

  for(real_t t = 0; t < 1; t+=dt)
  {
    frw.P1_step(dt);
    frw.P2_step(dt);
    frw.P3_step(dt);
    frw.RK_total_step(dt);
  }

  std::cout << "Final values: " << std::endl;
  std::cout << "  phi = " << frw.get_phi() << std::endl;
  std::cout << "  rho = " << frw.get_rho() << std::endl;
  std::cout << "  S = " << frw.get_S() << std::endl;
  std::cout << "  K^2 = " << PW2(frw.get_K()) << std::endl;
  std::cout << "  K_eq^2 = " << 75.3982237*frw.get_rho() << std::endl;
  std::cout << "  FRW equation residual: "
      << PW2(frw.get_K()) - 75.3982237*frw.get_rho() << std::endl;

  exit(EXIT_SUCCESS);
}
