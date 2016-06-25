// g++ --std=c++11 reference_frw.cc && ./a.out

#include "../utils/FRW.h"
#include <iostream>
#include <string>

#define PW2(x) ((x)*(x))

typedef float real_t;

int main()
{

  real_t H0 = 2.89440502;

  cosmo::FRW<real_t> frw (0.0, -3.0*H0);
  // equal parts matter & lambda (rho's must add to 1)
  frw.addFluid(0.5, 0.0);
  frw.addFluid(0.5, -1.0);

  real_t dt = 0.01;
  for(real_t t = 0; t < 1; t+=dt)
  {
    frw.P1_step(dt);
    frw.P2_step(dt);
    frw.P3_step(dt);
    frw.RK_total_step(dt);

    std::cout << "t = " << t << " | ";
    std::cout << "phi = " << frw.get_phi() << " | ";
    std::cout << "rho = " << frw.get_rho() << " | ";
    std::cout << "S = " << frw.get_S() << " | ";
    std::cout << "K^2 = " << PW2(frw.get_K()) << " | ";
    std::cout << "K_eq^2 = " << PW2(-3.0)*8.37758041*(
        frw.get_rho()
      ) << std::endl;
  }

  std::cout << "Testing." << std::endl;

  exit(EXIT_SUCCESS);
}
