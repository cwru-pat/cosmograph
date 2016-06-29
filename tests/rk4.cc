// g++ --std=c++11 rk4.cc -O0 && ./a.out

#include <cmath>
#include <iostream>
#include "../utils/RK4Register.h"

using namespace cosmo;

typedef double real_t;

real_t ev_x(real_t x)
{
  // difeq to integrate; x' = f(x)
  // f(x) is:
  return std::exp(-x);
}

real_t soln_t(real_t t)
{
  // solution, x(t)
  return std::log(t);
}

int main()
{
  real_t t0 = 1.0;
  real_t dt = 0.1;
  int max_steps = 10000;

  RK4Register<int, real_t> x;
  x.init(1, 1, 1, dt);
  x._array_p[0] = soln_t(t0); // "initial conditions"

  for(int step=0; step < max_steps; ++step)
  {
    x.stepInit();
    x._array_c[0] = ev_x(x._array_a[0]);
    x.K1Finalize();
    x._array_c[0] = ev_x(x._array_a[0]);
    x.K2Finalize();
    x._array_c[0] = ev_x(x._array_a[0]);
    x.K3Finalize();
    x._array_c[0] = ev_x(x._array_a[0]);
    x.K4Finalize();
  }

  std::cout << "Final answer is: " << x._array_p[0]
    << " (analytic soln is " << soln_t(1.0 + max_steps*dt) << ")"
    << std::endl;

  real_t residual = x._array_p[0] - soln_t(1.0 + max_steps*dt);
  std::cout << "Residual is: " << residual
    << " with N*dt^4 = " << dt*dt*dt*dt*max_steps << std::endl;

  if(std::abs(residual) > 10e-5)
  {
    std::cout << "Error: unusually large error detected in RK4 integrator!";
    throw -1;
  }

  exit(EXIT_SUCCESS);
}
