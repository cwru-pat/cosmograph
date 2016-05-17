// g++ --std=c++11 timer.cc ../utils/Timer.cc -lrt -O0 && ./a.out

#include "../utils/Timer.h"
#include <iostream>
#include <string>

int main()
{
  cosmo::TimerManager timer;

  timer["half sec pause"].start();
  system("sleep 0.5");
  timer["half sec pause"].stop();

  timer["quarter sec pause"].start();
  system("sleep 0.25");
  timer["quarter sec pause"].stop();

  // two methods of printing timing information
  std::cout << timer;
  std::cout << timer.getStateString();

  exit(EXIT_SUCCESS);
}
