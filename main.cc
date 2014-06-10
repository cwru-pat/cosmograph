
#include "cosmo.h"

using namespace std;
using namespace cosmo;

int main(int argc, char **argv)
{
  TimerManager _timer;
  _timer["MAIN"].start();

  ConfigParser _config("config.txt");

  const idx_t n = 1;

  BSSN simulation;

  _timer["loop"].start();
  for(idx_t i=0; i < n; ++i) {
    simulation.step();
  }
  _timer["loop"].stop();

  _timer["MAIN"].stop();

  cout << endl << _timer << endl;

  return EXIT_SUCCESS;
}
