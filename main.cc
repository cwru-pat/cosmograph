
#include "cosmo.h"

using namespace std;
using namespace cosmo;

int main(int argc, char **argv)
{
  TimerManager _timer;
  _timer["MAIN"].start();

  const idx_t n = 10;

  Wave wave;

  _timer["loop"].start();
  for(idx_t i=0; i < n; ++i) {
    wave.step();
    wave.dump_strip("phi", 1, N/2, N/2);
  }
  _timer["loop"].stop();

  _timer["MAIN"].stop();

  cout << endl << _timer << endl;

  return EXIT_SUCCESS;
}
