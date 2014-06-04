
#include "cosmo.h"

using namespace std;
using namespace cosmo;

int main(int argc, char **argv)
{
  _timer["MAIN"].start();

  const idx_t n = 10;

  WAVE wave;
  wave.init();

  _timer["loop"].start();
  for(idx_t i=0; i < n; ++i) {
    wave.step();
  }
  _timer["loop"].stop();

  _timer["MAIN"].stop();

  cout << endl << _timer << endl;

  return EXIT_SUCCESS;
}
