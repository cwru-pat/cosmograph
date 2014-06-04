
#include "cosmo.h"

using namespace std;
using namespace cosmo;

int main(int argc, char **argv)
{
  _timer["MAIN"].start();

  const idx_t N = 1000000;

  real_t v = 0;
  _timer["loop"].start();
  for(idx_t i=0; i < N; ++i) {
    v += i;
  }
  _timer["loop"].stop();
  cout << "v: " << v << endl;

  _timer["MAIN"].stop();

  cout << endl << _timer << endl;

  return EXIT_SUCCESS;
}
