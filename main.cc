
#include "cosmo.h"

using namespace std;
using namespace cosmo;

int main(int argc, char **argv)
{
  _timer["MAIN"].start();


  _timer["loop"].start();
  WAVE wave;
  wave.init();
  wave.step();
  _timer["loop"].stop();


  _timer["MAIN"].stop();

  cout << endl << _timer << endl;

  return EXIT_SUCCESS;
}
