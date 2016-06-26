// g++ --std=c++11 array.cc -O0 && ./a.out

#include <chrono>
#include <iostream>
#include "../utils/Array.h"
#include "../utils/math.h"
#include "../cosmo_globals.h"

using namespace cosmo;
real_t dx;

#define N_SMALL 10
#define N_LARGE 500
#define N_FN_CALLS 10000

// Different test functions that manipulate array values.
void performSomeComputation(arr_t & arr)
{
  for(int i=0; i<N_LARGE; ++i)
  {
    arr[i] = i/10;
  }
  for(int i=0; i<N_LARGE; ++i)
  {
    arr[0] += arr[(i+10)%N_LARGE];
  }
}
void performSomeComputation(real_t * arr)
{
  for(int i=0; i<N_LARGE; ++i)
  {
    arr[i] = i/10;
  }
  for(int i=0; i<N_LARGE; ++i)
  {
    arr[0] += arr[(i+10)%N_LARGE];
  }
}
void performSomeComputationWrapper(arr_t & arr)
{
  performSomeComputation(arr._array);
}
void performSomeComputationAlt(arr_t & arr)
{
  real_t * arr_ref = arr._array;
  for(int i=0; i<N_LARGE; ++i)
  {
    arr_ref[i] = i/10;
  }
  for(int i=0; i<N_LARGE; ++i)
  {
    arr_ref[0] += arr_ref[(i+10)%N_LARGE];
  }
}

int main()
{
  dx = 1.0;

  arr_t myArr1, myArr2;
  myArr1.init(1, 1, N_SMALL);
  myArr2.init(1, 1, N_SMALL);

  // check initialization
  for(int i=0; i<N_SMALL; ++i)
  {
    if(myArr1[i] != 0 || myArr2[i] != 0)
    {
      std::cout << "Error: initializing array failed.\n";
      throw -1;
    }
  }

  // check 'copy by value'
  myArr1[0] = 1;
  myArr2 = myArr1;
  if(myArr2[0] != 1 || myArr1[0] != 1)
  {
    std::cout << "Error: copying array failed (1).\n";
    throw -2;
  }
  myArr2[N_SMALL-1] = 2;
  if(myArr2[N_SMALL-1] != 2 || myArr1[N_SMALL-1] != 0)
  {
    throw -3;
    std::cout << "Error: copying array failed (2).\n";
  }
  myArr1[N_SMALL-2] = 3;
  if(myArr1[N_SMALL-2] != 3 || myArr2[N_SMALL-2] != 0)
  {
    throw -4;
    std::cout << "Error: copying array failed (3).\n";
  }

  // check swap
  cosmoArraySwap(myArr2, myArr1);
  if(myArr1[N_SMALL-1] != 2 || myArr2[N_SMALL-1] != 0)
  {
    throw -3;
    std::cout << "Error: copying array failed (4).\n";
  }
  myArr2[N_SMALL-3] = 2;
  if(myArr2[N_SMALL-3] != 2 || myArr1[N_SMALL-3] != 0)
  {
    throw -3;
    std::cout << "Error: copying array failed (5).\n";
  }
  myArr1[N_SMALL-4] = 3;
  if(myArr1[N_SMALL-4] != 3 || myArr2[N_SMALL-4] != 0)
  {
    throw -4;
    std::cout << "Error: copying array failed (6).\n";
  }


  // test speed vs. "normal" array
  std::chrono::steady_clock::time_point start, end;

  start = std::chrono::steady_clock::now();
  arr_t bigArr (N_LARGE); // N_LARGE^3 array
  real_t * big_arr = new real_t[N_LARGE*N_LARGE*N_LARGE];
  end = std::chrono::steady_clock::now();
  std::cout << "Allocation took "
            << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
            << "us.\n";

  start = std::chrono::steady_clock::now();
  for(int i=0; i<N_FN_CALLS; ++i)
    performSomeComputation(big_arr);
  end = std::chrono::steady_clock::now();
  std::cout << "Array function took "
            << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
            << "us.\n";

  start = std::chrono::steady_clock::now();
  for(int i=0; i<N_FN_CALLS; ++i)
    performSomeComputation(bigArr);
  end = std::chrono::steady_clock::now();
  std::cout << "Class function took "
            << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
            << "us.\n";

  start = std::chrono::steady_clock::now();
  for(int i=0; i<N_FN_CALLS; ++i)
    performSomeComputationWrapper(bigArr);
  end = std::chrono::steady_clock::now();
  std::cout << "Wrapper function took "
            << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
            << "us.\n";

  start = std::chrono::steady_clock::now();
  for(int i=0; i<N_FN_CALLS; ++i)
    performSomeComputationAlt(bigArr);
  end = std::chrono::steady_clock::now();
  std::cout << "Alternative Class function took "
            << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
            << "us.\n";


  arr_t cosmoArr (NX, NY, NZ); // Cosmo sized arrays
  real_t * cosmo_arr = new real_t[NX*NY*NZ];
  for(int i=0; i<POINTS; ++i)
  {
    cosmoArr[i] = ((double) i) / 100 / POINTS;
    cosmo_arr[i] = ((double) i) / 100 / POINTS;
  }

  exit(EXIT_SUCCESS);
}
