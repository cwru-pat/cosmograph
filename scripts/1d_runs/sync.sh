#!/bin/bash

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# up a directory should be the main codebase.
cd ../..
mkdir -p build/1d_dust_runs
cd build/1d_dust_runs
pwd

declare -a Ls=('0.1' '0.12865616' '0.16552408' '0.213' '0.274' '0.352' '0.456' '0.583' '0.751' '0.966' '1.24' '1.60' '2.06' '2.65' '3.40' '4.38' '5.63' '7.25' '9.33' '12.0')
declare -a As=('4.629e-03' '3.219e-03' '2.181e-03' '1.444e-03' '9.382e-04' '6.000e-04' '3.789e-04' '2.368e-04' '1.468e-04' '9.038e-05' '5.538e-05' '3.382e-05' '2.060e-05' '1.253e-05' '7.615e-06' '4.623e-06' '2.804e-06' '1.697e-06' '1.022e-06' '6.090e-07')

for ((i=0; i<20; i+=1));
do
  A=${As[$i]}
  L=${Ls[$i]}
  echo "Running for L = $L, A = $A..."
  ../../scripts/1d_runs/1d_runs.sh -C -N=0064 -A=$A -l=$L -g=Static
  ../../scripts/1d_runs/1d_runs.sh -C -N=0096 -A=$A -l=$L -g=Static
  ../../scripts/1d_runs/1d_runs.sh -C -N=0128 -A=$A -l=$L -g=Static
done

