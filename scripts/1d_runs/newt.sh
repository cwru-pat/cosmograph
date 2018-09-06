#!/bin/bash

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# up a directory should be the main codebase.
cd ../..
mkdir -p build/1d_runs_newt
cd build/1d_runs_newt

# zI = 500
# declare -a As=('0.0000523831' '0.0000363' '0.0000244' '0.0000158' '9.76e-6' '5.75-6' '3.25e-6' '1.78e-6' '9.58e-7' '5.08e-7' '2.68e-7' '1.41e-7' '7.37e-8' '3.86e-8' '2.02e-8' '1.06e-8' '5.51e-9' '2.87e-9' '1.47e-9' '7.27e-10')
# declare -a Ls=('0.1' '0.139' '0.192' '0.267' '0.370' '0.513' '0.712' '0.987' '1.37' '1.90' '2.63' '3.65' '5.07' '7.03' '9.74' '13.5' '18.7' '26.0' '36.1' '50.0')

# zI = 50
# declare -a Ls=('0.213' '0.274' '0.352' '0.456' '0.583' '0.751' '0.966' '1.24' '1.60' '2.06' '2.65' '3.40' '4.38' '5.63' '7.25' '9.33' '12.0')
# declare -a As=('1.444e-03' '9.382e-04' '6.000e-04' '3.789e-04' '2.368e-04' '1.468e-04' '9.038e-05' '5.538e-05' '3.382e-05' '2.060e-05' '1.253e-05' '7.615e-06' '4.623e-06' '2.804e-06' '1.697e-06' '1.022e-06' '6.090e-07')

declare -a Ls=('0.1' '0.12865616' '0.16552408')
declare -a As=('4.629e-03' '3.219e-03' '2.181e-03')

for ((i=0; i<3; i+=1));
do
  A=${As[$i]}
  L=${Ls[$i]}
  echo "Running for L = $L, A = $A..."
  ../../scripts/1d_runs.sh -C -N=0064 -A=$A -l=$L -g=GeneralizedNewton --GN_eta=0.01
  ../../scripts/1d_runs.sh -C -N=0096 -A=$A -l=$L -g=GeneralizedNewton --GN_eta=0.01
  ../../scripts/1d_runs.sh -C -N=0128 -A=$A -l=$L -g=GeneralizedNewton --GN_eta=0.01
done
