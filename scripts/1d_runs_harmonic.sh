#!/bin/bash

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# up a directory should be the main codebase.
cd ..
mkdir -p build/1d_runs
cd build/1d_runs

for L in 40 20 10 5 2.5 1.25 0.625 0.3125 0.15625 0.078125 # box length
do
  for A in 0.0008 0.0004 0.0002 0.0001 0.00005 0.000025 0.0000125 0.00000625 0.000003125 # mode amp.
  do
    ../../scripts/1d_runs.sh -C -N=0032 -c=4 -A=$A -l=$L -g=Harmonic
    ../../scripts/1d_runs.sh -C -N=0064 -c=8 -A=$A -l=$L -g=Harmonic
    ../../scripts/1d_runs.sh -C -N=0128 -c=16 -A=$A -l=$L -g=Harmonic
  done
done
