#!/bin/bash

RES_N=8
# read in options
for i in "$@"
do
  case $i in
      -h|--help)
      printf "Usage: ./particle_SVT_test.sh [(-N|--resolution-N)=128]\n"
      exit 0
      ;;
      -N=*|--resolution-N=*)
      RES_N="${i#*=}"
      shift # past argument=value
      ;;
      *)
        printf "Unrecognized option will not be used: ${i#*=}\n"
        # unknown option
      ;;
  esac
done

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# up a directory should be the main codebase.

cd ..

mkdir -p build
cd build

mkdir -p particle_SVT_test
cd particle_SVT_test

cp ../../config/particles_sinusoid.txt config.txt

do_runs () {
  H_FRAC=$1
  SIM_LEN_FRAC=$(bc -l <<< "1.0/$H_FRAC")
  STEPS=$((1000*$RES_N*$H_FRAC/8))
  OUT_INT=$((10*$RES_N*$H_FRAC/8))
  echo "Running with H_FRAC=$H_FRAC, SIM_LEN_FRAC=$SIM_LEN_FRAC, STEPS=$STEPS, OUT_INT=$OUT_INT."

  sed -i.bak -r "s/steps = [0-9]+/steps = ${STEPS}/" config.txt
  sed -i.bak -r "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = ${OUT_INT}/" config.txt
  sed -i.bak -r "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = ${OUT_INT}/" config.txt
  sed -i.bak -r "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = ${OUT_INT}/" config.txt
  cmake -DCOSMO_N=$RES_N -DCOSMO_NY=2 -DCOSMO_NZ=2 -DCOSMO_STENCIL_ORDER=2 -DCOSMO_H_LEN_FRAC=$SIM_LEN_FRAC ../.. && make -j32

  sed -i.bak -r "s/k_driver_coeff = [\.0-9]+/k_driver_coeff = 0.00/" config.txt
  sed -i.bak -r "s/output_dir = [A-Za-z0-9._-]+/output_dir = R-${RES_N}_Hfrac-${H_FRAC}_Kcoeff-0/" config.txt
  ./cosmo config.txt

  sed -i.bak -r "s/k_driver_coeff = [\.0-9]+/k_driver_coeff = 0.01/" config.txt
  sed -i.bak -r "s/output_dir = [A-Za-z0-9._-]+/output_dir = R-${RES_N}_Hfrac-${H_FRAC}_Kcoeff-0.01/" config.txt
  ./cosmo config.txt

  sed -i.bak -r "s/k_driver_coeff = [\.0-9]+/k_driver_coeff = 0.1/" config.txt
  sed -i.bak -r "s/output_dir = [A-Za-z0-9._-]+/output_dir = R-${RES_N}_Hfrac-${H_FRAC}_Kcoeff-0.1/" config.txt
  ./cosmo config.txt
}

do_runs 1
do_runs 2
do_runs 10
do_runs 20
do_runs 100
