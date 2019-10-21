#!/bin/bash

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# And up two directories should be the main codebase.
pwd

module load cmake

###
# Compile
###
mkdir -p ../build/dust
cd ../build/dust
TMP_CONFIG_FILE=dust_raytracing.txt.test2
cp ../../config/dust_raytracing.txt $TMP_CONFIG_FILE

do_runs () {
  RES=$1

  RES_INT=$(echo $RES | sed 's/^0*//')
  IO_INT=$((RES_INT*50))
  STEPS=$((RES_INT*160000))
  STEPSO2=$((RES_INT*160000/2))

  sed -i -E "s/steps = [0-9]+/steps = $STEPS/g" $TMP_CONFIG_FILE
  sed -i -E "s/ray_flip_step = [0-9]+/ray_flip_step = $STEPSO2/g" $TMP_CONFIG_FILE

  sed -i -E "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = $IO_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = $IO_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = $IO_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/IO_raytrace_interval = [0-9]+/IO_raytrace_interval = $IO_INT/g" $TMP_CONFIG_FILE

  DIR="dust_run-N_$RES_INT"
  sed -i -E "s,output_dir = [[:alnum:]_-\./]+,output_dir = $DIR,g" $TMP_CONFIG_FILE
  mkdir -p $DIR
  printf "Performing run with settings in dir $DIR.\n"

  # normal run
  cmake ../.. -DCOSMO_N=$RES_INT -DCOSMO_NY=1 -DCOSMO_NZ=1 -DCOSMO_STENCIL_ORDER=8 -DCOSMO_H_LEN_FRAC=0.15 -DCOSMO_USE_REFERENCE_FRW=0 && make -j24
  if [ $? -ne 0 ]; then
    echo "Error: compilation failed!"
  else
    ./cosmo $TMP_CONFIG_FILE
  fi

  # linear run
  # cmake .. -DCOSMO_N=$RES_INT -DCOSMO_NY=1 -DCOSMO_NZ=1 -DCOSMO_STENCIL_ORDER=8 -DCOSMO_H_LEN_FRAC=0.15 -DCOSMO_USE_REFERENCE_FRW=0 -DCOSMO_EXCLUDE_SECOND_ORDER_SMALL=1 && make -j24
  # if [ $? -ne 0 ]; then
  #   echo "Error: compilation failed!"
  # else
  #   ./cosmo $TMP_CONFIG_FILE
  # fi
}

do_runs_plain () {

  BOX_LENGTH=$1
  RES_STR=$2
  MODE_AMPLITUDE=$3

  RES_INT=$(echo $RES_STR | sed 's/^0*//')
  IO_INT=$((RES_INT*10))
  STEPS=$((RES_INT*10000))

  sed -i -E "s/steps = [0-9]+/steps = $STEPS/g" $TMP_CONFIG_FILE
  sed -i -E "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = $IO_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = $IO_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = $IO_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/SVT_constraint_interval = [0-9]+/SVT_constraint_interval = $IO_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/peak_amplitude_frac = [\.0-9]+/peak_amplitude_frac = $MODE_AMPLITUDE/g" $TMP_CONFIG_FILE
  sed -i -E "s/ray_flip_step = [0-9]+/ray_flip_step = $STEPS/g" $TMP_CONFIG_FILE

  DIR="dust_run-L_$BOX_LENGTH-r_$RES_STR-A_$MODE_AMPLITUDE"
  sed -i -E "s,output_dir = [[:alnum:]_-\./]+,output_dir = $DIR,g" $TMP_CONFIG_FILE
  mkdir -p $DIR
  printf "Performing run with settings in dir $DIR.\n"

  cmake ../.. -DCOSMO_N=$RES_INT -DCOSMO_NY=1 -DCOSMO_NZ=1 -DCOSMO_STENCIL_ORDER=8 -DCOSMO_USE_REFERENCE_FRW=0 -DCOSMO_H_LEN_FRAC=$BOX_LENGTH && make -j24
  if [ $? -ne 0 ]; then
    echo "Error: compilation failed!"
  else
    ./cosmo $TMP_CONFIG_FILE
  fi

}

do_runs_plain 0.02 "0064" 0.0005
# do_runs_plain 0.02 "0128" 0.0005
# do_runs_plain 0.02 "0256" 0.0005
# do_runs_plain 0.02 "0512" 0.0005

# do_runs "0016"
# do_runs "0032"
# do_runs "0064"
# do_runs "0128"

rm $TMP_CONFIG_FILE
