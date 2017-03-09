#!/bin/bash

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# And up two directories should be the main codebase.
pwd

###
# Compile
###
mkdir -p ../build
cd ../build
TMP_CONFIG_FILE=particles_vector.txt.test
cp ../config/particles_vector.txt $TMP_CONFIG_FILE


do_runs () {
  PEAK_AMPLITUDE=$1
  PARTICLES_PER_DY=$2
  SMOOTHING_RADIUS=$3
  INITIAL_SHIFT=$4
  SLICING=$5
  sed -i -E "s/peak_amplitude = [\.0-9]+/peak_amplitude = $PEAK_AMPLITUDE/g" $TMP_CONFIG_FILE
  sed -i -E "s/particles_per_dy = [0-9]+/particles_per_dy = $PARTICLES_PER_DY/g" $TMP_CONFIG_FILE
  sed -i -E "s/smoothing_radius = [\.0-9]+/smoothing_radius = $SMOOTHING_RADIUS/g" $TMP_CONFIG_FILE
  sed -i -E "s/use_initial_shift = [0-9]+/use_initial_shift = $INITIAL_SHIFT/g" $TMP_CONFIG_FILE

  if [ "$SLICING" -eq "0" ]; then
    sed -i -E "s/shift = [A-Za-z]+/shift = RedShift/g" $TMP_CONFIG_FILE
    sed -i -E "s/lapse = [A-Za-z]+/lapse = OnePlusLog/g" $TMP_CONFIG_FILE
    USE_GAMMA_DRIVER="true"
  else
    sed -i -E "s/shift = [A-Za-z]+/shift = RedShift/g" $TMP_CONFIG_FILE
    sed -i -E "s/lapse = [A-Za-z]+/lapse = Static/g" $TMP_CONFIG_FILE
    USE_GAMMA_DRIVER="false"
  fi

  DIR="vector_run-A_$PEAK_AMPLITUDE-ppdy_$PARTICLES_PER_DY-rs_$SMOOTHING_RADIUS-bi_$INITIAL_SHIFT-gauge_$SLICING"
  mkdir -p $DIR

  declare -a RESOLUTIONS=("0016" "0032" "0064" "0128" "1024")
  for r in "${RESOLUTIONS[@]}"
  do
    RES=$(echo $r | sed 's/^0*//')
    STEPS=$((RES*100))
    printf "Performing run with N = $RES\n"
    cmake -DCOSMO_N=$RES -DCOSMO_NX=1 -DCOSMO_NZ=1 -DCOSMO_STENCIL_ORDER=2 -DCOSMO_USE_GAMMA_DRIVER=$USE_GAMMA_DRIVER .. && make -j32
    sed -i -E "s/steps = [0-9]+/steps = $STEPS/g" $TMP_CONFIG_FILE
    sed -i -E "s,output_dir = [[:alnum:]_-\./]+,output_dir = $DIR/R$r,g" $TMP_CONFIG_FILE
    ./cosmo $TMP_CONFIG_FILE
  done
}

for B in 0.2 0.02 # amplitude
do
  for ppdy in 1 2 # particles per fy
  do
    for rs in 1.5 2.5 # softening radius
    do
      for bi in 0 1 # initial shift
      do
        for g in 0 1 # gauge
        do
          printf "Running: do_runs $B $ppdy $rs $bi $g\n"
          do_runs $B $ppdy $rs $bi $g
        done
      done
    done
  done
done

rm $TMP_CONFIG_FILE
