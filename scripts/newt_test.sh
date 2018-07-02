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
TMP_CONFIG_FILE=phase_space_sheets.txt.test
cp ../config/tests/phase_space_sheets.txt $TMP_CONFIG_FILE

do_runs () {
  BOX_LENGTH=$1
  RES_STR=$2
  CARRIERS_PER_DX=$3
  CARRIER_COUNT_SCHEME=$4
  DEPOSIT_SCHEME=$5
  METHOD_ORDER_RES=$6
  MODE_AMPLITUDE=$7
  GAUGE=$8

  RES_INT=$(echo $RES_STR | sed 's/^0*//')
  IO_INT=$((RES_INT/4))
  IO_INT_3D=$((RES_INT*50))
  STEPS=$((RES_INT*500))
  METHOD_ORDER=$((METHOD_ORDER_RES*2))
  IC_PPDX=$((RES_INT*20000))

  sed -i -E "s/ns1 = [\.0-9]+/ns1 = $RES_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/carriers_per_dx = [0-9]+/carriers_per_dx = $CARRIERS_PER_DX/g" $TMP_CONFIG_FILE
  sed -i -E "s/carrier_count_scheme = [\.0-9]+/carrier_count_scheme = $CARRIER_COUNT_SCHEME/g" $TMP_CONFIG_FILE
  sed -i -E "s/deposit_scheme = [0-9]+/deposit_scheme = $DEPOSIT_SCHEME/g" $TMP_CONFIG_FILE
  sed -i -E "s/integration_points_per_dx = [0-9]+/integration_points_per_dx = $IC_PPDX/g" $TMP_CONFIG_FILE
  sed -i -E "s/steps = [0-9]+/steps = $STEPS/g" $TMP_CONFIG_FILE
  sed -i -E "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = $RES_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/IO_3D_grid_interval = [0-9]+/IO_3D_grid_interval = $IO_INT_3D/g" $TMP_CONFIG_FILE
  sed -i -E "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = $IO_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = $IO_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/SVT_constraint_interval = [0-9]+/SVT_constraint_interval = $IO_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/peak_amplitude = [\.0-9]+/peak_amplitude = $MODE_AMPLITUDE/g" $TMP_CONFIG_FILE
  sed -i -E "s/lapse = [a-zA-Z]+/lapse = $GAUGE/g" $TMP_CONFIG_FILE
  sed -i -E "s/dt_frac = [\.0-9]+/dt_frac = 0.001/g" $TMP_CONFIG_FILE

  DIR="sheet_run-L_$BOX_LENGTH-r_$RES_STR-Odx_$METHOD_ORDER-cpdx_$CARRIERS_PER_DX-ccs_$CARRIER_COUNT_SCHEME-ds_$DEPOSIT_SCHEME-A_$MODE_AMPLITUDE_$GAUGE"
  sed -i -E "s,output_dir = [[:alnum:]_-\./]+,output_dir = $DIR,g" $TMP_CONFIG_FILE
  mkdir -p $DIR
  printf "Performing run with settings in dir $DIR.\n"

  cmake .. -DCOSMO_N=$RES_INT -DCOSMO_NY=1 -DCOSMO_NZ=1 -DCOSMO_USE_GENERALIZED_NEWTON=1\
     -DCOSMO_STENCIL_ORDER=$METHOD_ORDER -DCOSMO_USE_REFERENCE_FRW=0 -DCOSMO_H_LEN_FRAC=$BOX_LENGTH\
     -DCOSMO_USE_Z4c_DAMPING=1 && make -j24
  if [ $? -ne 0 ]; then
    echo "Error: compilation failed!"
  else
    ./cosmo $TMP_CONFIG_FILE
  fi
}

do_runs 0.5 "0016" 4 1 1 4 0.002 GeneralizedNewton
do_runs 0.5 "0032" 8 1 1 4 0.002 GeneralizedNewton
do_runs 0.5 "0064" 16 1 1 4 0.002 GeneralizedNewton
do_runs 0.5 "0128" 32 1 1 4 0.002 GeneralizedNewton

rm $TMP_CONFIG_FILE
