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
TMP_CONFIG_FILE=dust.txt.test
cp ../config/dust_SVT_test.txt $TMP_CONFIG_FILE

do_runs () {
  BOX_LENGTH=$1
  RES_STR=$2
  CARRIERS_PER_DX=$3
  CARRIER_COUNT_SCHEME=$4
  DEPOSIT_SCHEME=$5
  METHOD_ORDER_RES=$6
  MODE_AMPLITUDE=$7
  GAUGE=$8
  DT_FRAC=$9

  RES_INT=$(echo $RES_STR | sed 's/^0*//')
  IO_INT=$((RES_INT/4))
  IO_INT_3D=$((RES_INT*50))
  STEPS=$((RES_INT*1000))
  METHOD_ORDER=$((METHOD_ORDER_RES*2))
  IC_PPDX=$((RES_INT*20000))

  sed -i -E "s/ns1 = [\.0-9]+/ns1 = $RES_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/steps = [0-9]+/steps = $STEPS/g" $TMP_CONFIG_FILE
  sed -i -E "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = $RES_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/IO_3D_grid_interval = [0-9]+/IO_3D_grid_interval = $IO_INT_3D/g" $TMP_CONFIG_FILE
  sed -i -E "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = $IO_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = $IO_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/SVT_constraint_interval = [0-9]+/SVT_constraint_interval = $IO_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/peak_amplitude_frac = [\.0-9]+/peak_amplitude_frac = $MODE_AMPLITUDE/g" $TMP_CONFIG_FILE
  sed -i -E "s/dt_frac = [\.0-9]+/dt_frac = $DT_FRAC/g" $TMP_CONFIG_FILE  
  sed -i -E "s/lapse = [a-zA-Z]+/lapse = $GAUGE/g" $TMP_CONFIG_FILE
  
  DIR="dust_run-L_$BOX_LENGTH-r_$RES_STR-Odx_$METHOD_ORDER-A_$MODE_AMPLITUDE-GAUGE_$GAUGE"
  sed -i -E "s,output_dir = [[:alnum:]_-\./]+,output_dir = $DIR,g" $TMP_CONFIG_FILE
  mkdir -p $DIR
  printf "Performing run with settings in dir $DIR.\n"

  cmake .. -DCOSMO_N=$RES_INT -DCOSMO_NY=1 -DCOSMO_NZ=1 -DCOSMO_STENCIL_ORDER=$METHOD_ORDER -DCOSMO_USE_REFERENCE_FRW=0 -DCOSMO_H_LEN_FRAC=$BOX_LENGTH -DCOSMO_USE_Z4c_DAMPING=0 && make -j
  if [ $? -ne 0 ]; then
    echo "Error: compilation failed!"
  else
    ./cosmo $TMP_CONFIG_FILE
  fi
}

A=0.000019
L=0.01
Gauge=ExpansionSyncLapse
DT_FRAC=1

# r = 0.005
do_runs $L "0256" 32 1 1 4 $A $Gauge $DT_FRAC


#do_runs 0.5 "0064" 8 1 1 4 0.002 Static
#do_runs 0.5 "0128" 16 1 1 4 0.002 Static
#do_runs 0.5 "0256" 32 1 1 4 0.002 Static
#do_runs 0.5 "0512" 64 1 1 4 0.002 Static

# do_runs 0.05 "0064" 8 1 1 4 0.002 Harmonic
# do_runs 0.05 "0128" 16 1 1 4 0.002 Harmonic
# do_runs 0.05 "0256" 32 1 1 4 0.002 Harmonic
# do_runs 0.05 "0512" 64 1 1 4 0.002 Harmonic

# for L in 0.5 # box length
# do
#   for r in "0064" "0128" "0256" "0512" # resolution
#   do
#     for cpdx in 8 16 32 # carriers per dx
#     do
#       for ccs in 1 # carrier particle counting scheme (per dx, ds) # per dx is very unstable
#       do
#         for ds in 1 # deposit scheme (0=CIC, 1=PCS, 2=CINT)
#         do
#           for Odx in 4 # finite difference order
#           do
#             for A in 0.008 0.004 0.002 0.001 # 0.0005 0.00025
#             do
#               printf "Running: do_runs $L $r $cpdx $ccs $ds $Odx $A\n"
#               do_runs $L $r $cpdx $ccs $ds $Odx $A Harmonic
#             done
#           done
#         done
#       done
#     done
#   done
# done

# do_runs 0.5 8192 8 1 1 4 0.005
# do_runs 0.5 5120 8 1 0 1 0.005

rm $TMP_CONFIG_FILE


