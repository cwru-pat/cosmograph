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

sed -i -E "s/ns2 = [\.0-9]+/ns2 = 2/g" $TMP_CONFIG_FILE
sed -i -E "s/ns3 = [\.0-9]+/ns3 = 2/g" $TMP_CONFIG_FILE

do_runs () {
  BOX_LENGTH=$1
  RES_STR=$2
  CARRIERS_PER_DX=$3
  CARRIER_COUNT_SCHEME=$4
  DEPOSIT_SCHEME=$5

  RES_INT=$(echo $RES_STR | sed 's/^0*//')
  STEPS=$((RES_INT*100))

  sed -i -E "s/ns1 = [\.0-9]+/ns1 = $RES_INT/g" $TMP_CONFIG_FILE
  sed -i -E "s/carriers_per_dx = [0-9]+/carriers_per_dx = $CARRIERS_PER_DX/g" $TMP_CONFIG_FILE
  sed -i -E "s/carrier_count_scheme = [\.0-9]+/carrier_count_scheme = $CARRIER_COUNT_SCHEME/g" $TMP_CONFIG_FILE
  sed -i -E "s/deposit_scheme = [0-9]+/deposit_scheme = $DEPOSIT_SCHEME/g" $TMP_CONFIG_FILE
  sed -i -E "s/steps = [0-9]+/steps = $STEPS/g" $TMP_CONFIG_FILE

  DIR="sheet_run-L_$BOX_LENGTH-r_$RES_STR-cpdx_$CARRIERS_PER_DX-ccs_$CARRIER_COUNT_SCHEME-ds_$DEPOSIT_SCHEME"
  sed -i -E "s,output_dir = [[:alnum:]_-\./]+,output_dir = $DIR,g" $TMP_CONFIG_FILE
  mkdir -p $DIR
  printf "Performing run with settings in dir $DIR.\n"
  
  cmake .. -DCOSMO_N=$RES_INT -DCOSMO_NY=1 -DCOSMO_NZ=1 -DCOSMO_STENCIL_ORDER=2 -DCOSMO_USE_REFERENCE_FRW=0 -DCOSMO_H_LEN_FRAC=$BOX_LENGTH && make -j24
  if [ $? -ne 0 ]; then
    echo "Error: compilation failed!"
  else
    ./cosmo $TMP_CONFIG_FILE
  fi
}

# for L in 0.5 0.05 # box length
# do
#   for r in "0064" "0128" "0256" # resolution
#   do
#     for cpdx in 2 10 # carriers per dx
#     do
#       for ccs in 1 # carrier particle counting scheme (per dx, ds) # per dx is very unstable
#       do
#         for ds in 1 2 # deposit scheme (0=CIC, 1=PCS, 2=CINT)
#         do
#           printf "Running: do_runs $L $r $cpdx $ccs $ds\n"
#           do_runs $L $r $cpdx $ccs $ds
#         done
#       done
#     done
#   done
# done

# do_runs 0.05 64 10 1 1

do_runs 0.5 32768 100 1 1

rm $TMP_CONFIG_FILE
