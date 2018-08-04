#!/bin/bash

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# up a directory should be the main codebase.

# Some colors
RED='\033[31m'
GREEN='\033[32m'
YELLOW='\033[33m'
BLUE='\033[36m'
NC='\033[0m' # No Color


# default options
USE_CLUSTER=false
DRY_RUN=false
RES_STR='0064'
BOX_LENGTH=1
GAUGE='Harmonic'
CARRIERS_PER_DX=2
MODE_AMPLITUDE=0.002
# Add'l options, but not read in:
CARRIER_COUNT_SCHEME=1
DEPOSIT_SCHEME=1
METHOD_ORDER_RES=4

# read in options
for i in "$@"
do
  case $i in
      -h|--help)
      printf "Usage: ./1d_runs.sh\n"
      printf "         [-C|--cluster-run] [-d|--dry-run]\n"
      printf "         [(-R|--resolution-str)='0064'] [(-l|--l)=1] [(-A|--amplitude)=0.002]\n"
      printf "         [(-c|--carriers)=2] [(-g|--gauge)=Harmonic]\n"
      exit 0
      ;;
      -N=*|--resolution-N=*)
      RES_STR="${i#*=}"
      shift # past argument=value
      ;;
      -l=*|--l=*)
      BOX_LENGTH="${i#*=}"
      shift # past argument=value
      ;;
      -g=*|--gauge=*)
      GAUGE="${i#*=}"
      shift # past argument=value
      ;;
      -A=*|--amplitude=*)
      MODE_AMPLITUDE="${i#*=}"
      shift # past argument=value
      ;;
      -c=*|--carriers=*)
      CARRIERS_PER_DX="${i#*=}"
      shift # past argument
      ;;
      -C|--cluster-run)
      USE_CLUSTER=true
      shift # past argument
      ;;
      -d|--dry-run)
      DRY_RUN=true
      shift # past argument
      ;;
      *)
        printf "Unrecognized option will not be used: ${i#*=}\n"
        # unknown option
      ;;
  esac
done

# derived vars
RES_INT=$(echo $RES_STR | sed 's/^0*//')
STEPS=$(bc <<< "$RES_INT*8000/$BOX_LENGTH")
IO_INT=$((STEPS/1000))
IO_INT_1D=$((STEPS/10))
IO_INT_3D=$((STEPS/10))
METHOD_ORDER=$((METHOD_ORDER_RES*2))
IC_PPDX=$((RES_INT*20000))



mkdir -p ../build
cd ../build
DIR="1d_runs/sheet_run-L_$BOX_LENGTH-r_$RES_STR-Odx_$METHOD_ORDER-cpdx_$CARRIERS_PER_DX-ccs_$CARRIER_COUNT_SCHEME-ds_$DEPOSIT_SCHEME-A_${MODE_AMPLITUDE}_$GAUGE"

printf "${BLUE}Deploying runs:${NC}\n"
printf "  Will be using directory: $DIR\n"
if "$USE_CLUSTER"; then
  printf "  Using cluster (USE_CLUSTER=$USE_CLUSTER)\n"
else
  printf "  Running locally  (USE_CLUSTER=$USE_CLUSTER)\n"
fi
if "$DRY_RUN"; then
  printf "  ${YELLOW}Performing dry run (code will not be run)${NC}\n"
fi

read -r -t 3 -p "Continue to deploy? Will automatically proceed in 3 seconds... [Y/n]: " response
response=${response,,}    # tolower
if ! [[ $response =~ ^(|y|yes)$ ]] ; then
  printf "${RED}Aborting deployment.${NC}\n"
  exit 1
fi
printf "${GREEN}Deploying...${NC}\n"
printf "\n"

mkdir -p $DIR
cd $DIR

TMP_CONFIG_FILE=config.txt
cp ../../../config/tests/phase_space_sheets.txt $TMP_CONFIG_FILE


sed -i -E "s/ns1 = [\.0-9]+/ns1 = $RES_INT/g" $TMP_CONFIG_FILE
sed -i -E "s/carriers_per_dx = [0-9]+/carriers_per_dx = $CARRIERS_PER_DX/g" $TMP_CONFIG_FILE
sed -i -E "s/carrier_count_scheme = [\.0-9]+/carrier_count_scheme = $CARRIER_COUNT_SCHEME/g" $TMP_CONFIG_FILE
sed -i -E "s/deposit_scheme = [0-9]+/deposit_scheme = $DEPOSIT_SCHEME/g" $TMP_CONFIG_FILE
sed -i -E "s/integration_points_per_dx = [0-9]+/integration_points_per_dx = $IC_PPDX/g" $TMP_CONFIG_FILE
sed -i -E "s/steps = [0-9]+/steps = $STEPS/g" $TMP_CONFIG_FILE
sed -i -E "s/IO_1D_grid_interval = [0-9]+/IO_1D_grid_interval = $IO_INT_1D/g" $TMP_CONFIG_FILE
sed -i -E "s/IO_3D_grid_interval = [0-9]+/IO_3D_grid_interval = $IO_INT_3D/g" $TMP_CONFIG_FILE
sed -i -E "s/IO_bssnstats_interval = [0-9]+/IO_bssnstats_interval = $IO_INT/g" $TMP_CONFIG_FILE
sed -i -E "s/IO_constraint_interval = [0-9]+/IO_constraint_interval = $IO_INT/g" $TMP_CONFIG_FILE
sed -i -E "s/SVT_constraint_interval = [0-9]+/SVT_constraint_interval = $IO_INT/g" $TMP_CONFIG_FILE
sed -i -E "s/peak_amplitude = [\.0-9]+/peak_amplitude = $MODE_AMPLITUDE/g" $TMP_CONFIG_FILE
sed -i -E "s/lapse = [a-zA-Z]+/lapse = $GAUGE/g" $TMP_CONFIG_FILE
sed -i -E "s/dt_frac = [\.0-9]+/dt_frac = 0.1/g" $TMP_CONFIG_FILE
sed -i -E "s,output_dir = [[:alnum:]_-\./]+,output_dir = output,g" $TMP_CONFIG_FILE


# Load modules
if "$USE_CLUSTER"; then
  printf "Loading Modules...\n"
  module load cmake
fi

cmake ../../.. -DCOSMO_N=$RES_INT -DCOSMO_NY=1 -DCOSMO_NZ=1 -DCOSMO_USE_GENERALIZED_NEWTON=1\
   -DCOSMO_STENCIL_ORDER=$METHOD_ORDER -DCOSMO_USE_REFERENCE_FRW=0 -DCOSMO_H_LEN_FRAC=$BOX_LENGTH\
   && make -j24 # -DCOSMO_USE_Z4c_DAMPING=1 
if [ $? -ne 0 ]; then
  echo "Error: compilation failed!"
  exit 1
fi

if "$USE_CLUSTER"; then
  cp ../../../scripts/job_template.slurm job.slurm
  # Adjust job name
  sed -i.bak "s/JOBNAME/1d-$RES_STR-$BOX_LENGTH-$MODE_AMPLITUDE" job.slurm
fi

# Run job, go back up a dir
if [ "$DRY_RUN" = false ]; then
  if "$USE_CLUSTER"; then
    printf "${GREEN}Queueing job.${NC}\n"
    sbatch job.slurm
  else
    printf "${GREEN}Running job.${NC}\n"
    ./cosmo $TMP_CONFIG_FILE
  fi
fi

if "$USE_CLUSTER"; then
  printf "squeue -u jbm120\n"
  squeue -u jbm120
  printf "squeue -u jbm120 --start\n"
  squeue -u jbm120 --start
fi

printf "\n"
printf "All done!\n"
