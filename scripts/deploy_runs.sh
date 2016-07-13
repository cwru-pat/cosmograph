#!/bin/bash

RED='\033[31m'
GREEN='\033[32m'
YELLOW='\033[33m'
BLUE='\033[36m'
NC='\033[0m' # No Color

# default options
declare -a AMPLITUDES=("40.00")

USE_CLUSTER=false
DRY_RUN=false

RES=16
LAMBDA=false
LINEARIZE_FRW=false
LINEARIZE_SMALL=false

# read in options
for i in "$@"
do
  case $i in
      -h|--help)
      printf "Usage: ./deploy_runs\n"
      printf "         [-C|--cluster-run] [(-N|--resolution-N)=16]\n"
      printf "         [-L|--use-lambda] [-f|--linearize-frw]\n"
      printf "         [-s|--linearize-small] [-d|--dry-run]\n"
      printf "         [(-a|--amplitudes)=40.00[,05.00,..]]\n"
      printf "         [(-c|--config)=config/dust_fiducial.txt]\n"
      exit 0
      ;;
      -N=*|--resolution-N=*)
      RES="${i#*=}"
      shift # past argument=value
      ;;
      -C|--cluster-run)
      USE_CLUSTER=true
      shift # past argument=value
      ;;
      -L|--use-lambda)
      LAMBDA=true
      shift # past argument=value
      ;;
      -f|--linearize-frw)
      LINEARIZE_FRW=true
      shift # past argument=value
      ;;
      -s|--linearize-small)
      LINEARIZE_SMALL=true
      shift # past argument=value
      ;;
      -d|--dry-run)
      DRY_RUN=true
      shift # past argument=value
      ;;
      -a=*|--amplitudes=*)
      AMP_STR="${i#*=}"
      set -f
      AMPLITUDES=(${AMP_STR//,/ })
      shift # past argument=value
      ;;
      *)
        printf "Unrecognized option will not be used: ${i#*=}\n"
        # unknown option
      ;;
  esac
done

# Validate resolution
if (( "$RES" < 4 )) ; then
  printf "${RED}Error: resolution N should be larger than 4!${NC}\n" >&2; exit 1
fi

# set LINEARIZE_FRW CMakeflag
if "$LINEARIZE_FRW" ; then
  LINEARIZE_FRW_FLAG="-DCOSMO_EXCLUDE_SECOND_ORDER_FRW=1"
else
  LINEARIZE_FRW_FLAG=''
fi

# set LINEARIZE_SMALL CMake flag
if "$LINEARIZE_SMALL" ; then
  LINEARIZE_SMALL_FLAG="-DCOSMO_EXCLUDE_SECOND_ORDER_SMALL=1"
else
  LINEARIZE_SMALL_FLAG=''
fi

# simulation step information
STEPS=$((12000*$RES/128 + 1))
FLIP_STEP=$((6000*$RES/128))
IO3D=$((3000*$RES/128))
IO2D=$((1000*$RES/128))
IO1D=$((100*$RES/128))

# Job directory
JOBDIR="ampvary_$RES"
if "$LAMBDA" ; then
  JOBDIR="${JOBDIR}_lambda"
fi
if "$LINEARIZE_SMALL" ; then
  JOBDIR="${JOBDIR}_linearsmall"
fi
if "$LINEARIZE_FRW" ; then
  JOBDIR="${JOBDIR}_linearfrw"
fi

printf "${BLUE}Deploying runs:${NC}\n"
if "$USE_CLUSTER"; then
  printf "  Using cluster (USE_CLUSTER=$USE_CLUSTER)\n"
else
  printf "  Running locally  (USE_CLUSTER=$USE_CLUSTER)\n"
fi
if "$DRY_RUN"; then
  printf "  ${YELLOW}Performing dry run (code will not be run)${NC}\n"
fi
printf "  Res = $RES\n"
printf "  Lambda = $LAMBDA\n"
printf "  LINEARIZE_FRW = $LINEARIZE_FRW\n"
printf "  LINEARIZE_SMALL = $LINEARIZE_SMALL\n"
printf "  For amplitudes:\n"
printf "    ${AMPLITUDES[*]}\n"
printf "  Output will be in $JOBDIR\n"

if ! (( 6000*$RES % 128 == 0 )) ; then
  printf "${YELLOW}Warning: chosen resolution not directly comparable to 128^3${NC}\n" >&2;
fi

read -r -t 10 -p "Continue to deploy? Will automatically proceed in 10 seconds... [Y/n]: " response
response=${response,,}    # tolower
if ! [[ $response =~ ^(|y|yes)$ ]] ; then
  printf "${RED}Aborting deployment.${NC}\n"
  exit 1
fi
printf "${GREEN}Deploying...${NC}\n"
printf "\n"

# Load modules
 if "$USE_CLUSTER"; then
  printf "Loading Modules...\n"
  module load gcc/4.9.3
  module load git/2.4.8
  module load fftw/3.3.4
  module load hdf5/1.8.15
  module load depends
  module load cmake/3.2.2
fi

# create some dirs & build
printf "Creating Directories & Building Project...\n"
mkdir -p "../build/$JOBDIR"
cd ../build
rm -rf CMake* # clean up cmake
cmake -DCMAKE_CXX_COMPILER=g++ -DCOSMO_N=$RES -DCOSMO_USE_REFERENCE_FRW=1 "$LINEARIZE_SMALL_FLAG" "$LINEARIZE_FRW_FLAG" ..
if [ $? -ne 0 ]; then
  printf "  Unable to compile - Aborting.\n"
  exit 1
fi
make -j16
if [ $? -ne 0 ]; then
  printf "  Unable to compile - Aborting.\n"
  exit 1
fi
mv cosmo "$JOBDIR/."
cd "$JOBDIR"

for i in "${AMPLITUDES[@]}"
do
  printf "\n"
  printf "Deploying run with A = $i\n"
   
  mkdir "$i"
  cd "$i"
  cp ../cosmo cosmo
  cp ../../../config/dust_fiducial.txt config.txt

  if "$USE_CLUSTER"; then
    cp ../../../scripts/job_template.slurm job.slurm
    # Adjust job name
    sed -i.bak "s/JOBNAME/A_${i}_${JOBDIR}/" job.slurm
  fi

  # Adjust amplitude
  sed -i.bak "s/peak_amplitude_frac = 5\.0/peak_amplitude_frac = $i/" config.txt

  # Use Lambda if requested
  if "$LAMBDA" ; then
    sed -i.bak "s/dust/dust_lambda/g" config.txt
  fi

  # Adjust output/step parameters per resolution
  sed -i.bak "s/steps = 12001/steps = $STEPS/" config.txt
  sed -i.bak "s/ray_flip_step = 6000/ray_flip_step = $FLIP_STEP/" config.txt
  sed -i.bak "s/IO_3D_grid_interval = 3000/IO_3D_grid_interval = $IO3D/" config.txt
  sed -i.bak "s/IO_2D_grid_interval = 1000/IO_2D_grid_interval = $IO2D/" config.txt
  sed -i.bak "s/IO_1D_grid_interval = 100/IO_1D_grid_interval = $IO1D/" config.txt

  # Run job, go back up a dir
  if [ ! "$DRY_RUN" ]; then
    if "$USE_CLUSTER"; then
      printf "${GREEN}Queueing job A_${i}_${JOBDIR}.${NC}\n"
      #sbatch job.slurm
    else
      printf "${GREEN}Running job.${NC}\n"
      #./cosmo config.txt
    fi
  fi

  cd ..

done

if "$USE_CLUSTER"; then
  printf "squeue -u jbm120\n"
  squeue -u jbm120 --start
fi

printf "\n"
printf "All done!\n"
