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

RES=64
FLIP_DELAY=0
POWER_CUT=1
NSIDE=16
SEED=7

# read in options
for i in "$@"
do
  case $i in
      -h|--help)
      printf "Usage: ./dust_lensing.sh\n"
      printf "         [-C|--cluster-run] [(-N|--resolution-N)=64] [-d|--dry-run]\n"
      printf "         [(-f|--flip-delay)=0] [(-k|--k)=1] [(-n|-nside)=16] [(-s|-seed)=7]\n"
      exit 0
      ;;
      -N=*|--resolution-N=*)
      RES="${i#*=}"
      shift # past argument=value
      ;;
      -f=*|--flip-delay=*)
      FLIP_DELAY="${i#*=}"
      shift # past argument=value
      ;;
      -k=*|--k=*)
      POWER_CUT="${i#*=}"
      shift # past argument=value
      ;;
      -C|--cluster-run)
      USE_CLUSTER=true
      shift # past argument
      ;;
      -d|--dry-run)
      DRY_RUN=true
      shift # past argument
      ;;
      -s=*|--seed=*)
      SEED="${i#*=}"
      shift # past argument=value
      ;;
      -n=*|--nside=*)
      NSIDE="${i#*=}"
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

# simulation step information
STEPS=$((16000*$RES/64 + 1))
FLIP_STEP=$(((8000 + $FLIP_DELAY)*$RES/64))
IO3D=$((2000*$RES/64))

# Job directory
JOBDIR="proper_time_R-${RES}_F-${FLIP_STEP}_kcut-${POWER_CUT}_Nside-${NSIDE}_S-${SEED}"

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
printf "  Output will be in $JOBDIR\n"

if ! (( 16000*$RES % 64 == 0 )) ; then
  printf "${YELLOW}Warning: chosen resolution not directly comparable to 64^3${NC}\n" >&2;
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
  module load cmake
fi

# create some dirs & build
printf "Creating Directories & Building Project...\n"
mkdir -p "../build/$JOBDIR"
cd ../build
rm -rf CMake* # clean up cmake
cmake -DCOSMO_N=$RES .. # -DCOSMO_USE_Z4c_DAMPING=0
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

cp ../../config/proper_time.txt config.txt
cp ../../config/healpix_vecs/nside_$NSIDE.vecs nside_$NSIDE.vecs
sed -i.bak "s/healpix_vecs_file = ..\/config\/healpix_vecs\//healpix_vecs_file = /" config.txt

if "$USE_CLUSTER"; then
  cp ../../scripts/job_template.slurm job.slurm
  # Adjust job name
  sed -i.bak "s/JOBNAME/${JOBDIR}/" job.slurm
  JOBTIME=$((10*$RES/32*$RES/32*$RES/32*$RES/32))
  sed -i.bak "s/72:00:00/${JOBTIME}:00:00/" job.slurm
fi

sed -i.bak -r "s/mt19937_seed = [\.0-9]+/mt19937_seed = $SEED/" config.txt
sed -i.bak -r "s/ic_spec_cut = [\.0-9]+/ic_spec_cut = $POWER_CUT/" config.txt
sed -i.bak -r "s/healpix_vecs_file = nside_[0-9]+\.vecs/healpix_vecs_file = nside_$NSIDE\.vecs/" config.txt

# Adjust output/step parameters per resolution
sed -i.bak -r "s/steps = [\.0-9]+/steps = $STEPS/" config.txt
sed -i.bak -r "s/ray_flip_step = [\.0-9]+/ray_flip_step = $FLIP_STEP/" config.txt
sed -i.bak -r "s/IO_3D_grid_interval = [\.0-9]+/IO_3D_grid_interval = $IO3D/" config.txt

# Run job, go back up a dir
if [ "$DRY_RUN" = false ]; then
  if "$USE_CLUSTER"; then
    printf "${GREEN}Queueing job A_${i}_${JOBDIR}.${NC}\n"
    sbatch job.slurm
  else
    printf "${GREEN}Running job.${NC}\n"
    ./cosmo config.txt
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
