#!/bin/bash

# Check for correct syntax (4 arguments)
if [ $# -ne 4 ]; then
    echo "$0: usage: ./deploy_runs.sh RESOLUTION USE_LAMBDA LINEARIZE_FRW LINEARIZE_SMALL"
    exit 1
fi

# Amplitudes to run for:
declare -a AMPLITUDES=("40.00")

# Validate resolution
RES=$1
if (( "$RES" < 8 )) ; then
  echo "error: RESOLUTION should be between 8 and 512 " >&2; exit 1
fi
if (( "$RES" > 512 )) ; then
  echo "error: RESOLUTION should be between 8 and 512 " >&2; exit 1
fi

# Validate using Lambda
LAMBDA=$2
if (( $LAMBDA > 0 )) ; then
  LAMBDA=true
else
  LAMBDA=false
fi

# Validate LINEARIZE_FRW
LINEARIZE_FRW=$3
if (( $LINEARIZE_FRW > 0 )) ; then
  LINEARIZE_FRW=true
  LINEARIZE_FRW_FLAG="-DCOSMO_EXCLUDE_SECOND_ORDER_FRW=1"
else
  LINEARIZE_FRW=false
  LINEARIZE_FRW_FLAG=''
fi

# Validate LINEARIZE_SMALL
LINEARIZE_SMALL=$4
if (( $LINEARIZE_SMALL > 0 )) ; then
  LINEARIZE_SMALL=true
  LINEARIZE_SMALL_FLAG="-DCOSMO_EXCLUDE_SECOND_ORDER_SMALL=1"
else
  LINEARIZE_SMALL=false
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

echo "Deploying runs with:"
echo "  Res = $RES"
echo "  Lambda = $LAMBDA"
echo "  LINEARIZE_FRW = $LINEARIZE_FRW"
echo "  LINEARIZE_SMALL = $LINEARIZE_SMALL"
echo "  For amplitudes:"
echo "    ${AMPLITUDES[*]}"
echo "  Will output to $JOBDIR"

if ! (( 6000*$RES % 128 == 0 )) ; then
  echo "Warning: chosen resolution not directly comparable to 128^3 " >&2; exit 1
fi

read -r -t 10 -p "Continue to deploy? Will automatically proceed in 10 seconds... [Y/n]: " response
response=${response,,}    # tolower
if ! [[ $response =~ ^(|y|yes)$ ]] ; then
  echo "Aborting deployment."
  exit 1
fi
echo "Deploying..."
echo ""

# Load modules
echo "Loading Modules..."
module load gcc/4.9.3
module load git/2.4.8
module load fftw/3.3.4
module load hdf5/1.8.15
module load depends
module load cmake/3.2.2

# create some dirs & build
echo "Creating Directories & Building Project..."
mkdir -p "../build/$JOBDIR"
cd ../build
rm -rf CMake*
cmake -DCMAKE_CXX_COMPILER=g++ -DCOSMO_N=$RES -DCOSMO_USE_REFERENCE_FRW=1 "$LINEARIZE_SMALL_FLAG" "$LINEARIZE_FRW_FLAG" ..
if [ $? -ne 0 ]; then
  echo "  Unable to compile - Aborting."
  exit 1
fi
make -j16
if [ $? -ne 0 ]; then
  echo "  Unable to compile - Aborting."
  exit 1
fi
mv cosmo "$JOBDIR/."
cd "$JOBDIR"

for i in "${AMPLITUDES[@]}"
do
  echo ""
  echo "Deploying run with A = $i"
   
  mkdir "$i"
  cd "$i"
  cp ../cosmo cosmo
  cp ../../../config/dust_fiducial.txt config.txt
  cp ../../../scripts/job_template.slurm job.slurm

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

  # Adjust job name
  sed -i.bak "s/JOBNAME/A_${i}_${JOBDIR}/" job.slurm

  # Run job, go back up a dir
  sbatch job.slurm
  cd ..

done

echo ""
echo "All done!"
echo "squeue -u jbm120 is"
squeue -u jbm120 --start
