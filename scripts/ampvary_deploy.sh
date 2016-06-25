#!/bin/bash

# Check for correct syntax (2 arguments)
if [ $# -ne 2 ]; then
    echo $0: usage: ./test.sh RESOLUTION USE_LAMBDA
    exit 1
fi

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

# Amplitudes to run for:
declare -a AMPLITUDES=("40.00" "20.00" "10.00" "05.00" "02.50" "01.25")

# simulation step information
STEPS=$((12000*$RES/128 + 1))
FLIP_STEP=$((6000*$RES/128))
IO3D=$((3000*$RES/128))
IO2D=$((1000*$RES/128))
IO1D=$((100*$RES/128))

# Job directory
if "$LAMBDA" ; then
  JOBDIR="ampvary_lambda_$RES"
else
  JOBDIR="ampvary_$RES"
fi

echo "Deploying runs with:"
echo "  Res = $RES"
echo "  Lambda = $LAMBDA"
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
sed -i.bak "s/define N [0-9]\+/define N $RES/" ../cosmo_macros.h
rm -rf CMake*
cmake -DCMAKE_CXX_COMPILER=g++ ..
if [ $? -ne 0 ]; then
  echo "  Unable to compile - Aborting."
  rm ../cosmo_macros.h
  mv ../cosmo_macros.h.bak ../cosmo_macros.h
  exit 1
fi
make
if [ $? -ne 0 ]; then
  echo "  Unable to compile - Aborting."
  rm ../cosmo_macros.h
  mv ../cosmo_macros.h.bak ../cosmo_macros.h
  exit 1
fi
rm ../cosmo_macros.h
mv ../cosmo_macros.h.bak ../cosmo_macros.h
mv cosmo "$JOBDIR/."
cd "$JOBDIR"

for i in "${AMPLITUDES[@]}"
do
  echo ""
  echo "Deploying run with A = $i"
   
  mkdir "$i"
  cd "$i"
  cp ../cosmo cosmo
  cp ../../../config/fiducial_config.txt config.txt
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
squeue -u jbm120
