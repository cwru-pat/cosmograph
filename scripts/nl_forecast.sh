#!/bin/bash

# Switch to the directory containing this script,
cd "$(dirname "$0")"
# up a directory should be the main codebase.
module load gsl

mkdir -p ../build/nl_runs
cd ../build/nl_runs

do_run () {
  RES=$1
  P=$2
  RESCALE=$3
  NSIDE=$4
  MTSEED=$5
  LAPSE=$6
  LinRaysheet=$7
  NPIX=$((12*$NSIDE*$NSIDE))

  if [[ "$LAPSE" == "Harmonic" ]]; then
    STEPS=$((144*2*$RES/8))
    FLIP_STEP=$((144*$RES/8))
  elif [[ "$LAPSE" == "ConformalFLRW" ]]; then
    STEPS=$((3672*2*$RES/96))
    FLIP_STEP=$((3672*$RES/96))
  else
    echo "Use Harmonic or ConformalFLRW gauges."
    exit 1
  fi

  cmake -DCOSMO_N=${RES} -DCOSMO_H_LEN_FRAC=1.9 ../..
  make -j8

  RUN_DIR="v1_r${RES}_P${P}_M${MTSEED}_resc${RESCALE}_${LAPSE}"
  if [ "$LinHydroSrc" == "0" ] ; then
    RUN_DIR="${RUN_DIR}_nlsrc"
  fi
  mkdir -p $RUN_DIR
  cp cosmo "$RUN_DIR/."
  cp ../../config/dust_lensing.txt "$RUN_DIR/config.txt"
  cp ../../config/healpix_vecs/nside_${NSIDE}.vecs "$RUN_DIR/nside.vecs"
  cd $RUN_DIR

  sed -i.bak -r "s/P = [\.0-9]+/P = ${P}/" config.txt
  sed -i.bak -r "s/rescale_metric = [\.0-9]+/rescale_metric = ${RESCALE}/" config.txt
  sed -i.bak -r "s/rescale_sheet = [\.0-9]+/rescale_sheet = ${LinRaysheet}/" config.txt
  sed -i.bak -r "s/raysheet_flip_step = [0-9]+/raysheet_flip_step = ${FLIP_STEP}/" config.txt
  sed -i.bak -r "s/steps = [0-9]+/steps = ${STEPS}/" config.txt
  sed -i.bak -r "s/healpix_vecs_file = \.\.\/config\/healpix_vecs\/nside_16\.vecs/healpix_vecs_file = nside\.vecs/" config.txt
  sed -i.bak -r "s/ns1 = [0-9]+/ns1 = ${NPIX}/" config.txt
  sed -i.bak -r "s/nside = [0-9]+/nside = ${NSIDE}/" config.txt
  sed -i.bak -r "s/mt19937_seed = [0-9]+/mt19937_seed = ${MTSEED}/" config.txt
  sed -i.bak -r "s/lapse = [A-Za-z0-9]+/lapse = ${LAPSE}/" config.txt

  sbatch --exclusive -N 1 --mem=90g -n 24 --time=200:00:00 --wrap="./cosmo config.txt"

  cd ..
}

for ((i=1; i<20; i+=1));
do
  MTS=$i
  
  # do_run 96 0.13 1.0 32 "$MTS" Harmonic 0 # full nonlinear
  # do_run 96 0.13 1.0e-3 32 "$MTS" Harmonic 0 # linearize metric only
  # do_run 96 0.00000013 1.0 32 "$MTS" Harmonic 0 # linearize all
  do_run 96 0.13 1.0e-3 32 "$MTS" Harmonic 1.0e-3 # linearize metric & raytrace
  
  # do_run 96 0.13 1.0 32 "$MTS" ConformalFLRW 1
  # do_run 96 0.13 1.0e-3 32 "$MTS" ConformalFLRW 1
done

# do_run 128 0.13 1.0e-3 32 1 Harmonic 0
# do_run 128 0.13 1.0e-3 32 1 Harmonic 1
# do_run 128 0.13 1.0 32 1 Harmonic 1

# do_run 96 0.13 1.0 32 1 
# do_run 96 0.13 1.0e-3 32 1

# do_run 128 0.13 1.0 32 1
# do_run 128 0.13 1.0e-3 32 1

# do_run 144 0.13 1.0 32 1
# do_run 144 0.13 1.0e-3 32 1

# do_run 192 0.13 1.0 32 1
# do_run 192 0.13 1.0e-3 32 1

# do_run 256 0.13 1.0 32 1
# do_run 256 0.13 1.0e-3 32 1
