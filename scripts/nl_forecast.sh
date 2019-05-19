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
  NSIDE=$3
  MTSEED=$4
  NPIX=$((12*$NSIDE*$NSIDE))
  STEPS=$((3672*2*$RES/96))
  FLIP_STEP=$((3672*$RES/96))  

  cmake -DCOSMO_N=${RES} -DCOSMO_H_LEN_FRAC=1.9 ../..
  make -j8

  RUN_DIR="r${RES}_P${P}_M${MTSEED}"
  mkdir -p $RUN_DIR
  cp cosmo "$RUN_DIR/."
  cp ../../config/dust_lensing.txt "$RUN_DIR/config.txt"
  cp ../../config/healpix_vecs/nside_${NSIDE}.vecs "$RUN_DIR/nside.vecs"
  cd $RUN_DIR

  sed -i.bak -r "s/P = [\.0-9]+/P = ${P}/" config.txt
  sed -i.bak -r "s/raysheet_flip_step = [0-9]+/raysheet_flip_step = ${FLIP_STEP}/" config.txt
  sed -i.bak -r "s/steps = [0-9]+/steps = ${STEPS}/" config.txt
  sed -i.bak -r "s/healpix_vecs_file = \.\.\/config\/healpix_vecs\/nside_16\.vecs/healpix_vecs_file = nside\.vecs/" config.txt
  sed -i.bak -r "s/ns1 = [0-9]+/ns1 = ${NPIX}/" config.txt
  sed -i.bak -r "s/nside = [0-9]+/nside = ${NSIDE}/" config.txt
  sed -i.bak -r "s/mt19937_seed = [0-9]+/mt19937_seed = ${MTSEED}/" config.txt

  sbatch --exclusive -N 1 --mem=60g -n 24 --time=200:00:00 --wrap="./cosmo config.txt"

  cd ..
}

for ((i=1; i<10; i+=1));
do
  MTS=$i
  do_run 96 0.13 32 "$MTS"
  do_run 96 0.00000013 32 "$MTS"
done

do_run 128 0.13 32 1
do_run 128 0.00000013 32 1

do_run 144 0.13 32 1
do_run 144 0.00000013 32 1

do_run 192 0.13 32 1
do_run 192 0.00000013 32 1

# do_run 256 0.13 32 1
# do_run 256 0.00000013 32 1
