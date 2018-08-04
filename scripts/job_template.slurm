#!/bin/bash
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=72:00:00
#SBATCH --mem=60g
#SBATCH -A gds6
#SBATCH --mail-type=STAGE_OUT
#SBATCH -o job.out
#SBATCH -J JOBNAME

JOBDIR="job-$SLURM_JOBID"
JOBPATH="$TMPDIR/$JOBDIR"
mkdir -p $JOBPATH
cp -r cosmo $JOBPATH
cp -r config.txt $JOBPATH
cp -r nside_* $JOBPATH
cd $JOBPATH

rm job.out
echo "Working in: " . $JOBPATH

./cosmo config.txt

cd $TMPDIR
cp -r $JOBDIR $SLURM_SUBMIT_DIR
