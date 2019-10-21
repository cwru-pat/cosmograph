#!/bin/bash
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=72:00:00
#SBATCH --mem=60g
#SBATCH -A gds6
#SBATCH -J JOBNAME
#SBATCH --exclude=compt178

JOBDIR="job-$SLURM_JOBID"
JOBPATH="$JOBDIR"
mkdir -p $JOBPATH
cp -r cosmo $JOBPATH
cp -r config.txt $JOBPATH
cp -r nside_* $JOBPATH
cd $JOBPATH

rm job.out
echo "Working in: " . $JOBPATH

./cosmo config.txt
