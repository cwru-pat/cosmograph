#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=72:00:00
#SBATCH --mem=28g
#SBATCH -A gds6
#SBATCH --mail-type=STAGE_OUT
#SBATCH -o job.out
#SBATCH -J JOBNAME

cp -r * $TMPDIR
cd $TMPDIR
rm job.out
echo "Working in: " . $TMPDIR

./cosmo config.txt

cp -r * $SLURM_SUBMIT_DIR
