#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=300:00:00
#SBATCH --mem=20g
#SBATCH -A gds6
#SBATCH -o job.out
#SBATCH -J JOBNAME

cp -r * $PFSDIR
cd $PFSDIR
rm job.out
echo "Working in: " . $PFSDIR

./cosmo config.txt

cp -r * $SLURM_SUBMIT_DIR
