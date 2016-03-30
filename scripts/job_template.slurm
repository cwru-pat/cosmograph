#!/bin/bash
#SBATCH --mail-user=jbm120@case.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=100:00:00
#SBATCH --mem=20g
#SBATCH -o job.out
#SBATCH -J JOBNAME

cp -r * $PFSDIR
cd $PFSDIR
echo $PFSDIR

./cosmo config.txt

cp -r * $SLURM_SUBMIT_DIR
