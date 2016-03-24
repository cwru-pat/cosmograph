#!/bin/bash

# Load modules
module load gcc/4.9.3
module load git/2.4.8
module load fftw/3.3.4
module load hdf5/1.8.15
module load depends
module load cmake/3.2.2

# build & create some dirs
mkdir -p ../build/speccut
cd ../build
cmake -DCMAKE_CXX_COMPILER=g++ .. && make
mkdir -p speccut
cd speccut

# Set up individual runs

mkdir C7 &&  cp ../cosmo C7/cosmo
cp ../../config/fiducial_config.txt C7/config.txt
cp ../../scripts/job_template.slurm C7/job.slurm
sed -i.bak 's#ic_spec_cut = 10#ic_spec_cut = 7#g' C7/config.txt
sed -i.bak 's#JOBNAME#C_7#g' C7/job.slurm

mkdir C8 &&  cp ../cosmo C8/cosmo
cp ../../config/fiducial_config.txt C8/config.txt
cp ../../scripts/job_template.slurm C8/job.slurm
sed -i.bak 's#ic_spec_cut = 10#ic_spec_cut = 8#g' C8/config.txt
sed -i.bak 's#JOBNAME#C_8#g' C8/job.slurm

mkdir C9 &&  cp ../cosmo C9/cosmo
cp ../../config/fiducial_config.txt C9/config.txt
cp ../../scripts/job_template.slurm C9/job.slurm
sed -i.bak 's#ic_spec_cut = 10#ic_spec_cut = 9#g' C9/config.txt
sed -i.bak 's#JOBNAME#C_9#g' C9/job.slurm

mkdir C11 &&  cp ../cosmo C11/cosmo
cp ../../config/fiducial_config.txt C11/config.txt
cp ../../scripts/job_template.slurm C11/job.slurm
sed -i.bak 's#ic_spec_cut = 10#ic_spec_cut = 11#g' C11/config.txt
sed -i.bak 's#JOBNAME#C_11#g' C11/job.slurm

mkdir C12 &&  cp ../cosmo C12/cosmo
cp ../../config/fiducial_config.txt C12/config.txt
cp ../../scripts/job_template.slurm C12/job.slurm
sed -i.bak 's#ic_spec_cut = 10#ic_spec_cut = 12#g' C12/config.txt
sed -i.bak 's#JOBNAME#C_12#g' C12/job.slurm

mkdir C13 &&  cp ../cosmo C13/cosmo
cp ../../config/fiducial_config.txt C13/config.txt
cp ../../scripts/job_template.slurm C13/job.slurm
sed -i.bak 's#ic_spec_cut = 10#ic_spec_cut = 13#g' C13/config.txt
sed -i.bak 's#JOBNAME#C_13#g' C13/job.slurm

mkdir C14 &&  cp ../cosmo C14/cosmo
cp ../../config/fiducial_config.txt C14/config.txt
cp ../../scripts/job_template.slurm C14/job.slurm
sed -i.bak 's#ic_spec_cut = 10#ic_spec_cut = 14#g' C14/config.txt
sed -i.bak 's#JOBNAME#C_14#g' C14/job.slurm

mkdir C15 &&  cp ../cosmo C15/cosmo
cp ../../config/fiducial_config.txt C15/config.txt
cp ../../scripts/job_template.slurm C15/job.slurm
sed -i.bak 's#ic_spec_cut = 10#ic_spec_cut = 15#g' C15/config.txt
sed -i.bak 's#JOBNAME#C_15#g' C15/job.slurm

mkdir C20 &&  cp ../cosmo C20/cosmo
cp ../../config/fiducial_config.txt C20/config.txt
cp ../../scripts/job_template.slurm C20/job.slurm
sed -i.bak 's#ic_spec_cut = 10#ic_spec_cut = 20#g' C20/config.txt
sed -i.bak 's#JOBNAME#C_20#g' C20/job.slurm

mkdir C25 &&  cp ../cosmo C25/cosmo
cp ../../config/fiducial_config.txt C25/config.txt
cp ../../scripts/job_template.slurm C25/job.slurm
sed -i.bak 's#ic_spec_cut = 10#ic_spec_cut = 25#g' C25/config.txt
sed -i.bak 's#JOBNAME#C_25#g' C25/job.slurm

mkdir C30 &&  cp ../cosmo C30/cosmo
cp ../../config/fiducial_config.txt C30/config.txt
cp ../../scripts/job_template.slurm C30/job.slurm
sed -i.bak 's#ic_spec_cut = 10#ic_spec_cut = 30#g' C30/config.txt
sed -i.bak 's#JOBNAME#C_30#g' C30/job.slurm

mkdir C40 &&  cp ../cosmo C40/cosmo
cp ../../config/fiducial_config.txt C40/config.txt
cp ../../scripts/job_template.slurm C40/job.slurm
sed -i.bak 's#ic_spec_cut = 10#ic_spec_cut = 40#g' C40/config.txt
sed -i.bak 's#JOBNAME#C_40#g' C40/job.slurm

mkdir C50 &&  cp ../cosmo C50/cosmo
cp ../../config/fiducial_config.txt C50/config.txt
cp ../../scripts/job_template.slurm C50/job.slurm
sed -i.bak 's#ic_spec_cut = 10#ic_spec_cut = 50#g' C50/config.txt
sed -i.bak 's#JOBNAME#C_50#g' C50/job.slurm

mkdir C100 &&  cp ../cosmo C100/cosmo
cp ../../config/fiducial_config.txt C100/config.txt
cp ../../scripts/job_template.slurm C100/job.slurm
sed -i.bak 's#ic_spec_cut = 10#ic_spec_cut = 100#g' C100/config.txt
sed -i.bak 's#JOBNAME#C_100#g' C100/job.slurm

mkdir C200 &&  cp ../cosmo C200/cosmo
cp ../../config/fiducial_config.txt C200/config.txt
cp ../../scripts/job_template.slurm C200/job.slurm
sed -i.bak 's#ic_spec_cut = 10#ic_spec_cut = 200#g' C200/config.txt
sed -i.bak 's#JOBNAME#C_200#g' C200/job.slurm

(cd C7 && sbatch job.slurm)
(cd C8 && sbatch job.slurm)
(cd C9 && sbatch job.slurm)
(cd C11 && sbatch job.slurm)
(cd C12 && sbatch job.slurm)
(cd C13 && sbatch job.slurm)
(cd C14 && sbatch job.slurm)
(cd C15 && sbatch job.slurm)
(cd C20 && sbatch job.slurm)
(cd C25 && sbatch job.slurm)
(cd C30 && sbatch job.slurm)
(cd C40 && sbatch job.slurm)
(cd C50 && sbatch job.slurm)
(cd C100 && sbatch job.slurm)
(cd C200 && sbatch job.slurm)
