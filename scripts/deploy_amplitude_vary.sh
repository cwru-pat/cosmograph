#!/bin/bash

# Load modules
module load gcc/4.9.3
module load git/2.4.8
module load fftw/3.3.4
module load hdf5/1.8.15
module load depends
module load cmake/3.2.2

# build & create some dirs
mkdir -p ../build/ampvary
cd ../build
cmake -DCMAKE_CXX_COMPILER=g++ .. && make
mkdir -p ampvary
cd ampvary

# Set up individual runs

mkdir A40.00 &&  cp ../cosmo A40.00/cosmo
cp ../../config/fiducial_config.txt A40.00/config.txt
cp ../../scripts/job_template.slurm A40.00/job.slurm
sed -i.bak 's#peak_amplitude_frac = 5\.0#peak_amplitude_frac = 40\.00#g' A40.00/config.txt
sed -i.bak 's#JOBNAME#A_40\.00#g' A40.00/job.slurm

mkdir A20.00 &&  cp ../cosmo A20.00/cosmo
cp ../../config/fiducial_config.txt A20.00/config.txt
cp ../../scripts/job_template.slurm A20.00/job.slurm
sed -i.bak 's#peak_amplitude_frac = 5\.0#peak_amplitude_frac = 20\.00#g' A20.00/config.txt
sed -i.bak 's#JOBNAME#A_20\.00#g' A20.00/job.slurm

mkdir A10.00 &&  cp ../cosmo A10.00/cosmo
cp ../../config/fiducial_config.txt A10.00/config.txt
cp ../../scripts/job_template.slurm A10.00/job.slurm
sed -i.bak 's#peak_amplitude_frac = 5\.0#peak_amplitude_frac = 10\.00#g' A10.00/config.txt
sed -i.bak 's#JOBNAME#A_10\.00#g' A10.00/job.slurm

mkdir A02.50 &&  cp ../cosmo A02.50/cosmo
cp ../../config/fiducial_config.txt A02.50/config.txt
cp ../../scripts/job_template.slurm A02.50/job.slurm
sed -i.bak 's#peak_amplitude_frac = 5\.0#peak_amplitude_frac = 02\.50#g' A02.50/config.txt
sed -i.bak 's#JOBNAME#A_02\.50#g' A02.50/job.slurm

mkdir A01.25 &&  cp ../cosmo A01.25/cosmo
cp ../../config/fiducial_config.txt A01.25/config.txt
cp ../../scripts/job_template.slurm A01.25/job.slurm
sed -i.bak 's#peak_amplitude_frac = 5\.0#peak_amplitude_frac = 01\.25#g' A01.25/config.txt
sed -i.bak 's#JOBNAME#A_01\.25#g' A01.25/job.slurm

