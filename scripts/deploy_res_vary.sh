#!/bin/bash

# Load modules
module load gcc/4.9.3
module load git/2.4.8
module load fftw/3.3.4
module load hdf5/1.8.15
module load depends
module load cmake/3.2.2

# build & create some dirs
cd ../build
mkdir -p resvary


## -- 64^3 = (2/4 * 128)^3

# Make dir for run
mkdir -p resvary/R64
# Change resolution in cosmo_macros.h; compile; move executable
sed -i.bak 's#N 128#N 64#g' ../cosmo_macros.h
rm -rf CMake* && cmake -DCMAKE_CXX_COMPILER=g++ .. && make
mv cosmo resvary/R64
# "undo" changes
rm ../cosmo_macros.h
mv ../cosmo_macros.h.bak ../cosmo_macros.h
# Set up config file, job script
cp ../config/fiducial_config.txt resvary/R64/config.txt
sed -i.bak 's#steps = 12001#steps = 6001#g' resvary/R64/config.txt
sed -i.bak 's#ray_flip_step = 6000#ray_flip_step = 3000#g' resvary/R64/config.txt
sed -i.bak 's#IO_3D_grid_interval = 3000#IO_3D_grid_interval = 1500#g' resvary/R64/config.txt
sed -i.bak 's#IO_2D_grid_interval = 1000#IO_2D_grid_interval = 500#g' resvary/R64/config.txt
sed -i.bak 's#IO_1D_grid_interval = 100#IO_1D_grid_interval = 50#g' resvary/R64/config.txt
cp ../scripts/job_template.slurm resvary/R64/job.slurm
sed -i.bak 's#JOBNAME#R64#g' resvary/R64/job.slurm


## -- 96^3 = (3/4 * 128)^3

# Make dir for run
mkdir -p resvary/R96
# Change resolution in cosmo_macros.h; compile; move executable
sed -i.bak 's#N 128#N 96#g' ../cosmo_macros.h
rm -rf CMake* && cmake -DCMAKE_CXX_COMPILER=g++ .. && make
mv cosmo resvary/R96
# "undo" changes
rm ../cosmo_macros.h
mv ../cosmo_macros.h.bak ../cosmo_macros.h
# Set up config file, job script
cp ../config/fiducial_config.txt resvary/R96/config.txt
sed -i.bak 's#steps = 12001#steps = 9001#g' resvary/R96/config.txt
sed -i.bak 's#ray_flip_step = 6000#ray_flip_step = 4500#g' resvary/R96/config.txt
sed -i.bak 's#IO_3D_grid_interval = 3000#IO_3D_grid_interval = 2250#g' resvary/R96/config.txt
sed -i.bak 's#IO_2D_grid_interval = 1000#IO_2D_grid_interval = 750#g' resvary/R96/config.txt
sed -i.bak 's#IO_1D_grid_interval = 100#IO_1D_grid_interval = 75#g' resvary/R96/config.txt
cp ../scripts/job_template.slurm resvary/R96/job.slurm
sed -i.bak 's#JOBNAME#R96#g' resvary/R96/job.slurm


## -- 128^3

# Make dir for run
mkdir -p resvary/R128
# Change resolution in cosmo_macros.h; compile; move executable
sed -i.bak 's#N 128#N 128#g' ../cosmo_macros.h
rm -rf CMake* && cmake -DCMAKE_CXX_COMPILER=g++ .. && make
mv cosmo resvary/R128
# "undo" changes
rm ../cosmo_macros.h
mv ../cosmo_macros.h.bak ../cosmo_macros.h
# Set up config file, job script
cp ../config/fiducial_config.txt resvary/R128/config.txt
sed -i.bak 's#steps = 12001#steps = 12001#g' resvary/R128/config.txt
sed -i.bak 's#ray_flip_step = 6000#ray_flip_step = 6000#g' resvary/R128/config.txt
sed -i.bak 's#IO_3D_grid_interval = 3000#IO_3D_grid_interval = 3000#g' resvary/R128/config.txt
sed -i.bak 's#IO_2D_grid_interval = 1000#IO_2D_grid_interval = 1000#g' resvary/R128/config.txt
sed -i.bak 's#IO_1D_grid_interval = 100#IO_1D_grid_interval = 100#g' resvary/R128/config.txt
cp ../scripts/job_template.slurm resvary/R128/job.slurm
sed -i.bak 's#JOBNAME#R128#g' resvary/R128/job.slurm


## -- 160^3 = (5/4 * 128)^3

# Make dir for run
mkdir -p resvary/R160
# Change resolution in cosmo_macros.h; compile; move executable
sed -i.bak 's#N 128#N 160#g' ../cosmo_macros.h
rm -rf CMake* && cmake -DCMAKE_CXX_COMPILER=g++ .. && make
mv cosmo resvary/R160
# "undo" changes
rm ../cosmo_macros.h
mv ../cosmo_macros.h.bak ../cosmo_macros.h
# Set up config file, job script
cp ../config/fiducial_config.txt resvary/R160/config.txt
sed -i.bak 's#steps = 12001#steps = 15001#g' resvary/R160/config.txt
sed -i.bak 's#ray_flip_step = 6000#ray_flip_step = 7500#g' resvary/R160/config.txt
sed -i.bak 's#IO_3D_grid_interval = 3000#IO_3D_grid_interval = 3750#g' resvary/R160/config.txt
sed -i.bak 's#IO_2D_grid_interval = 1000#IO_2D_grid_interval = 1250#g' resvary/R160/config.txt
sed -i.bak 's#IO_1D_grid_interval = 100#IO_1D_grid_interval = 125#g' resvary/R160/config.txt
cp ../scripts/job_template.slurm resvary/R160/job.slurm
sed -i.bak 's#JOBNAME#R160#g' resvary/R160/job.slurm
sed -i.bak 's#100:00:00#200:00:00#g' resvary/R160/job.slurm


# Execute
(cd resvary/R64 && sbatch job.slurm)
(cd resvary/R96 && sbatch job.slurm)
(cd resvary/R128 && sbatch job.slurm)
(cd resvary/R160 && sbatch job.slurm)

