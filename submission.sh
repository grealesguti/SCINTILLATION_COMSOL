#!/bin/sh
#SBATCH --account=innovation
#SBATCH --job-name="matlab_demo"
#SBATCH --partition=compute
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --output=OUT/hello_world_mpi.%j.out
#SBATCH --error=ERR/hello_world_mpi.%j.err

srun scomsoltest.sh
