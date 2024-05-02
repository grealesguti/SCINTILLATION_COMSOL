#!/bin/sh
#SBATCH --account=innovation
#SBATCH --job-name="matlab_demo"
#SBATCH --partition=compute
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --output=hello_world_mpi.%j.out
#SBATCH --error=hello_world_mpi.%j.err

srun scomsoltest.sh
