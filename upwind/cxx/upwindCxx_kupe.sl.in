#!/bin/bash

#SBATCH -J upwindCxx
#SBATCH --partition=NeSI
#SBATCH --account=nesi99999
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --nodes=1

export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

exe="@CMAKE_BINARY_DIR@/upwind/cxx/upwindCxx"
time srun --hint=nomultithread $exe -numCells 800 -numSteps 10



