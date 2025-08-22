#!/bin/bash
#SBATCH --job-name=fidibench                                                                                
#SBATCH --time=00:05:00       # Walltime (HH:MM:SS)                                                                                                     
#SBATCH --hint=nomultithread                                                                                                                            
#SBATCH --ntasks=8           # number of tasks (e.g. MPI)                                                                                              
#SBATCH --cpus-per-task=1     # number of cores per task (e.g. OpenMP)                                                                                  
#SBATCH --nodes=2                                                                                                                                       
module purge
export I_MPI_ROOT=/opt/nesi/CS400_centos7_bdw/impi/2021.5.1-intel-compilers-2022.0.2/mpi/2021.5.1
export PATH=$I_MPI_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$I_MPI_ROOT/lib:$LD_LIBRARY_PATH
mpiexec -np ${SLURM_NTASKS} --bind-to none --map-by slot \
	apptainer exec --bind $I_MPI_ROOT:$I_MPI_ROOT fidibench_intel.aif /software/fidibench/bin/upwindMpiCxx

