#!/bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --partition=Debug

#export PAT_RT_PERFCNTR=2

# to see the comm pattern between ranks
export PAT_RT_SUMMARY=0

srun ./upwindMpiCxx+apa -numCells 800 -numSteps 10 


