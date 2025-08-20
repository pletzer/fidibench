#!/bin/bash -e
#SBATCH --job-name=fidibench_build
#SBATCH --time=0-02:00:00
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=2
module purge
# build the container
apptainer build --force --fakeroot fidibench_intel.sif fidibench_intel.def
