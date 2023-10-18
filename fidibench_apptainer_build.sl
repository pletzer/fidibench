#!/bin/bash -e
#SBATCH --job-name=fidibench_build
#SBATCH --partition=milan
#SBATCH --time=0-02:00:00
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=2
# load environment module
module purge
module load Apptainer/1.2.2
# recent Apptainer modules set APPTAINER_BIND, which typically breaks
# container builds, so unset it here
unset APPTAINER_BIND
# create a build and cache directory on nobackup storage
export APPTAINER_CACHEDIR="/nesi/nobackup//pletzera/apptainer_cache"
export APPTAINER_TMPDIR="/nesi/nobackup//pletzera/apptainer_tmpdir"
mkdir -p $APPTAINER_CACHEDIR $APPTAINER_TMPDIR 
setfacl -b $APPTAINER_TMPDIR 
# build the container
apptainer build --force --fakeroot fidibench.sif fidibench.def
