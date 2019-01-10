#!/bin/bash

#SBATCH --partition=west
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --output=output/output.out
rm output/*

. /etc/profile.d/modules.sh
. /etc/profile.d/wr-spack.sh
spack load --dependencies mpi

mpiexec ./partdiff 1 1 64 2 2 120 >> output/thread.out
