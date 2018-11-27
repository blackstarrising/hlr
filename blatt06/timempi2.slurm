#!/bin/bash

#SBATCH --partition=west
#SBATCH --nodes=3
#SBATCH --ntasks=24
#SBATCH --output=timempi2.out
. /etc/profile.d/modules.sh
. /etc/profile.d/wr-spack.sh
spack load --dependencies mpi

mpiexec ./timempi2
