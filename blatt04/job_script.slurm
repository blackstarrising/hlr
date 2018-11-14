#!/bin/bash

#SBATCH --partition=amd
#SBATCH -N 1
#SBATCH --output=output.out

export OMP_NUM_THREADS=$1

for i in `seq 1 5`;
do
	srun partdiff 1 2 512 2 2 500 >> threads$1.out
done
