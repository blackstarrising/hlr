#!/bin/bash

#SBATCH --partition=amd
#SBATCH -N 1
#SBATCH --output=output.out

export OMP_NUM_THREADS=12

for i in `seq 1 5`;
do
	srun partdiff 1 2 $1  2 2 3162 >> interlines$1.out
done
