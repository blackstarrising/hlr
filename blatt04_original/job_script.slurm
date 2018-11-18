#!/bin/bash

#SBATCH --partition=west
#SBATCH -N 1
#SBATCH -w west10
#SBATCH --output=output.out

srun partdiff 1 2 512 2 2 500 >> thread.out
