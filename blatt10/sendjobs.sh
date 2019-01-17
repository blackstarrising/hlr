#!/bin/bash
for slurmfile in slurm/*; do
	sbatch $slurmfile
done
