#!/bin/bash

#SBATCH --exclusive
#SBATCH --output=slurm.out
#SBATCH -p RT_study
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --comment "ionwake"

mpirun -np 2 ./a.out
