#!/bin/bash

# Job name:
#SBATCH --job-name=vec_kmeans
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1  
#SBATCH --cpus-per-task=1 
#SBATCH --time=00:05:00
#SBATCH --output=test2.out
#SBATCH --constraint=amd
#SBATCH --reservation=fri
#SBATCH --propagate=STACK


FILE=dwt
#FILE=example



gcc -O2 -lm dwt.c -o dwt

srun dwt
