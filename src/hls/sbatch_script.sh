#!/bin/bash
#SBATCG --mem=16G
#SBATCH --job-name=hw_synthesis 
#SBATCH --output=output_hwSynthesis.txt
#SBATCH --error=error_hwSynthesis.txt
#SBATCH --time=01:00:00
#SBATCH --reservation=fri 
#SBATCH --constraint=fpga

srun bash run_hls.sh
