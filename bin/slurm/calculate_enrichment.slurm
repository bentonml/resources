#!/bin/bash
#SBATCH --mem=2G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --tasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --output=calculate_enrichment.out

# load modules
module load Anaconda3
source activate enh_gain-loss  # can clone from ../resources/conda_envs/enh_gain-loss.yml

source /accre/usr/bin/setup_accre_runtime_dir  # IMPORTANT to remove .tmp files even after failure 

python ./calculate_enrichment.py  # potentially need to change this path

