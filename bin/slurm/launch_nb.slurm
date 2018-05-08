#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=2G               # 2Gb memory
#SBATCH --time=0-04:00:00      # 4 hours
#SBATCH --output=notebook.out

# modified from: https://github.com/accre/Python/blob/master/notebooks/tunnel.slurm

PORT=7777
echo To open a tunnel from local machine, 
echo Execute in a new terminal window:
echo "ssh -L 9999:$HOSTNAME:$PORT $(whoami)@login.accre.vanderbilt.edu"
printf '=%.0s' {1..80}
echo

module load GCC Anaconda3   # load python 
source activate capra_py3   # load virtual environment

jt -t solarizedd -f anonymous -ofs 11 -cellw 90%  # optionally set theme

jupyter notebook --no-browser --ip='*' --port=${PORT}

