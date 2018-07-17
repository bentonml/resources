#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=8G               # Gb memory
#SBATCH --time=0-08:00:00      # runtime 
#SBATCH --mail-type=BEGIN      # email when started
#SBATCH --mail-user=username@vanderbilt.edu
#SBATCH --output=notebook.out

# when job launches, open tunnel from your local machine

PORT=7777
echo To open a tunnel from local machine, 
echo Execute in a new terminal window:
echo "ssh -L 8888:$HOSTNAME:$PORT $(whoami)@login.accre.vanderbilt.edu"
printf '=%.0s' {1..80}
echo

module load Anaconda3       # load python 
source activate sample_env  # load virtual environment

jt -t solarizedd -f anonymous -fs 12 -nfs 14 -ofs 12 -cellw 90%  # optionally set theme

jupyter notebook --no-browser --ip='*' --port=${PORT}

