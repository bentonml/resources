#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=00-05:00:00
#SBATCH --mem=4GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marylauren.benton@vanderbilt.edu
#SBATCH --output=results/slurm/slurm_%j.out

###
# set environment
###
module load GCC Anaconda3  #todo lmod appropriate libraries
source activate proj_env   #todo change to correct virtual env


###
# variables
###
bin_dir='./bin'
log_dir='./logs'
res_dir='./results'  #todo add analysis subdir if necessary

# make directory for day's results; note the $() expansion
today=$(date +%F)

# string concatenation with += works here too; separate commands with ';'
res_dir+="/$today"; mkdir -p $res_dir 


###
# run scripts 
###
echo "pwd"  # write commands twice: 1 to echo, 2 to execute
pwd

echo "time python $bin_dir/sample_script.py './data/sample_file.bed' >> $res_dir/result.out 2>> $log_dir/result.err"
time python $bin_dir/sample_script.py './data/sample_file.bed' >> $res_dir/result.out 2>> $log_dir/result.err


###
# print job details
###
echo -e "\n$SLURM_JOBID | $today"

