#!/bin/bash
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --tasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --array=1-20
#SBATCH --output=out/array_job_%A_%a.out

# load modules
module load Anaconda3
source activate enh_gain-loss

rloopinput=$1

DAT='/path/to/data/here'
OUTDIR="path/to/results/here/$(date '+%Y-%m-%d')"  # date stamped results directory

mkdir -p $OUTDIR  # generate directory if it does not exist


# pull filenames to use from directory
# .. data path from variable at top of script
# .. uses `ls` to list one filename per line 
# .. alt approach: use a text file that contains each filename
# .. $SLURM_ARRAY_TASK_ID pulls each individual line
shuffle=`ls $DAT/* | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

# run enrichments
SCRIPT="/dors/capra_lab/users/bentonml/resources/bin/calculate_enrichment.py"
python $SCRIPT -i 1000 -s dm3 $shuffle $DAT/not_shuffled_file.bed > $OUTDIR/outfile_${SLURM_ARRAY_TASK_ID}.out

