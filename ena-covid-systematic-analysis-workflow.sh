#!/bin/bash

# SLURM job parameters
#SBATCH --time=7-00:00:00   # walltime
#SBATCH --mem=500G        # memory per node
#SBATCH -J "s_a_w_i&n"   # job name
#SBATCH -e "systematic-analysis-workflow.err"   # job error file
#SBATCH --mail-user=khadim@ebi.ac.uk   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=48
#SBATCH --cpus-per-task=2

rm -rf slurm*

# Load Conda and activate the environment

source /hps/nobackup/cochrane/ena/users/analyser/miniconda3/etc/profile.d/conda.sh
export PATH="/hps/nobackup/cochrane/ena/users/analyser/miniconda3/bin:$PATH"

conda activate analysis_env



/hps/nobackup/cochrane/ena/users/analyser/miniconda3/bin/python fetch_new_data.py -p 'illumina' &
/hps/nobackup/cochrane/ena/users/analyser/miniconda3/bin/python fetch_new_data.py -p 'oxford_nanopore' &

wait

# Run illumina.py and nanopore.py in parallel
/hps/nobackup/cochrane/ena/users/analyser/miniconda3/bin/python illumina.py &
/hps/nobackup/cochrane/ena/users/analyser/miniconda3/bin/python nanopore.py &

wait



