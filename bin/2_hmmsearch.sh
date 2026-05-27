#!/bin/bash

#SBATCH --time=18:10:00   # walltime --> 15 hrs for both Archaea, Bacteria
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 4GB   # memory per CPU core
#SBATCH -J hmmsearch.e%j   # job name

# Notify at the beginning, end of job and on failure.
#SBATCH --mail-user=zshivji@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

## /SBATCH -p general # partition (queue)
## /SBATCH -o slurm.%N.%j.out # STDOUT
## /SBATCH -e slurm.%N.%j.err # STDERR

eval "$(conda shell.bash hook)"
conda activate /resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "======================================================"
echo ""

#HMM search on GTDB archaea
for file in ../../diazoDB-HPC/protein_faa_reps/archaea/*; do
        f="${file##*/}"
        hmmsearch ../HMMs/methanol_dehydrogenase_pqq_xoxF_mxaF.hmm "$file" >> "../results/hmmsearch_results/archaea/${f%.faa}_MDH.out"
done

chgrp hpc_enviromics ../results/hmmsearch_results/archaea/*

echo "Archaea completed"

#HMM search on GTDB bacteria
for file in ../../diazoDB-HPC/protein_faa_reps/bacteria/*; do
        f="${file##*/}"
        hmmsearch ../HMMs/methanol_dehydrogenase_pqq_xoxF_mxaF.hmm "$file" >> "../results/hmmsearch_results/bacteria/${f%.faa}_nif.out"
done

for file in ../results/hmmsearch_results/bacteria/*; do
        chgrp hpc_enviromics ../results/hmmsearch_results/bacteria/"$file"
done

echo "Bacteria completed"

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""

