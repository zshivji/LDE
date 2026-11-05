#!/bin/bash

#SBATCH --time=30:10:00   # walltime --> 15 hrs for both Archaea, Bacteria
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 4GB   # memory per CPU core
#SBATCH --job-name=hmmsearch   # job name
#SBATCH -o logs/%x-%j.out # STDOUT

echo "====================================================="
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

GTDB_REP_DIR="/resnick/groups/enviromics/zahra/diazoDB-HPC/protein_faa_reps_232"

eval "$(conda shell.bash hook)"
conda activate /resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm

#HMM search on GTDB archaea
rm -R ../results/hmmsearch_results/archaea/
mkdir -p ../results/hmmsearch_results/archaea/

for file in "$GTDB_REP_DIR"/archaea/*; do
        f="${file##*/}"
        hmmsearch \
        --domtblout "../results/hmmsearch_results/archaea/${f%.faa}_ADH.domtblout" \
        -o "../results/hmmsearch_results/archaea/${f%.faa}_ADH.out" \
        ../HMMs/PQQ_ADH.hmm "$file"
done

chgrp hpc_enviromics ../results/hmmsearch_results/archaea/*

echo "Archaea completed"

#HMM search on GTDB bacteria
rm -R ../results/hmmsearch_results/bacteria/
mkdir -p ../results/hmmsearch_results/bacteria/

for file in "$GTDB_REP_DIR"/bacteria/*; do
        f="${file##*/}"
        hmmsearch \
        --domtblout "../results/hmmsearch_results/bacteria/${f%.faa}_ADH.domtblout" \
        -o "../results/hmmsearch_results/bacteria/${f%.faa}_ADH.out" \
        ../HMMs/PQQ_ADH.hmm "$file"
done

for file in ../results/hmmsearch_results/bacteria/*; do
        chgrp hpc_enviromics "$file"
done

echo "Bacteria completed"

# HMM search on AFDB
#hmmsearch --domtblout ../results/hmmsearch_results/AFDB/AFDB_ADH.domtblout -o ../results/hmmsearch_results/AFDB/AFDB_ADH.out ../HMMs/PQQ_ADH.hmm ../AFDB_ungapped.fasta

echo "AFDB completed"

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""

