#!/bin/bash

#SBATCH --time=12:04:00   # walltime #8hrs?
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 40GB   # memory per CPU core
#SBATCH --job-name=conserved-res   # job name
#SBATCH -o logs/%x-%j.out # STDOUT
##SBATCH --dependency=afterok:65687517

eval "$(conda shell.bash hook)"
conda activate /resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "======================================================"
echo ""

#python conserved_res.py --align_fasta --test

python conserved_res.py \
 --min_genes 1 \
 --proteins_dir /resnick/groups/enviromics/zahra/diazoDB-HPC/protein_faa_reps_232 \
 --test \
 --out_prefix GTDB
  # --reload_fasta \
  # --align_fasta \

python conserved_res.py \
  --hits_file /resnick/groups/enviromics/zahra/LDE/results/hmmsearch_results/AFDB_hits.csv \
  --reload_fasta \
  --min_genes 1 \
  --proteins_dir /resnick/groups/enviromics/zahra/LDE/AFDB_ungapped.fasta \
  --test \
  --out_prefix AFDB

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""

