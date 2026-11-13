#!/bin/bash

#SBATCH --time=3:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 4GB   # memory per CPU core
#SBATCH --job-name=parse-hmm   # job name
#SBATCH -o logs/%x-%j.out # STDOUT
##SBATCH --dependency=afterok:6232361

eval "$(conda shell.bash hook)"
conda activate /resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm/

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "======================================================"
echo ""

DIR="/resnick/groups/enviromics/zahra/LDE/results/hmmsearch_results"

#python parse_hmm.py \
#    --hits "$DIR"/archaea/ \
#    --outdir "$DIR"/archaea \
#    --min_genes 1 \
#    --gene_range 15 \
#    --refine

#python parse_hmm.py \
#	--hits "$DIR"/bacteria/ \
#	--outdir "$DIR"/bacteria \
#    --min_genes 1 \
#    --gene_range 15 \
#    --refine

python parse_hmm.py \
    --hits "$DIR"/AFDB/ \
    --outdir "$DIR"/AFDB \
    --min_genes 1 \
    --gene_range 15 \
# don't use refine for AFDB since data is not organized into contigs
echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""

