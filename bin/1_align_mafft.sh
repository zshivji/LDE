#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --time=00:12:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem 4GB   # memory per CPU core
#SBATCH --job-name=hmm_seed   # job name
#SBATCH -o logs/%x-%j.out # STDOUT

echo "====================================================="
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

eval "$(conda shell.bash hook)"
conda activate /resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm

module load mafft/7.505-gcc-13.2.0-nklkvtc

# combine seed sequences
cat ../HMM_seeds/ExaF.faa ../HMM_seeds/XoxF.faa > ../HMM_seeds/combined_ExaF_XoxF.faa

# align sequences
mafft --auto --thread 4 ../HMM_seeds/combined_ExaF_XoxF.faa > ../alignments/combined_ExaF_XoxF.aln

# build hmm profiles
hmmbuild ../HMMs/combined_ExaF_XoxF.hmm ../alignments/combined_ExaF_XoxF.aln

# combine hmm profiles
#cat ../HMMs/*.hmm > ../HMMs/PQQ_ADH.hmm

# # convert Banfield HMM to HMMER3 format (only need to do this once)
# hmmconvert ../HMMs/methanol_dehydrogenase_pqq_xoxF_mxaF.hmm \
#   > /tmp/methanol_dehydrogenase.hmm

# cat ../HMMs/combined_ExaF_XoxF.hmm \
#     /tmp/methanol_dehydrogenase.hmm \
#     > ../HMMs/PQQ_ADH.hmm

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
