#!/usr/bin/env bash

# Submit with runner/runner.py, not by hand.
# Required env vars:
#   LANTREE_JOB_ID   - API job UUID
#   LANTREE_INPUT    - input FASTA path in the job workspace
#   LANTREE_INTERMEDIATE   - output directory for intermediate files
#   LANTREE_OUTPUT   - final CSV/HTML/PDF path to post back to the API

#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem=4GB
#SBATCH -J lantree_classify
#SBATCH -o /resnick/scratch/zshivji/LanTree/logs/%x-%j.out

set -euo pipefail

REPO_ROOT="/resnick/groups/enviromics/zahra/LDE"
SCRIPT_DIR="$REPO_ROOT/bin"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

# Prefer env vars injected by runner; fall back to positional args if provided.
LANTREE_JOB_ID="${LANTREE_JOB_ID:-${1:-}}"
INPUT_FASTA="${LANTREE_INPUT:-${2:-}}"
INTERMEDIATE="${LANTREE_INTERMEDIATE:-${3:-}}"
FINAL_DIR="${LANTREE_OUTPUT:-${4:-}}"

if [[ -z "$LANTREE_JOB_ID" || -z "$INPUT_FASTA" || -z "$INTERMEDIATE" || -z "$FINAL_DIR" ]]; then
  echo "Missing required inputs. Provide LANTREE_JOB_ID, LANTREE_INPUT, LANTREE_INTERMEDIATE, LANTREE_OUTPUT (or positional args)." >&2
  exit 2
fi

mkdir -p "$INTERMEDIATE" "$FINAL_DIR"

HMMSEARCH_BIN="${LANTREE_HMMSEARCH_BIN:-hmmsearch}"
PRODIGAL_BIN="${LANTREE_PRODIGAL_BIN:-prodigal}"
USE_PRODIGAL="${LANTREE_USE_PRODIGAL:-false}"
LANTREE_CONDA_ENV="${LANTREE_CONDA_ENV:-/resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm}"
LANTREE_PRODIGAL_ENV="${LANTREE_PRODIGAL_ENV:-/resnick/groups/enviromics/zahra/miniconda3/envs/prodigal}"
CONDA_BIN="${LANTREE_CONDA_BIN:-/resnick/groups/enviromics/zahra/miniconda3/bin/conda}"
MICROBE_ENV="${LANTREE_MICROBE_ENV:-/resnick/groups/enviromics/zahra/miniconda3/envs/microbeannotator}"
MICROBE_DB="${LANTREE_MICROBE_DB:-/resnick/groups/enviromics/databases/microbeannotator-db}"
MICROBE_THREADS="${LANTREE_MICROBE_THREADS:-1}"

QUERY_FASTA="$INPUT_FASTA"
HMM_DB="$REPO_ROOT/HMMs/PQQ_ADH.hmm"
JOB_HMM_DIR="$INTERMEDIATE/hmmsearch_results"
JOB_PARSED_HITS="$JOB_HMM_DIR/hits.csv"
TMP="$INTERMEDIATE/tmp"

mkdir -p "$JOB_HMM_DIR" "$TMP"

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID      : $LANTREE_JOB_ID"
echo "Input       : $INPUT_FASTA"
echo "Output      : $FINAL_DIR"
echo "Slurm ID    : ${SLURM_JOBID:-local}"
echo "====================================================="

if [[ ! -x "$CONDA_BIN" ]]; then
  echo "Conda executable not found or not executable: $CONDA_BIN" >&2
  exit 127
fi
eval "$("$CONDA_BIN" shell.bash hook)"
  conda activate "$LANTREE_CONDA_ENV"

if [[ "$USE_PRODIGAL" == "true" ]]; then
  QUERY_FASTA="$OUTDIR/predicted_proteins.faa"
  "$PRODIGAL_BIN" \
    -i "$INPUT_FASTA" \
    -a "$QUERY_FASTA" \
    -p meta \
    -q
fi

INPUT_STEM="$(basename "$INPUT_FASTA")"
INPUT_STEM="${INPUT_STEM%.*}"
SAFE_INPUT_STEM="$(printf '%s' "$INPUT_STEM" | tr -c 'A-Za-z0-9._-' '_')"
SAFE_JOB_ID="$(printf '%s' "$LANTREE_JOB_ID" | tr -c 'A-Za-z0-9._-' '_')"
SAFE_INPUT_STEM="${SAFE_INPUT_STEM:-input}"
HMM_PREFIX="${SAFE_INPUT_STEM}__job_${SAFE_JOB_ID}_nif"
JOB_DOMTBLOUT="$JOB_HMM_DIR/${HMM_PREFIX}.domtblout"
JOB_HMM_OUT="$JOB_HMM_DIR/${HMM_PREFIX}.out"

# Match bin/2_hmmsearch.sh: parse_hmm.py consumes HMMER domtblout files.
"$HMMSEARCH_BIN" \
  --domtblout "$JOB_DOMTBLOUT" \
  -o "$JOB_HMM_OUT" \
  "$HMM_DB" \
  "$QUERY_FASTA"

# run jobs from LanTree/analysis/bin to ensure relative paths to config files are correct
cd "$SCRIPT_DIR"

# Preserve a parser result containing only this job's HMM output.
python parse_hmm.py \
  --hits "$JOB_HMM_DIR" \
  --outdir "$JOB_HMM_DIR" \
  --output_file "$JOB_PARSED_HITS" \
  --refine \
  --min_genes 1

# Run conserved-residue classification entirely within this job's directories.
if python conserved_res.py \
  --hits_file "$JOB_PARSED_HITS" \
  --reload_fasta \
  --min_genes 1 \
  --proteins_dir "$(dirname "$QUERY_FASTA")" \
  --tmp "$TMP" \
  --final_dir "$FINAL_DIR" \
  --skip_metadata \
  --out_prefix pqq \
  --external
then
  :
else
  CLASSIFIER_STATUS=$?
  if [[ "$CLASSIFIER_STATUS" -eq 3 ]]; then
    log "Job completed, no results found"
  fi
  exit "$CLASSIFIER_STATUS"
fi

# Create operon organization diagrams in this job's isolated workspace.
# The classifier has already produced pqq_clusters.csv and pqq_final.csv in
# FINAL_DIR.  The metadata helper uses those files to select neighborhoods,
# while MicrobeAnnotator supplies annotations for surrounding genes.
OPERON_DIR="$INTERMEDIATE/operon-org"
OPERON_INPUT_DIR="$OPERON_DIR/input-fastas"
OPERON_ANNOT_DIR="$OPERON_DIR/microbeannotator"
OPERON_METADATA="$FINAL_DIR/operon_metadata.json"
OPERON_CLUSTERS="$FINAL_DIR/pqq_operons.csv"
FINAL_CSV="$FINAL_DIR/pqq_final.csv"

log "Preparing operon FASTA inputs in $OPERON_INPUT_DIR"
mkdir -p "$OPERON_INPUT_DIR" "$OPERON_ANNOT_DIR"
if [[ -s "$OPERON_CLUSTERS" && -s "$FINAL_CSV" ]]; then
  (
    cd "$SCRIPT_DIR"
    "$CONDA_BIN" run -p "$LANTREE_CONDA_ENV" python metadata.py \
      --prepare \
      --proteins_dir "$(dirname "$QUERY_FASTA")" \
      --clusters_file "$OPERON_CLUSTERS" \
      --operon_dir "$OPERON_DIR"
  )
  
  mapfile -t OPERON_INPUTS < <(find "$OPERON_INPUT_DIR" -maxdepth 1 -type f -name '*.fasta' | sort)
  if ((${#OPERON_INPUTS[@]})); then
    if command -v module >/dev/null 2>&1; then
      module load diamond/2.1.7-gcc-13.2.0-cfkl5pd
    fi
    conda activate "$MICROBE_ENV"
    if [[ ! -d "$MICROBE_DB" ]]; then
      echo "MicrobeAnnotator database does not exist: $MICROBE_DB" >&2
      exit 2
    fi

    log "Annotating ${#OPERON_INPUTS[@]} operon neighborhoods"
    microbeannotator \
      --input "${OPERON_INPUTS[@]}" \
      --outdir "$OPERON_ANNOT_DIR" \
      --method diamond \
      --database "$MICROBE_DB" \
      -t "$MICROBE_THREADS" \
      --refine \
      --light \
      --no_plot

    (
      cd "$SCRIPT_DIR"
      "$CONDA_BIN" run -p "$LANTREE_CONDA_ENV" python metadata.py \
        --data \
        --clusters_file "$OPERON_CLUSTERS" \
        --final_file "$FINAL_CSV" \
        --operon_dir "$OPERON_DIR" \
        --metadata_file "$OPERON_METADATA"
    )
  else
    log "No candidate operons found; skipping operon diagram"
  fi
else
  log "No nif result tables found; skipping operon diagram"
fi

# show identified nif genes on nifH reference tree
#module load blast/2.15.0
#makeblastdb
#blastp

# Emit final result for the runner
if [[ ! -s "$FINAL_CSV" ]]; then
  echo "Expected result was not created: $FINAL_CSV" >&2
  exit 1
fi


echo "====================================================="
echo "End Time    : $(date)"
echo "====================================================="
