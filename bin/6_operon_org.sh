#!/bin/bash
# Run, for example:
#   bash bin/6_operon_org.sh --prefix AFDB --proteins_dir /path/to/proteins.fasta

#SBATCH --time=4:10:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=150GB
#SBATCH --job-name=operon-org
#SBATCH -o logs/%x-%j.out # STDOUT
##SBATCH --dependency=afterok:66482956

echo "====================================================="
echo "Start Time  : $(date)"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "======================================================"
echo ""

# run with
# bin/6_operon_org.sh \
#   --prefix AFDB \
#   --proteins_dir /resnick/groups/enviromics/zahra/LDE/results/tmp/AFDB_hits.fasta

set -euo pipefail

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

usage() {
    cat >&2 <<'EOF'
Usage: 6_operon_org.sh --prefix PREFIX --proteins_dir PATH [submit|run|postprocess]

  --prefix PREFIX          Prefix used for the results files, e.g. AFDB
  --proteins_dir PATH      Protein FASTA file or directory passed to metadata.py
  submit                   Prepare batches and submit the workflow (default)
  run                      Run one Slurm array task
  postprocess              Process completed MicrobeAnnotator results
EOF
}

PREFIX=""
PROTEINS_DIR=""
COMMAND="submit"

while (($# > 0)); do
    case "$1" in
        --prefix)
            [[ $# -ge 2 ]] || { log "--prefix requires a value"; usage; exit 2; }
            PREFIX="$2"
            shift 2
            ;;
        --proteins_dir)
            [[ $# -ge 2 ]] || { log "--proteins_dir requires a value"; usage; exit 2; }
            PROTEINS_DIR="$2"
            shift 2
            ;;
        submit|run|postprocess)
            COMMAND="$1"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            log "Unknown argument: $1"
            usage
            exit 2
            ;;
    esac
done

if [[ -z "${PREFIX}" || -z "${PROTEINS_DIR}" ]]; then
    log "--prefix and --proteins_dir must be specified"
    usage
    exit 2
fi

# Accept either AFDB or AFDB_; result filenames use one separator.
FILE_PREFIX="${PREFIX%_}_"
ROOT="/resnick/groups/enviromics/zahra/LDE"
SCRIPT="${ROOT}/bin/6_operon_org.sh"
DIR="${ROOT}/operon-org/${PREFIX%_}"
INPUT_DIR="${DIR}/input-fastas"
BATCH_DIR="${DIR}/microbeannotator-batches"
OUTPUT_DIR="${DIR}/microbeannotator"
DB="/resnick/groups/enviromics/databases/microbeannotator-db"
MICROBE_ENV="/resnick/groups/enviromics/zahra/miniconda3/envs/microbeannotator"
PARSE_ENV="/resnick/groups/enviromics/zahra/miniconda3/envs/parse_hmm"
BATCH_SIZE=50

FINAL_FILE="${ROOT}/results/final/${FILE_PREFIX}final_20.csv"

prepare_batches() {
    log "Creating operon FASTA inputs"
    mkdir -p "$BATCH_DIR"
    (
        cd "${ROOT}/bin"
        conda run -p "${PARSE_ENV}" python metadata.py \
		--prepare \
		--clusters_file "${FINAL_FILE}" \
        --operon_dir "${DIR}" \
		--final_file "${FINAL_FILE}" \
		--proteins_dir "${PROTEINS_DIR}"
    )
    find "${BATCH_DIR}" -maxdepth 1 -type f -name 'batch_*' -delete
    find "${INPUT_DIR}" -maxdepth 1 -type f -name '*.fasta' | sort \
        > "${BATCH_DIR}/inputs.txt"
    [[ -s "${BATCH_DIR}/inputs.txt" ]] || {
        log "No operon FASTA files found in ${INPUT_DIR}"
        exit 1
    }
    split -l "${BATCH_SIZE}" -d -a 5 \
        "${BATCH_DIR}/inputs.txt" "${BATCH_DIR}/batch_"
}

run_batch() {
    local batch output
    batch=$(printf '%s/batch_%05d' "${BATCH_DIR}" "${SLURM_ARRAY_TASK_ID}")
    output=$(printf '%s/batch_%05d' "${OUTPUT_DIR}" "${SLURM_ARRAY_TASK_ID}")
    [[ -s "${batch}" ]] || {
        log "Missing or empty batch file: ${batch}"
        exit 1
    }

    mapfile -t inputs < "${batch}"
    eval "$(conda shell.bash hook)"
    conda activate "${MICROBE_ENV}"
    module load diamond/2.1.7-gcc-13.2.0-cfkl5pd

    log "Annotating ${#inputs[@]} operon FASTA files"
    microbeannotator \
        --input "${inputs[@]}" \
        --outdir "${output}" \
        --method diamond \
        --database "${DB}" \
        -p 8 \
        -t 4 \
        --refine \
        --no_plot
}

postprocess() {
    local annotation annotation_dir="${OUTPUT_DIR}/annotation_results" annotation_count=0
    mkdir -p "${annotation_dir}"
    find "${annotation_dir}" -maxdepth 1 -type l -delete

    while IFS= read -r annotation; do
        ln -sf "${annotation}" "${annotation_dir}/$(basename "${annotation}")"
        ((annotation_count += 1))
    done < <(
        find "${OUTPUT_DIR}"/batch_* \
            -type f -path '*/annotation_results/*.annot' | sort
    )

    [[ "${annotation_count}" -gt 0 ]] || {
        log "No annotation files found under ${OUTPUT_DIR}/batch_*"
        exit 1
    }

    log "Collected ${annotation_count} annotation files; calling get_plot_data()"
    (
        cd "${ROOT}/bin"
        conda run -p "${PARSE_ENV}" python metadata.py \
            --data \
            --clusters_file "${FINAL_FILE}" \
            --final_file "${FINAL_FILE}" \
            --operon_dir "${DIR}"
    )
}

submit_workflow() {
    # local batch_count array_job post_job
    local post_job

    # prepare_batches
    # batch_count=$(find "${BATCH_DIR}" -maxdepth 1 \
    #     -type f -name 'batch_[0-9][0-9][0-9][0-9][0-9]' | wc -l)

    # array_job=$(sbatch --parsable \
    #     --output="logs/%x-%A_%a.out" \
    #     --array="0-$((batch_count - 1))" \
    #     "${SCRIPT}" run \
    #     --prefix "${PREFIX}" \
    #     --proteins_dir "${PROTEINS_DIR}")
    post_job=$(sbatch --parsable \
        --output="logs/%x-%j.out" \
        "${SCRIPT}" postprocess \
        --prefix "${PREFIX}" \
        --proteins_dir "${PROTEINS_DIR}")

            # --dependency="afterok:${array_job}" \

   # log "Submitted MicrobeAnnotator array ${array_job} (${batch_count} batches)"
   # log "Submitted get_plot_data job ${post_job} after ${array_job}"
}

case "${COMMAND}" in
    submit) submit_workflow ;;
    run) run_batch ;;
    postprocess) postprocess ;;
esac

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
echo ""
