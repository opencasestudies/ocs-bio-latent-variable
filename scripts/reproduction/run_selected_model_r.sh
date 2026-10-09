#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
IMAGE="${COGAPS_RUNTIME_IMAGE:-othomas2/pycogaps-runtime-guide:0.3.1}"
OUTDIR="data/processed/selected_model_k6/r"
TAG="K6_seed2_iter2000"
FULL_RESULT="${ROOT}/${OUTDIR}/cogaps_${TAG}.rds"
FORCE_RERUN="${FORCE_RERUN:-0}"

DOCKER_RUN_ARGS=(--rm)
if [[ -t 1 ]]; then
  DOCKER_RUN_ARGS+=(-t)
fi

cd "${ROOT}"

if [[ -e "${FULL_RESULT}" && "${FORCE_RERUN}" != "1" ]]; then
  echo "Full R model already exists: ${FULL_RESULT}" >&2
  echo "Set FORCE_RERUN=1 to overwrite it in a visible Terminal session." >&2
  exit 2
fi

if ! command -v docker >/dev/null 2>&1; then
  echo "Docker is not available in this terminal." >&2
  echo "Run this command in your computer's Terminal, not in the RStudio Terminal. The wrapper starts its own Docker container." >&2
  echo "If you are already in the Host Terminal, check that Docker is installed and available there." >&2
  exit 1
fi

if ! docker info >/dev/null; then
  echo "The Docker engine could not be reached; no model was started." >&2
  echo "In the Host Terminal, open Docker Desktop, wait for the engine to start, and retry docker info." >&2
  echo "If it still fails, check Docker's connection settings and the error above before rerunning this wrapper." >&2
  exit 1
fi

mkdir -p "${OUTDIR}"

docker run "${DOCKER_RUN_ARGS[@]}" --platform linux/amd64 \
  -v "${ROOT}:/workspace/case-study" \
  -w /workspace/case-study \
  "${IMAGE}" \
  bash -lc 'set -euo pipefail
    export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    Rscript -e "cat(\"R CoGAPS OpenMP support:\", CoGAPS::compiledWithOpenMPSupport(), \"\n\")"
    Rscript scripts/cogaps_run_one_singleprocess_r.R \
      --preprocessed-h5ad data/processed/input/preprocessed_cells_hvg3000.h5ad \
      --outdir data/processed/selected_model_k6/r \
      --k 6 \
      --seed 2 \
      --n-iter 2000 \
      --top-genes 50 \
      --stim-label stim \
      --cogaps-threads 4 \
      --async-updates \
      --use-sparse-opt \
      --output-frequency 100 \
      --checkpoint-interval 500 \
      --n-snapshots 10 \
      --snapshot-phase all \
      --take-pump-samples \
      --force-rerun'
