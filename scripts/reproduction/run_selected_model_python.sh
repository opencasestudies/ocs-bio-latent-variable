#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
IMAGE="${COGAPS_RUNTIME_IMAGE:-othomas2/pycogaps-runtime-guide:0.3.0}"
OUTDIR="data/processed/selected_model_k6/python"
FULL_RESULT="${ROOT}/${OUTDIR}/cogaps_K6_seed2_iter2000.h5ad"
INPUT_H5AD="${ROOT}/data/processed/input/cogaps_input_genesxcells_hvg3000_float64.h5ad"
FORCE_RERUN="${FORCE_RERUN:-0}"

DOCKER_TTY=()
if [[ -t 1 ]]; then
  DOCKER_TTY=(-t)
fi

cd "${ROOT}"
mkdir -p "${OUTDIR}"

if [[ ! -f "${INPUT_H5AD}" ]]; then
  echo "Missing dense Python CoGAPS input: ${INPUT_H5AD}" >&2
  echo "Create it from the full reproduction guide or obtain it from the external artifact bundle." >&2
  exit 1
fi

if [[ -e "${FULL_RESULT}" && "${FORCE_RERUN}" != "1" ]]; then
  echo "Full Python model already exists: ${FULL_RESULT}" >&2
  echo "Set FORCE_RERUN=1 to overwrite it in a visible Terminal session." >&2
  exit 2
fi

docker run --rm "${DOCKER_TTY[@]}" --platform linux/amd64 \
  -v "${ROOT}:/workspace/case-study" \
  -w /workspace/case-study \
  "${IMAGE}" \
  bash -lc 'set -euo pipefail
    export PYTHONWARNINGS=ignore::FutureWarning
    python scripts/cogaps_run_one_singleprocess.py \
      --use-sparse-opt \
      --outdir data/processed/selected_model_k6/python \
      --cogaps-input-h5ad data/processed/input/cogaps_input_genesxcells_hvg3000_float64.h5ad \
      --preprocessed-h5ad data/processed/input/preprocessed_cells_hvg3000.h5ad \
      --k 6 \
      --seed 2 \
      --n-iter 2000 \
      --top-genes 50 \
      --stim-label stim \
      --blas-threads 1 \
      --cogaps-threads 1 \
      --output-frequency 1000 \
      --checkpoint-interval 0 \
      --n-snapshots 0 \
      --snapshot-phase sampling'
