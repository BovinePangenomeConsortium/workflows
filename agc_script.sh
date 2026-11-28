#!/usr/bin/env bash
set -euo pipefail

# Preferred workflow:
#   snakemake -s Snakefile --configfile config.yaml data/<graph>.agc
#
# Standalone helper to build an AGC archive from PanSN-formatted assemblies
# listed in a BPC manifest (TSV or CSV).

if [ "$#" -lt 3 ] || [ "$#" -gt 5 ]; then
  echo "Usage: $0 <manifest> <assemblies_dir> <out.agc> [graph=All] [threads=8]" >&2
  exit 1
fi

MANIFEST="$1"
ASSEMBLIES_DIR="$2"
OUT_AGC="$3"
GRAPH="${4:-All}"
THREADS="${5:-8}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLES=$(python3 "${SCRIPT_DIR}/scripts/manifest.py" list-graphs --manifest "${MANIFEST}" | awk -F'\t' -v g="${GRAPH}" '$1==g {print $3}')

if [ -z "${SAMPLES}" ]; then
  echo "No samples resolved for graph '${GRAPH}' in '${MANIFEST}'" >&2
  exit 1
fi

INPUT_FASTAS=()
IFS=',' read -r -a SAMPLE_ARR <<< "${SAMPLES}"
for S in "${SAMPLE_ARR[@]}"; do
  F="${ASSEMBLIES_DIR}/${S}.fa.gz"
  if [ ! -f "${F}" ]; then
    echo "Missing formatted assembly: ${F}" >&2
    exit 1
  fi
  INPUT_FASTAS+=("${F}")
done

mkdir -p "$(dirname "${OUT_AGC}")"
agc create -d -t "${THREADS}" -o "${OUT_AGC}" "${INPUT_FASTAS[@]}"
