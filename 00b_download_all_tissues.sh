#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
PEAK_DIR="${PROJECT_ROOT}/data/raw/chip_peaks"
LOG_FILE="${PROJECT_ROOT}/logs/00b_download_all_tissues.log"

GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182952/suppl"

TISSUES=("heart" "kidney" "liver" "lung" "smallintestine" "testis" "thyroid")
MARKS=("ctcf" "h3k27ac" "h3k4me3")

mkdir -p "${PEAK_DIR}" "$(dirname "${LOG_FILE}")"

log() {
    printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$1" | tee -a "${LOG_FILE}"
}

log "=== downloading peak files for ${#TISSUES[@]} tissues ==="

for tissue in "${TISSUES[@]}"; do
    for mark in "${MARKS[@]}"; do
        fname="GSE182952_${tissue}-${mark}_FinalReplicatedPeaks.narrowPeak.gz"
        dest="${PEAK_DIR}/${fname}"

        existing=$(find "${PEAK_DIR}" -maxdepth 1 -name "*${tissue}-${mark}*narrowPeak*" | head -n 1)
        if [[ -n "${existing}" ]]; then
            log "SKIP ${tissue}/${mark} present as ${existing##*/}"
            continue
        fi

        log "GET ${fname}"
        curl -fL --retry 5 --retry-delay 5 -o "${dest}.part" "${GEO_BASE}/${fname}"
        mv "${dest}.part" "${dest}"
        gzip -t "${dest}"
        log "OK ${fname} ($(du -h "${dest}" | cut -f1))"
    done
done

log "=== inventory ==="
for tissue in "${TISSUES[@]}"; do
    for mark in "${MARKS[@]}"; do
        f=$(find "${PEAK_DIR}" -maxdepth 1 -name "*${tissue}-${mark}*narrowPeak*" | head -n 1)
        if [[ -n "${f}" ]]; then
            n=$(gzip -dc "${f}" | wc -l | tr -d ' ')
            log "${tissue}/${mark}: ${n} peaks"
        else
            log "${tissue}/${mark}: MISSING"
        fi
    done
done
