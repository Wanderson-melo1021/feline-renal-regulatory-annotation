#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
PEAK_DIR="${PROJECT_ROOT}/data/raw/chip_peaks"
ANNOT_DIR="${PROJECT_ROOT}/data/raw/annotation"
LOG_DIR="${PROJECT_ROOT}/logs"
LOG_FILE="${LOG_DIR}/00_download_chromatin.log"

GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182952/suppl"
ENSEMBL_BASE="https://ftp.ensembl.org/pub/release-113"
ASSEMBLY="Felis_catus_9.0"
RELEASE="113"

MARKS=("ctcf" "h3k27ac" "h3k4me3")
TISSUE="kidney"

mkdir -p "${PEAK_DIR}" "${ANNOT_DIR}" "${LOG_DIR}"

log() {
    printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$1" | tee -a "${LOG_FILE}"
}

fetch() {
    local url="$1"
    local dest="$2"
    if [[ -s "${dest}" ]]; then
        log "SKIP ${dest##*/} already present"
        return 0
    fi
    log "GET ${url}"
    curl -fL --retry 5 --retry-delay 5 --connect-timeout 30 -o "${dest}.part" "${url}"
    mv "${dest}.part" "${dest}"
    log "OK ${dest##*/} ($(du -h "${dest}" | cut -f1))"
}

verify_gzip() {
    local target="$1"
    if gzip -t "${target}" 2>/dev/null; then
        log "VALID ${target##*/}"
    else
        log "CORRUPT ${target##*/}"
        return 1
    fi
}

log "=== module 00 start ==="
log "project root ${PROJECT_ROOT}"

for mark in "${MARKS[@]}"; do
    fname="GSE182952_${TISSUE}-${mark}_FinalReplicatedPeaks.narrowPeak.gz"
    fetch "${GEO_BASE}/${fname}" "${PEAK_DIR}/${fname}"
    verify_gzip "${PEAK_DIR}/${fname}"
done

GTF="${ASSEMBLY}.${RELEASE}.gtf.gz"
fetch "${ENSEMBL_BASE}/gtf/felis_catus/Felis_catus.${GTF}" "${ANNOT_DIR}/Felis_catus.${GTF}"
verify_gzip "${ANNOT_DIR}/Felis_catus.${GTF}"

log "=== inventory ==="
for mark in "${MARKS[@]}"; do
    fname="${PEAK_DIR}/GSE182952_${TISSUE}-${mark}_FinalReplicatedPeaks.narrowPeak.gz"
    n_peaks=$(gzip -dc "${fname}" | wc -l | tr -d ' ')
    first_seqname=$(gzip -dc "${fname}" | head -n 1 | cut -f1)
    n_seqnames=$(gzip -dc "${fname}" | cut -f1 | sort -u | wc -l | tr -d ' ')
    log "${mark}: ${n_peaks} peaks | ${n_seqnames} seqnames | example '${first_seqname}'"
done

GTF_PATH="${ANNOT_DIR}/Felis_catus.${GTF}"
gtf_seqname=$(gzip -dc "${GTF_PATH}" | grep -v '^#' | head -n 1 | cut -f1)
gtf_genes=$(gzip -dc "${GTF_PATH}" | awk -F'\t' '$3=="gene"' | wc -l | tr -d ' ')
log "gtf: ${gtf_genes} gene features | example seqname '${gtf_seqname}'"

log "=== module 00 end ==="
