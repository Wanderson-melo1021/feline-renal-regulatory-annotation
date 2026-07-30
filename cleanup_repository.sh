#!/usr/bin/env bash
set -euo pipefail

ARCHIVE="_archive"
mkdir -p "${ARCHIVE}/scripts" "${ARCHIVE}/results"

archive() {
    local target="$1"
    local dest="$2"
    if [[ -e "${target}" ]]; then
        mv "${target}" "${dest}/"
        echo "archived ${target}"
    fi
}

echo "=== archiving superseded scripts ==="
for f in 02_error.R 02_error_2.R 02_error_diagnostics.R \
         04_motif_enrichment_correct.R 04_motif_enrichment_third_try.R \
         03_coexpression_modules.R 03b_module_characterization.R \
         03c_composition_adjustment.R \
         07_tissue_specificity_2.R 07_tissue_specificity_3.R \
         08_figures_tables_2.R 08_figures_tables_3.R \
         08_figures_tables_4.R 08_figures_tables_5.R; do
    archive "$f" "${ARCHIVE}/scripts"
done

echo
echo "=== promoting final result directories ==="
if [[ -d results/07_tissue_specificity_3 ]]; then
    archive results/07_tissue_specificity "${ARCHIVE}/results"
    archive results/07_tissue_specificity_2 "${ARCHIVE}/results"
    mv results/07_tissue_specificity_3 results/07_tissue_specificity
    echo "promoted 07_tissue_specificity_3 -> 07_tissue_specificity"
fi

if [[ -d results/08_figures_5 ]]; then
    for d in results/08_figures results/08_figures_2 results/08_figures_3 \
             results/08_figures_4; do
        archive "$d" "${ARCHIVE}/results"
    done
    mv results/08_figures_5 results/08_figures
    echo "promoted 08_figures_5 -> 08_figures"
fi

echo
echo "=== archiving superseded analyses ==="
archive results/02b_diagnostics "${ARCHIVE}/results"
archive results/03_coexpression_modules "${ARCHIVE}/results"
archive results/03b_module_characterization "${ARCHIVE}/results"
archive results/03c_composition_adjustment "${ARCHIVE}/results"
archive results/04_motif_enrichment "${ARCHIVE}/results"
archive results/04b_gc_diagnostics "${ARCHIVE}/results"

echo
echo "=== normalising peak file names ==="
cd data/raw/chip_peaks 2>/dev/null && {
    for f in *_narrowPeak.gz; do
        [[ -e "$f" ]] || continue
        mv "$f" "${f/_narrowPeak.gz/.narrowPeak.gz}"
        echo "renamed $f"
    done
    cd - > /dev/null
}

echo
echo "=== files larger than 50 MB outside ignored paths ==="
find . -type f -size +50M \
    -not -path "./data/*" -not -path "./_archive/*" -not -path "./.git/*" \
    -exec ls -lh {} \; | awk '{print $5, $9}'

echo
echo "=== remaining scripts ==="
ls -1 *.sh *.R 2>/dev/null
