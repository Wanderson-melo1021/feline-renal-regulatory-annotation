#!/usr/bin/env bash

declare -a CHECKS=(
"01_differential_expression.R|ADJUST_FOR_AGE|age in the model"
"02_regulatory_annotation.R|PEAK_WIDTH_QUANTILE|peak width filter"
"03_module_pipeline.R|commandArgs|command-line parameters"
"03_module_pipeline.R|direction_pc1|per-direction principal component"
"04_motif_enrichment.R|load_summit_peaks|peak summits"
"04_motif_enrichment.R|name_sequences|unique sequence names"
"04_motif_enrichment.R|GC_BREAKS|extended GC bins"
"04_motif_enrichment.R|MOTIF_WORKERS|parallelism and resumability"
"04_motif_enrichment.R|direction_coherent|directional coherence"
"05_err_target_characterization.R|CONTROL_MOTIFS|CTCF and SP1 controls"
"06_architecture_controls.R|control2_interaction_model|interaction model"
"07_tissue_specificity.R|PEAK_WIDTH_QUANTILE|peak width filter"
"07_tissue_specificity.R|mean_fold_change|signal-adjusted models"
"08_figures_tables.R|AGE_DIR|cohort table and age figure"
"09_cohort_and_age.R|age_sensitivity_genes|per-gene age sensitivity"
)

fail=0
printf "%-36s %-30s %s\n" "SCRIPT" "FEATURE" "STATUS"
printf '%.0s-' {1..80}; echo

for entry in "${CHECKS[@]}"; do
    IFS='|' read -r file pattern label <<< "${entry}"
    if [[ ! -f "${file}" ]]; then
        printf "%-36s %-30s %s\n" "${file}" "${label}" "MISSING FILE"
        fail=1
    elif grep -q "${pattern}" "${file}"; then
        printf "%-36s %-30s %s\n" "${file}" "${label}" "ok"
    else
        printf "%-36s %-30s %s\n" "${file}" "${label}" "OUTDATED"
        fail=1
    fi
done

echo
if [[ ${fail} -eq 0 ]]; then
    echo "All modules are at the expected version."
else
    echo "One or more modules are outdated. Replace them before running the pipeline."
fi
exit ${fail}
