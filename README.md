# Cis-regulatory annotation of the feline renal cortex

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21710845.svg)](https://doi.org/10.5281/zenodo.21710845)

Analysis code for the manuscript *A cis-regulatory annotation of the feline renal cortex and the tissue-specificity shift in chronic kidney disease*.

The study annotates promoter and enhancer architecture for genes expressed in the renal cortex of the domestic cat, combining replicated ChIP-seq peak sets with the RefSeq annotation, and asks whether the transcriptomic signature of chronic kidney disease has recoverable coregulatory structure. It does not: none of 69 candidate coexpression modules survives resampling of animals, and the dominant axis of each direction tracks cell-type composition instead. Genes that lose expression are enriched for enhancers active in the kidney and nowhere else, and genes that gain expression for enhancers shared across organs.

This repository contains code and derived result tables. Primary data are publicly deposited and are not versioned here.

## Data sources

| Data | Source | Notes |
|---|---|---|
| Renal expression matrices | Supplementary Data 3 and 4 of Li et al., *Communications Biology* 2025, doi:10.1038/s42003-025-09164-8 | GEO accession GSE303653; log2 CPM, TMM-normalised |
| Chromatin peak sets | GEO GSE182952, BioProject PRJNA758414 | seven tissues; deposited without an associated publication |
| Genome annotation and sequence | NCBI RefSeq GCF_018350175.1, F.catus_Fca126_mat1.0 | Annotation Release 105 |

## Setting up

The `data/` directory is not versioned. Create it before running anything:

```
mkdir -p data/raw/expression data/raw/annotation data/raw/chip_peaks data/processed
```

One file must be downloaded by hand, because it is distributed as journal supplementary material rather than through a public API: the Supplementary Data workbook of Li et al., file `42003_2025_9164_MOESM3_ESM.xlsx`, available from https://www.nature.com/articles/s42003-025-09164-8. Place it in `data/raw/expression/` under that exact filename, which is what module 01 expects. Sheets `S3 RNA cortex` and `S4 RNA medulla` carry the expression matrices and the sample metadata.

Peak sets, genome annotation and genomic sequence are fetched by the scripts in step 0 and do not need to be retrieved manually. The genomic FASTA must be decompressed before module 04, which indexes it on first use.

## Assembly caveat

The two primary datasets are aligned to different assemblies that share chromosome names. Harmonising them by chromosome name retains over 99% of intervals and produces no error, while yielding promoter overlap at the level of chance. Module 02 annotates on the assembly of the chromatin data and aborts if promoter H3K4me3 enrichment falls below fivefold over a permuted null. Anyone combining feline datasets deposited before and after 2021 should apply an equivalent check.

## Pipeline

Run in order. Each module writes to its own directory under `results/` and reads only from earlier outputs.

```
00_download_chromatin.sh          kidney peak sets and genome annotation
00b_download_all_tissues.sh       peak sets for the remaining six tissues
01_differential_expression.R      limma-trend on the deposited log2 CPM matrices
02_regulatory_annotation.R        coordinate harmonisation, promoter and domain assignment
03_module_pipeline.R              coexpression modules, stability, principal components
04_motif_enrichment.R             binned motif enrichment against a matched background
05_err_target_characterization.R  functional characterisation of ERR-motif genes
06_architecture_controls.R        controls separating architecture from regulation
07_tissue_specificity.R           tissue breadth of regulatory elements
08_figures_tables.R               manuscript figures and tables
09_cohort_and_age.R               cohort table and age sensitivity
```

Module 03 takes its parameters on the command line and is run once per compartment:

```
Rscript 03_module_pipeline.R cortex late_vs_control 0.10 0
Rscript 03_module_pipeline.R medulla early_vs_control 0.10 0
```

The four arguments are compartment, contrast, false discovery threshold and absolute log2 fold-change threshold. The same thresholds are used in both compartments.

Module 01 includes age in the differential expression model by default; set `ADJUST_FOR_AGE` to `FALSE` to reproduce the specification used by the source study, and `MIN_AGE` to a positive value to exclude younger animals as a sensitivity analysis.

## Environment

R 4.6.0 with:

| Package | Version | Package | Version |
|---|---|---|---|
| GenomicRanges | 1.64.0 | monaLisa | 1.18.0 |
| rtracklayer | 1.72.0 | TFBSTools | 1.50.0 |
| Rsamtools | 2.28.0 | JASPAR2020 | 0.99.10 |
| Biostrings | 2.80.0 | gprofiler2 | 0.2.4 |
| limma | 3.68.3 | org.Hs.eg.db | 3.23.1 |
| dynamicTreeCut | 1.63-1 | readxl | 1.5.0 |
| ggplot2 | 4.0.3 | patchwork | 1.3.2 |

Motif enrichment is the only computationally heavy step. Set the `MOTIF_WORKERS` environment variable to control parallelism; on a machine with limited memory, use `MOTIF_WORKERS=1`. That module writes each result as it completes and skips combinations already present, so an interrupted run can be resumed.

## Reproducibility notes

This was an exploratory analysis without a prespecified plan. Decisions made in response to intermediate results are enumerated in the manuscript Methods. Superseded script versions from development are not included; the commit history is the record of the analytical path.

Result tables under `results/` correspond to the code as committed. Anyone re-running the pipeline should obtain the same values.

## Licensing

Three layers, licensed separately.

**Code** — all `.R` and `.sh` files, under the MIT License. See `LICENSE`.

**Derived result tables** — the contents of `results/`, under Creative Commons Attribution 4.0 International. See `results/LICENSE`.

**Primary data** — not covered by either licence and not redistributed here. The expression matrices originate from an open-access article under the terms set by its publisher; the chromatin peak sets and the genome annotation originate from NCBI. Anyone reusing those data should attribute them to their original sources and observe the terms attached there.

## Artificial intelligence assistance

See `ai_use_statement.md`.

## Citation

```
Melo WGG. Cis-regulatory annotation of the feline renal cortex and the
tissue-specificity shift in chronic kidney disease: analysis code.
Zenodo. doi:10.5281/zenodo.21710845
```

[manuscript citation, to be added on acceptance]

The DOI above is the concept identifier and always resolves to the most recent version.
