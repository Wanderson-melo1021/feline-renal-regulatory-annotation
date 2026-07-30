# Cis-regulatory annotation of the feline renal cortex

Analysis code for the manuscript *A cis-regulatory annotation of the feline renal cortex and the tissue-specificity shift in chronic kidney disease*.

The study annotates promoter and enhancer architecture for genes expressed in the renal cortex of the domestic cat, combining replicated ChIP-seq peak sets with the RefSeq annotation, and asks whether the transcriptomic signature of chronic kidney disease has recoverable coregulatory structure. It does not: the only reproducible axis tracks tubular content. Genes that lose expression are instead enriched for enhancers active in the kidney and nowhere else, and genes that gain expression for enhancers shared across organs.

This repository contains code and derived result tables. Primary data are publicly deposited and are downloaded by the scripts; they are not versioned here.

## Data sources

| Data | Source | Notes |
|---|---|---|
| Renal expression matrices | Supplementary Data 3 and 4 of Li et al., *Communications Biology* 2025, doi:10.1038/s42003-025-09164-8 | GEO accession GSE303653; log2 CPM, TMM-normalised |
| Chromatin peak sets | GEO GSE182952, BioProject PRJNA758414 | seven tissues; deposited without an associated publication |
| Genome annotation and sequence | NCBI RefSeq GCF_018350175.1, F.catus_Fca126_mat1.0 | Annotation Release 105 |

Download the expression supplement manually into `data/raw/expression/`. Everything else is fetched by the scripts in step 0.

## Assembly caveat

The two primary datasets are aligned to different assemblies that share chromosome names. Harmonising them by chromosome name retains over 99% of intervals and produces no error, while yielding promoter overlap at the level of chance. Module 02 annotates on the assembly of the chromatin data and aborts if promoter H3K4me3 enrichment falls below fivefold over a permuted null. Anyone combining feline datasets deposited before and after 2021 should apply an equivalent check.

## Pipeline

Run in order. Each module writes to its own directory under `results/` and reads only from earlier outputs.

```
00_download_chromatin.sh          kidney peak sets and genome annotation
00b_download_all_tissues.sh       peak sets for the remaining six tissues
01_differential_expression.R      limma-trend on the deposited log2 CPM matrices
02_regulatory_annotation.R        coordinate harmonisation, promoter and domain assignment
03_module_pipeline.R              coexpression modules with stability assessment
04_motif_enrichment.R             binned motif enrichment against a matched background
05_err_target_characterization.R  functional characterisation of ERR-motif genes
06_architecture_controls.R        controls separating architecture from regulation
07_tissue_specificity.R           tissue breadth of regulatory elements
08_figures_tables.R               manuscript figures and tables
09_cohort_and_age.R               cohort table and age sensitivity
```

Module 03 is parameterised by compartment and contrast in its first two lines.

## Requirements

R with limma, dynamicTreeCut, rtracklayer, GenomicRanges, Rsamtools, Biostrings, monaLisa, TFBSTools, JASPAR2020, gprofiler2, org.Hs.eg.db, readxl, ggplot2 and patchwork.

Motif enrichment is the only computationally heavy step. Set the `MOTIF_WORKERS` environment variable to control parallelism; on a machine with limited memory, use `MOTIF_WORKERS=1`. That module writes each result as it completes and skips combinations already present, so an interrupted run can be resumed.

## Reproducibility notes

This was an exploratory analysis without a prespecified plan. Decisions made in response to intermediate results are enumerated in the manuscript Methods. Superseded script versions from development are not included; the commit history is the record of the analytical path.

## Licensing

Three layers, licensed separately.

**Code** — all `.R` and `.sh` files, under the MIT License. See `LICENSE`.

**Derived result tables** — the contents of `results/`, under Creative Commons Attribution 4.0 International. See `results/LICENSE`.

**Primary data** — not covered by either licence and not redistributed here. The expression matrices originate from an open-access article under the terms set by its publisher; the chromatin peak sets and the genome annotation originate from NCBI. Anyone reusing those data should attribute them to their original sources and observe the terms attached there.

## Citation

If you use this code or the derived tables, please cite the manuscript and this repository:

```
[manuscript citation, to be completed on acceptance]

[author list]. Cis-regulatory annotation of the feline renal cortex.
Zenodo. doi:[concept DOI]
```

The Zenodo concept DOI resolves to the most recent version and should be preferred over a version-specific DOI.
