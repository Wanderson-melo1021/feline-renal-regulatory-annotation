library(GenomicRanges)
library(Rsamtools)
library(Biostrings)
library(monaLisa)
library(TFBSTools)
library(JASPAR2020)

FASTA_PATH <- "data/raw/annotation/GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna"
REG_OBJ <- "data/processed/regulatory_objects.rds"
EXPR_OBJ <- "data/processed/expression_objects.rds"
DE_DIR <- "results/01_differential_expression"
OUT_ROOT <- "results/04_motif_enrichment"

SETS <- list(
  cortex_late = list(compartment = "cortex", contrast = "late_vs_control"),
  medulla_early = list(compartment = "medulla", contrast = "early_vs_control")
)

FDR_CUTOFF <- 0.10
REGION_WIDTH <- 300
MAX_BACKGROUND_RATIO <- 5
MIN_REGIONS <- 100
MOTIF_FDR <- 0.05
SEED <- 1

CELL_TYPE_MARKERS <- c(
  "PTPRC", "CD68", "CSF1R", "ITGAM", "AIF1", "LYZ", "C1QA", "C1QB",
  "TYROBP", "FCER1G", "COL1A1", "COL1A2", "COL3A1", "ACTA2", "PDGFRB",
  "FN1", "LUM", "DCN"
)

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

ensure_fasta_index <- function(path) {
  if (!file.exists(path)) stop("fasta not found: ", path)
  if (!file.exists(paste0(path, ".fai"))) {
    message("indexing fasta, this runs once")
    indexFa(path)
  }
  FaFile(path)
}

load_jaspar <- function() {
  pfm <- getMatrixSet(JASPAR2020,
                      list(species = 9606, collection = "CORE",
                           matrixtype = "PFM"))
  toPWM(pfm)
}

standardize_regions <- function(gr, width) {
  gr <- resize(gr, width = 1, fix = "center")
  gr <- resize(gr, width = width, fix = "center")
  unique(gr)
}

regions_for_genes <- function(genes, annotation, peaks, promoters_gr,
                              domains_gr, region_type) {
  idx <- which(annotation$gene_symbol %in% genes)
  if (!length(idx)) return(GRanges())

  if (region_type == "promoter") {
    windows <- promoters_gr[idx]
    hits <- subsetByOverlaps(peaks$h3k4me3, windows)
  } else {
    windows <- domains_gr[idx]
    candidates <- subsetByOverlaps(peaks$h3k27ac, windows)
    hits <- candidates[!overlapsAny(candidates, promoters_gr)]
  }
  standardize_regions(hits, REGION_WIDTH)
}

extract_sequences <- function(fa, gr) {
  seqs <- getSeq(fa, gr)
  keep <- !grepl("N", as.character(seqs), fixed = TRUE)
  seqs[keep]
}

run_enrichment <- function(fg_seqs, bg_seqs, pwms) {
  seqs <- c(fg_seqs, bg_seqs)
  bins <- factor(c(rep("foreground", length(fg_seqs)),
                   rep("background", length(bg_seqs))),
                 levels = c("background", "foreground"))
  res <- calcBinnedMotifEnrR(seqs = seqs, bins = bins, pwmL = pwms,
                             background = "otherBins", verbose = FALSE)
  data.frame(
    motif_id = rownames(res),
    motif_name = SummarizedExperiment::rowData(res)$motif.name,
    log2_enrichment = SummarizedExperiment::assay(res, "log2enr")[, "foreground"],
    neg_log10_padj = SummarizedExperiment::assay(res, "negLog10Padj")[, "foreground"],
    row.names = NULL, stringsAsFactors = FALSE
  )
}

motif_to_symbols <- function(motif_name) {
  parts <- unlist(strsplit(motif_name, "::|\\+"))
  toupper(trimws(parts))
}

fa <- ensure_fasta_index(FASTA_PATH)
pwms <- load_jaspar()
reg <- readRDS(REG_OBJ)
expr_objects <- readRDS(EXPR_OBJ)

annotation <- reg$annotation
promoters_gr <- promoters(reg$genes, upstream = 2000, downstream = 500)
promoters_gr <- promoters_gr[match(annotation$gene_symbol,
                                  mcols(reg$genes)$gene_symbol)]
domains_gr <- reg$domains[match(annotation$gene_symbol,
                                mcols(reg$domains)$gene_symbol)]

set.seed(SEED)
all_results <- list()

for (set_name in names(SETS)) {
  cfg <- SETS[[set_name]]
  expr <- expr_objects[[cfg$compartment]]$expr
  de <- read.csv(file.path(DE_DIR, sprintf("de_%s_%s.csv", cfg$compartment,
                                           cfg$contrast)),
                 stringsAsFactors = FALSE)

  coding <- annotation$gene_symbol[!is.na(annotation$gene_biotype) &
                                     annotation$gene_biotype == "protein_coding"]
  universe <- intersect(intersect(rownames(expr), coding), de$gene_symbol)
  de <- de[de$gene_symbol %in% universe, ]

  for (direction in c("down", "up")) {
    fg_genes <- de$gene_symbol[de$fdr < FDR_CUTOFF &
                                 (if (direction == "up") de$log_fc > 0
                                  else de$log_fc < 0)]
    bg_genes <- setdiff(universe, de$gene_symbol[de$fdr < FDR_CUTOFF])

    for (exclude_markers in c(FALSE, TRUE)) {
      fg_use <- if (exclude_markers) setdiff(fg_genes, CELL_TYPE_MARKERS)
                else fg_genes

      for (region_type in c("promoter", "enhancer")) {
        fg_gr <- regions_for_genes(fg_use, annotation, reg$peaks,
                                   promoters_gr, domains_gr, region_type)
        bg_gr <- regions_for_genes(bg_genes, annotation, reg$peaks,
                                   promoters_gr, domains_gr, region_type)
        bg_gr <- bg_gr[!overlapsAny(bg_gr, fg_gr)]

        if (length(bg_gr) > MAX_BACKGROUND_RATIO * length(fg_gr)) {
          bg_gr <- bg_gr[sample(length(bg_gr),
                                MAX_BACKGROUND_RATIO * length(fg_gr))]
        }

        label <- sprintf("%s_%s_%s%s", set_name, direction, region_type,
                         if (exclude_markers) "_nomarkers" else "")

        if (length(fg_gr) < MIN_REGIONS) {
          message(sprintf("%s: only %d foreground regions, skipped",
                          label, length(fg_gr)))
          next
        }

        fg_seqs <- extract_sequences(fa, fg_gr)
        bg_seqs <- extract_sequences(fa, bg_gr)

        message(sprintf("%s: %d foreground and %d background sequences",
                        label, length(fg_seqs), length(bg_seqs)))

        res <- try(run_enrichment(fg_seqs, bg_seqs, pwms), silent = TRUE)
        if (inherits(res, "try-error")) {
          message(sprintf("%s: enrichment failed", label))
          next
        }

        res$set <- set_name
        res$direction <- direction
        res$region_type <- region_type
        res$markers_excluded <- exclude_markers
        res$padj <- 10^(-res$neg_log10_padj)

        tf_expression <- do.call(rbind, lapply(seq_len(nrow(res)), function(i) {
          symbols <- motif_to_symbols(res$motif_name[i])
          matched <- de[de$gene_symbol %in% symbols, ]
          if (!nrow(matched)) {
            data.frame(tf_log_fc = NA_real_, tf_fdr = NA_real_,
                       tf_matched = NA_character_, stringsAsFactors = FALSE)
          } else {
            best <- matched[which.min(matched$fdr), ]
            data.frame(tf_log_fc = best$log_fc, tf_fdr = best$fdr,
                       tf_matched = best$gene_symbol,
                       stringsAsFactors = FALSE)
          }
        }))
        res <- cbind(res, tf_expression)

        res <- res[order(res$padj), ]
        write.csv(res, file.path(OUT_ROOT, sprintf("motifs_%s.csv", label)),
                  row.names = FALSE)
        all_results[[label]] <- res

        significant <- res[!is.na(res$padj) & res$padj < MOTIF_FDR, ]
        message(sprintf("%s: %d motifs at padj < %.2f", label,
                        nrow(significant), MOTIF_FDR))
      }
    }
  }
}

combined <- do.call(rbind, all_results)
write.csv(combined, file.path(OUT_ROOT, "motifs_all.csv"), row.names = FALSE)

converging <- combined[!is.na(combined$padj) & combined$padj < MOTIF_FDR &
                         !is.na(combined$tf_fdr) & combined$tf_fdr < FDR_CUTOFF, ]
converging <- converging[order(converging$set, converging$direction,
                               converging$padj), ]
write.csv(converging, file.path(OUT_ROOT, "converging_candidates.csv"),
          row.names = FALSE)

print(head(converging[, c("set", "direction", "region_type", "motif_name",
                          "log2_enrichment", "padj", "tf_matched",
                          "tf_log_fc", "tf_fdr")], 30))

saveRDS(all_results, file.path(OUT_ROOT, "motif_objects.rds"))
