library(GenomicRanges)
library(Rsamtools)
library(Biostrings)
library(monaLisa)
library(TFBSTools)
library(JASPAR2020)
library(BiocParallel)

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
MAX_BACKGROUND_RATIO <- 3
MIN_REGIONS <- 100
MOTIF_FDR <- 0.05
GC_BREAKS <- c(0, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 1)
N_WORKERS <- min(as.integer(Sys.getenv("MOTIF_WORKERS", "2")),
                 max(parallel::detectCores() - 1L, 1L))
BP <- if (N_WORKERS > 1L) MulticoreParam(workers = N_WORKERS,
                                         progressbar = TRUE) else SerialParam()
MAX_FOREGROUND <- 1500
MAX_BACKGROUND <- 4500
SEED <- 1

PEAK_DIR <- "data/raw/chip_peaks"
MARK_PATTERNS <- c(h3k4me3 = "kidney-h3k4me3.*narrowPeak",
                   h3k27ac = "kidney-h3k27ac.*narrowPeak")

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

load_summit_peaks <- function(dir, patterns, seq_map, valid_seqlevels) {
  available <- list.files(dir)
  out <- lapply(names(patterns), function(mark) {
    hit <- grep(patterns[[mark]], available, value = TRUE, ignore.case = TRUE)
    if (length(hit) != 1L) stop("expected one file for ", mark)
    tab <- read.delim(file.path(dir, hit), header = FALSE,
                      stringsAsFactors = FALSE)
    if (ncol(tab) < 10) stop("narrowPeak lacks the summit column for ", mark)
    accession <- unname(seq_map[sub("^chr", "", tab$V1)])
    ok <- !is.na(accession) & accession %in% valid_seqlevels & tab$V10 >= 0
    summit <- tab$V2[ok] + tab$V10[ok] + 1L
    GRanges(seqnames = accession[ok],
            ranges = IRanges(start = summit, width = 1L))
  })
  names(out) <- names(patterns)
  out
}

filter_pwms_expressed <- function(pwms, expressed) {
  keep <- vapply(seq_along(pwms), function(i) {
    any(motif_to_symbols(name(pwms[[i]])) %in% expressed)
  }, logical(1))
  pwms[keep]
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
  gr <- unique(gr)
  sort(gr)
}

name_sequences <- function(seqs, gr) {
  names(seqs) <- sprintf("%s:%d-%d", as.character(seqnames(gr)),
                         start(gr), end(gr))
  seqs
}

regions_for_genes <- function(genes, annotation, peaks, promoters_gr,
                              domains_gr, region_type) {
  stopifnot(all(width(peaks$h3k4me3) == 1L))
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
  seqs <- name_sequences(getSeq(fa, gr), gr)
  keep <- !grepl("N", as.character(seqs), fixed = TRUE)
  seqs <- seqs[keep]
  seqs[!duplicated(names(seqs))]
}

run_enrichment <- function(fg_seqs, bg_seqs, pwms) {
  seqs <- c(fg_seqs, bg_seqs)
  bins <- factor(c(rep("foreground", length(fg_seqs)),
                   rep("background", length(bg_seqs))),
                 levels = c("background", "foreground"))
  stopifnot(!any(duplicated(names(seqs))))
  res <- calcBinnedMotifEnrR(seqs = seqs, bins = bins, pwmL = pwms,
                             background = "otherBins",
                             GCbreaks = GC_BREAKS, BPPARAM = BP,
                             verbose = TRUE)
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
pwms_all <- load_jaspar()
reg <- readRDS(REG_OBJ)
expr_objects <- readRDS(EXPR_OBJ)

annotation <- reg$annotation
promoters_gr <- promoters(reg$genes, upstream = 2000, downstream = 500)
promoters_gr <- promoters_gr[match(annotation$gene_symbol,
                                  mcols(reg$genes)$gene_symbol)]
domains_gr <- reg$domains[match(annotation$gene_symbol,
                                mcols(reg$domains)$gene_symbol)]

summit_peaks <- load_summit_peaks(PEAK_DIR, MARK_PATTERNS, reg$seq_map,
                                  unique(as.character(seqnames(reg$genes))))
message(sprintf("summit peaks: %s",
                paste(sprintf("%s=%d", names(summit_peaks),
                              lengths(summit_peaks)), collapse = ", ")))

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

  pwms <- filter_pwms_expressed(pwms_all, rownames(expr))
  message(sprintf("%s: %d of %d motifs map to an expressed gene",
                  set_name, length(pwms), length(pwms_all)))

  for (direction in c("down", "up")) {
    fg_genes <- de$gene_symbol[de$fdr < FDR_CUTOFF &
                                 (if (direction == "up") de$log_fc > 0
                                  else de$log_fc < 0)]
    bg_genes <- setdiff(universe, de$gene_symbol[de$fdr < FDR_CUTOFF])

    for (exclude_markers in c(FALSE, TRUE)) {
      fg_use <- if (exclude_markers) setdiff(fg_genes, CELL_TYPE_MARKERS)
                else fg_genes

      if (exclude_markers && identical(sort(fg_use), sort(fg_genes))) {
        message(sprintf("%s / %s: no cell-type marker in the set, sensitivity run skipped",
                        set_name, direction))
        next
      }

      for (region_type in c("promoter", "enhancer")) {
        fg_gr <- regions_for_genes(fg_use, annotation, summit_peaks,
                                   promoters_gr, domains_gr, region_type)
        bg_gr <- regions_for_genes(bg_genes, annotation, summit_peaks,
                                   promoters_gr, domains_gr, region_type)
        bg_gr <- bg_gr[!overlapsAny(bg_gr, fg_gr)]

        if (length(fg_gr) > MAX_FOREGROUND) {
          fg_gr <- fg_gr[sort(sample(length(fg_gr), MAX_FOREGROUND))]
        }
        target_bg <- min(MAX_BACKGROUND, MAX_BACKGROUND_RATIO * length(fg_gr))
        if (length(bg_gr) > target_bg) {
          bg_gr <- bg_gr[sort(sample(length(bg_gr), target_bg))]
        }

        label <- sprintf("%s_%s_%s%s", set_name, direction, region_type,
                         if (exclude_markers) "_nomarkers" else "")
        out_file <- file.path(OUT_ROOT, sprintf("motifs_%s.csv", label))

        if (file.exists(out_file)) {
          message(sprintf("%s: already done, skipped", label))
          all_results[[label]] <- read.csv(out_file, stringsAsFactors = FALSE)
          next
        }

        if (length(fg_gr) < MIN_REGIONS) {
          message(sprintf("%s: only %d foreground regions, skipped",
                          label, length(fg_gr)))
          next
        }

        fg_seqs <- extract_sequences(fa, fg_gr)
        bg_seqs <- extract_sequences(fa, bg_gr)

        message(sprintf("%s: %d foreground and %d background sequences",
                        label, length(fg_seqs), length(bg_seqs)))

        elapsed <- system.time(
          res <- run_enrichment(fg_seqs, bg_seqs, pwms)
        )
        message(sprintf("%s: enrichment took %.1f min", label,
                        elapsed[["elapsed"]] / 60))

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
        write.csv(res, out_file, row.names = FALSE)
        all_results[[label]] <- res
        rm(fg_seqs, bg_seqs, fg_gr, bg_gr)
        gc(verbose = FALSE)

        significant <- res[!is.na(res$padj) & res$padj < MOTIF_FDR, ]
        message(sprintf("%s: %d motifs at padj < %.2f", label,
                        nrow(significant), MOTIF_FDR))
      }
    }
  }
}

if (!length(all_results)) {
  stop("no enrichment result was produced")
}

combined <- do.call(rbind, all_results)
write.csv(combined, file.path(OUT_ROOT, "motifs_all.csv"), row.names = FALSE)

combined$direction_coherent <- with(combined,
  (direction == "up" & tf_log_fc > 0) | (direction == "down" & tf_log_fc < 0))

converging <- combined[!is.na(combined$padj) & combined$padj < MOTIF_FDR &
                         !is.na(combined$tf_fdr) & combined$tf_fdr < FDR_CUTOFF &
                         !is.na(combined$direction_coherent) &
                         combined$direction_coherent, ]
converging <- converging[order(converging$set, converging$direction,
                               converging$padj), ]
write.csv(converging, file.path(OUT_ROOT, "converging_candidates.csv"),
          row.names = FALSE)

print(head(converging[, c("set", "direction", "region_type", "motif_name",
                          "log2_enrichment", "padj", "tf_matched",
                          "tf_log_fc", "tf_fdr")], 30))

saveRDS(all_results, file.path(OUT_ROOT, "motif_objects.rds"))
