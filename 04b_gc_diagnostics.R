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
OUT_DIR <- "results/04b_gc_diagnostics"

COMPARTMENT <- "cortex"
CONTRAST <- "late_vs_control"
FDR_CUTOFF <- 0.10
REGION_WIDTH <- 300
DEFAULT_BREAKS <- c(0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8)
EXTENDED_BREAKS <- c(0, DEFAULT_BREAKS, 1)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

fa <- FaFile(FASTA_PATH)
reg <- readRDS(REG_OBJ)
expr_objects <- readRDS(EXPR_OBJ)

annotation <- reg$annotation
promoters_gr <- promoters(reg$genes, upstream = 2000, downstream = 500)
promoters_gr <- promoters_gr[match(annotation$gene_symbol,
                                  mcols(reg$genes)$gene_symbol)]

expr <- expr_objects[[COMPARTMENT]]$expr
de <- read.csv(file.path(DE_DIR, sprintf("de_%s_%s.csv", COMPARTMENT, CONTRAST)),
               stringsAsFactors = FALSE)
coding <- annotation$gene_symbol[!is.na(annotation$gene_biotype) &
                                   annotation$gene_biotype == "protein_coding"]
universe <- intersect(intersect(rownames(expr), coding), de$gene_symbol)
de <- de[de$gene_symbol %in% universe, ]

fg_genes <- de$gene_symbol[de$fdr < FDR_CUTOFF & de$log_fc < 0]
bg_genes <- setdiff(universe, de$gene_symbol[de$fdr < FDR_CUTOFF])

standardize <- function(gr) {
  gr <- resize(resize(gr, width = 1, fix = "center"),
               width = REGION_WIDTH, fix = "center")
  unique(gr)
}

promoter_regions <- function(genes) {
  idx <- which(annotation$gene_symbol %in% genes)
  standardize(subsetByOverlaps(reg$peaks$h3k4me3, promoters_gr[idx]))
}

fg_gr <- promoter_regions(fg_genes)
bg_gr <- promoter_regions(bg_genes)
bg_gr <- bg_gr[!overlapsAny(bg_gr, fg_gr)]

get_clean_seqs <- function(gr) {
  s <- getSeq(fa, gr)
  s[!grepl("N", as.character(s), fixed = TRUE)]
}

fg_seqs <- get_clean_seqs(fg_gr)
bg_seqs <- get_clean_seqs(bg_gr)

gc_fraction <- function(s) {
  f <- letterFrequency(s, letters = c("G", "C"), as.prob = FALSE)
  rowSums(f) / width(s)
}

fg_gc <- gc_fraction(fg_seqs)
bg_gc <- gc_fraction(bg_seqs)

gc_summary <- data.frame(
  group = c("foreground", "background"),
  n = c(length(fg_gc), length(bg_gc)),
  min = round(c(min(fg_gc), min(bg_gc)), 3),
  q1 = round(c(quantile(fg_gc, 0.25), quantile(bg_gc, 0.25)), 3),
  median = round(c(median(fg_gc), median(bg_gc)), 3),
  q3 = round(c(quantile(fg_gc, 0.75), quantile(bg_gc, 0.75)), 3),
  max = round(c(max(fg_gc), max(bg_gc)), 3),
  frac_above_0.8 = round(c(mean(fg_gc > 0.8), mean(bg_gc > 0.8)), 4),
  frac_below_0.2 = round(c(mean(fg_gc < 0.2), mean(bg_gc < 0.2)), 4),
  stringsAsFactors = FALSE
)
write.csv(gc_summary, file.path(OUT_DIR, "gc_summary.csv"), row.names = FALSE)
print(gc_summary)

bin_table <- function(breaks) {
  fg_bin <- cut(fg_gc, breaks = breaks, include.lowest = TRUE)
  bg_bin <- cut(bg_gc, breaks = breaks, include.lowest = TRUE)
  tab <- data.frame(
    bin = levels(fg_bin),
    n_foreground = as.integer(table(fg_bin)),
    n_background = as.integer(table(bg_bin)),
    stringsAsFactors = FALSE
  )
  tab$usable <- tab$n_foreground > 0 & tab$n_background > 0
  tab$dropped_foreground <- ifelse(tab$usable, 0L, tab$n_foreground)
  tab
}

default_bins <- bin_table(DEFAULT_BREAKS)
extended_bins <- bin_table(EXTENDED_BREAKS)

cat("\ndefault GCbreaks:\n")
print(default_bins)
cat(sprintf("foreground outside break range: %d of %d\n",
            length(fg_gc) - sum(default_bins$n_foreground), length(fg_gc)))

cat("\nextended GCbreaks:\n")
print(extended_bins)
cat(sprintf("foreground outside break range: %d of %d\n",
            length(fg_gc) - sum(extended_bins$n_foreground), length(fg_gc)))

write.csv(default_bins, file.path(OUT_DIR, "bins_default.csv"), row.names = FALSE)
write.csv(extended_bins, file.path(OUT_DIR, "bins_extended.csv"), row.names = FALSE)

pwms <- toPWM(getMatrixSet(JASPAR2020,
                           list(species = 9606, collection = "CORE",
                                matrixtype = "PFM")))

seqs <- c(fg_seqs, bg_seqs)
bins <- factor(c(rep("foreground", length(fg_seqs)),
                 rep("background", length(bg_seqs))),
               levels = c("background", "foreground"))

cat("\nattempting enrichment with extended GCbreaks\n")
res <- calcBinnedMotifEnrR(seqs = seqs, bins = bins, pwmL = pwms,
                           background = "otherBins",
                           GCbreaks = EXTENDED_BREAKS,
                           verbose = TRUE)

print(class(res))
print(dim(res))
print(SummarizedExperiment::assayNames(res))
print(colnames(SummarizedExperiment::rowData(res)))
print(colnames(res))
print(head(SummarizedExperiment::assay(res, 1)))
