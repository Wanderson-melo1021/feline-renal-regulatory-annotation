library(rtracklayer)
library(GenomicRanges)

GFF_PATH <- "data/raw/annotation/GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.gff.gz"
REPORT_PATH <- "data/raw/annotation/GCF_018350175.1_F.catus_Fca126_mat1.0_assembly_report.txt"
PEAK_DIR <- "data/raw/chip_peaks"
EXPR_OBJ <- "data/processed/expression_objects.rds"
OUT_DIR <- "results/02_regulatory_annotation"
OBJ_DIR <- "data/processed"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJ_DIR, recursive = TRUE, showWarnings = FALSE)

PROMOTER_UPSTREAM <- 2000
PROMOTER_DOWNSTREAM <- 500
MAX_DOMAIN_HALF_WIDTH <- 250000
PEAK_WIDTH_QUANTILE <- 0.995
MIN_CHROM_LENGTH <- 1e7
N_PERMUTATIONS <- 100
MIN_H3K4ME3_FOLD <- 5
SEED <- 1

MARK_PATTERNS <- c(h3k4me3 = "kidney-h3k4me3.*narrowPeak",
                   h3k27ac = "kidney-h3k27ac.*narrowPeak",
                   ctcf    = "kidney-ctcf.*narrowPeak")

discover_peak_files <- function(dir, patterns) {
  available <- list.files(dir)
  if (!length(available)) stop("no files found in ", dir)
  paths <- vapply(names(patterns), function(mark) {
    hit <- grep(patterns[[mark]], available, value = TRUE, ignore.case = TRUE)
    if (length(hit) != 1L) {
      stop(sprintf("expected exactly one file for %s, found %d in %s:\n%s",
                   mark, length(hit), dir, paste(available, collapse = "\n")))
    }
    file.path(dir, hit)
  }, character(1))
  message("peak files:\n", paste(paths, collapse = "\n"))
  paths
}

read_assembly_map <- function(path) {
  lines <- readLines(path)
  body <- lines[!startsWith(lines, "#") & nzchar(lines)]
  tab <- read.delim(text = paste(body, collapse = "\n"), header = FALSE,
                    stringsAsFactors = FALSE, quote = "", comment.char = "")

  refseq <- as.character(tab[[7]])
  name_cols <- intersect(c(1, 3, 10), seq_len(ncol(tab)))

  pairs <- do.call(rbind, lapply(name_cols, function(i) {
    data.frame(key = sub("^chr", "", as.character(tab[[i]])),
               value = refseq, stringsAsFactors = FALSE)
  }))

  pairs <- pairs[!is.na(pairs$key) & nzchar(pairs$key) &
                   pairs$key != "na" & pairs$value != "na", ]
  pairs <- pairs[!duplicated(pairs$key), ]

  message(sprintf("assembly map: %d keys, e.g. %s", nrow(pairs),
                  paste(sprintf("%s=%s", head(pairs$key, 3),
                                head(pairs$value, 3)), collapse = ", ")))

  setNames(pairs$value, pairs$key)
}

load_gene_annotation <- function(path) {
  gff <- import(path, format = "gff", feature.type = "gene")
  meta <- mcols(gff)
  symbol <- if ("gene" %in% colnames(meta)) as.character(meta$gene) else
    as.character(meta$Name)
  biotype <- if ("gene_biotype" %in% colnames(meta))
    as.character(meta$gene_biotype) else rep(NA_character_, length(gff))
  keep <- !is.na(symbol) & nzchar(symbol)
  out <- gff[keep]
  mcols(out) <- data.frame(gene_symbol = symbol[keep],
                           gene_biotype = biotype[keep],
                           stringsAsFactors = FALSE)
  out
}

load_peaks <- function(path, seq_map, valid_seqlevels, mark) {
  tab <- read.delim(path, header = FALSE, stringsAsFactors = FALSE)
  n_total <- nrow(tab)

  chrom_name <- sub("^chr", "", tab$V1)
  accession <- unname(seq_map[chrom_name])
  ok <- !is.na(accession) & accession %in% valid_seqlevels

  if (!any(ok)) {
    stop(sprintf("no peak mapped to annotation for %s.\npeak names: %s\nmap keys: %s",
                 mark, paste(head(sort(unique(chrom_name)), 8), collapse = ", "),
                 paste(head(names(seq_map), 8), collapse = ", ")))
  }

  tab <- tab[ok, , drop = FALSE]
  gr <- GRanges(seqnames = accession[ok],
                ranges = IRanges(start = tab$V2 + 1L, end = tab$V3),
                score = tab$V5,
                fold_change = tab$V7)
  n_mapped <- length(gr)
  cutoff <- quantile(width(gr), PEAK_WIDTH_QUANTILE)
  gr <- gr[width(gr) <= cutoff]
  attr(gr, "counts") <- c(total = n_total, mapped = n_mapped,
                          retained = length(gr))
  sort(gr)
}

build_ctcf_domains <- function(genes, ctcf) {
  tss <- start(resize(genes, width = 1, fix = "start"))
  chrom <- as.character(seqnames(genes))
  mid <- (start(ctcf) + end(ctcf)) %/% 2L
  sites_by_chrom <- lapply(split(mid, as.character(seqnames(ctcf))), sort)

  lower <- numeric(length(genes))
  upper <- numeric(length(genes))

  for (chr in unique(chrom)) {
    idx <- which(chrom == chr)
    sites <- sites_by_chrom[[chr]]
    if (is.null(sites) || !length(sites)) {
      lower[idx] <- tss[idx] - MAX_DOMAIN_HALF_WIDTH
      upper[idx] <- tss[idx] + MAX_DOMAIN_HALF_WIDTH
      next
    }
    pos <- findInterval(tss[idx], sites)
    left <- rep(NA_real_, length(idx))
    has_left <- pos >= 1L
    left[has_left] <- sites[pos[has_left]]
    right <- rep(NA_real_, length(idx))
    has_right <- (pos + 1L) <= length(sites)
    right[has_right] <- sites[pos[has_right] + 1L]
    lower[idx] <- ifelse(is.na(left), tss[idx] - MAX_DOMAIN_HALF_WIDTH, left)
    upper[idx] <- ifelse(is.na(right), tss[idx] + MAX_DOMAIN_HALF_WIDTH, right)
  }

  lower <- pmax(pmax(lower, tss - MAX_DOMAIN_HALF_WIDTH), 1)
  upper <- pmin(upper, tss + MAX_DOMAIN_HALF_WIDTH)
  upper <- pmax(upper, lower + 1)

  gr <- GRanges(seqnames = chrom,
                ranges = IRanges(start = as.integer(lower),
                                 end = as.integer(upper)))
  mcols(gr) <- mcols(genes)
  gr
}

shuffle_within_chromosome <- function(gr, lengths) {
  chrom <- as.character(seqnames(gr))
  w <- width(gr)
  limits <- pmax(lengths[chrom] - w, 1)
  GRanges(seqnames = chrom,
          ranges = IRanges(start = floor(runif(length(gr)) * limits) + 1,
                           width = w))
}

seq_map <- read_assembly_map(REPORT_PATH)
peak_paths <- discover_peak_files(PEAK_DIR, MARK_PATTERNS)

genes <- load_gene_annotation(GFF_PATH)
chrom_length <- tapply(end(genes), as.character(seqnames(genes)), max)
main_seqs <- names(chrom_length)[chrom_length > MIN_CHROM_LENGTH]
genes <- genes[as.character(seqnames(genes)) %in% main_seqs]
chrom_length <- chrom_length[main_seqs]

message(sprintf("annotation: %d genes on %d chromosome-level sequences",
                length(genes), length(main_seqs)))

peaks <- lapply(names(peak_paths), function(m) {
  load_peaks(peak_paths[[m]], seq_map, main_seqs, m)
})
names(peaks) <- names(peak_paths)

peak_report <- do.call(rbind, lapply(names(peaks), function(m) {
  cnt <- attr(peaks[[m]], "counts")
  data.frame(mark = m, n_total = cnt[["total"]], n_mapped = cnt[["mapped"]],
             n_retained = cnt[["retained"]],
             pct_retained = round(100 * cnt[["retained"]] / cnt[["total"]], 2),
             stringsAsFactors = FALSE)
}))
write.csv(peak_report, file.path(OUT_DIR, "peak_harmonization_report.csv"),
          row.names = FALSE)
print(peak_report)

prom <- promoters(genes, upstream = PROMOTER_UPSTREAM,
                  downstream = PROMOTER_DOWNSTREAM)
domains <- build_ctcf_domains(genes, peaks$ctcf)

set.seed(SEED)
enrichment <- do.call(rbind, lapply(names(peaks), function(m) {
  observed <- sum(countOverlaps(prom, peaks[[m]]) > 0)
  null <- vapply(seq_len(N_PERMUTATIONS), function(i) {
    sum(countOverlaps(prom, shuffle_within_chromosome(peaks[[m]],
                                                      chrom_length)) > 0)
  }, numeric(1))
  data.frame(mark = m, observed_genes_hit = observed,
             observed_fraction = round(observed / length(prom), 4),
             expected_genes_hit = round(mean(null), 1),
             fold_enrichment = round(observed / mean(null), 2),
             stringsAsFactors = FALSE)
}))
write.csv(enrichment, file.path(OUT_DIR, "promoter_enrichment_check.csv"),
          row.names = FALSE)
print(enrichment)

h3k4me3_fold <- enrichment$fold_enrichment[enrichment$mark == "h3k4me3"]
if (h3k4me3_fold < MIN_H3K4ME3_FOLD) {
  stop(sprintf("H3K4me3 promoter enrichment is %.2f-fold, below the %d-fold floor: coordinates still inconsistent",
               h3k4me3_fold, MIN_H3K4ME3_FOLD))
}

annotation <- data.frame(
  gene_symbol = mcols(genes)$gene_symbol,
  gene_biotype = mcols(genes)$gene_biotype,
  seqnames = as.character(seqnames(genes)),
  start = start(genes),
  end = end(genes),
  strand = as.character(strand(genes)),
  domain_start = start(domains),
  domain_end = end(domains),
  domain_width = width(domains),
  stringsAsFactors = FALSE
)

for (mark in names(peaks)) {
  annotation[[paste0("promoter_", mark)]] <- countOverlaps(prom, peaks[[mark]])
  annotation[[paste0("domain_", mark)]] <- countOverlaps(domains, peaks[[mark]])
}

annotation$distal_h3k27ac <- pmax(annotation$domain_h3k27ac -
                                    annotation$promoter_h3k27ac, 0L)

annotation$regulatory_class <- with(annotation, ifelse(
  promoter_h3k4me3 > 0 & promoter_h3k27ac > 0, "active_promoter",
  ifelse(promoter_h3k4me3 > 0, "primed_promoter",
  ifelse(distal_h3k27ac > 0, "distal_enhancer_only", "no_active_mark"))))

expr_objects <- readRDS(EXPR_OBJ)
mapping_report <- do.call(rbind, lapply(names(expr_objects), function(cp) {
  symbols <- rownames(expr_objects[[cp]]$expr)
  matched <- symbols %in% annotation$gene_symbol
  data.frame(compartment = cp, n_expressed = length(symbols),
             n_mapped = sum(matched),
             pct_mapped = round(100 * mean(matched), 2),
             stringsAsFactors = FALSE)
}))
write.csv(mapping_report, file.path(OUT_DIR, "symbol_mapping_report.csv"),
          row.names = FALSE)
print(mapping_report)

writeLines(sort(setdiff(rownames(expr_objects$cortex$expr),
                        annotation$gene_symbol)),
           file.path(OUT_DIR, "unmapped_symbols_cortex.txt"))

annotation <- annotation[!duplicated(annotation$gene_symbol), ]
write.csv(annotation, file.path(OUT_DIR, "gene_regulatory_annotation.csv"),
          row.names = FALSE)

class_table <- as.data.frame(table(annotation$regulatory_class),
                             stringsAsFactors = FALSE)
colnames(class_table) <- c("regulatory_class", "n_genes")
write.csv(class_table, file.path(OUT_DIR, "regulatory_class_counts.csv"),
          row.names = FALSE)
print(class_table)

cat("\ndomain width distribution:\n")
print(summary(annotation$domain_width))

saveRDS(list(annotation = annotation, peaks = peaks, genes = genes,
             domains = domains, seq_map = seq_map),
        file.path(OBJ_DIR, "regulatory_objects.rds"))
