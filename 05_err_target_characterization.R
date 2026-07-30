library(GenomicRanges)
library(Rsamtools)
library(Biostrings)
library(monaLisa)
library(TFBSTools)
library(JASPAR2020)
library(gprofiler2)

FASTA_PATH <- "data/raw/annotation/GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna"
REG_OBJ <- "data/processed/regulatory_objects.rds"
EXPR_OBJ <- "data/processed/expression_objects.rds"
DE_DIR <- "results/01_differential_expression"
PEAK_DIR <- "data/raw/chip_peaks"
OUT_DIR <- "results/05_err_target_characterization"

COMPARTMENT <- "cortex"
CONTRAST <- "late_vs_control"
FDR_CUTOFF <- 0.10
REGION_WIDTH <- 300
MIN_SCORE <- 10
TARGET_MOTIFS <- c("ESRRA", "ESRRB", "ESRRG")
CONTROL_MOTIFS <- c("CTCF", "SP1")
PROMOTER_UPSTREAM <- 2000
PROMOTER_DOWNSTREAM <- 500

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

reg <- readRDS(REG_OBJ)
expr_objects <- readRDS(EXPR_OBJ)
fa <- FaFile(FASTA_PATH)

annotation <- reg$annotation
genes_gr <- reg$genes[match(annotation$gene_symbol,
                            mcols(reg$genes)$gene_symbol)]
promoters_gr <- promoters(genes_gr, upstream = PROMOTER_UPSTREAM,
                          downstream = PROMOTER_DOWNSTREAM)
mcols(promoters_gr) <- data.frame(gene_symbol = annotation$gene_symbol,
                                  stringsAsFactors = FALSE)

load_summits <- function(dir, pattern, seq_map, valid_seqlevels) {
  hit <- grep(pattern, list.files(dir), value = TRUE, ignore.case = TRUE)
  stopifnot(length(hit) == 1L)
  tab <- read.delim(file.path(dir, hit), header = FALSE,
                    stringsAsFactors = FALSE)
  accession <- unname(seq_map[sub("^chr", "", tab$V1)])
  ok <- !is.na(accession) & accession %in% valid_seqlevels & tab$V10 >= 0
  GRanges(seqnames = accession[ok],
          ranges = IRanges(start = tab$V2[ok] + tab$V10[ok] + 1L, width = 1L))
}

summits <- load_summits(PEAK_DIR, "kidney-h3k4me3.*narrowPeak", reg$seq_map,
                        unique(as.character(seqnames(reg$genes))))
regions <- resize(summits, width = REGION_WIDTH, fix = "center")

pairs <- findOverlaps(regions, promoters_gr)
region_gene <- data.frame(
  region_index = queryHits(pairs),
  gene_symbol = mcols(promoters_gr)$gene_symbol[subjectHits(pairs)],
  stringsAsFactors = FALSE
)

used_regions <- sort(unique(region_gene$region_index))
region_gr <- regions[used_regions]
names(region_gr) <- sprintf("%s:%d-%d", as.character(seqnames(region_gr)),
                            start(region_gr), end(region_gr))
region_key_lookup <- setNames(names(region_gr), used_regions)
region_gene$region_key <- region_key_lookup[as.character(region_gene$region_index)]

seqs <- getSeq(fa, region_gr)
names(seqs) <- names(region_gr)
clean <- !grepl("N", as.character(seqs), fixed = TRUE)
seqs <- seqs[clean]
message(sprintf("%d promoter regions scanned", length(seqs)))

pfm <- getMatrixSet(JASPAR2020, list(species = 9606, collection = "CORE",
                                     matrixtype = "PFM"))
pwms <- toPWM(pfm)
motif_names <- vapply(seq_along(pwms), function(i) name(pwms[[i]]), character(1))

select_motifs <- function(targets) {
  idx <- which(toupper(motif_names) %in% toupper(targets))
  if (!length(idx)) stop("no motif matched: ", paste(targets, collapse = ", "))
  message(sprintf("matched motifs: %s",
                  paste(motif_names[idx], collapse = ", ")))
  pwms[idx]
}

scan_motifs <- function(pwm_subset, seqs) {
  hits <- findMotifHits(query = pwm_subset, subject = seqs,
                        min.score = MIN_SCORE, method = "matchPWM")
  message("findMotifHits columns: ",
          paste(colnames(mcols(hits)), collapse = ", "))
  data.frame(region_key = as.character(seqnames(hits)),
             motif = as.character(mcols(hits)$pwmname),
             score = mcols(hits)$score,
             stringsAsFactors = FALSE)
}

err_hits <- scan_motifs(select_motifs(TARGET_MOTIFS), seqs)
control_hits <- scan_motifs(select_motifs(CONTROL_MOTIFS), seqs)

write.csv(err_hits, file.path(OUT_DIR, "err_motif_hits.csv"), row.names = FALSE)

de <- read.csv(file.path(DE_DIR, sprintf("de_%s_%s.csv", COMPARTMENT, CONTRAST)),
               stringsAsFactors = FALSE)
expr <- expr_objects[[COMPARTMENT]]$expr
coding <- annotation$gene_symbol[!is.na(annotation$gene_biotype) &
                                   annotation$gene_biotype == "protein_coding"]
universe <- intersect(intersect(rownames(expr), coding), de$gene_symbol)
de <- de[de$gene_symbol %in% universe, ]

down_genes <- de$gene_symbol[de$fdr < FDR_CUTOFF & de$log_fc < 0]

genes_with_motif <- function(hits) {
  unique(region_gene$gene_symbol[region_gene$region_key %in%
                                   unique(hits$region_key)])
}

err_genes <- genes_with_motif(err_hits)
control_genes <- genes_with_motif(control_hits)

classified <- data.frame(
  gene_symbol = down_genes,
  err_positive = down_genes %in% err_genes,
  control_positive = down_genes %in% control_genes,
  stringsAsFactors = FALSE
)
classified <- merge(classified, de[, c("gene_symbol", "log_fc", "fdr",
                                       "ave_expr")], by = "gene_symbol")
write.csv(classified, file.path(OUT_DIR, "repressed_genes_by_motif.csv"),
          row.names = FALSE)

message(sprintf("repressed genes: %d total, %d with ERR site, %d without",
                nrow(classified), sum(classified$err_positive),
                sum(!classified$err_positive)))

effect_test <- wilcox.test(log_fc ~ err_positive, data = classified)
effect_summary <- aggregate(log_fc ~ err_positive, data = classified,
                            FUN = function(v) c(median = median(v), n = length(v)))
print(effect_summary)
message(sprintf("effect size difference, Wilcoxon p = %.4g",
                effect_test$p.value))

run_gost <- function(query, background, label) {
  if (length(query) < 10) {
    message(sprintf("%s: only %d genes, enrichment skipped", label,
                    length(query)))
    return(NULL)
  }
  res <- gost(query = toupper(query), organism = "hsapiens",
              custom_bg = toupper(background), correction_method = "gSCS",
              sources = c("GO:BP", "GO:CC", "KEGG", "REAC"),
              evcodes = FALSE)
  if (is.null(res)) {
    message(sprintf("%s: no enriched term", label))
    return(NULL)
  }
  out <- res$result[, c("source", "term_id", "term_name", "p_value",
                        "term_size", "query_size", "intersection_size")]
  out$comparison <- label
  out[order(out$p_value), ]
}

err_vs_rest <- run_gost(classified$gene_symbol[classified$err_positive],
                        classified$gene_symbol, "err_positive_vs_repressed")
control_vs_rest <- run_gost(classified$gene_symbol[classified$control_positive],
                            classified$gene_symbol,
                            "control_positive_vs_repressed")

enrichment <- do.call(rbind, Filter(Negate(is.null),
                                    list(err_vs_rest, control_vs_rest)))
if (!is.null(enrichment)) {
  write.csv(enrichment, file.path(OUT_DIR, "functional_enrichment.csv"),
            row.names = FALSE)
  print(head(err_vs_rest[, c("source", "term_name", "p_value",
                             "intersection_size")], 25))
}

saveRDS(list(classified = classified, err_hits = err_hits,
             enrichment = enrichment, region_gene = region_gene),
        file.path(OUT_DIR, "err_objects.rds"))
