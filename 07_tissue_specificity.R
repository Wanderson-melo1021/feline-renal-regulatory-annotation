library(GenomicRanges)

REG_OBJ <- "data/processed/regulatory_objects.rds"
EXPR_OBJ <- "data/processed/expression_objects.rds"
DE_DIR <- "results/01_differential_expression"
PEAK_DIR <- "data/raw/chip_peaks"
OUT_DIR <- "results/07_tissue_specificity"

COMPARTMENT <- "cortex"
CONTRAST <- "late_vs_control"
FDR_CUTOFF <- 0.10
PROMOTER_UPSTREAM <- 2000
PROMOTER_DOWNSTREAM <- 500
REFERENCE_TISSUE <- "kidney"
TISSUES <- c("heart", "kidney", "liver", "lung", "smallintestine",
             "testis", "thyroid")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

reg <- readRDS(REG_OBJ)
expr_objects <- readRDS(EXPR_OBJ)
annotation <- reg$annotation
seq_map <- reg$seq_map
valid_seqlevels <- unique(as.character(seqnames(reg$genes)))

load_peaks <- function(tissue, mark) {
  pattern <- sprintf("%s-%s.*narrowPeak", tissue, mark)
  hit <- grep(pattern, list.files(PEAK_DIR), value = TRUE, ignore.case = TRUE)
  if (length(hit) != 1L) {
    stop(sprintf("expected one file for %s/%s, found %d", tissue, mark,
                 length(hit)))
  }
  tab <- read.delim(file.path(PEAK_DIR, hit), header = FALSE,
                    stringsAsFactors = FALSE)
  accession <- unname(seq_map[sub("^chr", "", tab$V1)])
  ok <- !is.na(accession) & accession %in% valid_seqlevels
  sort(GRanges(seqnames = accession[ok],
               ranges = IRanges(start = tab$V2[ok] + 1L, end = tab$V3[ok])))
}

score_specificity <- function(mark) {
  peaks <- lapply(TISSUES, load_peaks, mark = mark)
  names(peaks) <- TISSUES

  reference <- peaks[[REFERENCE_TISSUE]]
  presence <- vapply(TISSUES, function(t) {
    as.integer(overlapsAny(reference, peaks[[t]]))
  }, integer(length(reference)))

  mcols(reference)$n_tissues <- rowSums(presence)
  mcols(reference)$specificity <- cut(
    mcols(reference)$n_tissues,
    breaks = c(0, 1, 3, 5, 7),
    labels = c("kidney_specific", "restricted", "shared", "ubiquitous"))
  reference
}

message("scoring H3K27ac")
enhancers <- score_specificity("h3k27ac")
message("scoring H3K4me3")
promoter_marks <- score_specificity("h3k4me3")

specificity_summary <- do.call(rbind, lapply(
  list(h3k27ac = enhancers, h3k4me3 = promoter_marks), function(gr) {
    tab <- table(mcols(gr)$specificity)
    data.frame(class = names(tab), n = as.integer(tab),
               pct = round(100 * as.integer(tab) / length(gr), 2),
               stringsAsFactors = FALSE)
  }))
specificity_summary$mark <- rep(c("h3k27ac", "h3k4me3"), each = 4)
write.csv(specificity_summary, file.path(OUT_DIR, "element_specificity.csv"),
          row.names = FALSE)
print(specificity_summary)

genes_gr <- reg$genes[match(annotation$gene_symbol,
                            mcols(reg$genes)$gene_symbol)]
promoters_gr <- promoters(genes_gr, upstream = PROMOTER_UPSTREAM,
                          downstream = PROMOTER_DOWNSTREAM)
domains_gr <- reg$domains[match(annotation$gene_symbol,
                                mcols(reg$domains)$gene_symbol)]

distal <- enhancers[!overlapsAny(enhancers, promoters_gr)]

gene_pairs <- findOverlaps(domains_gr, distal)
gene_enhancers <- data.frame(
  gene_symbol = annotation$gene_symbol[queryHits(gene_pairs)],
  n_tissues = mcols(distal)$n_tissues[subjectHits(gene_pairs)],
  kidney_specific = mcols(distal)$n_tissues[subjectHits(gene_pairs)] == 1L,
  stringsAsFactors = FALSE
)

gene_level <- aggregate(cbind(n_tissues, kidney_specific) ~ gene_symbol,
                        data = gene_enhancers,
                        FUN = function(v) c(mean = mean(v), n = length(v)))
gene_level <- data.frame(
  gene_symbol = gene_level$gene_symbol,
  mean_tissues = round(gene_level$n_tissues[, "mean"], 3),
  frac_kidney_specific = round(gene_level$kidney_specific[, "mean"], 3),
  n_enhancers = as.integer(gene_level$n_tissues[, "n"]),
  stringsAsFactors = FALSE
)

expr <- expr_objects[[COMPARTMENT]]$expr
de <- read.csv(file.path(DE_DIR, sprintf("de_%s_%s.csv", COMPARTMENT, CONTRAST)),
               stringsAsFactors = FALSE)
coding <- annotation$gene_symbol[!is.na(annotation$gene_biotype) &
                                   annotation$gene_biotype == "protein_coding"]
universe <- intersect(intersect(rownames(expr), coding), de$gene_symbol)
de <- de[de$gene_symbol %in% universe, ]

dat <- merge(de[, c("gene_symbol", "log_fc", "fdr")], gene_level,
             by = "gene_symbol")
dat$state <- with(dat, ifelse(fdr < FDR_CUTOFF & log_fc < 0, "repressed",
                       ifelse(fdr < FDR_CUTOFF & log_fc > 0, "induced",
                              "unchanged")))
dat$state <- factor(dat$state, levels = c("unchanged", "repressed", "induced"))
write.csv(dat, file.path(OUT_DIR, "gene_enhancer_specificity.csv"),
          row.names = FALSE)

state_summary <- aggregate(cbind(mean_tissues, frac_kidney_specific,
                                 n_enhancers) ~ state, data = dat,
                           FUN = median)
state_summary$n_genes <- as.integer(table(dat$state)[state_summary$state])
write.csv(state_summary, file.path(OUT_DIR, "specificity_by_state.csv"),
          row.names = FALSE)
cat("\nmedian enhancer specificity by expression state\n")
print(state_summary)

kruskal_tissues <- kruskal.test(mean_tissues ~ state, data = dat)
kruskal_specific <- kruskal.test(frac_kidney_specific ~ state, data = dat)
cat(sprintf("\nmean_tissues across states: Kruskal p = %.3g\n",
            kruskal_tissues$p.value))
cat(sprintf("frac_kidney_specific across states: Kruskal p = %.3g\n",
            kruskal_specific$p.value))

pairwise <- do.call(rbind, lapply(c("repressed", "induced"), function(s) {
  sub <- dat[dat$state %in% c("unchanged", s), ]
  wt <- wilcox.test(mean_tissues ~ droplevels(sub$state), data = sub)
  ws <- wilcox.test(frac_kidney_specific ~ droplevels(sub$state), data = sub)
  data.frame(comparison = sprintf("%s vs unchanged", s),
             median_tissues = median(dat$mean_tissues[dat$state == s]),
             median_tissues_unchanged = median(
               dat$mean_tissues[dat$state == "unchanged"]),
             p_tissues = signif(wt$p.value, 3),
             median_specific = median(dat$frac_kidney_specific[dat$state == s]),
             median_specific_unchanged = median(
               dat$frac_kidney_specific[dat$state == "unchanged"]),
             p_specific = signif(ws$p.value, 3),
             stringsAsFactors = FALSE)
}))
write.csv(pairwise, file.path(OUT_DIR, "pairwise_tests.csv"), row.names = FALSE)
print(pairwise)

model_repressed <- glm(I(state == "repressed") ~ frac_kidney_specific +
                         log(n_enhancers + 1),
                       data = dat[dat$state != "induced", ], family = binomial())
model_induced <- glm(I(state == "induced") ~ frac_kidney_specific +
                       log(n_enhancers + 1),
                     data = dat[dat$state != "repressed", ], family = binomial())

extract_coef <- function(model, label) {
  out <- as.data.frame(coef(summary(model)))
  out$term <- rownames(out)
  out$odds_ratio <- round(exp(out$Estimate), 3)
  out$outcome <- label
  rownames(out) <- NULL
  out[, c("outcome", "term", "Estimate", "odds_ratio", "Std. Error",
          "Pr(>|z|)")]
}

models <- rbind(extract_coef(model_repressed, "repressed"),
                extract_coef(model_induced, "induced"))
write.csv(models, file.path(OUT_DIR, "logistic_models.csv"), row.names = FALSE)
cat("\nassociation between enhancer tissue specificity and expression state\n")
print(models)

saveRDS(list(enhancers = enhancers, promoter_marks = promoter_marks,
             gene_level = gene_level, data = dat, models = models),
        file.path(OUT_DIR, "specificity_objects.rds"))
