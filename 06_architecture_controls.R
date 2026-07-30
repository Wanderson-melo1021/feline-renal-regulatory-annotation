library(org.Hs.eg.db)
library(AnnotationDbi)
library(gprofiler2)

ERR_OBJ <- "results/05_err_target_characterization/err_objects.rds"
REG_OBJ <- "data/processed/regulatory_objects.rds"
EXPR_OBJ <- "data/processed/expression_objects.rds"
DE_DIR <- "results/01_differential_expression"
OUT_DIR <- "results/06_architecture_controls"

COMPARTMENT <- "cortex"
CONTRAST <- "late_vs_control"
FDR_CUTOFF <- 0.10
MITO_GO <- "GO:0005739"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

err <- readRDS(ERR_OBJ)
reg <- readRDS(REG_OBJ)
expr_objects <- readRDS(EXPR_OBJ)

annotation <- reg$annotation
expr <- expr_objects[[COMPARTMENT]]$expr
de <- read.csv(file.path(DE_DIR, sprintf("de_%s_%s.csv", COMPARTMENT, CONTRAST)),
               stringsAsFactors = FALSE)

coding <- annotation$gene_symbol[!is.na(annotation$gene_biotype) &
                                   annotation$gene_biotype == "protein_coding"]
universe <- intersect(intersect(rownames(expr), coding), de$gene_symbol)
de <- de[de$gene_symbol %in% universe, ]

region_gene <- err$region_gene
err_regions <- unique(err$err_hits$region_key)

regions_per_gene <- as.data.frame(table(region_gene$gene_symbol),
                                  stringsAsFactors = FALSE)
colnames(regions_per_gene) <- c("gene_symbol", "n_regions")

err_genes <- unique(region_gene$gene_symbol[region_gene$region_key %in%
                                              err_regions])

mito <- AnnotationDbi::select(org.Hs.eg.db, keys = MITO_GO,
                              keytype = "GOALL", columns = "SYMBOL")
mito_symbols <- unique(toupper(mito$SYMBOL))
message(sprintf("mitochondrial reference set: %d symbols", length(mito_symbols)))

dat <- data.frame(gene_symbol = universe, stringsAsFactors = FALSE)
dat <- merge(dat, de[, c("gene_symbol", "log_fc", "fdr")], by = "gene_symbol")
dat <- merge(dat, regions_per_gene, by = "gene_symbol", all.x = TRUE)
dat$n_regions[is.na(dat$n_regions)] <- 0L
dat$err_site <- dat$gene_symbol %in% err_genes
dat$mitochondrial <- toupper(dat$gene_symbol) %in% mito_symbols
dat$repressed <- dat$fdr < FDR_CUTOFF & dat$log_fc < 0
dat$scanned <- dat$n_regions > 0

analysis <- dat[dat$scanned, ]
message(sprintf("%d genes with at least one scanned promoter region", nrow(analysis)))

odds_ratio_table <- function(df, label) {
  tab <- table(err_site = df$err_site, mitochondrial = df$mitochondrial)
  if (any(dim(tab) < 2)) return(NULL)
  ft <- fisher.test(tab)
  data.frame(stratum = label,
             n = nrow(df),
             n_err_site = sum(df$err_site),
             pct_mito_err_pos = round(100 * mean(df$mitochondrial[df$err_site]), 2),
             pct_mito_err_neg = round(100 * mean(df$mitochondrial[!df$err_site]), 2),
             odds_ratio = round(unname(ft$estimate), 3),
             ci_low = round(ft$conf.int[1], 3),
             ci_high = round(ft$conf.int[2], 3),
             p_value = signif(ft$p.value, 3),
             stringsAsFactors = FALSE)
}

strata <- rbind(
  odds_ratio_table(analysis[analysis$repressed, ], "repressed"),
  odds_ratio_table(analysis[!analysis$repressed, ], "not_repressed")
)
write.csv(strata, file.path(OUT_DIR, "control1_mito_odds_by_stratum.csv"),
          row.names = FALSE)
cat("\ncontrol 1: association between ERR site and mitochondrial annotation\n")
print(strata)

interaction_model <- glm(mitochondrial ~ err_site * repressed,
                         data = analysis, family = binomial())
interaction_summary <- as.data.frame(coef(summary(interaction_model)))
interaction_summary$term <- rownames(interaction_summary)
rownames(interaction_summary) <- NULL
write.csv(interaction_summary,
          file.path(OUT_DIR, "control1_interaction_model.csv"),
          row.names = FALSE)
cat("\ncontrol 1: does the ERR-mitochondria association differ between strata\n")
print(interaction_summary[, c("term", "Estimate", "Std. Error", "Pr(>|z|)")])

model_null <- glm(repressed ~ mitochondrial + log(n_regions + 1),
                  data = analysis, family = binomial())
model_full <- glm(repressed ~ err_site + mitochondrial + log(n_regions + 1),
                  data = analysis, family = binomial())
model_inter <- glm(repressed ~ err_site * mitochondrial + log(n_regions + 1),
                   data = analysis, family = binomial())

lrt_err <- anova(model_null, model_full, test = "LRT")
lrt_inter <- anova(model_full, model_inter, test = "LRT")

coef_full <- as.data.frame(coef(summary(model_full)))
coef_full$term <- rownames(coef_full)
coef_full$odds_ratio <- round(exp(coef_full$Estimate), 3)
rownames(coef_full) <- NULL

coef_inter <- as.data.frame(coef(summary(model_inter)))
coef_inter$term <- rownames(coef_inter)
coef_inter$odds_ratio <- round(exp(coef_inter$Estimate), 3)
rownames(coef_inter) <- NULL

write.csv(coef_full, file.path(OUT_DIR, "control2_additive_model.csv"),
          row.names = FALSE)
write.csv(coef_inter, file.path(OUT_DIR, "control2_interaction_model.csv"),
          row.names = FALSE)

cat("\ncontrol 2: repression as a function of ERR site, conditioning on mitochondrial annotation\n")
print(coef_full[, c("term", "Estimate", "odds_ratio", "Std. Error", "Pr(>|z|)")])
cat("\nlikelihood ratio test for the ERR term:\n")
print(lrt_err)
cat("\nlikelihood ratio test for the interaction:\n")
print(lrt_inter)
print(coef_inter[, c("term", "Estimate", "odds_ratio", "Pr(>|z|)")])

non_repressed <- analysis[!analysis$repressed, ]
gost_non_repressed <- gost(
  query = toupper(non_repressed$gene_symbol[non_repressed$err_site]),
  organism = "hsapiens",
  custom_bg = toupper(non_repressed$gene_symbol),
  correction_method = "gSCS",
  sources = c("GO:BP", "GO:CC", "KEGG", "REAC"))

if (!is.null(gost_non_repressed)) {
  out <- gost_non_repressed$result[, c("source", "term_id", "term_name",
                                       "p_value", "term_size",
                                       "intersection_size")]
  out <- out[order(out$p_value), ]
  write.csv(out, file.path(OUT_DIR, "control1_enrichment_non_repressed.csv"),
            row.names = FALSE)
  cat("\ncontrol 1: enrichment of ERR-positive genes among the non-repressed\n")
  print(head(out, 20))
} else {
  message("control 1: no enriched term among non-repressed ERR-positive genes")
}

saveRDS(list(data = analysis, strata = strata, coef_full = coef_full,
             coef_inter = coef_inter, lrt_err = lrt_err,
             lrt_inter = lrt_inter),
        file.path(OUT_DIR, "control_objects.rds"))
