library(readxl)
library(limma)

XLSX_PATH <- "data/raw/expression/42003_2025_9164_MOESM3_ESM.xlsx"
EXPR_OBJ <- "data/processed/expression_objects.rds"
CHAR_DIR <- "results/03_modules/cortex_late_vs_control"
OUT_DIR <- "results/09_cohort_and_age"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SHEETS <- c(cortex = "S3 RNA cortex", medulla = "S4 RNA medulla")
META_ROWS <- c(id = 2, group = 3, sex = 4, age = 5, breed = 6)
GROUP_LEVELS <- c("Control", "CKD 1/2", "CKD 3/4")
GROUP_RECODE <- c("Control" = "Control", "CKD1/2" = "CKD 1/2", "CKD3/4" = "CKD 3/4")
CONTRASTS <- c(cortex = "CKD34 - Control", medulla = "CKD12 - Control")
FDR_CUTOFF <- 0.10

to_numeric <- function(x) as.numeric(gsub(",", ".", as.character(x)))

read_metadata <- function(path, sheet) {
  raw <- as.data.frame(read_excel(path, sheet = sheet, col_names = FALSE,
                                  .name_repair = "minimal", n_max = 6,
                                  progress = FALSE), stringsAsFactors = FALSE)
  ids <- as.character(unlist(raw[META_ROWS["id"], -1]))
  keep <- !is.na(ids) & nzchar(ids)
  data.frame(
    sample_id = ids[keep],
    group = factor(unname(GROUP_RECODE[as.character(unlist(raw[META_ROWS["group"], -1]))[keep]]),
                   levels = GROUP_LEVELS),
    sex = factor(as.character(unlist(raw[META_ROWS["sex"], -1]))[keep]),
    age = to_numeric(unlist(raw[META_ROWS["age"], -1]))[keep],
    breed = factor(as.character(unlist(raw[META_ROWS["breed"], -1]))[keep]),
    stringsAsFactors = FALSE
  )
}

meta <- lapply(names(SHEETS), function(cp) {
  m <- read_metadata(XLSX_PATH, SHEETS[[cp]])
  m$compartment <- cp
  m
})
names(meta) <- names(SHEETS)

summarise_cohort <- function(m, compartment) {
  rows <- lapply(levels(m$group), function(g) {
    s <- m[m$group == g, ]
    data.frame(
      compartment = compartment,
      group = g,
      n = nrow(s),
      age_mean = round(mean(s$age), 1),
      age_sd = round(sd(s$age), 1),
      age_median = round(median(s$age), 1),
      age_min = round(min(s$age), 1),
      age_max = round(max(s$age), 1),
      sex_female = sum(s$sex == "Female"),
      sex_male = sum(s$sex == "Male"),
      breed_dsh = sum(s$breed == "DSH"),
      breed_dlh = sum(s$breed == "DLH"),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$age_p <- signif(kruskal.test(age ~ group, data = m)$p.value, 3)
  out$sex_p <- signif(fisher.test(table(m$sex, m$group))$p.value, 3)
  out$breed_p <- signif(fisher.test(table(m$breed, m$group))$p.value, 3)
  out
}

cohort <- do.call(rbind, lapply(names(meta), function(cp)
  summarise_cohort(meta[[cp]], cp)))
write.csv(cohort, file.path(OUT_DIR, "table1_cohort.csv"), row.names = FALSE)
print(cohort)

all_ids <- lapply(meta, function(m) m$sample_id)
pairing <- data.frame(
  cat_id = sort(union(all_ids$cortex, all_ids$medulla)),
  stringsAsFactors = FALSE)
pairing$cortex <- pairing$cat_id %in% all_ids$cortex
pairing$medulla <- pairing$cat_id %in% all_ids$medulla
pairing$both <- pairing$cortex & pairing$medulla
write.csv(pairing, file.path(OUT_DIR, "compartment_pairing.csv"), row.names = FALSE)
cat(sprintf("\n%d cats total, %d with both compartments, %d cortex only, %d medulla only\n",
            nrow(pairing), sum(pairing$both),
            sum(pairing$cortex & !pairing$medulla),
            sum(pairing$medulla & !pairing$cortex)))

expr_objects <- readRDS(EXPR_OBJ)

fit_contrast <- function(expr, samples, formula_rhs, contrast_string) {
  design <- model.matrix(formula_rhs, data = samples)
  colnames(design) <- make.names(sub("^group", "", colnames(design)))
  cm <- makeContrasts(contrasts = contrast_string, levels = design)
  colnames(cm) <- "target"
  fit <- eBayes(contrasts.fit(lmFit(expr, design), cm),
                trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "target", number = Inf, sort.by = "none")
  data.frame(gene_symbol = rownames(tt), log_fc = tt$logFC,
             fdr = tt$adj.P.Val, row.names = NULL, stringsAsFactors = FALSE)
}

concordance <- do.call(rbind, lapply(names(SHEETS), function(cp) {
  expr <- expr_objects[[cp]]$expr
  s <- expr_objects[[cp]]$samples
  m <- meta[[cp]]
  s$age <- m$age[match(s$sample_id, m$sample_id)]
  s$group <- factor(make.names(sub("/", "", as.character(s$group))))
  stopifnot(!any(is.na(s$age)))

  without_age <- fit_contrast(expr, s, ~ 0 + group + sex, CONTRASTS[[cp]])
  with_age <- fit_contrast(expr, s, ~ 0 + group + sex + age, CONTRASTS[[cp]])
  merged <- merge(with_age, without_age, by = "gene_symbol",
                  suffixes = c("_with_age", "_without_age"))

  sig_age <- merged$fdr_with_age < FDR_CUTOFF
  sig_base <- merged$fdr_without_age < FDR_CUTOFF

  write.csv(merged, file.path(OUT_DIR,
            sprintf("age_sensitivity_genes_%s.csv", cp)), row.names = FALSE)

  data.frame(
    compartment = cp,
    contrast = CONTRASTS[[cp]],
    n_primary_with_age = sum(sig_age),
    n_sensitivity_without_age = sum(sig_base),
    n_shared = sum(sig_base & sig_age),
    jaccard = round(sum(sig_base & sig_age) / sum(sig_base | sig_age), 3),
    logfc_spearman = round(cor(merged$log_fc_with_age,
                               merged$log_fc_without_age,
                               method = "spearman"), 4),
    stringsAsFactors = FALSE)
}))

write.csv(concordance, file.path(OUT_DIR, "age_sensitivity_concordance.csv"),
          row.names = FALSE)
cat("\nconcordance between the primary model (with age) and the sensitivity model (without age)\n")
print(concordance)

composition <- read.csv(file.path(CHAR_DIR, "composition_and_eigengenes.csv"),
                        stringsAsFactors = FALSE)
idx <- match(composition$sample_id, meta$cortex$sample_id)
composition$age <- meta$cortex$age[idx]
composition$group <- droplevels(meta$cortex$group[idx])
stopifnot(!any(is.na(composition$age)), !any(is.na(composition$group)),
          nlevels(composition$group) >= 2)

MARKER_COLUMNS <- intersect(c("epithelial", "tubular", "immune", "fibroblast"),
                            colnames(composition))
stopifnot(length(MARKER_COLUMNS) > 0)

age_composition <- do.call(rbind, lapply(MARKER_COLUMNS,
                                         function(cell) {
  ct <- suppressWarnings(cor.test(composition$age, composition[[cell]],
                                  method = "spearman"))
  full <- lm(composition[[cell]] ~ composition$age + composition$group)
  reduced <- lm(composition[[cell]] ~ composition$group)
  lrt <- anova(reduced, full)
  data.frame(marker_set = cell,
             spearman_rho = round(unname(ct$estimate), 3),
             spearman_p = signif(ct$p.value, 3),
             age_p_given_group = signif(lrt[["Pr(>F)"]][2], 3),
             stringsAsFactors = FALSE)
}))

write.csv(age_composition, file.path(OUT_DIR, "age_vs_composition.csv"),
          row.names = FALSE)
cat("\nage against composition scores, unadjusted and conditional on disease group\n")
print(age_composition)

cat("\nage by group, cortex\n")
print(aggregate(age ~ group, data = meta$cortex,
                FUN = function(v) c(n = length(v), mean = round(mean(v), 2),
                                    min = min(v), max = max(v))))
