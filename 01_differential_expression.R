library(readxl)
library(limma)

XLSX_PATH <- "data/raw/expression/42003_2025_9164_MOESM3_ESM.xlsx"
OUT_DIR <- "results/01_differential_expression"
OBJ_DIR <- "data/processed"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJ_DIR, recursive = TRUE, showWarnings = FALSE)

SHEETS <- c(cortex = "S3 RNA cortex", medulla = "S4 RNA medulla")
META_ROWS <- c(id = 2, group = 3, sex = 4)
FIRST_GENE_ROW <- 7

GROUP_LEVELS <- c("Control", "CKD12", "CKD34")
GROUP_RECODE <- c("Control" = "Control", "CKD1/2" = "CKD12", "CKD3/4" = "CKD34")

CONTRASTS <- c(
  early_vs_control = "CKD12 - Control",
  late_vs_control  = "CKD34 - Control",
  late_vs_early    = "CKD34 - CKD12"
)

read_compartment <- function(path, sheet) {
  raw <- read_excel(path, sheet = sheet, col_names = FALSE,
                    .name_repair = "minimal", progress = FALSE)
  raw <- as.data.frame(raw, stringsAsFactors = FALSE)

  sample_ids <- as.character(unlist(raw[META_ROWS["id"], -1]))
  keep <- !is.na(sample_ids) & nzchar(sample_ids)
  sample_ids <- sample_ids[keep]

  group_raw <- as.character(unlist(raw[META_ROWS["group"], -1]))[keep]
  sex_raw <- as.character(unlist(raw[META_ROWS["sex"], -1]))[keep]

  gene_block <- raw[FIRST_GENE_ROW:nrow(raw), , drop = FALSE]
  genes <- as.character(gene_block[[1]])
  valid <- !is.na(genes) & nzchar(genes)
  gene_block <- gene_block[valid, , drop = FALSE]
  genes <- genes[valid]

  expr <- as.matrix(gene_block[, -1, drop = FALSE][, keep, drop = FALSE])
  expr <- matrix(as.numeric(expr), nrow = nrow(expr),
                 dimnames = list(genes, sample_ids))

  samples <- data.frame(
    sample_id = sample_ids,
    group = factor(unname(GROUP_RECODE[group_raw]), levels = GROUP_LEVELS),
    sex = factor(sex_raw),
    stringsAsFactors = FALSE
  )

  stopifnot(!any(is.na(expr)), !any(is.na(samples$group)))
  list(expr = expr, samples = samples)
}

fit_compartment <- function(expr, samples) {
  design <- model.matrix(~ 0 + group + sex, data = samples)
  colnames(design) <- make.names(sub("^group", "", colnames(design)))

  cm <- makeContrasts(contrasts = unname(CONTRASTS), levels = design)
  colnames(cm) <- names(CONTRASTS)

  fit <- lmFit(expr, design)
  fit <- contrasts.fit(fit, cm)
  eBayes(fit, trend = TRUE, robust = TRUE)
}

extract_table <- function(fit, coef_name, compartment) {
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  data.frame(
    gene_symbol = rownames(tt),
    compartment = compartment,
    contrast = coef_name,
    log_fc = tt$logFC,
    ave_expr = tt$AveExpr,
    t_stat = tt$t,
    p_value = tt$P.Value,
    fdr = tt$adj.P.Val,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

results <- list()
objects <- list()

for (compartment in names(SHEETS)) {
  dat <- read_compartment(XLSX_PATH, SHEETS[[compartment]])
  fit <- fit_compartment(dat$expr, dat$samples)

  tabs <- lapply(names(CONTRASTS), extract_table, fit = fit,
                 compartment = compartment)
  names(tabs) <- names(CONTRASTS)

  for (nm in names(tabs)) {
    write.csv(tabs[[nm]],
              file.path(OUT_DIR, sprintf("de_%s_%s.csv", compartment, nm)),
              row.names = FALSE)
  }

  results[[compartment]] <- do.call(rbind, tabs)
  objects[[compartment]] <- list(expr = dat$expr, samples = dat$samples, fit = fit)

  message(sprintf("%s: %d genes, %d samples (%s)",
                  compartment, nrow(dat$expr), ncol(dat$expr),
                  paste(levels(dat$samples$group),
                        table(dat$samples$group), sep = "=", collapse = ", ")))
}

saveRDS(objects, file.path(OBJ_DIR, "expression_objects.rds"))
write.csv(do.call(rbind, results),
          file.path(OUT_DIR, "de_all_contrasts.csv"), row.names = FALSE)

summary_counts <- do.call(rbind, lapply(names(results), function(cp) {
  df <- results[[cp]]
  do.call(rbind, lapply(names(CONTRASTS), function(ct) {
    sub <- df[df$contrast == ct, ]
    data.frame(
      compartment = cp,
      contrast = ct,
      n_fdr_010 = sum(sub$fdr < 0.10),
      n_fdr_005 = sum(sub$fdr < 0.05),
      n_fdr_005_lfc1 = sum(sub$fdr < 0.05 & abs(sub$log_fc) > 1),
      stringsAsFactors = FALSE
    )
  }))
}))

write.csv(summary_counts, file.path(OUT_DIR, "de_summary_counts.csv"),
          row.names = FALSE)
print(summary_counts)
