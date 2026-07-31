library(readxl)
library(limma)

XLSX_PATH <- "data/raw/expression/42003_2025_9164_MOESM3_ESM.xlsx"
OUT_DIR <- "results/01_differential_expression"
OBJ_PATH <- "data/processed/expression_objects.rds"

MIN_AGE <- 0
ADJUST_FOR_AGE <- TRUE

SHEETS <- c(cortex = "S3 RNA cortex", medulla = "S4 RNA medulla")
META_ROWS <- c(id = 2, group = 3, sex = 4, age = 5, breed = 6)
FIRST_GENE_ROW <- 7

GROUP_LEVELS <- c("Control", "CKD12", "CKD34")
GROUP_RECODE <- c("Control" = "Control", "CKD1/2" = "CKD12", "CKD3/4" = "CKD34")

CONTRASTS <- c(
  early_vs_control = "CKD12 - Control",
  late_vs_control  = "CKD34 - Control",
  late_vs_early    = "CKD34 - CKD12"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OBJ_PATH), recursive = TRUE, showWarnings = FALSE)

to_numeric <- function(x) as.numeric(gsub(",", ".", as.character(x)))

read_compartment <- function(path, sheet) {
  raw <- as.data.frame(read_excel(path, sheet = sheet, col_names = FALSE,
                                  .name_repair = "minimal", progress = FALSE),
                       stringsAsFactors = FALSE)

  sample_ids <- as.character(unlist(raw[META_ROWS["id"], -1]))
  keep <- !is.na(sample_ids) & nzchar(sample_ids)
  sample_ids <- sample_ids[keep]

  samples <- data.frame(
    sample_id = sample_ids,
    group = factor(unname(GROUP_RECODE[
      as.character(unlist(raw[META_ROWS["group"], -1]))[keep]]),
      levels = GROUP_LEVELS),
    sex = factor(as.character(unlist(raw[META_ROWS["sex"], -1]))[keep]),
    age = to_numeric(unlist(raw[META_ROWS["age"], -1]))[keep],
    breed = factor(as.character(unlist(raw[META_ROWS["breed"], -1]))[keep]),
    stringsAsFactors = FALSE
  )

  gene_block <- raw[FIRST_GENE_ROW:nrow(raw), , drop = FALSE]
  genes <- as.character(gene_block[[1]])
  valid <- !is.na(genes) & nzchar(genes)
  gene_block <- gene_block[valid, , drop = FALSE]
  genes <- genes[valid]

  expr <- as.matrix(gene_block[, -1, drop = FALSE][, keep, drop = FALSE])
  expr <- matrix(as.numeric(expr), nrow = nrow(expr),
                 dimnames = list(genes, sample_ids))

  stopifnot(!any(is.na(expr)), !any(is.na(samples$group)),
            !any(is.na(samples$age)))
  list(expr = expr, samples = samples)
}

apply_age_filter <- function(dat, compartment) {
  if (MIN_AGE <= 0) return(dat)
  dropped <- dat$samples[dat$samples$age < MIN_AGE, ]
  if (nrow(dropped)) {
    message(sprintf("%s: excluding %d animal(s) below %g years: %s",
                    compartment, nrow(dropped), MIN_AGE,
                    paste(sprintf("%s (%s, %.1f y)", dropped$sample_id,
                                  dropped$group, dropped$age),
                          collapse = "; ")))
  }
  keep <- dat$samples$age >= MIN_AGE
  list(expr = dat$expr[, keep, drop = FALSE],
       samples = droplevels(dat$samples[keep, , drop = FALSE]))
}

fit_compartment <- function(expr, samples) {
  rhs <- if (ADJUST_FOR_AGE) ~ 0 + group + sex + age else ~ 0 + group + sex
  design <- model.matrix(rhs, data = samples)
  colnames(design) <- make.names(sub("^group", "", colnames(design)))

  contrast_terms <- function(cs) {
    parts <- unlist(strsplit(cs, "[-+*/() ]+"))
    parts[nzchar(parts) & is.na(suppressWarnings(as.numeric(parts)))]
  }

  usable <- names(CONTRASTS)[vapply(CONTRASTS, function(cs) {
    all(contrast_terms(cs) %in% colnames(design))
  }, logical(1))]

  cm <- makeContrasts(contrasts = unname(CONTRASTS[usable]), levels = design)
  colnames(cm) <- usable

  fit <- eBayes(contrasts.fit(lmFit(expr, design), cm),
                trend = TRUE, robust = TRUE)
  attr(fit, "contrasts_fitted") <- usable
  fit
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
    row.names = NULL, stringsAsFactors = FALSE
  )
}

results <- list()
objects <- list()

for (compartment in names(SHEETS)) {
  dat <- apply_age_filter(read_compartment(XLSX_PATH, SHEETS[[compartment]]),
                          compartment)
  fit <- fit_compartment(dat$expr, dat$samples)
  fitted_contrasts <- attr(fit, "contrasts_fitted")

  tabs <- lapply(fitted_contrasts, extract_table, fit = fit,
                 compartment = compartment)
  names(tabs) <- fitted_contrasts

  for (nm in names(tabs)) {
    write.csv(tabs[[nm]],
              file.path(OUT_DIR, sprintf("de_%s_%s.csv", compartment, nm)),
              row.names = FALSE)
  }

  results[[compartment]] <- do.call(rbind, tabs)
  objects[[compartment]] <- list(expr = dat$expr, samples = dat$samples,
                                 fit = fit)

  message(sprintf("%s: %d genes, %d samples (%s); age %s",
                  compartment, nrow(dat$expr), ncol(dat$expr),
                  paste(levels(dat$samples$group),
                        table(dat$samples$group), sep = "=", collapse = ", "),
                  if (ADJUST_FOR_AGE) "in model" else "not in model"))
}

saveRDS(objects, OBJ_PATH)
write.csv(do.call(rbind, results),
          file.path(OUT_DIR, "de_all_contrasts.csv"), row.names = FALSE)

summary_counts <- do.call(rbind, lapply(names(results), function(cp) {
  df <- results[[cp]]
  do.call(rbind, lapply(unique(df$contrast), function(ct) {
    sub <- df[df$contrast == ct, ]
    data.frame(compartment = cp, contrast = ct,
               n_fdr_010 = sum(sub$fdr < 0.10),
               n_fdr_005 = sum(sub$fdr < 0.05),
               n_fdr_005_lfc1 = sum(sub$fdr < 0.05 & abs(sub$log_fc) > 1),
               stringsAsFactors = FALSE)
  }))
}))

write.csv(summary_counts, file.path(OUT_DIR, "de_summary_counts.csv"),
          row.names = FALSE)
print(summary_counts)
