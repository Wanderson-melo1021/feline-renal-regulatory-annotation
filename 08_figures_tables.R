library(ggplot2)
library(patchwork)

OUT_DIR <- "results/08_figures"
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")

REG_DIR <- "results/02_regulatory_annotation"
DE_DIR <- "results/01_differential_expression"
CHAR_DIR <- "results/03b_module_characterization"
MOD_ROOT <- "results/03_modules"
MOD_CORTEX <- file.path(MOD_ROOT, "cortex_late_vs_control")
MOD_MEDULLA <- file.path(MOD_ROOT, "medulla_early_vs_control")
SPEC_DIR <- "results/07_tissue_specificity"
AGE_DIR <- "results/09_cohort_and_age"
MOTIF_DIR <- "results/04_motif_enrichment"
PEAK_DIR <- "data/raw/chip_peaks"

MARK_LABELS <- c(h3k4me3 = "H3K4me3", h3k27ac = "H3K27ac", ctcf = "CTCF")
GROUP_LABELS <- c(Control = "Control", CKD12 = "CKD 1/2", CKD34 = "CKD 3/4")
DOMAIN_CAP_KB <- 500

FIG_WIDTH <- 180
FIG_DPI <- 600

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

theme_manuscript <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2),
          strip.background = element_rect(fill = "grey92", colour = NA),
          legend.key.size = unit(0.4, "cm"),
          plot.tag = element_text(face = "bold", size = base_size + 2))
}

theme_set(theme_manuscript())

read_if_exists <- function(path) {
  if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else NULL
}

clean_names <- function(df) {
  colnames(df) <- gsub("_", " ", colnames(df))
  colnames(df) <- gsub("\\.+", " ", colnames(df))
  colnames(df) <- trimws(colnames(df))
  df
}

round_numeric <- function(df, digits = 3) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(v) {
    ifelse(abs(v) < 1e-3 & v != 0, signif(v, 3), round(v, digits))
  })
  df
}

format_table <- function(df, name, caption) {
  if (is.null(df)) return(invisible(NULL))
  df <- clean_names(round_numeric(df))
  write.csv(df, file.path(TAB_DIR, sprintf("%s.csv", name)), row.names = FALSE)
  writeLines(c(caption, "", paste(colnames(df), collapse = "\t"),
               apply(df, 1, paste, collapse = "\t")),
             file.path(TAB_DIR, sprintf("%s.tsv", name)))
  message("saved table ", name)
}


save_figure <- function(plot, name, height_mm) {
  for (ext in c("pdf", "png")) {
    ggsave(file.path(FIG_DIR, sprintf("%s.%s", name, ext)), plot,
           width = FIG_WIDTH, height = height_mm, units = "mm",
           dpi = FIG_DPI, limitsize = FALSE)
  }
  message("saved ", name)
}

enrichment <- read_if_exists(file.path(REG_DIR, "promoter_enrichment_check.csv"))
classes <- read_if_exists(file.path(REG_DIR, "regulatory_class_counts.csv"))
annotation <- read_if_exists(file.path(REG_DIR, "gene_regulatory_annotation.csv"))

if (!is.null(enrichment) && !is.null(classes) && !is.null(annotation)) {
  enrichment$mark_label <- factor(unname(MARK_LABELS[enrichment$mark]),
                                  levels = unname(MARK_LABELS))
  enr_long <- data.frame(
    mark = rep(enrichment$mark_label, 2),
    type = rep(c("Observed", "Expected"), each = nrow(enrichment)),
    genes = c(enrichment$observed_genes_hit, enrichment$expected_genes_hit))
  enr_long$type <- factor(enr_long$type, levels = c("Expected", "Observed"))

  p1a <- ggplot(enr_long, aes(mark, genes, fill = type)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(data = enrichment,
              aes(x = mark_label, y = observed_genes_hit,
                  label = sprintf("%.1f-fold", fold_enrichment)),
              inherit.aes = FALSE, vjust = -0.5, size = 2.6) +
    scale_fill_manual(values = c(Expected = "grey70", Observed = "#2c6fbb")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Genes with a promoter peak", fill = NULL)

  classes$regulatory_class <- factor(
    classes$regulatory_class,
    levels = c("active_promoter", "primed_promoter", "distal_enhancer_only",
               "no_active_mark"))
  levels(classes$regulatory_class) <- c("Active promoter", "Primed promoter",
                                        "Distal enhancer only", "No active mark")

  p1b <- ggplot(classes, aes(regulatory_class, n_genes)) +
    geom_col(fill = "#2c6fbb", width = 0.65) +
    geom_text(aes(label = format(n_genes, big.mark = ",")), vjust = -0.4,
              size = 2.6) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Genes") +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  p1c <- ggplot(annotation, aes(domain_width / 1000)) +
    geom_histogram(bins = 60, fill = "#2c6fbb", colour = NA) +
    geom_vline(xintercept = DOMAIN_CAP_KB, linetype = "dashed",
               linewidth = 0.3, colour = "grey30") +
    annotate("text", x = DOMAIN_CAP_KB, y = Inf, label = "method cap",
             hjust = 1.1, vjust = 1.6, size = 2.4, colour = "grey30") +
    scale_x_log10() +
    labs(x = "CTCF-delimited domain width (kb, log scale)", y = "Genes")

  fig1 <- (p1a | p1b) / p1c + plot_annotation(tag_levels = "A")
  save_figure(fig1, "figure1_regulatory_annotation", 150)
}

de_summary <- read_if_exists(file.path(DE_DIR, "de_summary_counts.csv"))
de_cortex <- read_if_exists(file.path(DE_DIR, "de_cortex_late_vs_control.csv"))

if (!is.null(de_summary)) {
  de_summary$contrast <- factor(de_summary$contrast,
                                levels = c("early_vs_control",
                                           "late_vs_control", "late_vs_early"))
  levels(de_summary$contrast) <- c("CKD 1/2 vs control", "CKD 3/4 vs control",
                                   "CKD 3/4 vs CKD 1/2")
  de_summary$compartment <- factor(de_summary$compartment,
                                   levels = c("cortex", "medulla"),
                                   labels = c("Cortex", "Medulla"))

  p2a <- ggplot(de_summary, aes(contrast, n_fdr_010, fill = compartment)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = n_fdr_010),
              position = position_dodge(width = 0.7), vjust = -0.4, size = 2.5) +
    scale_fill_manual(values = c(Cortex = "#2c6fbb", Medulla = "#d1793a")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Differentially expressed genes (FDR < 0.10)",
         fill = NULL) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  if (!is.null(de_cortex)) {
    de_cortex$state <- with(de_cortex, ifelse(
      fdr < 0.10 & log_fc < 0, "Repressed",
      ifelse(fdr < 0.10 & log_fc > 0, "Induced", "Unchanged")))
    de_cortex$state <- factor(de_cortex$state,
                              levels = c("Unchanged", "Repressed", "Induced"))

    p2b <- ggplot(de_cortex, aes(log_fc, -log10(p_value), colour = state)) +
      geom_point(size = 0.35, alpha = 0.5) +
      scale_colour_manual(values = c(Unchanged = "grey75",
                                     Repressed = "#2c6fbb",
                                     Induced = "#c0392b")) +
      guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
      labs(x = expression(log[2]~fold~change), y = expression(-log[10]~P),
           colour = NULL)

    fig2 <- (p2a | p2b) + plot_annotation(tag_levels = "A")
  } else {
    fig2 <- p2a + plot_annotation(tag_levels = "A")
  }
  save_figure(fig2, "figure2_differential_expression", 80)
}

load_compartment <- function(dir, label) {
  st <- read_if_exists(file.path(dir, "module_stability.csv"))
  if (!is.null(st)) st$compartment <- label
  cm <- read_if_exists(file.path(dir, "composition_and_eigengenes.csv"))
  if (!is.null(cm)) cm$compartment <- label
  pc <- read_if_exists(file.path(dir, "direction_pc1_composition.csv"))
  if (!is.null(pc)) pc$compartment <- label
  list(stability = st, composition = cm, pc1 = pc)
}

cortex_mod <- load_compartment(MOD_CORTEX, "Cortex")
medulla_mod <- load_compartment(MOD_MEDULLA, "Medulla")

stability <- rbind(cortex_mod$stability, medulla_mod$stability)
pc1_all <- rbind(cortex_mod$pc1, medulla_mod$pc1)

if (!is.null(stability)) {
  stability$compartment <- factor(stability$compartment,
                                  levels = c("Cortex", "Medulla"))
  stability$direction <- factor(stability$direction,
                                levels = c("down", "up"),
                                labels = c("Repressed", "Induced"))
  stability <- stability[order(stability$compartment, stability$mean_jaccard), ]
  stability$rank <- ave(stability$mean_jaccard, stability$compartment,
                        FUN = seq_along)

  p3a <- ggplot(stability, aes(rank, mean_jaccard, colour = direction)) +
    geom_hline(yintercept = 0.6, linetype = "dashed", linewidth = 0.3) +
    geom_linerange(aes(ymin = pmax(mean_jaccard - sd_jaccard, 0),
                       ymax = pmin(mean_jaccard + sd_jaccard, 1)),
                   linewidth = 0.3, alpha = 0.6) +
    geom_point(size = 1.3) +
    facet_wrap(~ compartment, scales = "free_x") +
    scale_colour_manual(values = c(Repressed = "#2c6fbb", Induced = "#c0392b")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Candidate modules, ordered by stability",
         y = "Mean Jaccard index under subsampling", colour = NULL)
}

if (!is.null(cortex_mod$composition) && "down_PC1" %in% colnames(cortex_mod$composition)) {
  cx <- cortex_mod$composition
  cx$group <- factor(unname(GROUP_LABELS[cx$group]), levels = unname(GROUP_LABELS))
  rho_cx <- cor(cx$down_PC1, cx$epithelial, method = "spearman")

  p3b <- ggplot(cx, aes(epithelial, down_PC1, colour = group)) +
    geom_smooth(method = "lm", se = FALSE, colour = "grey40", linewidth = 0.4,
                inherit.aes = FALSE, aes(epithelial, down_PC1)) +
    geom_point(size = 1.6) +
    scale_colour_manual(values = c("Control" = "#4a7fb5", "CKD 1/2" = "#e0a458",
                                   "CKD 3/4" = "#b5453a")) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 2.8,
             label = sprintf("rho == %.2f", rho_cx), parse = TRUE) +
    labs(x = "Proximal tubule marker score",
         y = "First component, repressed genes", colour = NULL)
}

if (!is.null(pc1_all)) {
  pc1_all$compartment <- factor(pc1_all$compartment,
                                levels = c("Cortex", "Medulla"))
  pc1_all$direction <- factor(pc1_all$direction, levels = c("down", "up"),
                              labels = c("Repressed", "Induced"))
  pc1_all$marker_set <- factor(pc1_all$marker_set,
                               levels = c("epithelial", "immune", "fibroblast"),
                               labels = c("Epithelial", "Immune", "Fibroblast"))

  p3c <- ggplot(pc1_all, aes(marker_set, direction, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 2.6) +
    facet_wrap(~ compartment) +
    scale_fill_gradient2(low = "#c0392b", mid = "white", high = "#2c6fbb",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(x = NULL, y = NULL, fill = "Spearman rho")
}

if (exists("p3a") && exists("p3b") && exists("p3c")) {
  fig3 <- (p3a | p3b) / p3c + plot_layout(heights = c(1.4, 1)) +
    plot_annotation(tag_levels = "A")
  save_figure(fig3, "figure3_composition_dominance", 150)
}

motifs <- read_if_exists(file.path(MOTIF_DIR, "motifs_all.csv"))

if (!is.null(motifs)) {
  motifs$region_type <- factor(motifs$region_type,
                               levels = c("promoter", "enhancer"),
                               labels = c("Promoters", "Enhancers"))
  motifs$compartment <- ifelse(grepl("^cortex", motifs$set), "Cortex", "Medulla")
  motifs$direction_label <- factor(motifs$direction, levels = c("down", "up"),
                                   labels = c("Repressed", "Induced"))
  motifs$significant <- !is.na(motifs$padj) & motifs$padj < 0.05
  motifs$coherent <- with(motifs,
    (direction == "up" & tf_log_fc > 0) | (direction == "down" & tf_log_fc < 0))

  n_missing <- sum(is.na(motifs$log2_enrichment) | is.na(motifs$neg_log10_padj))
  if (n_missing) {
    message(sprintf("motif panel: %d of %d tests dropped for missing values",
                    n_missing, nrow(motifs)))
  }
  motifs_plot <- motifs[!is.na(motifs$log2_enrichment) &
                          !is.na(motifs$neg_log10_padj), ]

  p4a <- ggplot(motifs_plot, aes(log2_enrichment, neg_log10_padj,
                                 colour = significant)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               linewidth = 0.3) +
    geom_point(size = 0.5, alpha = 0.5) +
    facet_wrap(~ region_type) +
    scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "#c0392b"),
                        guide = "none") +
    labs(x = expression(log[2]~motif~enrichment),
         y = expression(-log[10]~adjusted~P))

  hits <- motifs[motifs$significant & !is.na(motifs$tf_fdr) &
                   motifs$tf_fdr < 0.10 & !is.na(motifs$coherent) &
                   motifs$coherent, ]

  if (nrow(hits)) {
    hits <- hits[!duplicated(paste(hits$compartment, hits$motif_name)), ]
    hits$label <- sprintf("%s (%s)", hits$motif_name, hits$compartment)
    hits <- hits[order(hits$log2_enrichment), ]
    hits$label <- factor(hits$label, levels = hits$label)

    long <- rbind(
      data.frame(label = hits$label, compartment = hits$compartment,
                 panel = "Motif enrichment", value = hits$log2_enrichment,
                 stringsAsFactors = FALSE),
      data.frame(label = hits$label, compartment = hits$compartment,
                 panel = "Factor expression", value = hits$tf_log_fc,
                 stringsAsFactors = FALSE))
    long$panel <- factor(long$panel,
                         levels = c("Motif enrichment", "Factor expression"))

    p4b <- ggplot(long, aes(value, label, colour = compartment)) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
      geom_segment(aes(x = 0, xend = value, yend = label), linewidth = 0.4) +
      geom_point(size = 2.2) +
      facet_wrap(~ panel, scales = "free_x") +
      scale_colour_manual(values = c(Cortex = "#2c6fbb", Medulla = "#d1793a")) +
      labs(x = expression(log[2]~scale), y = NULL, colour = NULL)

    fig4 <- (p4a / p4b) + plot_layout(heights = c(1, 1.1)) +
      plot_annotation(tag_levels = "A")
    save_figure(fig4, "figure4_motif_enrichment", 135)

    write.csv(hits[, c("compartment", "direction", "region_type", "motif_name",
                       "log2_enrichment", "padj", "tf_matched", "tf_log_fc",
                       "tf_fdr")],
              file.path(TAB_DIR, "table6_converging_motifs.csv"),
              row.names = FALSE)
  }

  motif_summary <- as.data.frame(table(motifs$compartment, motifs$region_type,
                                       motifs$significant))
  colnames(motif_summary) <- c("compartment", "region_type", "significant", "n")
  format_table(motif_summary, "tableS5_motif_summary",
               "Supplementary Table S5. Number of motif tests reaching a corrected P below 0.05, by compartment and region type.")
}

element_spec <- read_if_exists(file.path(SPEC_DIR, "element_specificity.csv"))
gene_spec <- read_if_exists(file.path(SPEC_DIR,
                                      "gene_enhancer_specificity.csv"))
models <- read_if_exists(file.path(SPEC_DIR, "logistic_models.csv"))

if (!is.null(element_spec)) {
  element_spec$class <- factor(element_spec$class,
                               levels = c("kidney_specific", "restricted",
                                          "shared", "ubiquitous"))
  levels(element_spec$class) <- c("Kidney-specific", "Restricted (2-3)",
                                  "Shared (4-5)", "Ubiquitous (6-7)")
  element_spec$mark <- factor(toupper(element_spec$mark),
                              levels = c("H3K4ME3", "H3K27AC"),
                              labels = c("H3K4me3", "H3K27ac"))

  p5a <- ggplot(element_spec, aes(mark, pct, fill = class)) +
    geom_col(width = 0.6) +
    scale_fill_brewer(palette = "Blues", direction = -1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Percentage of kidney elements", fill = "Tissue breadth")
}

if (!is.null(gene_spec)) {
  gene_spec$state <- factor(gene_spec$state,
                            levels = c("unchanged", "repressed", "induced"),
                            labels = c("Unchanged", "Repressed", "Induced"))

  p5b <- ggplot(gene_spec, aes(state, mean_tissues, fill = state)) +
    geom_violin(scale = "width", linewidth = 0.25, colour = "grey40") +
    geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.3,
                 fill = "white") +
    scale_fill_manual(values = c(Unchanged = "grey75", Repressed = "#2c6fbb",
                                 Induced = "#c0392b"), guide = "none") +
    labs(x = NULL, y = "Mean tissue breadth of gene enhancers")
}

if (!is.null(models)) {
  forest <- models[models$term == "frac_kidney_specific", ]
  forest$ci_low <- exp(forest$Estimate - 1.96 * forest$Std..Error)
  forest$ci_high <- exp(forest$Estimate + 1.96 * forest$Std..Error)
  forest$outcome <- factor(forest$outcome, levels = c("induced", "repressed"),
                           labels = c("Induced", "Repressed"))
  forest$or_per_0.1 <- round(exp(forest$Estimate * 0.1), 3)

  p5c <- ggplot(forest, aes(odds_ratio, outcome)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                  width = 0.12, linewidth = 0.4) +
    geom_point(size = 2.4, colour = "#2c6fbb") +
    scale_x_log10() +
    labs(x = "Odds ratio, no kidney-specific enhancers vs all kidney-specific",
         y = NULL)

  fig4_spec <- (p5a | p5b) / p5c + plot_layout(heights = c(1.6, 1)) +
    plot_annotation(tag_levels = "A")
  save_figure(fig4_spec, "figure5_tissue_specificity", 150)
}


cohort <- read_if_exists(file.path(AGE_DIR, "table1_cohort.csv"))
if (!is.null(cohort)) {
  cohort_table <- data.frame(
    Compartment = tools::toTitleCase(cohort$compartment),
    Group = cohort$group,
    n = cohort$n,
    `Age, mean (SD)` = sprintf("%.1f (%.1f)", cohort$age_mean, cohort$age_sd),
    `Age, median [range]` = sprintf("%.1f [%.1f-%.1f]", cohort$age_median,
                                    cohort$age_min, cohort$age_max),
    `Sex, F/M` = sprintf("%d/%d", cohort$sex_female, cohort$sex_male),
    `Breed, DSH/DLH` = sprintf("%d/%d", cohort$breed_dsh, cohort$breed_dlh),
    check.names = FALSE, stringsAsFactors = FALSE)
  format_table(cohort_table, "table1_cohort",
               sprintf("Table 1. Cohort characteristics. Homogeneity across disease groups: age P = %s (cortex) and %s (medulla); sex P = %s and %s; breed P = %s and %s.",
                       cohort$age_p[1], cohort$age_p[4], cohort$sex_p[1],
                       cohort$sex_p[4], cohort$breed_p[1], cohort$breed_p[4]))
}

format_table(de_summary, "table2_de_counts",
             "Table 2. Differentially expressed genes by compartment and contrast.")
format_table(enrichment, "table3_promoter_enrichment",
             "Table 3. Permutation validation of promoter mark enrichment.")
format_table(element_spec, "table4_element_specificity",
             "Table 4. Tissue breadth of renal cortical regulatory elements.")
if (!is.null(models)) {
  models$ci_low <- round(exp(models$Estimate - 1.96 * models$Std..Error), 3)
  models$ci_high <- round(exp(models$Estimate + 1.96 * models$Std..Error), 3)
  models <- models[, c("outcome", "term", "Estimate", "odds_ratio",
                       "ci_low", "ci_high", "Pr...z..")]
  colnames(models)[colnames(models) == "Pr...z.."] <- "p_value"
}

format_table(models, "table5_specificity_models",
             "Table 5. Association between enhancer tissue specificity and expression state.")
if (!is.null(stability)) {
  stability_out <- stability[, c("compartment", "direction", "module",
                                 "n_genes", "mean_jaccard", "sd_jaccard",
                                 "retained")]
  stability_out <- stability_out[order(stability_out$compartment,
                                       -stability_out$mean_jaccard), ]
} else {
  stability_out <- NULL
}

format_table(stability_out, "tableS1_module_stability",
             "Supplementary Table S1. Stability of all candidate coexpression modules under subsampling, both compartments.")
format_table(pc1_all, "tableS4_pc1_composition",
             "Supplementary Table S4. Association between the first principal component of each direction and cell-type marker scores.")
format_table(read_if_exists(file.path(AGE_DIR, "age_sensitivity_concordance.csv")),
             "tableS3_age_sensitivity",
             "Supplementary Table S3. Concordance between the primary differential expression model, which includes age, and a sensitivity model without it.")

tissues <- c("heart", "kidney", "liver", "lung", "smallintestine",
             "testis", "thyroid")
marks <- c("h3k4me3", "h3k27ac", "ctcf")

peak_counts <- do.call(rbind, lapply(tissues, function(t) {
  do.call(rbind, lapply(marks, function(m) {
    hit <- grep(sprintf("%s-%s.*narrowPeak", t, m), list.files(PEAK_DIR),
                value = TRUE, ignore.case = TRUE)
    if (length(hit) != 1L) return(NULL)
    n <- length(readLines(gzfile(file.path(PEAK_DIR, hit[1]))))
    data.frame(tissue = t, mark = unname(MARK_LABELS[m]), n_peaks = n,
               stringsAsFactors = FALSE)
  }))
}))

if (!is.null(peak_counts)) {
  peak_counts$tissue <- factor(peak_counts$tissue, levels = tissues)
  peak_counts$mark <- factor(peak_counts$mark, levels = unname(MARK_LABELS))
  peak_counts$reference <- peak_counts$tissue == "kidney"

  figS1 <- ggplot(peak_counts, aes(tissue, n_peaks, fill = reference)) +
    geom_col(width = 0.65) +
    facet_wrap(~ mark, scales = "free_y") +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#2c6fbb"),
                      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Replicated peaks") +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))

  save_figure(figS1, "figureS1_peak_counts_by_tissue", 70)
  format_table(peak_counts[, c("tissue", "mark", "n_peaks")],
               "tableS2_peak_counts_by_tissue",
               "Supplementary Table S2. Replicated peak counts per tissue and mark.")
}

age_genes <- do.call(rbind, lapply(c("cortex", "medulla"), function(cp) {
  d <- read_if_exists(file.path(AGE_DIR,
                                sprintf("age_sensitivity_genes_%s.csv", cp)))
  if (is.null(d)) return(NULL)
  d$compartment <- cp
  d
}))

if (!is.null(age_genes)) {
  age_genes$compartment <- factor(age_genes$compartment,
                                  levels = c("cortex", "medulla"),
                                  labels = c("Cortex", "Medulla"))
  labs_df <- do.call(rbind, lapply(split(age_genes, age_genes$compartment),
                                   function(d) data.frame(
    compartment = d$compartment[1],
    label = sprintf("rho == %.3f", cor(d$log_fc_with_age,
                                       d$log_fc_without_age,
                                       method = "spearman")),
    stringsAsFactors = FALSE)))

  figS2 <- ggplot(age_genes, aes(log_fc_without_age, log_fc_with_age)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.3, colour = "grey50") +
    geom_point(size = 0.3, alpha = 0.3, colour = "#2c6fbb") +
    geom_text(data = labs_df, aes(x = -Inf, y = Inf, label = label),
              parse = TRUE, hjust = -0.15, vjust = 1.6, size = 2.8,
              inherit.aes = FALSE) +
    facet_wrap(~ compartment) +
    labs(x = "log2 fold change, sensitivity model without age",
         y = "log2 fold change, primary model with age")

  save_figure(figS2, "figureS2_age_sensitivity", 80)
}

message("figures in ", FIG_DIR)
message("tables in ", TAB_DIR)
