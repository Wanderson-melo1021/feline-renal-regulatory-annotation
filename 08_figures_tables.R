library(ggplot2)
library(patchwork)

OUT_DIR <- "results/08_figures"
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")

REG_DIR <- "results/02_regulatory_annotation"
DE_DIR <- "results/01_differential_expression"
CHAR_DIR <- "results/03b_module_characterization"
MOD_DIR <- "results/03_coexpression_modules"
SPEC_DIR <- "results/07_tissue_specificity"
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

stability <- read_if_exists(file.path(MOD_DIR, "module_stability.csv"))
composition <- read_if_exists(file.path(CHAR_DIR,
                                        "composition_and_eigengenes.csv"))

if (!is.null(stability)) {
  stability$label <- sprintf("%s M%s", stability$direction, stability$module)
  stability <- stability[order(stability$mean_jaccard), ]
  stability$label <- factor(stability$label, levels = stability$label)

  p3a <- ggplot(stability, aes(label, mean_jaccard,
                               fill = mean_jaccard >= 0.6)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = pmax(mean_jaccard - sd_jaccard, 0),
                      ymax = pmin(mean_jaccard + sd_jaccard, 1)),
                  width = 0.25, linewidth = 0.3) +
    geom_hline(yintercept = 0.6, linetype = "dashed", linewidth = 0.3) +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#2c6fbb"),
                      guide = "none") +
    coord_flip() +
    labs(x = NULL, y = "Mean Jaccard index under subsampling")
}

if (!is.null(composition) && "down_M1" %in% colnames(composition)) {
  rho <- cor(composition$down_M1, composition$tubular, method = "spearman")
  composition$group <- factor(unname(GROUP_LABELS[composition$group]),
                              levels = unname(GROUP_LABELS))
  p3b <- ggplot(composition, aes(tubular, down_M1, colour = group)) +
    scale_colour_manual(values = c("Control" = "#4a7fb5",
                                   "CKD 1/2" = "#e0a458",
                                   "CKD 3/4" = "#b5453a")) +
    geom_point(size = 1.6) +
    geom_smooth(method = "lm", se = FALSE, colour = "grey40",
                linewidth = 0.4, inherit.aes = FALSE,
                aes(tubular, down_M1)) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 2.8,
             label = sprintf("rho == %.2f", rho), parse = TRUE) +
    labs(x = "Proximal tubule marker score",
         y = "Module eigengene (down_M1)", colour = NULL)

  fig3 <- (p3a | p3b) + plot_annotation(tag_levels = "A")
  save_figure(fig3, "figure3_composition_dominance", 85)
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

  p4a <- ggplot(element_spec, aes(mark, pct, fill = class)) +
    geom_col(width = 0.6) +
    scale_fill_brewer(palette = "Blues", direction = -1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Percentage of kidney elements", fill = "Tissue breadth")
}

if (!is.null(gene_spec)) {
  gene_spec$state <- factor(gene_spec$state,
                            levels = c("unchanged", "repressed", "induced"),
                            labels = c("Unchanged", "Repressed", "Induced"))

  p4b <- ggplot(gene_spec, aes(state, mean_tissues, fill = state)) +
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

  p4c <- ggplot(forest, aes(odds_ratio, outcome)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.12,
                   linewidth = 0.4) +
    geom_point(size = 2.4, colour = "#2c6fbb") +
    scale_x_log10() +
    labs(x = "Odds ratio, no kidney-specific enhancers vs all kidney-specific",
         y = NULL)

  fig4 <- (p4a | p4b) / p4c + plot_layout(heights = c(1.6, 1)) +
    plot_annotation(tag_levels = "A")
  save_figure(fig4, "figure4_tissue_specificity", 150)
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

format_table(de_summary, "table1_de_counts",
             "Table 1. Differentially expressed genes by compartment and contrast.")
format_table(enrichment, "table2_promoter_enrichment",
             "Table 2. Permutation validation of promoter mark enrichment.")
format_table(element_spec, "table3_element_specificity",
             "Table 3. Tissue breadth of renal cortical regulatory elements.")
if (!is.null(models)) {
  models$ci_low <- round(exp(models$Estimate - 1.96 * models$Std..Error), 3)
  models$ci_high <- round(exp(models$Estimate + 1.96 * models$Std..Error), 3)
  models <- models[, c("outcome", "term", "Estimate", "odds_ratio",
                       "ci_low", "ci_high", "Pr...z..")]
  colnames(models)[colnames(models) == "Pr...z.."] <- "p_value"
}

format_table(models, "table4_specificity_models",
             "Table 4. Association between enhancer tissue specificity and expression state.")
if (!is.null(stability)) {
  stability <- stability[, c("direction", "module", "n_genes", "mean_jaccard",
                             "sd_jaccard", "retained")]
}

format_table(stability, "tableS1_module_stability",
             "Supplementary Table S1. Coexpression module stability under subsampling.")

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

message("figures in ", FIG_DIR)
message("tables in ", TAB_DIR)
