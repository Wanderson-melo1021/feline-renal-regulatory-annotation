library(dynamicTreeCut)

ARGS <- commandArgs(trailingOnly = TRUE)

COMPARTMENT <- if (length(ARGS) >= 1) ARGS[1] else "cortex"
CONTRAST <- if (length(ARGS) >= 2) ARGS[2] else "late_vs_control"

EXPR_OBJ <- "data/processed/expression_objects.rds"
REG_OBJ <- "data/processed/regulatory_objects.rds"
DE_DIR <- "results/01_differential_expression"
OUT_ROOT <- "results/03_modules"

FDR_CUTOFF <- if (length(ARGS) >= 3) as.numeric(ARGS[3]) else 0.10
LFC_CUTOFF <- if (length(ARGS) >= 4) as.numeric(ARGS[4]) else 0
MIN_MODULE_SIZE <- 30
DEEP_SPLIT <- 2
N_BOOTSTRAP <- 200
SUBSAMPLE_FRACTION <- 0.8
MIN_JACCARD <- 0.6
CONFOUND_CEILING <- 0.6
SEED <- 1

MARKER_SETS <- list(
  cortex = list(
    epithelial = c("LRP2", "CUBN", "SLC34A1", "SLC5A2", "SLC5A12", "MIOX",
                   "GGT1", "SLC22A6", "SLC22A8", "SLC3A1", "SLC7A9", "HNF4A"),
    immune = c("PTPRC", "CD68", "CSF1R", "ITGAM", "AIF1", "LYZ", "C1QA",
               "C1QB", "TYROBP", "FCER1G"),
    fibroblast = c("COL1A1", "COL1A2", "COL3A1", "ACTA2", "PDGFRB", "FN1",
                   "LUM", "DCN")
  ),
  medulla = list(
    epithelial = c("UMOD", "SLC12A1", "CLDN16", "CLDN10", "KCNJ1", "ESRRB",
                   "AQP2", "AQP3", "SCNN1G", "SCNN1B", "GATA3", "ELF5",
                   "ATP6V1B1", "ATP6V0D2", "SLC4A1", "FOXI1", "AQP1"),
    immune = c("PTPRC", "CD68", "CSF1R", "ITGAM", "AIF1", "LYZ", "C1QA",
               "C1QB", "TYROBP", "FCER1G"),
    fibroblast = c("COL1A1", "COL1A2", "COL3A1", "ACTA2", "PDGFRB", "FN1",
                   "LUM", "DCN")
  )
)

OUT_DIR <- file.path(OUT_ROOT, sprintf("%s_%s", COMPARTMENT, CONTRAST))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message(sprintf("compartment %s, contrast %s, FDR < %g, |log2FC| > %g",
                COMPARTMENT, CONTRAST, FDR_CUTOFF, LFC_CUTOFF))

build_distance <- function(expr) {
  as.dist(1 - cor(t(expr), method = "spearman"))
}

cluster_genes <- function(expr) {
  d <- build_distance(expr)
  tree <- hclust(d, method = "average")
  labels <- cutreeDynamic(dendro = tree, distM = as.matrix(d),
                          deepSplit = DEEP_SPLIT,
                          minClusterSize = MIN_MODULE_SIZE,
                          pamRespectsDendro = FALSE)
  setNames(as.integer(labels), rownames(expr))
}

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

module_stability <- function(expr, reference) {
  ref <- split(names(reference), reference)
  ref <- ref[names(ref) != "0"]
  if (!length(ref)) return(NULL)
  n_keep <- max(round(ncol(expr) * SUBSAMPLE_FRACTION), 5)

  scores <- matrix(NA_real_, nrow = length(ref), ncol = N_BOOTSTRAP,
                   dimnames = list(names(ref), NULL))

  for (b in seq_len(N_BOOTSTRAP)) {
    sub <- try(cluster_genes(expr[, sample(ncol(expr), n_keep), drop = FALSE]),
               silent = TRUE)
    if (inherits(sub, "try-error")) next
    sub_mod <- split(names(sub), sub)
    sub_mod <- sub_mod[names(sub_mod) != "0"]
    if (!length(sub_mod)) next
    for (m in names(ref)) {
      scores[m, b] <- max(vapply(sub_mod, jaccard, numeric(1), a = ref[[m]]))
    }
  }

  data.frame(module = names(ref),
             n_genes = vapply(ref, length, integer(1)),
             mean_jaccard = round(rowMeans(scores, na.rm = TRUE), 3),
             sd_jaccard = round(apply(scores, 1, sd, na.rm = TRUE), 3),
             row.names = NULL, stringsAsFactors = FALSE)
}

module_eigengene <- function(expr, genes) {
  x <- t(scale(t(expr[genes, , drop = FALSE])))
  pc <- prcomp(t(x), center = FALSE, scale. = FALSE)
  eig <- pc$x[, 1]
  if (sign(mean(pc$rotation[, 1])) < 0) eig <- -eig
  list(eigengene = eig,
       var_explained = round(pc$sdev[1]^2 / sum(pc$sdev^2), 3))
}

score_marker_set <- function(expr, markers) {
  present <- intersect(markers, rownames(expr))
  attr_present <- length(present)
  out <- if (attr_present) {
    colMeans(t(scale(t(expr[present, , drop = FALSE]))))
  } else {
    rep(NA_real_, ncol(expr))
  }
  attr(out, "n_markers") <- attr_present
  out
}

expr_objects <- readRDS(EXPR_OBJ)
reg <- readRDS(REG_OBJ)

expr <- expr_objects[[COMPARTMENT]]$expr
samples <- expr_objects[[COMPARTMENT]]$samples

de <- read.csv(file.path(DE_DIR, sprintf("de_%s_%s.csv", COMPARTMENT, CONTRAST)),
               stringsAsFactors = FALSE)

coding <- reg$annotation$gene_symbol[
  !is.na(reg$annotation$gene_biotype) &
    reg$annotation$gene_biotype == "protein_coding"]

selected <- de[de$fdr < FDR_CUTOFF & abs(de$log_fc) > LFC_CUTOFF &
                 de$gene_symbol %in% coding &
                 de$gene_symbol %in% rownames(expr), ]

message(sprintf("%s / %s: %d genes selected (%d up, %d down) from %d samples",
                COMPARTMENT, CONTRAST, nrow(selected),
                sum(selected$log_fc > 0), sum(selected$log_fc < 0),
                ncol(expr)))

markers <- MARKER_SETS[[COMPARTMENT]]
composition <- data.frame(sample_id = samples$sample_id,
                          group = samples$group,
                          stringsAsFactors = FALSE)
marker_coverage <- list()
for (cell in names(markers)) {
  s <- score_marker_set(expr, markers[[cell]])
  composition[[cell]] <- as.numeric(s)
  marker_coverage[[cell]] <- data.frame(
    marker_set = cell,
    n_requested = length(markers[[cell]]),
    n_found = attr(s, "n_markers"),
    stringsAsFactors = FALSE)
}
marker_coverage <- do.call(rbind, marker_coverage)
write.csv(marker_coverage, file.path(OUT_DIR, "marker_coverage.csv"),
          row.names = FALSE)
print(marker_coverage)

identifiability <- do.call(rbind, lapply(names(markers), function(cell) {
  if (all(is.na(composition[[cell]]))) return(NULL)
  fit <- lm(composition[[cell]] ~ samples$group)
  data.frame(marker_set = cell,
             r_squared_on_group = round(summary(fit)$r.squared, 3),
             anova_p = signif(anova(fit)[["Pr(>F)"]][1], 3),
             stringsAsFactors = FALSE)
}))
write.csv(identifiability, file.path(OUT_DIR, "identifiability_report.csv"),
          row.names = FALSE)
print(identifiability)

set.seed(SEED)
results <- list()

for (direction in c("up", "down")) {
  genes <- selected$gene_symbol[
    if (direction == "up") selected$log_fc > 0 else selected$log_fc < 0]

  if (length(genes) < MIN_MODULE_SIZE) {
    message(sprintf("%s: only %d genes, skipped", direction, length(genes)))
    next
  }

  sub_expr <- expr[genes, , drop = FALSE]
  labels <- cluster_genes(sub_expr)
  stability <- module_stability(sub_expr, labels)
  if (is.null(stability)) next

  stability$direction <- direction
  stability$retained <- stability$mean_jaccard >= MIN_JACCARD

  assignment <- data.frame(gene_symbol = names(labels), direction = direction,
                           module = as.character(labels),
                           stringsAsFactors = FALSE)
  assignment$retained <- assignment$module %in%
    stability$module[stability$retained]

  kept <- stability$module[stability$retained]
  if (length(kept)) {
    eig <- lapply(kept, function(m) {
      module_eigengene(sub_expr, assignment$gene_symbol[assignment$module == m])
    })
    names(eig) <- kept
    eig_matrix <- do.call(cbind, lapply(eig, `[[`, "eigengene"))
    colnames(eig_matrix) <- sprintf("%s_M%s", direction, kept)
    var_explained <- vapply(eig, `[[`, numeric(1), "var_explained")
  } else {
    eig_matrix <- NULL
    var_explained <- NULL
  }

  results[[direction]] <- list(assignment = assignment, stability = stability,
                               eigengenes = eig_matrix,
                               var_explained = var_explained)
}

direction_pc1 <- do.call(rbind, lapply(c("up", "down"), function(direction) {
  genes <- selected$gene_symbol[
    if (direction == "up") selected$log_fc > 0 else selected$log_fc < 0]
  genes <- intersect(genes, rownames(expr))
  if (length(genes) < MIN_MODULE_SIZE) return(NULL)

  pc <- module_eigengene(expr, genes)
  composition[[paste0(direction, "_PC1")]] <<- pc$eigengene

  do.call(rbind, lapply(names(markers), function(cell) {
    if (all(is.na(composition[[cell]]))) return(NULL)
    ct <- suppressWarnings(cor.test(pc$eigengene, composition[[cell]],
                                    method = "spearman"))
    data.frame(direction = direction, n_genes = length(genes),
               var_explained = pc$var_explained,
               marker_set = cell,
               rho = round(unname(ct$estimate), 3),
               p_value = signif(ct$p.value, 3),
               stringsAsFactors = FALSE)
  }))
}))

if (!is.null(direction_pc1)) {
  direction_pc1$fdr <- signif(p.adjust(direction_pc1$p_value, "BH"), 3)
  write.csv(direction_pc1, file.path(OUT_DIR, "direction_pc1_composition.csv"),
            row.names = FALSE)
  cat("\nfirst principal component of each direction against composition\n")
  print(direction_pc1)

  group_pc1 <- do.call(rbind, lapply(c("up", "down"), function(direction) {
    key <- paste0(direction, "_PC1")
    if (is.null(composition[[key]])) return(NULL)
    data.frame(direction = direction,
               kruskal_p = signif(kruskal.test(composition[[key]] ~
                                                 samples$group)$p.value, 3),
               stringsAsFactors = FALSE)
  }))
  write.csv(group_pc1, file.path(OUT_DIR, "direction_pc1_group.csv"),
            row.names = FALSE)
  print(group_pc1)

  write.csv(composition, file.path(OUT_DIR, "composition_and_eigengenes.csv"),
            row.names = FALSE)
}

if (!length(results)) stop("no direction yielded modules")

stability_all <- do.call(rbind, lapply(results, `[[`, "stability"))
assignment_all <- do.call(rbind, lapply(results, `[[`, "assignment"))
eigengenes <- do.call(cbind, Filter(Negate(is.null),
                                    lapply(results, `[[`, "eigengenes")))

write.csv(stability_all, file.path(OUT_DIR, "module_stability.csv"),
          row.names = FALSE)
write.csv(assignment_all, file.path(OUT_DIR, "module_assignment.csv"),
          row.names = FALSE)
print(stability_all)

abundance <- merge(assignment_all, de[, c("gene_symbol", "ave_expr")],
                   by = "gene_symbol")
module_abundance <- aggregate(ave_expr ~ direction + module, data = abundance,
                              FUN = median)
stability_abundance <- merge(stability_all, module_abundance,
                             by = c("direction", "module"))
write.csv(stability_abundance, file.path(OUT_DIR, "stability_vs_abundance.csv"),
          row.names = FALSE)

if (is.null(eigengenes)) {
  message("no module passed the stability threshold; nothing to test for confounding")
} else {
  composition <- cbind(composition, as.data.frame(eigengenes))
  write.csv(composition, file.path(OUT_DIR, "composition_and_eigengenes.csv"),
            row.names = FALSE)

  confound <- do.call(rbind, lapply(colnames(eigengenes), function(m) {
    do.call(rbind, lapply(names(markers), function(cell) {
      if (all(is.na(composition[[cell]]))) return(NULL)
      ct <- suppressWarnings(cor.test(eigengenes[, m], composition[[cell]],
                                     method = "spearman"))
      data.frame(module = m, marker_set = cell,
                 rho = round(unname(ct$estimate), 3),
                 p_value = signif(ct$p.value, 3),
                 stringsAsFactors = FALSE)
    }))
  }))
  confound$fdr <- signif(p.adjust(confound$p_value, "BH"), 3)
  write.csv(confound, file.path(OUT_DIR, "eigengene_composition_confounding.csv"),
            row.names = FALSE)
  print(confound)

  group_test <- do.call(rbind, lapply(colnames(eigengenes), function(m) {
    data.frame(module = m,
               var_explained = unname(unlist(lapply(results, `[[`,
                                                    "var_explained"))[
                 match(m, colnames(eigengenes))]),
               kruskal_p = signif(kruskal.test(eigengenes[, m] ~
                                                 samples$group)$p.value, 3),
               stringsAsFactors = FALSE)
  }))
  group_test$kruskal_fdr <- signif(p.adjust(group_test$kruskal_p, "BH"), 3)
  write.csv(group_test, file.path(OUT_DIR, "module_group_association.csv"),
            row.names = FALSE)
  print(group_test)

  worst <- confound[which.max(abs(confound$rho)), ]
  message(sprintf("strongest confounding: %s vs %s, rho = %.3f (ceiling %.1f)",
                  worst$module, worst$marker_set, worst$rho, CONFOUND_CEILING))

  flagged <- unique(confound$module[abs(confound$rho) > CONFOUND_CEILING])
  if (length(flagged)) {
    message("modules exceeding the ceiling on any marker set: ",
            paste(flagged, collapse = ", "))
  } else {
    message("no module exceeds the ceiling on any marker set")
  }

  clean <- setdiff(colnames(eigengenes), flagged)
  writeLines(clean, file.path(OUT_DIR, "modules_passing_confounding.txt"))
  message(sprintf("%d of %d modules pass on every marker set",
                  length(clean), ncol(eigengenes)))
}

retained <- merge(assignment_all[assignment_all$retained, ],
                  de[, c("gene_symbol", "log_fc", "fdr")], by = "gene_symbol")
write.csv(retained, file.path(OUT_DIR, "retained_module_genes.csv"),
          row.names = FALSE)

for (m in unique(paste(retained$direction, retained$module, sep = "_M"))) {
  parts <- strsplit(m, "_M")[[1]]
  genes <- retained$gene_symbol[retained$direction == parts[1] &
                                  retained$module == parts[2]]
  writeLines(sort(genes), file.path(OUT_DIR, sprintf("genes_%s.txt", m)))
  message(sprintf("%s: %d genes", m, length(genes)))
}

saveRDS(list(results = results, composition = composition,
             identifiability = identifiability,
             marker_coverage = marker_coverage),
        file.path(OUT_DIR, "module_objects.rds"))
