#!/usr/bin/env Rscript

# Complete reproduction of published Figure 5 (panels A-D)
# CD96-related T-cell immunity and immune-checkpoint markers in glioma
#
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")

# Required input files
#   1. CGGA.mRNAseq_325.RSEM-genes.20190909.txt  (gene symbols x 325 samples)
#   2. easy_input_expr.csv                       (TCGA Ensembl IDs x samples)
#   3. easy_input_GO.txt                         (seven GO IDs)
#   4. corrgram_CGGA.txt                         (checkpoint expression)
#   5. corrgram_TCGA.txt                         (checkpoint expression)
#
# Run
#   Rscript Figure5_CD96_correlation_analysis.R \
#     CGGA.mRNAseq_325.RSEM-genes.20190909.txt \
#     easy_input_expr.csv easy_input_GO.txt \
#     corrgram_CGGA.txt corrgram_TCGA.txt Figure5_outputs

required <- c(
  "GSVA", "GSEABase", "AnnotationDbi", "org.Hs.eg.db",
  "GO.db", "corrplot", "RColorBrewer"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install required packages: ", paste(missing, collapse = ", "))
}

args <- c(
  "CGGA_expression.txt",
  "TCGA_expression.csv",
  "easy_input_GO.txt",
  "corrgram_CGGA.txt",
  "corrgram_TCGA.txt",
  "Figure5_outputs"
)
cgga_file <- args[[1L]]
tcga_file <- args[[2L]]
go_file <- args[[3L]]
cgga_checkpoint_file <- args[[4L]]
tcga_checkpoint_file <- args[[5L]]
output_dir <- if (length(args) >= 6L) args[[6L]] else "Figure5_outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

go_order <- scan(go_file, what = character(), quiet = TRUE)
expected_go <- c(
  "GO:0042088", "GO:0042092", "GO:2000556", "GO:2000553",
  "GO:0002860", "GO:0002842", "GO:0002852"
)
if (!identical(go_order, expected_go)) {
  warning("GO list/order differs from the seven pathways in the published figure.")
}

read_numeric_matrix <- function(path, sep) {
  x <- utils::read.table(
    path, sep = sep, header = TRUE, row.names = 1,
    check.names = FALSE, stringsAsFactors = FALSE, quote = "",
    comment.char = "", fill = FALSE
  )
  x[] <- lapply(x, function(z) as.numeric(as.character(z)))
  as.matrix(x)
}

collapse_duplicate_genes <- function(x) {
  if (!anyDuplicated(rownames(x))) return(x)
  rowsum(x, group = rownames(x), reorder = FALSE) /
    as.vector(table(factor(rownames(x), levels = unique(rownames(x)))))
}

read_cgga <- function(path) {
  x <- read_numeric_matrix(path, sep = "\t")
  rownames(x) <- trimws(rownames(x))
  collapse_duplicate_genes(x)
}

read_tcga <- function(path) {
  x <- read_numeric_matrix(path, sep = ",")
  ensembl <- sub("\\..*$", "", rownames(x))
  map <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(ensembl), keytype = "ENSEMBL", columns = "SYMBOL"
  )
  map <- map[!is.na(map$SYMBOL) & nzchar(map$SYMBOL), ]
  map <- map[!duplicated(map$ENSEMBL), ]
  keep <- match(map$ENSEMBL, ensembl)
  x <- x[keep, , drop = FALSE]
  rownames(x) <- map$SYMBOL

  # Retain primary tumors only when TCGA sample types are encoded in barcodes.
  sample_type <- substr(colnames(x), 14L, 15L)
  if (any(sample_type == "01")) x <- x[, sample_type == "01", drop = FALSE]
  collapse_duplicate_genes(x)
}

make_go_gene_sets <- function(go_ids) {
  go_map <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = go_ids, keytype = "GOALL", columns = c("SYMBOL", "ONTOLOGYALL")
  )
  go_map <- go_map[
    !is.na(go_map$SYMBOL) & go_map$ONTOLOGYALL == "BP",
    c("GOALL", "SYMBOL")
  ]
  sets <- lapply(go_ids, function(id) unique(go_map$SYMBOL[go_map$GOALL == id]))
  names(sets) <- go_ids

  # Historical replication rule used in the supplied original code.
  # GO:0002852 is obsolete/not returned by current annotation resources.
  if ("GO:0002852" %in% go_ids && length(sets[["GO:0002852"]]) == 0L) {
    sets[["GO:0002852"]] <- "MUC4"
  }
  sets
}

run_gsva_compatible <- function(expr, gene_sets) {
  if (utils::packageVersion("GSVA") >= "1.50.0") {
    param <- GSVA::gsvaParam(
      exprData = expr,
      geneSets = gene_sets,
      kcdf = "Gaussian",
      minSize = 1L,
      maxSize = Inf
    )
    GSVA::gsva(param, verbose = FALSE)
  } else {
    GSVA::gsva(
      expr = expr, gset.idx.list = gene_sets,
      method = "gsva", kcdf = "Gaussian", min.sz = 1L,
      verbose = FALSE
    )
  }
}

prepare_pathway_analysis <- function(expr, cohort, go_sets) {
  if (!"CD96" %in% rownames(expr)) stop("CD96 is absent from ", cohort, " expression.")

  usable_sets <- lapply(go_sets, intersect, y = rownames(expr))
  sizes <- lengths(usable_sets)
  if (any(sizes == 0L)) {
    stop(cohort, " has no mapped genes for: ", paste(names(sizes)[sizes == 0L], collapse = ", "))
  }

  # TCGA input is already log2(FPKM+1). CGGA RSEM values are non-negative;
  # log2(x+1) stabilizes their scale before the Gaussian-kernel GSVA step.
  if (max(expr, na.rm = TRUE) > 50) expr <- log2(expr + 1)

  cd96 <- expr["CD96", ]
  ord <- order(cd96, na.last = NA)
  expr <- expr[, ord, drop = FALSE]
  cd96 <- cd96[ord]
  scores <- run_gsva_compatible(expr, usable_sets)
  scores <- scores[names(go_sets), , drop = FALSE]

  tests <- lapply(seq_len(nrow(scores)), function(i) {
    stats::cor.test(scores[i, ], cd96, method = "spearman", exact = FALSE)
  })
  go_info <- AnnotationDbi::select(
    GO.db::GO.db,
    keys = rownames(scores),
    columns = "TERM",
    keytype = "GOID"
  )
  go_info <- go_info[!duplicated(go_info$GOID), , drop = FALSE]

  stats_table <- data.frame(
    cohort = rep(cohort, nrow(scores)),
    GO_ID = rownames(scores),
    GO_term = go_info$TERM[match(rownames(scores), go_info$GOID)],
    genes_in_annotation = lengths(go_sets)[rownames(scores)],
    genes_in_expression = sizes[rownames(scores)],
    spearman_rho = vapply(tests, function(z) unname(z$estimate), numeric(1)),
    p_value = vapply(tests, function(z) z$p.value, numeric(1)),
    stringsAsFactors = FALSE
  )
  stats_table$FDR <- p.adjust(stats_table$p_value, method = "BH")
  list(scores = scores, cd96 = cd96, stats = stats_table)
}

sig_stars <- function(p) {
  ifelse(p < 1e-4, "****", ifelse(p < 1e-3, "***",
    ifelse(p < 1e-2, "**", ifelse(p < 0.05, "*", "ns"))))
}

draw_sector <- function(cx, cy, radius_in, fraction, positive, color) {
  fraction <- min(1, max(0, fraction))

  # Convert a physical radius (inches) to separate x/y data-unit radii.
  # This keeps each glyph circular even though the heatmap axes have very
  # different numerical ranges (hundreds of columns versus seven rows).
  usr <- par("usr")
  pin <- par("pin")
  radius_x <- radius_in * diff(usr[1:2]) / pin[1]
  radius_y <- radius_in * diff(usr[3:4]) / pin[2]

  full_theta <- seq(0, 2 * pi, length.out = 121L)
  polygon(
    cx + radius_x * cos(full_theta),
    cy + radius_y * sin(full_theta),
    col = "white", border = NA
  )

  direction <- if (positive) -1 else 1
  theta <- seq(pi / 2, pi / 2 + direction * 2 * pi * fraction,
               length.out = max(3L, ceiling(120 * fraction)))
  if (fraction > 0) {
    polygon(
      c(cx, cx + radius_x * cos(theta)),
      c(cy, cy + radius_y * sin(theta)),
      col = color, border = NA
    )
  }

  lines(
    cx + radius_x * cos(full_theta),
    cy + radius_y * sin(full_theta),
    col = "grey35", lwd = 0.8
  )

  invisible(c(radius_x = radius_x, radius_y = radius_y))
}

draw_gsva_panel <- function(result, cohort, panel_letter) {
  z <- as.matrix(result$scores)
  storage.mode(z) <- "numeric"
  if (length(dim(z)) != 2L) {
    stop(cohort, ": result$scores is not two-dimensional. Class: ",
         paste(class(result$scores), collapse = ", "))
  }

  nr <- nrow(z)
  nc <- ncol(z)
  if (nr == 0L || nc == 0L) stop(cohort, ": GSVA score matrix is empty.")

  pal <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(201)
  lim <- max(abs(z), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1

  # Subassignment preserves matrix dimensions; pmin/pmax with a scalar may not.
  scaled <- z
  scaled[scaled > lim] <- lim
  scaled[scaled < -lim] <- -lim
  idx <- round((scaled + lim) / (2 * lim) * 200) + 1L
  idx <- matrix(as.integer(idx), nrow = nr, ncol = nc, dimnames = dimnames(z))
  idx[!is.finite(z)] <- 101L
  idx[idx < 1L] <- 1L
  idx[idx > 201L] <- 201L

  par(mar = c(2.2, 7.5, 2.5, 2), xaxs = "i", yaxs = "i")
  plot(c(0, nc + max(10, nc * 0.07)), c(0, nr), type = "n",
       axes = FALSE, xlab = "", ylab = "")
  for (j in seq_len(nc)) for (i in seq_len(nr)) {
    rect(j - 1, nr - i, j, nr - i + 1, col = pal[idx[i, j]], border = NA)
  }
  text(-0.01 * nc, nr - seq_len(nr) + 0.5, rownames(z), pos = 2, cex = 0.75, xpd = NA)
  title(main = paste0(panel_letter, "   ", cohort, " — CD96-related T-cell immunity"),
        adj = 0, font.main = 2, cex.main = 1.05)

  pie_x <- nc + max(5, nc * 0.03)
  radius_in <- 0.12
  stats_table <- result$stats[
    match(rownames(z), result$stats$GO_ID), , drop = FALSE
  ]
  for (i in seq_len(nr)) {
    rho <- stats_table$spearman_rho[i]
    p_value <- stats_table$p_value[i]
    if (!is.finite(rho)) rho <- 0
    col <- if (rho >= 0) grDevices::adjustcolor("#E84A5F", alpha.f = 0.82) else
      grDevices::adjustcolor("#70AD47", alpha.f = 0.82)
    glyph_size <- draw_sector(
      pie_x, nr - i + 0.5, radius_in,
      abs(rho), rho >= 0, col
    )
    text(pie_x + glyph_size[["radius_x"]] * 1.25, nr - i + 0.5,
         sig_stars(p_value), pos = 4, cex = 0.55)
  }
}

checkpoint_result <- function(path, genes, cohort) {
  x <- utils::read.delim(path, row.names = 1, check.names = FALSE)
  missing <- setdiff(genes, names(x))
  if (length(missing)) stop(cohort, " checkpoint table missing: ", paste(missing, collapse = ", "))
  x <- x[, genes, drop = FALSE]
  r <- stats::cor(x, use = "pairwise.complete.obs", method = "pearson")
  list(data = x, r = r)
}

draw_checkpoint_panel <- function(r, cohort, panel_letter) {
  corrplot::corrplot.mixed(
    r, lower = "circle", upper = "number", tl.pos = "lt", diag = "u",
    tl.col = "#C51B29", number.cex = 0.65,
    title = paste0(panel_letter, "   ", cohort, " — immune checkpoints"),
    mar = c(0, 0, 2.2, 0)
  )
}

message("Reading full expression matrices...")
cgga_expr <- read_cgga(cgga_file)
tcga_expr <- read_tcga(tcga_file)
go_sets <- make_go_gene_sets(go_order)

message("Running GSVA and Spearman correlations...")
cgga_gsva <- prepare_pathway_analysis(cgga_expr, "CGGA", go_sets)
tcga_gsva <- prepare_pathway_analysis(tcga_expr, "TCGA", go_sets)

cgga_cp_genes <- c("CD96", "PD-1", "CTLA4", "TIGIT", "TIM-3", "NR2F6", "GITR")
tcga_cp_genes <- c("CD96", "PD-1", "CTLA-4", "TIGIT", "TIM-3", "NR2F6", "GITR")
cgga_cp <- checkpoint_result(cgga_checkpoint_file, cgga_cp_genes, "CGGA")
tcga_cp <- checkpoint_result(tcga_checkpoint_file, tcga_cp_genes, "TCGA")

utils::write.csv(cgga_gsva$scores, file.path(output_dir, "Figure5A_CGGA_GSVA_scores.csv"))
utils::write.csv(tcga_gsva$scores, file.path(output_dir, "Figure5B_TCGA_GSVA_scores.csv"))
utils::write.csv(cgga_gsva$stats, file.path(output_dir, "Figure5A_CGGA_pathway_correlations.csv"), row.names = FALSE)
utils::write.csv(tcga_gsva$stats, file.path(output_dir, "Figure5B_TCGA_pathway_correlations.csv"), row.names = FALSE)
utils::write.csv(cgga_cp$r, file.path(output_dir, "Figure5C_CGGA_checkpoint_correlations.csv"))
utils::write.csv(tcga_cp$r, file.path(output_dir, "Figure5D_TCGA_checkpoint_correlations.csv"))

grDevices::pdf(file.path(output_dir, "Figure5_complete_reproduction.pdf"), width = 13, height = 9)
layout(matrix(1:4, nrow = 2, byrow = TRUE), heights = c(1, 1.05))
draw_gsva_panel(cgga_gsva, "CGGA", "A")
draw_gsva_panel(tcga_gsva, "TCGA", "B")
draw_checkpoint_panel(cgga_cp$r, "CGGA", "C")
draw_checkpoint_panel(tcga_cp$r, "TCGA", "D")
grDevices::dev.off()

message("Completed all Figure 5 panels. Output directory: ", normalizePath(output_dir))
