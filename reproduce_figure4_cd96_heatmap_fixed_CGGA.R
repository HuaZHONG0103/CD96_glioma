# ============================================================
# Reproduce Figure 4: CD96-related immune-gene heatmap
# CGGA and TCGA glioma cohorts
# ============================================================
# Concept:
#   1) Order samples from low to high CD96 expression.
#   2) For every immune-response gene, calculate Spearman correlation
#      with CD96 expression across samples.
#   3) Keep genes significantly correlated with CD96.
#   4) Row-standardize expression into z-scores and draw a heatmap.
# ============================================================

suppressPackageStartupMessages({
  library(pheatmap)
  library(RColorBrewer)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# -----------------------------
# User settings
# -----------------------------
workdir <- "."   # change to the folder containing the input files
setwd(workdir)

# Input files supplied for this analysis
cgga_expr_file <- "easy_input_CGGA.csv"       # genes x samples; gene symbols in column 1
cgga_clin_file <- "CGGA_CD96.csv"             # clinical annotation; sample ID column = ID

tcga_expr_file <- "TCGA_533_exp_trans.csv"    # samples x Entrez genes; sample IDs in column 1
tcga_clin_file <- "TCGA_533_clin.csv"         # clinical annotation; sample ID column = Sample ID

gene_list_file <- "easy_input_genes.txt"      # immune-response gene symbols

# -----------------------------
# Colors matching the manuscript style
# -----------------------------
heatmap_BlWtRd <- c("#6699CC", "white", "#FF3C38")
red    <- "#EA6767"
green  <- "#70C17A"
blue   <- "#445EAD"
grey   <- "#7A7A7A"
orange <- "#F7C07C"
yellow <- "#CEDC7C"

ann_colors <- list(
  CD96 = colorRampPalette(heatmap_BlWtRd)(128),
  IDH = c("Wildtype" = red, "Mutant" = green, "NA" = grey),
  Subtype = c(
    "Classical" = yellow,
    "Mesenchymal" = red,
    "Neural" = orange,
    "Proneural" = blue,
    "G-CIMP" = green,
    "NA" = grey
  )
)

# -----------------------------
# Helper functions
# -----------------------------
# Repair gene symbols that were automatically converted to dates by Excel.
# In this CGGA file, for example, MARC1 and MARCH1 both became "1-Mar",
# and MARC2 and MARCH2 both became "2-Mar".
repair_excel_gene_symbols <- function(genes) {
  genes <- trimws(as.character(genes))

  genes[genes == "1-Dec"] <- "DEC1"
  genes[genes == "15-Sep"] <- "SEPT15"

  # The rows are alphabetically ordered in the supplied expression matrix:
  # first 1-Mar/2-Mar = MARC1/MARC2; second 1-Mar/2-Mar = MARCH1/MARCH2.
  idx_1mar <- which(genes == "1-Mar")
  if (length(idx_1mar) >= 1) genes[idx_1mar[1]] <- "MARC1"
  if (length(idx_1mar) >= 2) genes[idx_1mar[-1]] <- "MARCH1"

  idx_2mar <- which(genes == "2-Mar")
  if (length(idx_2mar) >= 1) genes[idx_2mar[1]] <- "MARC2"
  if (length(idx_2mar) >= 2) genes[idx_2mar[-1]] <- "MARCH2"

  for (n in c(3:11)) {
    genes[genes == paste0(n, "-Mar")] <- paste0("MARCH", n)
  }
  for (n in c(2:12)) {
    genes[genes == paste0(n, "-Sep")] <- paste0("SEPT", n)
  }

  genes
}

read_cgga_expression <- function(file) {
  # Do NOT use row.names = 1 here: the supplied file contains duplicated
  # Excel-converted labels (1-Mar and 2-Mar), which causes read.csv() to stop.
  raw <- read.csv(file, header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE, row.names = NULL)

  genes <- repair_excel_gene_symbols(raw[[1]])
  raw <- raw[, -1, drop = FALSE]

  # Convert expression columns safely to numeric.
  raw[] <- lapply(raw, function(x) as.numeric(as.character(x)))
  expr <- as.matrix(raw)
  rownames(expr) <- genes

  # Final safeguard for any other accidental duplicated symbols.
  if (anyDuplicated(rownames(expr))) {
    warning("Duplicated gene symbols remain after repair; keeping the first occurrence.")
    expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  }

  expr
}

read_tcga_expression <- function(file) {
  raw <- read.csv(file, header = TRUE, row.names = 1,
                  check.names = FALSE, stringsAsFactors = FALSE)
  raw <- as.matrix(raw)
  mode(raw) <- "numeric"

  # TCGA file has samples as rows and Entrez IDs as columns.
  # Convert Entrez IDs to gene symbols and then transpose to genes x samples.
  entrez_ids <- colnames(raw)
  map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = entrez_ids,
    keytype = "ENTREZID",
    columns = c("SYMBOL")
  )
  map <- map[!is.na(map$SYMBOL), ]
  map <- map[!duplicated(map$SYMBOL), ]
  map <- map[map$ENTREZID %in% colnames(raw), ]

  expr <- t(raw[, map$ENTREZID, drop = FALSE])
  rownames(expr) <- map$SYMBOL
  expr
}

read_gene_list <- function(file) {
  genes <- read.table(file, header = TRUE, sep = "\t",
                      check.names = FALSE, stringsAsFactors = FALSE)[[1]]
  unique(genes[!is.na(genes) & genes != ""])
}

standardize_idh <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "NA"
  x <- ifelse(x %in% c("WT", "Wildtype", "wildtype", "IDHwt"), "Wildtype",
              ifelse(x == "NA", "NA", "Mutant"))
  x
}

zscore_by_gene <- function(mat, clip = 3) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z[z > clip] <- clip
  z[z < -clip] <- -clip
  z
}

reproduce_cd96_heatmap <- function(expr,
                                   clin_file,
                                   gene_list,
                                   sample_id_col,
                                   idh_col,
                                   subtype_col,
                                   output_prefix,
                                   rho_cutoff = 0.30,
                                   p_cutoff = 0.05) {
  clin <- read.csv(clin_file, header = TRUE, check.names = FALSE,
                   stringsAsFactors = FALSE)
  rownames(clin) <- clin[[sample_id_col]]

  # Keep only samples present in both expression and clinical annotation.
  common_samples <- intersect(colnames(expr), rownames(clin))
  expr <- expr[, common_samples, drop = FALSE]
  clin <- clin[common_samples, , drop = FALSE]

  if (!"CD96" %in% rownames(expr)) {
    stop("CD96 is not present in the expression matrix after gene-symbol processing.")
  }

  # Order samples by increasing CD96 expression.
  cd96 <- as.numeric(expr["CD96", ])
  names(cd96) <- colnames(expr)
  sample_order <- names(sort(cd96, decreasing = FALSE))

  # Restrict to immune-response genes available in the expression matrix.
  immune_genes <- intersect(rownames(expr), gene_list)
  immune_expr <- expr[immune_genes, sample_order, drop = FALSE]

  # Spearman correlation: each immune gene versus CD96.
  cor_res <- do.call(rbind, lapply(rownames(immune_expr), function(gene) {
    test <- suppressWarnings(cor.test(
      as.numeric(immune_expr[gene, ]),
      as.numeric(cd96[sample_order]),
      method = "spearman",
      exact = FALSE
    ))
    data.frame(gene = gene,
               rho = unname(test$estimate),
               p = test$p.value,
               stringsAsFactors = FALSE)
  }))

  write.csv(cor_res, paste0(output_prefix, "_CD96_correlations.csv"), row.names = FALSE)

  pos_genes <- cor_res[cor_res$rho > rho_cutoff & cor_res$p < p_cutoff, c("gene", "rho")]
  neg_genes <- cor_res[cor_res$rho < -rho_cutoff & cor_res$p < p_cutoff, c("gene", "rho")]

  # Sorting by rho helps create the low-to-high wedge pattern seen in Figure 4.
  pos_genes <- pos_genes[order(pos_genes$rho, decreasing = FALSE), ]
  neg_genes <- neg_genes[order(neg_genes$rho, decreasing = FALSE), ]
  plot_genes <- c(pos_genes$gene, neg_genes$gene)

  if (length(plot_genes) == 0) {
    stop("No genes passed the correlation cutoff. Try lowering rho_cutoff or check input data.")
  }

  ann_col <- data.frame(
    CD96 = as.numeric(cd96[sample_order]),
    Subtype = clin[sample_order, subtype_col],
    IDH = standardize_idh(clin[sample_order, idh_col]),
    row.names = sample_order,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  ann_col$Subtype[is.na(ann_col$Subtype) | ann_col$Subtype == ""] <- "NA"

  heat_mat <- zscore_by_gene(immune_expr[plot_genes, sample_order, drop = FALSE], clip = 3)

  pdf_file <- paste0(output_prefix, "_correlationHeatmap_reproduced.pdf")
  pheatmap(
    mat = heat_mat,
    scale = "none",
    border_color = NA,
    fontsize_row = 1,
    color = colorRampPalette(heatmap_BlWtRd)(128),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    gaps_row = if (nrow(neg_genes) > 0) nrow(pos_genes) else NULL,
    filename = pdf_file,
    width = 8,
    height = 8
  )

  message("Saved: ", pdf_file)
  message("Genes plotted: ", length(plot_genes),
          " | positive: ", nrow(pos_genes),
          " | negative: ", nrow(neg_genes),
          " | samples: ", length(sample_order))

  invisible(list(
    correlations = cor_res,
    positive_genes = pos_genes,
    negative_genes = neg_genes,
    sample_order = sample_order,
    heatmap_matrix = heat_mat
  ))
}

# -----------------------------
# Run CGGA heatmap
# -----------------------------
gene_list <- read_gene_list(gene_list_file)

cgga_expr <- read_cgga_expression(cgga_expr_file)
cgga_result <- reproduce_cd96_heatmap(
  expr = cgga_expr,
  clin_file = cgga_clin_file,
  gene_list = gene_list,
  sample_id_col = "ID",
  idh_col = "IDH_mutation_status",
  subtype_col = "Subtype",
  output_prefix = "CGGA"
)

# -----------------------------
# Run TCGA heatmap
# -----------------------------
tcga_expr <- read_tcga_expression(tcga_expr_file)
tcga_result <- reproduce_cd96_heatmap(
  expr = tcga_expr,
  clin_file = tcga_clin_file,
  gene_list = gene_list,
  sample_id_col = "Sample ID",
  idh_col = "IDH_mutation_status",
  subtype_col = "Subtype",
  output_prefix = "TCGA"
)

sessionInfo()
