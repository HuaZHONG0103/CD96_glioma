# =============================================================================
# Reproduce the glioma immune-expression landscape in RStudio
# =============================================================================
#
# How to run:
#   1. Open this file in RStudio.
#   2. Click Source, or select all lines and click Run.
#   3. Input files are read from the input folder.
#   4. Results are saved automatically in the output folder.
#
# Important interpretation:
#   The observations are 1,027 bulk tumor samples (patients), not single cells.
#   Monocle is used here to construct a low-dimensional expression landscape.
#   The resulting graph should not be interpreted as a temporal cell lineage.
# =============================================================================

# Install the required packages once if needed:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("monocle")
# install.packages(c("ggplot2", "rstudioapi"))

suppressPackageStartupMessages({
  library(monocle)
  library(DDRTree)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 1. Project settings
# -----------------------------------------------------------------------------

# Identify the folder containing this script. This works when the file is
# opened and sourced from RStudio. Outside RStudio, the current working
# directory is used.
if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getSourceEditorContext()$path
  project_dir <- if (nzchar(script_path)) dirname(script_path) else getwd()
} else {
  project_dir <- getwd()
}

input_dir <- file.path(project_dir, "input")
output_dir <- file.path(project_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Use the supplied 2,000 ordering genes for direct reproduction.
# Change this to FALSE to recalculate them using differentialGeneTest().
use_existing_ordering_genes <- TRUE

# Number of CPU cores used only when ordering genes are recalculated.
n_cores <- 1

# Set a root state only if a scientifically justified root is known.
# Otherwise keep NA and let orderCells() use its default.
root_state <- NA_integer_

random_seed <- 2026
set.seed(random_seed)

input_file <- function(x) file.path(input_dir, x)
output_file <- function(x) file.path(output_dir, x)

required_files <- c(
  "data_preprocess.txt",
  "IS_information.txt",
  "IS.csv",
  "order_gene.txt"
)

missing_files <- required_files[!file.exists(file.path(input_dir, required_files))]
if (length(missing_files) > 0) {
  stop("Missing input file(s): ", paste(missing_files, collapse = ", "))
}

# -----------------------------------------------------------------------------
# 2. Read the expression matrix
# -----------------------------------------------------------------------------

raw <- read.delim(
  input_file("data_preprocess.txt"),
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(all(c("geneID", "symbol") %in% colnames(raw)))
stopifnot(!anyDuplicated(raw$geneID))

gene_annotation <- data.frame(
  gene_short_name = raw$symbol,
  row.names = raw$geneID,
  stringsAsFactors = FALSE
)

# Keep every column after geneID and symbol. The original script removed the
# first patient because row.names = 1 had already removed geneID before it
# deleted columns 1 and 2.
patient_columns <- setdiff(colnames(raw), c("geneID", "symbol"))
expr_matrix <- as.matrix(raw[, patient_columns, drop = FALSE])
storage.mode(expr_matrix) <- "numeric"
rownames(expr_matrix) <- raw$geneID

# -----------------------------------------------------------------------------
# 3. Read and align patient metadata
# -----------------------------------------------------------------------------

sample_sheet <- read.delim(
  input_file("IS_information.txt"),
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(all(c("UniqueCell_ID", "majorCluster") %in% colnames(sample_sheet)))
stopifnot(!anyDuplicated(sample_sheet$UniqueCell_ID))
rownames(sample_sheet) <- sample_sheet$UniqueCell_ID

missing_ids <- setdiff(rownames(sample_sheet), colnames(expr_matrix))
if (length(missing_ids) > 0) {
  stop("These patient IDs are missing from the expression matrix: ",
       paste(head(missing_ids, 10), collapse = ", "))
}

# Reorder the expression columns to match the metadata exactly.
expr_matrix <- expr_matrix[, rownames(sample_sheet), drop = FALSE]
stopifnot(identical(colnames(expr_matrix), rownames(sample_sheet)))

cat("Genes:", nrow(expr_matrix), "\n")
cat("Patients:", ncol(expr_matrix), "\n")
print(table(sample_sheet$majorCluster))

# -----------------------------------------------------------------------------
# 4. PCA
# -----------------------------------------------------------------------------

# Rows must be patients and columns must be genes for patient-level PCA.
patient_by_gene <- t(expr_matrix)
pca_results <- prcomp(patient_by_gene, center = TRUE, scale. = FALSE)

# Patient coordinates and gene loadings are different outputs.
write.csv(
  pca_results$x,
  output_file("patient_PCA_coordinates.csv"),
  quote = FALSE
)
write.csv(
  pca_results$rotation,
  output_file("gene_PCA_loadings.csv"),
  quote = FALSE
)

# -----------------------------------------------------------------------------
# 5. Construct the Monocle 2 CellDataSet
# -----------------------------------------------------------------------------

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

cds <- newCellDataSet(
  expr_matrix,
  phenoData = pd,
  featureData = fd,
  expressionFamily = uninormal()
)

cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- rownames(
  subset(fData(cds), num_cells_expressed >= 10)
)

# -----------------------------------------------------------------------------
# 6. Select ordering genes
# -----------------------------------------------------------------------------

if (use_existing_ordering_genes) {
  ordering_gene_table <- read.table(
    input_file("order_gene.txt"),
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  ordering_genes <- as.character(ordering_gene_table[[1]])
  ordering_genes <- unique(ordering_genes)
  ordering_genes <- intersect(ordering_genes, rownames(cds))
  cat("Using", length(ordering_genes), "supplied ordering genes.\n")
} else {
  de_results <- differentialGeneTest(
    cds[expressed_genes, ],
    fullModelFormulaStr = "~ majorCluster",
    cores = n_cores
  )
  de_results <- de_results[
    order(de_results$qval, de_results$pval, na.last = TRUE),
    , drop = FALSE
  ]
  ordering_genes <- head(rownames(de_results), 2000)
  write.csv(
    de_results,
    output_file("subtype_differential_gene_test.csv"),
    quote = FALSE
  )
}

if (length(ordering_genes) < 2) stop("Too few ordering genes are available.")

write.table(
  ordering_genes,
  output_file("ordering_genes_used.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cds <- setOrderingFilter(cds, ordering_genes)

# -----------------------------------------------------------------------------
# 7. DDRTree embedding and graph states
# -----------------------------------------------------------------------------

cds <- reduceDimension(
  cds,
  max_components = 2,
  norm_method = "none",
  method = "DDRTree"
)

if (is.na(root_state)) {
  cds <- orderCells(cds)
} else {
  cds <- orderCells(cds, root_state = root_state)
}

saveRDS(cds, output_file("CDS_DDRTree_ordered.rds"))

# -----------------------------------------------------------------------------
# 8. Plot the landscape by major immune subtype, State and pseudotime
# -----------------------------------------------------------------------------

major_colors <- c(
  IS1 = "#E6AB02",
  IS2 = "#1F78B4",
  IS3 = "#E31A1C"
)

state_colors <- c(
  "1" = "#CFD600",
  "2" = "#9E9F9F",
  "3" = "#E52622",
  "4" = "#7570B3",
  "5" = "#EABF00",
  "6" = "#1CAC55",
  "7" = "#FFAA00",
  "8" = "#673327",
  "9" = "#006DBB"
)

p_major <- plot_cell_trajectory(
  cds,
  color_by = "majorCluster",
  cell_size = 1.25
) +
  scale_color_manual(values = major_colors, drop = FALSE) +
  theme_classic(base_size = 11)

ggsave(
  output_file("monocle_majorCluster_reproduced.pdf"),
  plot = p_major,
  width = 5.5,
  height = 5
)

p_state <- plot_cell_trajectory(
  cds,
  color_by = "State",
  cell_size = 1.25
) +
  scale_color_manual(values = state_colors, drop = FALSE) +
  theme_classic(base_size = 11)

ggsave(
  output_file("monocle_State_reproduced.pdf"),
  plot = p_state,
  width = 5.5,
  height = 5
)

p_pseudotime <- plot_cell_trajectory(
  cds,
  color_by = "Pseudotime",
  cell_size = 1.25
) + theme_classic(base_size = 11)

ggsave(
  output_file("monocle_pseudotime_reproduced.pdf"),
  plot = p_pseudotime,
  width = 5.5,
  height = 5
)

# Display the three plots in the RStudio Plots panel.
print(p_major)
print(p_state)
print(p_pseudotime)

# -----------------------------------------------------------------------------
# 9. Add manual IS branch labels without overwriting Monocle State
# -----------------------------------------------------------------------------

manual_labels <- read.csv(
  input_file("IS.csv"),
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(all(c("ID", "IS") %in% colnames(manual_labels)))
stopifnot(!anyDuplicated(manual_labels$ID))

label_index <- match(pData(cds)$UniqueCell_ID, manual_labels$ID)
if (anyNA(label_index)) {
  stop("Manual labels are missing for ", sum(is.na(label_index)), " patients.")
}

# Create a new variable. Do not replace cds$State.
pData(cds)$Sub_IS <- factor(
  manual_labels$IS[label_index],
  levels = c("IS1", "IS", "IS2_4", "IS2_6", "IS2_7", "IS3_1", "IS3_4")
)

sub_is_colors <- c(
  IS1 = "#1CAC55",
  IS = "#9E9F9F",
  IS2_4 = "#7570B3",
  IS2_6 = "#006DBB",
  IS2_7 = "#FFAA00",
  IS3_1 = "#E52622",
  IS3_4 = "#673327"
)

p_sub_is <- plot_cell_trajectory(
  cds,
  color_by = "Sub_IS",
  cell_size = 1.25
) +
  scale_color_manual(values = sub_is_colors, drop = FALSE) +
  theme_classic(base_size = 11)

print(p_sub_is)

ggsave(
  output_file("monocle_Sub_IS_reproduced.pdf"),
  plot = p_sub_is,
  width = 5.8,
  height = 5
)

# -----------------------------------------------------------------------------
# 10. Export patient coordinates and labels
# -----------------------------------------------------------------------------

ddr_coordinates <- t(reducedDimS(cds))

patient_results <- data.frame(
  ID = pData(cds)$UniqueCell_ID,
  majorCluster = pData(cds)$majorCluster,
  State = pData(cds)$State,
  Pseudotime = pData(cds)$Pseudotime,
  Sub_IS = pData(cds)$Sub_IS,
  Component_1 = ddr_coordinates[, 1],
  Component_2 = ddr_coordinates[, 2],
  num_genes_expressed = pData(cds)$num_genes_expressed,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write.csv(
  patient_results,
  output_file("monocle_patient_coordinates_and_labels.csv"),
  row.names = FALSE,
  quote = FALSE
)

saveRDS(cds, output_file("CDS_DDRTree_with_Sub_IS.rds"))
writeLines(capture.output(sessionInfo()), output_file("sessionInfo.txt"))

cat("\nAnalysis completed.\n")
cat("Results were saved in:\n", normalizePath(output_dir), "\n")
