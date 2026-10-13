[200~## =========================================================
## Fig. 1 simplified script
## CD96 expression across clinical and molecular groups
## Input:
##   CGGA_CD96.txt
##   TCGA_CD96.txt
## Output:
##   Fig1_CD96_clinical_expression.pdf
## =========================================================

library(ggplot2)
library(ggpubr)
library(gridExtra)

## ---------- 1. Read and clean data ----------

read_cd96 <- function(file, cohort_name) {
  dat <- read.table(file, sep = "\t", header = TRUE, check.names = FALSE)

  ## remove extra spaces in column names, e.g. "Age "
  names(dat) <- trimws(names(dat))

  dat$Cohort <- cohort_name
  dat$CD96_expression <- log2(as.numeric(dat$profile_values) + 0.001)

  ## clean common missing values
  dat[dat == "NA"] <- NA
  dat[dat == ""] <- NA

  ## harmonize 1p/19q spelling if needed
  if ("h1p19q_codeletion_status" %in% names(dat)) {
    dat$h1p19q_codeletion_status <- gsub("Non-code$", "Non-codel",
                                          dat$h1p19q_codeletion_status)
  }

  return(dat)
}

cgga <- read_cd96("CGGA_CD96.txt", "CGGA")
tcga <- read_cd96("TCGA_CD96.txt", "TCGA")

## ---------- 2. Common plotting function ----------

box_plot <- function(dat, group_col, xlab, title, group_order = NULL) {

  plot_dat <- dat[, c(group_col, "CD96_expression")]
  colnames(plot_dat) <- c("Group", "CD96_expression")

  plot_dat <- plot_dat[!is.na(plot_dat$Group) &
                         !is.na(plot_dat$CD96_expression), ]

  if (!is.null(group_order)) {
    plot_dat <- plot_dat[plot_dat$Group %in% group_order, ]
    plot_dat$Group <- factor(plot_dat$Group, levels = group_order)
  } else {
    plot_dat$Group <- factor(plot_dat$Group)
  }

  ## If 2 groups: Wilcoxon test; if >2 groups: Kruskal-Wallis test
  test_method <- ifelse(length(unique(plot_dat$Group)) == 2,
                        "wilcox.test", "kruskal.test")

  ggplot(plot_dat, aes(x = Group, y = CD96_expression, fill = Group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.75,
                 color = "black", linewidth = 0.4) +
    geom_jitter(width = 0.18, size = 0.8, alpha = 0.55) +
    stat_compare_means(method = test_method, label = "p.signif") +
    labs(
      title = title,
      x = xlab,
      y = "CD96 Gene Expression"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      panel.grid.minor = element_blank()
    )
}

## ---------- 3. Define group orders ----------

grade_order <- c("WHO II", "WHO III", "WHO IV")

histology_order <- c(
  "astrocytoma",
  "oligoastrocytoma",
  "oligodendroglioma",
  "glioblastoma"
)

idh_order <- c("Mutant", "Wildtype")

subtype_order <- c(
  "Classical",
  "Mesenchymal",
  "Neural",
  "Proneural"
)

## ---------- 4. Make CGGA panels ----------

p1 <- box_plot(cgga, "Grade", "WHO Grade",
               "CGGA: WHO Grade", grade_order)

p2 <- box_plot(cgga, "Histology", "Histology",
               "CGGA: Histology", histology_order)

p3 <- box_plot(cgga, "IDH_mutation_status", "IDH Mutation Status",
               "CGGA: IDH Status", idh_order)

p4 <- box_plot(cgga, "subtype", "Molecular Subtype",
               "CGGA: Subtype", subtype_order)

## ---------- 5. Make TCGA panels ----------

p5 <- box_plot(tcga, "Grade", "WHO Grade",
               "TCGA: WHO Grade", grade_order)

p6 <- box_plot(tcga, "Histology", "Histology",
               "TCGA: Histology", histology_order)

p7 <- box_plot(tcga, "IDH_mutation_status", "IDH Mutation Status",
               "TCGA: IDH Status", idh_order)

p8 <- box_plot(tcga, "subtype", "Molecular Subtype",
               "TCGA: Subtype", subtype_order)

## ---------- 6. Save Fig. 1 ----------

pdf("Fig1_CD96_clinical_expression.pdf", width = 13, height = 8)
grid.arrange(
  p1, p2, p3, p4,
  p5, p6, p7, p8,
  ncol = 4,
  top = "Figure 1. CD96 expression across clinical and molecular glioma phenotypes"
)
dev.off()

cat("Done: Fig1_CD96_clinical_expression.pdf\n")

