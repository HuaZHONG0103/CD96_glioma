# ============================================================================
# Reproduce Figure 7: Kaplan-Meier survival analysis of CD96 in glioma
#
# Panels:
#   A-D: CGGA, all grades / WHO IV / WHO II / WHO III
#   E-H: TCGA, all grades / WHO IV / WHO II / WHO III
#
# Input files:
#   CD96_CGGA.txt  (UTF-16LE, tab-delimited)
#   TCGA_CD96.txt  (tab-delimited)
#
# Main outputs:
#   Figure7_CD96_Kaplan_Meier.pdf
#   Figure7_CD96_Kaplan_Meier.png
#   Figure7_KM_statistics.csv
#   Figure7_sample_counts.csv
# ============================================================================

# -----------------------------
# 1. User settings
# -----------------------------

# Change this path if the files are stored elsewhere.
data_dir <- "/Users/apple/Desktop/zunyi/CD96/survival"

tcga_file <- file.path(data_dir, "TCGA_CD96.txt")
cgga_file <- file.path(data_dir, "CD96_CGGA.txt")

output_dir <- file.path(data_dir, "Figure7_output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Both uploaded datasets appear to use 1 = death/event and 0 = censored.
# If CGGA documentation says Censor uses the opposite definition, set this FALSE.
cgga_one_means_event <- TRUE

# TRUE reproduces the nearly equal high/low sample sizes shown in the paper:
# the CD96 median is recalculated separately within every plotted subgroup.
use_subgroup_specific_median <- TRUE

# Colors similar to the published figure.
high_color <- "#E76F51"
low_color  <- "#45B8B0"


# -----------------------------
# 2. Install/load packages
# -----------------------------

required_packages <- c("survival", "survminer", "ggplot2")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org",
    dependencies = TRUE
  )
}

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
})


# -----------------------------
# 3. Check and read input files
# -----------------------------

if (!file.exists(tcga_file)) {
  stop("TCGA file not found: ", tcga_file)
}

if (!file.exists(cgga_file)) {
  stop("CGGA file not found: ", cgga_file)
}

TCGA <- read.delim(
  tcga_file,
  header = TRUE,
  sep = "\t",
  na.strings = c("NA", "", "N/A"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

CGGA <- read.delim(
  cgga_file,
  header = TRUE,
  sep = "\t",
  fileEncoding = "UTF-16LE",
  na.strings = c("NA", "", "N/A"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Remove leading/trailing spaces in column names.
names(TCGA) <- trimws(names(TCGA))
names(CGGA) <- trimws(names(CGGA))

required_tcga <- c("ID", "Grade", "OS", "Status", "profile_values")
required_cgga <- c("ID", "Grade", "OS", "Censor", "profile_values")

missing_tcga_columns <- setdiff(required_tcga, names(TCGA))
missing_cgga_columns <- setdiff(required_cgga, names(CGGA))

if (length(missing_tcga_columns) > 0) {
  stop("Missing TCGA columns: ", paste(missing_tcga_columns, collapse = ", "))
}

if (length(missing_cgga_columns) > 0) {
  stop("Missing CGGA columns: ", paste(missing_cgga_columns, collapse = ", "))
}


# -----------------------------
# 4. Clean variables
# -----------------------------

clean_grade <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[[:space:]_-]", "", x)
  x <- gsub("GRADE", "WHO", x)
  x[x %in% c("2", "II", "WHO2")] <- "WHOII"
  x[x %in% c("3", "III", "WHO3")] <- "WHOIII"
  x[x %in% c("4", "IV", "WHO4", "GBM")] <- "WHOIV"
  x
}

TCGA$OS             <- suppressWarnings(as.numeric(TCGA$OS))
TCGA$Status         <- suppressWarnings(as.numeric(TCGA$Status))
TCGA$profile_values <- suppressWarnings(as.numeric(TCGA$profile_values))
TCGA$Grade          <- clean_grade(TCGA$Grade)

CGGA$OS             <- suppressWarnings(as.numeric(CGGA$OS))
CGGA$Censor         <- suppressWarnings(as.numeric(CGGA$Censor))
CGGA$profile_values <- suppressWarnings(as.numeric(CGGA$profile_values))
CGGA$Grade          <- clean_grade(CGGA$Grade)

if (cgga_one_means_event) {
  CGGA$Event <- CGGA$Censor
} else {
  CGGA$Event <- 1 - CGGA$Censor
}

cat("\nTCGA dimensions:", nrow(TCGA), "rows x", ncol(TCGA), "columns\n")
cat("CGGA dimensions:", nrow(CGGA), "rows x", ncol(CGGA), "columns\n")

cat("\nTCGA grade counts:\n")
print(table(TCGA$Grade, useNA = "ifany"))
cat("\nTCGA event counts:\n")
print(table(TCGA$Status, useNA = "ifany"))

cat("\nCGGA grade counts:\n")
print(table(CGGA$Grade, useNA = "ifany"))
cat("\nCGGA event counts used in Surv():\n")
print(table(CGGA$Event, useNA = "ifany"))


# -----------------------------
# 5. Define cohort-level medians
# -----------------------------

tcga_cohort_median <- median(TCGA$profile_values, na.rm = TRUE)
cgga_cohort_median <- median(CGGA$profile_values, na.rm = TRUE)

cat("\nTCGA overall CD96 median:", tcga_cohort_median, "\n")
cat("CGGA overall CD96 median:", cgga_cohort_median, "\n")


# -----------------------------
# 6. Kaplan-Meier plotting function
# -----------------------------

make_km_panel <- function(
    data,
    cohort,
    time_var,
    event_var,
    panel_label,
    panel_title,
    fixed_cutoff = NULL
) {

  needed <- c(time_var, event_var, "profile_values")
  dat <- data[complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]

  dat <- dat[
    is.finite(dat[[time_var]]) &
      is.finite(dat$profile_values) &
      dat[[time_var]] >= 0 &
      dat[[event_var]] %in% c(0, 1),
    , drop = FALSE
  ]

  if (nrow(dat) < 2) {
    stop(cohort, " ", panel_title, ": fewer than two complete observations.")
  }

  cutoff <- if (use_subgroup_specific_median || is.null(fixed_cutoff)) {
    median(dat$profile_values, na.rm = TRUE)
  } else {
    fixed_cutoff
  }

  # With an odd sample size, <= median is assigned to Low.
  dat$CD96_group <- ifelse(
    dat$profile_values > cutoff,
    "High expression",
    "Low expression"
  )

  dat$CD96_group <- factor(
    dat$CD96_group,
    levels = c("High expression", "Low expression")
  )

  group_counts <- table(dat$CD96_group)

  if (any(group_counts == 0)) {
    stop(cohort, " ", panel_title, ": one expression group has zero samples.")
  }

  # Copy the selected outcome columns to fixed internal names. This allows
  # survfit() to store a literal formula in its call. Some versions of
  # survminer fail when the stored formula is only the symbol
  # "survival_formula" and report:
  #   Error in x$formula: object of type 'symbol' is not subsettable
  dat$surv_time <- dat[[time_var]]
  dat$surv_event <- dat[[event_var]]

  survival_formula <- Surv(surv_time, surv_event) ~ CD96_group

  km_fit <- survfit(
    Surv(surv_time, surv_event) ~ CD96_group,
    data = dat
  )

  logrank_fit <- survdiff(
    Surv(surv_time, surv_event) ~ CD96_group,
    data = dat
  )
  logrank_p <- pchisq(logrank_fit$chisq, df = 1, lower.tail = FALSE)

  cox_fit <- coxph(
    Surv(surv_time, surv_event) ~ CD96_group,
    data = dat
  )
  cox_summary <- summary(cox_fit)

  # HR compares Low expression with High expression because High is reference.
  hr <- unname(cox_summary$conf.int[1, "exp(coef)"])
  hr_lower <- unname(cox_summary$conf.int[1, "lower .95"])
  hr_upper <- unname(cox_summary$conf.int[1, "upper .95"])
  cox_p <- unname(cox_summary$coefficients[1, "Pr(>|z|)"])

  maximum_time <- max(dat[[time_var]], na.rm = TRUE)
  time_break <- max(100, signif(maximum_time / 4, 1))

  full_title <- paste0(panel_label, "  ", panel_title)

  km_plot <- ggsurvplot(
    km_fit,
    data = dat,
    pval = TRUE,
    pval.method = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.27,
    risk.table.col = "strata",
    risk.table.y.text = FALSE,
    surv.median.line = "hv",
    conf.int = FALSE,
    censor = TRUE,
    censor.shape = 3,
    censor.size = 1.8,
    palette = c(high_color, low_color),
    legend.title = "Strata",
    legend.labs = c(
      paste0("High expression (", unname(group_counts["High expression"]), ")"),
      paste0("Low expression (", unname(group_counts["Low expression"]), ")")
    ),
    title = full_title,
    xlab = "Time (days)",
    ylab = "Survival probability",
    break.time.by = time_break,
    ggtheme = theme_bw(base_size = 9),
    tables.theme = theme_cleantable(base_size = 8),
    font.title = c(10, "bold"),
    font.x = c(9),
    font.y = c(9),
    font.tickslab = c(8),
    font.legend = c(7)
  )

  km_plot$plot <- km_plot$plot +
    theme(
      plot.title = element_text(hjust = 0, size = 10, face = "bold"),
      legend.position = "top",
      legend.justification = "left",
      legend.margin = margin(0, 0, 0, 0),
      panel.grid.minor = element_blank()
    )

  km_plot$table <- km_plot$table +
    theme(
      axis.title.x = element_blank(),
      plot.margin = margin(0, 5, 2, 5)
    )

  statistics <- data.frame(
    Panel = panel_label,
    Cohort = cohort,
    Analysis = panel_title,
    N = nrow(dat),
    Events = sum(dat[[event_var]] == 1),
    CD96_median_cutoff = cutoff,
    High_N = unname(group_counts["High expression"]),
    Low_N = unname(group_counts["Low expression"]),
    Logrank_P = logrank_p,
    HR_Low_vs_High = hr,
    HR_lower_95CI = hr_lower,
    HR_upper_95CI = hr_upper,
    Cox_P = cox_p,
    stringsAsFactors = FALSE
  )

  sample_counts <- data.frame(
    Panel = panel_label,
    Cohort = cohort,
    Analysis = panel_title,
    CD96_group = names(group_counts),
    N = as.integer(group_counts),
    stringsAsFactors = FALSE
  )

  cat("\n", full_title, "\n", sep = "")
  cat("N =", nrow(dat), "; events =", sum(dat[[event_var]] == 1), "\n")
  cat("CD96 cutoff =", cutoff, "\n")
  print(group_counts)
  cat("Log-rank P =", format.pval(logrank_p, digits = 4), "\n")

  list(
    plot = km_plot,
    statistics = statistics,
    sample_counts = sample_counts,
    data = dat,
    fit = km_fit
  )
}


# -----------------------------
# 7. Create eight Figure 7 panels
# -----------------------------

A <- make_km_panel(
  data = CGGA,
  cohort = "CGGA",
  time_var = "OS",
  event_var = "Event",
  panel_label = "A",
  panel_title = "All WHO Grade Survival (Glioma)",
  fixed_cutoff = cgga_cohort_median
)

B <- make_km_panel(
  data = subset(CGGA, Grade == "WHOIV"),
  cohort = "CGGA",
  time_var = "OS",
  event_var = "Event",
  panel_label = "B",
  panel_title = "WHO Grade IV Survival (GBM)",
  fixed_cutoff = cgga_cohort_median
)

C <- make_km_panel(
  data = subset(CGGA, Grade == "WHOII"),
  cohort = "CGGA",
  time_var = "OS",
  event_var = "Event",
  panel_label = "C",
  panel_title = "WHO Grade II Survival",
  fixed_cutoff = cgga_cohort_median
)

D <- make_km_panel(
  data = subset(CGGA, Grade == "WHOIII"),
  cohort = "CGGA",
  time_var = "OS",
  event_var = "Event",
  panel_label = "D",
  panel_title = "WHO Grade III Survival",
  fixed_cutoff = cgga_cohort_median
)

E <- make_km_panel(
  data = TCGA,
  cohort = "TCGA",
  time_var = "OS",
  event_var = "Status",
  panel_label = "E",
  panel_title = "All WHO Grade Survival (Glioma)",
  fixed_cutoff = tcga_cohort_median
)

F <- make_km_panel(
  data = subset(TCGA, Grade == "WHOIV"),
  cohort = "TCGA",
  time_var = "OS",
  event_var = "Status",
  panel_label = "F",
  panel_title = "WHO Grade IV Survival (GBM)",
  fixed_cutoff = tcga_cohort_median
)

G <- make_km_panel(
  data = subset(TCGA, Grade == "WHOII"),
  cohort = "TCGA",
  time_var = "OS",
  event_var = "Status",
  panel_label = "G",
  panel_title = "WHO Grade II Survival",
  fixed_cutoff = tcga_cohort_median
)

H <- make_km_panel(
  data = subset(TCGA, Grade == "WHOIII"),
  cohort = "TCGA",
  time_var = "OS",
  event_var = "Status",
  panel_label = "H",
  panel_title = "WHO Grade III Survival",
  fixed_cutoff = tcga_cohort_median
)


# -----------------------------
# 8. Combine panels like published Figure 7
# -----------------------------

# Layout:
#   A  B  E  F
#   C  D  G  H
figure7 <- arrange_ggsurvplots(
  list(
    A$plot, B$plot, E$plot, F$plot,
    C$plot, D$plot, G$plot, H$plot
  ),
  print = FALSE,
  ncol = 4,
  nrow = 2,
  risk.table.height = 0.27
)

print(figure7)


# -----------------------------
# 9. Export figure and statistics
# -----------------------------

pdf_file <- file.path(output_dir, "Figure7_CD96_Kaplan_Meier.pdf")
png_file <- file.path(output_dir, "Figure7_CD96_Kaplan_Meier.png")

pdf(
  pdf_file,
  width = 18,
  height = 10,
  onefile = TRUE,
  useDingbats = FALSE
)
print(figure7)
dev.off()

png(
  png_file,
  width = 6000,
  height = 3400,
  res = 300,
  type = ifelse(capabilities("aqua"), "quartz", "cairo")
)
print(figure7)
dev.off()

all_statistics <- do.call(
  rbind,
  lapply(list(A, B, C, D, E, F, G, H), function(x) x$statistics)
)

all_sample_counts <- do.call(
  rbind,
  lapply(list(A, B, C, D, E, F, G, H), function(x) x$sample_counts)
)

write.csv(
  all_statistics,
  file.path(output_dir, "Figure7_KM_statistics.csv"),
  row.names = FALSE
)

write.csv(
  all_sample_counts,
  file.path(output_dir, "Figure7_sample_counts.csv"),
  row.names = FALSE
)

# Save cleaned data with the event variables actually used.
write.csv(
  TCGA,
  file.path(output_dir, "TCGA_cleaned_for_survival.csv"),
  row.names = FALSE
)

write.csv(
  CGGA,
  file.path(output_dir, "CGGA_cleaned_for_survival.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("Figure 7 analysis completed successfully.\n")
cat("Output directory:\n", output_dir, "\n", sep = "")
cat("\nCreated files:\n")
cat(" - ", pdf_file, "\n", sep = "")
cat(" - ", png_file, "\n", sep = "")
cat(" - Figure7_KM_statistics.csv\n")
cat(" - Figure7_sample_counts.csv\n")
cat(" - TCGA_cleaned_for_survival.csv\n")
cat(" - CGGA_cleaned_for_survival.csv\n")
cat("============================================================\n")

sessionInfo()
