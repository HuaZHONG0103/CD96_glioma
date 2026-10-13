## =========================================================
## Figure 2. ROC analysis of CD96 in glioma
## Article-style version
##
## Input:
##   CGGA_CD96.txt
##   TCGA_CD96.txt
##
## Output:
##   Fig2_CD96_ROC_article_style.pdf
##   Fig2_CD96_ROC_AUC_summary.tsv
## =========================================================

library(pROC)

## -------------------------
## 1. Read data
## -------------------------
read_cd96_data <- function(file, cohort_name) {
  dat <- read.table(file, sep = "\t", header = TRUE, check.names = FALSE, quote = "")
  names(dat) <- trimws(names(dat))
  
  ## unify expression column
  if ("profile_values" %in% names(dat)) {
    dat$CD96_expr <- suppressWarnings(as.numeric(dat$profile_values))
  } else if ("CD96" %in% names(dat)) {
    dat$CD96_expr <- suppressWarnings(as.numeric(dat$CD96))
  } else {
    stop(paste0("No expression column found in ", file,
                ". Need either 'profile_values' or 'CD96'."))
  }
  
  dat$Cohort <- cohort_name
  
  ## clean common missing text
  dat[dat == ""] <- NA
  dat[dat == "NA"] <- NA
  dat[dat == "N/A"] <- NA
  
  return(dat)
}

cgga <- read_cd96_data("CGGA_CD96.txt", "CGGA")
tcga <- read_cd96_data("TCGA_CD96.txt", "TCGA")

## -------------------------
## 2. Helper functions
## -------------------------
fmt_pct1 <- function(x) sprintf("%.1f%%", x)
fmt_auc3 <- function(x) sprintf("%.3f", x)

## A/C panels: Mesenchymal subtype ROC with CI band
plot_mes_roc <- function(dat, panel_letter = "A") {
  d <- dat[!is.na(dat$subtype) & !is.na(dat$CD96_expr), ]
  
  ## keep standard TCGA-like subtypes
  d <- d[d$subtype %in% c("Classical", "Mesenchymal", "Neural", "Proneural"), ]
  
  ## positive class = Mesenchymal
  d$response <- ifelse(d$subtype == "Mesenchymal", 1, 0)
  
  roc_obj <- roc(
    response = d$response,
    predictor = d$CD96_expr,
    levels = c(0, 1),
    direction = "<",
    percent = TRUE,
    quiet = TRUE
  )
  
  ci_auc <- ci.auc(roc_obj)
  ci_band <- ci.se(roc_obj, specificities = seq(0, 100, 1))
  
  ## empty plot
  plot(
    NA,
    xlim = c(100, 0),
    ylim = c(0, 100),
    xlab = "Specificity (%)",
    ylab = "Sensitivity (%)",
    main = "",
    type = "n",
    axes = TRUE
  )
  
  ## CI polygon
  plot(ci_band,
       type = "shape",
       col = rgb(52, 93, 167, maxColorValue = 255, alpha = 90),
       border = NA,
       xlim = c(100, 0),
       ylim = c(0, 100),
       add = TRUE)
  
  ## diagonal
  segments(100, 0, 0, 100, col = "grey75", lwd = 1)
  
  ## ROC line
  plot(
    roc_obj,
    add = TRUE,
    legacy.axes = FALSE,
    percent = TRUE,
    col = "black",
    lwd = 1.2,
    print.auc = FALSE,
    xlab = "",
    ylab = ""
  )
  
  ## annotations
  usr <- par("usr")
  
  text(usr[1] - 1.5, usr[4] + 8, panel_letter, font = 2, xpd = NA)
  
  text(25, 90, "Confidence intervals", cex = 0.75)
  text(32, 28, "Mes", cex = 2)
  
  auc_txt <- paste0(
    "AUC: ",
    fmt_pct1(as.numeric(auc(roc_obj))),
    " (",
    fmt_pct1(as.numeric(ci_auc[1])),
    "...",
    fmt_pct1(as.numeric(ci_auc[3])),
    ")"
  )
  text(50, 20, auc_txt, font = 2, cex = 0.9)
  
  invisible(data.frame(
    analysis = paste0(unique(dat$Cohort), "_Mesenchymal"),
    n = nrow(d),
    n_positive = sum(d$response == 1),
    n_negative = sum(d$response == 0),
    auc = as.numeric(auc(roc_obj)),
    ci_low = as.numeric(ci_auc[1]),
    ci_high = as.numeric(ci_auc[3])
  ))
}

## B/D panels: WHO II / III / IV ROC overlay
plot_grade_rocs <- function(dat, panel_letter = "B") {
  d <- dat[!is.na(dat$Grade) & !is.na(dat$CD96_expr), ]
  d <- d[d$Grade %in% c("WHO II", "WHO III", "WHO IV"), ]
  
  ## Use fixed grade order
  grades <- c("WHO II", "WHO III", "WHO IV")
  colors <- c("red", "goldenrod2", "dodgerblue")
  
  roc_list <- list()
  auc_vals <- numeric(length(grades))
  
  for (i in seq_along(grades)) {
    g <- grades[i]
    d$response <- ifelse(d$Grade == g, 1, 0)
    
    roc_list[[g]] <- roc(
      response = d$response,
      predictor = d$CD96_expr,
      levels = c(0, 1),
      direction = "auto",
      quiet = TRUE
    )
    auc_vals[i] <- as.numeric(auc(roc_list[[g]]))
  }
  
  ## base plot using first curve
  plot(
    roc_list[[grades[1]]],
    legacy.axes = TRUE,
    col = colors[1],
    lwd = 1.2,
    main = "",
    xlab = "1 - Specificity",
    ylab = "Sensitivity",
    print.auc = FALSE
  )
  
  abline(0, 1, col = "grey75", lwd = 1)
  
  if (length(grades) > 1) {
    for (i in 2:length(grades)) {
      plot(roc_list[[grades[i]]], add = TRUE, col = colors[i], lwd = 1.2)
    }
  }
  
  usr <- par("usr")
  text(usr[1] - 0.10, usr[4] + 0.12, panel_letter, font = 2, xpd = NA)
  
  legend(
    "bottomright",
    legend = c(
      paste0("WHO II AUC ", fmt_auc3(auc_vals[1])),
      paste0("WHO III AUC ", fmt_auc3(auc_vals[2])),
      paste0("WHO IV AUC ", fmt_auc3(auc_vals[3]))
    ),
    col = colors,
    lwd = 1.2,
    bty = "n",
    cex = 0.85
  )
  
  out <- data.frame(
    analysis = paste0(unique(dat$Cohort), "_", grades),
    n = nrow(d),
    n_positive = c(
      sum(d$Grade == "WHO II"),
      sum(d$Grade == "WHO III"),
      sum(d$Grade == "WHO IV")
    ),
    n_negative = c(
      sum(d$Grade != "WHO II"),
      sum(d$Grade != "WHO III"),
      sum(d$Grade != "WHO IV")
    ),
    auc = auc_vals,
    ci_low = NA,
    ci_high = NA
  )
  
  invisible(out)
}

## -------------------------
## 3. Draw article-style figure
## -------------------------
pdf("Fig2_CD96_ROC_article_style.pdf", width = 8.6, height = 9.2)

layout(
  matrix(c(
    1, 1,
    2, 3,
    4, 4,
    5, 6
  ), nrow = 4, byrow = TRUE),
  heights = c(0.12, 1, 0.12, 1)
)

## ---- Title row 1 ----
par(mar = c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, "ROC Curve in CGGA Dataset", cex = 1.5)

## ---- Panel A ----
par(mar = c(4, 4, 1.5, 1.5))
res_A <- plot_mes_roc(cgga, "A")

## ---- Panel B ----
par(mar = c(4, 4, 1.5, 1.5))
res_B <- plot_grade_rocs(cgga, "B")

## ---- Title row 2 ----
par(mar = c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, "ROC Curve in TCGA Dataset", cex = 1.5)

## ---- Panel C ----
par(mar = c(4, 4, 1.5, 1.5))
res_C <- plot_mes_roc(tcga, "C")

## ---- Panel D ----
par(mar = c(4, 4, 1.5, 1.5))
res_D <- plot_grade_rocs(tcga, "D")

dev.off()

## -------------------------
## 4. Save summary table
## -------------------------
auc_summary <- rbind(
  res_A,
  res_B,
  res_C,
  res_D
)

write.table(
  auc_summary,
  "Fig2_CD96_ROC_AUC_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Done:\n")
cat("  Fig2_CD96_ROC_article_style.pdf\n")
cat("  Fig2_CD96_ROC_AUC_summary.tsv\n")

