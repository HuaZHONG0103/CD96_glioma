#!/usr/bin/env Rscript

# ============================================================
# Reproduce CD96 glioma Figure 3-style GO plots
# Fixed version v6
#
# Based on the previous v5 style, with two requested changes:
#   1. Remove GOCluster output completely.
#   2. Keep GOplot::GOBubble, but make it less crowded by:
#        - plotting only the top GO terms
#        - hiding the large table legend
#        - showing fewer GO ID labels
#
# Required input files in working directory:
#   CGGA_top50.csv
#   TCGA_top50.csv
#
# Required columns:
#   SYMBOL, log2fc
# ============================================================

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

# ---------------------------
# 0. Package setup
# ---------------------------
required_pkgs <- c(
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi",
  "enrichplot", "GOplot", "ggplot2", "dplyr", "readr",
  "stringr", "tibble", "forcats", "tidyr", "grid"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  message("Missing packages: ", paste(missing_pkgs, collapse = ", "))
  message("Install required packages first, for example:")
  message("if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')")
  message("BiocManager::install(c('clusterProfiler','org.Hs.eg.db','enrichplot','GOplot'))")
  message("install.packages(c('ggplot2','dplyr','readr','stringr','tibble','forcats','tidyr'))")
  stop("Please install missing packages and rerun.")
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(GOplot)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(forcats)
  library(tidyr)
  library(grid)
})

# ---------------------------
# 1. User settings
# ---------------------------
input_dir  <- "."
out_dir    <- "CD96_GO_reproduce_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ont_use <- "BP"
p_cut   <- 0.05
q_cut   <- 0.20
simplify_cutoff <- 0.70

show_terms_fig3       <- 6
show_terms_bubble     <- 12
show_terms_gobubble   <- 12   # keep the original-style table, but only show top 12 terms
gobubble_label_number <- 4    # label only top 4 bubbles to reduce overlap

make_GOplot_GOBubble <- TRUE
make_clusterProfiler_dotplot <- TRUE
make_clusterProfiler_barplot <- TRUE

# ---------------------------
# 2. Helper functions
# ---------------------------
read_gene_fc <- function(file) {
  if (!file.exists(file)) {
    stop("Input file not found: ", file)
  }

  x <- readr::read_csv(file, show_col_types = FALSE)

  if (!all(c("SYMBOL", "log2fc") %in% colnames(x))) {
    stop("Input file must contain columns: SYMBOL and log2fc. File: ", file)
  }

  x %>%
    dplyr::mutate(
      SYMBOL = as.character(SYMBOL),
      log2fc = as.numeric(log2fc)
    ) %>%
    dplyr::filter(!is.na(SYMBOL), !is.na(log2fc)) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)
}

map_symbol_to_entrez <- function(gene_fc, prefix = "dataset") {
  id_map <- suppressMessages(
    clusterProfiler::bitr(
      gene_fc$SYMBOL,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db::org.Hs.eg.db
    )
  )

  unmapped <- setdiff(gene_fc$SYMBOL, id_map$SYMBOL)
  if (length(unmapped) > 0) {
    readr::write_lines(unmapped, file.path(out_dir, paste0(prefix, "_unmapped_symbols.txt")))
    message(prefix, ": unmapped symbols saved to ", file.path(out_dir, paste0(prefix, "_unmapped_symbols.txt")))
  }

  gene_fc %>%
    dplyr::left_join(id_map, by = "SYMBOL") %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
    dplyr::select(SYMBOL, ENTREZID, log2fc)
}

safe_simplify_go <- function(ego, prefix) {
  message(prefix, ": simplifying enriched GO terms with clusterProfiler::simplify()")

  ego_simplified <- tryCatch({
    clusterProfiler::simplify(
      ego,
      cutoff     = simplify_cutoff,
      by         = "p.adjust",
      select_fun = min
    )
  }, error = function(e) {
    message(prefix, ": clusterProfiler::simplify() failed.")
    message("Reason: ", e$message)
    message(prefix, ": using unsimplified enrichGO results instead.")
    ego
  })

  ego_simplified
}

run_go_enrichment <- function(entrez_fc, prefix) {
  message(prefix, ": running GO BP enrichment")

  ego <- clusterProfiler::enrichGO(
    gene          = entrez_fc$ENTREZID,
    OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = ont_use,
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cut,
    qvalueCutoff  = q_cut,
    readable      = FALSE
  )

  ego_df <- as.data.frame(ego)
  readr::write_csv(ego_df, file.path(out_dir, paste0(prefix, "_enrichGO_BP.csv")))

  if (nrow(ego_df) == 0) {
    stop(prefix, ": no significant GO BP terms found. Try increasing p_cut/q_cut or checking input genes.")
  }

  ego_simplified <- safe_simplify_go(ego, prefix)

  readr::write_csv(
    as.data.frame(ego_simplified),
    file.path(out_dir, paste0(prefix, "_enrichGO_simplify_BP.csv"))
  )

  list(ego = ego, ego_simplified = ego_simplified)
}

make_gene_list <- function(entrez_fc) {
  gene_list <- entrez_fc$log2fc
  names(gene_list) <- entrez_fc$ENTREZID
  gene_list
}

make_goplot_data <- function(ego_df, entrez_fc) {
  go <- data.frame(
    Category = "BP",
    ID       = ego_df$ID,
    Term     = ego_df$Description,
    Genes    = gsub("/", ", ", ego_df$geneID),
    adj_pval = ego_df$p.adjust,
    stringsAsFactors = FALSE
  )

  genelist <- data.frame(
    ID    = as.character(entrez_fc$ENTREZID),
    logFC = entrez_fc$log2fc,
    stringsAsFactors = FALSE
  )

  circ <- GOplot::circle_dat(go, genelist)

  id_gsym <- suppressMessages(
    clusterProfiler::bitr(
      unique(circ$genes),
      fromType = "ENTREZID",
      toType   = "SYMBOL",
      OrgDb    = org.Hs.eg.db::org.Hs.eg.db
    )
  )
  id_gsym <- id_gsym[!duplicated(id_gsym$ENTREZID), ]
  rownames(id_gsym) <- id_gsym$ENTREZID

  circ_symbol <- circ
  circ_symbol$genes <- ifelse(
    circ$genes %in% rownames(id_gsym),
    id_gsym[circ$genes, "SYMBOL"],
    circ$genes
  )

  list(go = go, genelist = genelist, circ = circ, circ_symbol = circ_symbol)
}

safe_pdf <- function(filename, width = 8, height = 6, plot_expr) {
  grDevices::pdf(filename, width = width, height = height, useDingbats = FALSE)
  tryCatch({
    force(plot_expr)
  }, error = function(e) {
    message("Failed while writing ", filename, ": ", e$message)
  }, finally = {
    grDevices::dev.off()
  })
}

safe_ggsave <- function(filename, plot, width = 9, height = 6) {
  ggplot2::ggsave(filename, plot = plot, width = width, height = height, device = cairo_pdf)
}

safe_barplot_enrich <- function(egox, showCategory) {
  # clusterProfiler::barplot is available in many older versions.
  # enrichplot::barplot may not be exported in older Bioconductor versions.
  if ("barplot" %in% getNamespaceExports("clusterProfiler")) {
    return(clusterProfiler::barplot(egox, showCategory = showCategory))
  }

  return(barplot(egox, showCategory = showCategory))
}

# ---------------------------
# 3. clusterProfiler / enrichplot figures
# ---------------------------
plot_clusterprofiler <- function(ego_obj, gene_list, prefix) {
  egox <- clusterProfiler::setReadable(
    ego_obj,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID"
  )

  if (make_clusterProfiler_dotplot) {
    tryCatch({
      p_dot <- enrichplot::dotplot(egox, showCategory = show_terms_bubble) +
        ggplot2::ggtitle(paste0(prefix, " GO BP dot plot")) +
        ggplot2::theme_bw(base_size = 12)
      safe_ggsave(file.path(out_dir, paste0(prefix, "_clusterProfiler_dotplot.pdf")), p_dot, 9, 6)
    }, error = function(e) message(prefix, ": dotplot failed: ", e$message))
  }

  if (make_clusterProfiler_barplot) {
    tryCatch({
      p_bar <- safe_barplot_enrich(egox, showCategory = show_terms_bubble) +
        ggplot2::ggtitle(paste0(prefix, " GO BP bar plot")) +
        ggplot2::theme_bw(base_size = 12)
      safe_ggsave(file.path(out_dir, paste0(prefix, "_clusterProfiler_barplot.pdf")), p_bar, 9, 6)
    }, error = function(e) message(prefix, ": barplot failed: ", e$message))
  }

  tryCatch({
    p_cnet_circle <- enrichplot::cnetplot(
      egox,
      foldChange   = gene_list,
      circular     = TRUE,
      colorEdge    = TRUE,
      showCategory = show_terms_fig3
    ) + ggplot2::ggtitle(paste0(prefix, " circular gene-concept network"))
    safe_ggsave(file.path(out_dir, paste0(prefix, "_clusterProfiler_cnet_circular.pdf")), p_cnet_circle, 10, 8)
  }, error = function(e) message(prefix, ": circular cnetplot failed: ", e$message))

  tryCatch({
    p_cnet <- enrichplot::cnetplot(
      egox,
      foldChange   = gene_list,
      circular     = FALSE,
      colorEdge    = TRUE,
      showCategory = show_terms_bubble
    ) + ggplot2::ggtitle(paste0(prefix, " gene-concept network"))
    safe_ggsave(file.path(out_dir, paste0(prefix, "_clusterProfiler_cnet_network.pdf")), p_cnet, 10, 8)
  }, error = function(e) message(prefix, ": cnetplot failed: ", e$message))
}

# ---------------------------
# 4. Clean ggplot bubble and lollipop
# ---------------------------
plot_custom_gg <- function(ego_df, prefix) {
  df <- ego_df %>%
    dplyr::mutate(
      neglog10 = -log10(p.adjust),
      GeneRatioNum = as.numeric(sub("/.*", "", GeneRatio)) /
        as.numeric(sub(".*/", "", GeneRatio)),
      Description_wrapped = stringr::str_wrap(Description, 45)
    ) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = show_terms_bubble)

  p_bubble <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = GeneRatioNum, y = forcats::fct_reorder(Description_wrapped, GeneRatioNum))
  ) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = neglog10), alpha = 0.85) +
    ggplot2::scale_color_viridis_c(option = "C") +
    ggplot2::labs(
      x = "Gene ratio", y = NULL,
      color = "-log10(adj. p)", size = "Gene count",
      title = paste0(prefix, " GO BP bubble plot")
    ) +
    ggplot2::theme_bw(base_size = 12)
  safe_ggsave(file.path(out_dir, paste0(prefix, "_ggplot_bubble.pdf")), p_bubble, 9, 6)

  p_lollipop <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = GeneRatioNum, y = forcats::fct_reorder(Description_wrapped, GeneRatioNum))
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = GeneRatioNum, yend = Description_wrapped),
      linewidth = 0.6, color = "grey60"
    ) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = neglog10), alpha = 0.9) +
    ggplot2::scale_color_viridis_c(option = "C") +
    ggplot2::labs(
      x = "Gene ratio", y = NULL,
      color = "-log10(adj. p)", size = "Gene count",
      title = paste0(prefix, " GO BP lollipop plot")
    ) +
    ggplot2::theme_bw(base_size = 12)
  safe_ggsave(file.path(out_dir, paste0(prefix, "_ggplot_lollipop.pdf")), p_lollipop, 9, 6)
}

# ---------------------------
# 5. GOplot-style figures
# ---------------------------
plot_gobubble_clean <- function(goplot_data, prefix) {
  if (!make_GOplot_GOBubble) {
    return(invisible(NULL))
  }

  circ <- goplot_data$circ
  go <- goplot_data$go

  # Use top enriched terms only, but keep the original GOplot style:
  # left bubble plot + right ID/Description table.
  n_terms <- min(show_terms_gobubble, nrow(go))
  selected_terms <- go$Term[seq_len(n_terms)]
  circ_small <- circ[circ$term %in% selected_terms, , drop = FALSE]

  # This version is closer to the previous/original style:
  #   display = "single" removes the grey BP facet strip
  #   table.legend = TRUE brings back the green ID/Description table
  #   labels = 4 labels only the top bubbles, avoiding label crowding
  #   wider PDF gives the table enough room
  tryCatch({
    safe_pdf(file.path(out_dir, paste0(prefix, "_GOplot_GOBubble_previous_style.pdf")), width = 14.5, height = 7.0, {
      print(GOplot::GOBubble(
        circ_small,
        title = paste0(prefix, " GO BP bubble plot"),
        display = "single",
        labels = gobubble_label_number,
        ID = TRUE,
        table.legend = TRUE
      ))
    })

    # Also save using the old expected filename.
    safe_pdf(file.path(out_dir, paste0(prefix, "_GOplot_GOBubble_clean.pdf")), width = 14.5, height = 7.0, {
      print(GOplot::GOBubble(
        circ_small,
        title = paste0(prefix, " GO BP bubble plot"),
        display = "single",
        labels = gobubble_label_number,
        ID = TRUE,
        table.legend = TRUE
      ))
    })
  }, error = function(e) {
    message(prefix, ": GOplot GOBubble previous-style with ID=TRUE failed: ", e$message)
    message(prefix, ": retrying with ID=FALSE.")
    tryCatch({
      safe_pdf(file.path(out_dir, paste0(prefix, "_GOplot_GOBubble_previous_style.pdf")), width = 14.5, height = 7.0, {
        print(GOplot::GOBubble(
          circ_small,
          title = paste0(prefix, " GO BP bubble plot"),
          display = "single",
          labels = gobubble_label_number,
          ID = FALSE,
          table.legend = TRUE
        ))
      })
      safe_pdf(file.path(out_dir, paste0(prefix, "_GOplot_GOBubble_clean.pdf")), width = 14.5, height = 7.0, {
        print(GOplot::GOBubble(
          circ_small,
          title = paste0(prefix, " GO BP bubble plot"),
          display = "single",
          labels = gobubble_label_number,
          ID = FALSE,
          table.legend = TRUE
        ))
      })
    }, error = function(e2) {
      message(prefix, ": GOplot GOBubble failed again: ", e2$message)
    })
  })
}

plot_gocircle <- function(goplot_data, prefix) {
  circ <- goplot_data$circ

  tryCatch({
    safe_pdf(file.path(out_dir, paste0(prefix, "_GOplot_GOCircle.pdf")), width = 12, height = 12, {
      print(GOplot::GOCircle(circ))
    })
  }, error = function(e) message(prefix, ": GOCircle failed: ", e$message))
}

plot_gochord <- function(goplot_data, prefix) {
  go <- goplot_data$go
  genelist <- goplot_data$genelist
  circ <- goplot_data$circ

  n_terms <- min(show_terms_fig3, nrow(go))
  if (n_terms < 1) {
    warning(prefix, ": no GO terms available for GOChord.")
    return(invisible(NULL))
  }

  selected_terms <- go$Term[seq_len(n_terms)]

  tryCatch({
    chord <- GOplot::chord_dat(circ, genelist, selected_terms)

    id_gsym <- suppressMessages(
      clusterProfiler::bitr(
        row.names(chord),
        fromType = "ENTREZID",
        toType   = "SYMBOL",
        OrgDb    = org.Hs.eg.db::org.Hs.eg.db
      )
    )
    id_gsym <- id_gsym[!duplicated(id_gsym$ENTREZID), ]
    rownames(id_gsym) <- id_gsym$ENTREZID

    chord_symbol <- chord
    row.names(chord_symbol) <- ifelse(
      row.names(chord) %in% rownames(id_gsym),
      id_gsym[row.names(chord), "SYMBOL"],
      row.names(chord)
    )

    safe_pdf(file.path(out_dir, paste0(prefix, "_GOplot_GOChord_Fig3_style.pdf")), width = 12, height = 13, {
      print(GOplot::GOChord(
        chord_symbol,
        space         = 0.02,
        gene.order    = "logFC",
        lfc.col       = c("darkgoldenrod1", "black", "cyan1"),
        gene.space    = 0.25,
        gene.size     = 5,
        border.size   = 0.1,
        process.label = 8
      ))
    })
  }, error = function(e) message(prefix, ": GOChord failed: ", e$message))
}

plot_goplot_family <- function(goplot_data, prefix) {
  plot_gobubble_clean(goplot_data, prefix)
  plot_gocircle(goplot_data, prefix)
  plot_gochord(goplot_data, prefix)
  # GOCluster intentionally removed per user request.
}

# ---------------------------
# 6. Main function
# ---------------------------
run_one_dataset <- function(prefix, top50_file) {
  message("\n============================")
  message("Processing ", prefix)
  message("============================")

  gene_fc <- read_gene_fc(top50_file)
  entrez_fc <- map_symbol_to_entrez(gene_fc, prefix)

  if (nrow(entrez_fc) < 5) {
    stop(prefix, ": too few mapped ENTREZ genes. Check gene symbols in input file.")
  }

  readr::write_csv(entrez_fc, file.path(out_dir, paste0(prefix, "_very_easy_input_ENTREZ_logFC.csv")))

  enrich_res <- run_go_enrichment(entrez_fc, prefix)
  ego2 <- enrich_res$ego_simplified
  ego2_df <- as.data.frame(ego2)

  if (nrow(ego2_df) == 0) {
    warning(prefix, ": no enriched GO terms after simplification.")
    return(invisible(NULL))
  }

  gene_list <- make_gene_list(entrez_fc)

  message(prefix, ": drawing clusterProfiler/enrichplot figures")
  plot_clusterprofiler(ego2, gene_list, prefix)

  message(prefix, ": drawing ggplot bubble/lollipop figures")
  plot_custom_gg(ego2_df, prefix)

  message(prefix, ": drawing GOplot figures")
  goplot_data <- make_goplot_data(ego2_df, entrez_fc)
  plot_goplot_family(goplot_data, prefix)

  message(prefix, " finished. Outputs saved to: ", normalizePath(out_dir))
  invisible(TRUE)
}

# ---------------------------
# 7. Run CGGA and TCGA
# ---------------------------
run_one_dataset("CGGA", file.path(input_dir, "CGGA_top50.csv"))
run_one_dataset("TCGA", file.path(input_dir, "TCGA_top50.csv"))

message("\nAll done.")
message("Main output files:")
message("  *_GOplot_GOChord_Fig3_style.pdf       # Figure 3-style chord plot")
message("  *_GOplot_GOBubble_previous_style.pdf  # original-style GOplot bubble with table
  *_GOplot_GOBubble_clean.pdf           # same plot, saved with previous filename")
message("  *_ggplot_bubble.pdf                   # clean ggplot bubble")
message("  *_clusterProfiler_dotplot.pdf         # dot plot")
message("  *_ggplot_lollipop.pdf                 # lollipop plot")
message("  *_GOplot_GOCircle.pdf                 # GO circle plot")
message("GOCluster output is intentionally removed.")
