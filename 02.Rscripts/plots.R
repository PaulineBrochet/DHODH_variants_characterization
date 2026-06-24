# ============================================================
# Generic plotting helpers
# ============================================================
#
# This file contains reusable plotting functions that are not tied to a
# specific omics layer. The goal is to keep QC and differential-analysis
# plotting logic in one place, so transcriptomics, proteomics, metabolomics,
# and multi-omics scripts can call the same helpers.
#
# Design principles:
# - Functions return plot objects when possible.
# - Optional output_file arguments save plots without forcing file output.
# - Defaults preserve the visual style used in the legacy manuscript scripts.

project_library_file <- if (file.exists("02.Rscripts/library.R")) {
  "02.Rscripts/library.R"
} else if (file.exists("library.R")) {
  "library.R"
} else {
  NA_character_
}
if (!is.na(project_library_file)) {
  source(project_library_file)
  load_project_libraries()
}
rm(project_library_file)

# ------------------------------------------------------------
# PCA plot
# ------------------------------------------------------------
# Build a PCA scatter plot from a feature-by-sample matrix.
#
# Args:
#   mat: numeric matrix with features in rows and samples in columns.
#   metadata: data frame containing sample annotations.
#   sample_col: column in metadata matching colnames(mat).
#   color_col: optional metadata column used for point color.
#   shape_col: optional metadata column used for point shape.
#   label_col: optional metadata column used for point labels.
#   title: plot title.
#   pc_x, pc_y: principal components to display.
#   scale_features: passed to prcomp(scale. = ...).
#   center_features: passed to prcomp(center = ...).
#   output_file: optional path where the plot should be saved.
#
# Returns:
#   A ggplot object.
plot_pca_qc <- function(mat,
                        metadata,
                        sample_col,
                        color_col = NULL,
                        shape_col = NULL,
                        label_col = NULL,
                        title = "PCA",
                        pc_x = 1,
                        pc_y = 2,
                        scale_features = FALSE,
                        center_features = TRUE,
                        point_size = 3,
                        label_size = 2.5,
                        max_overlaps = 60,
                        output_file = NULL,
                        width = 7,
                        height = 5,
                        dpi = 200) {
  stopifnot(is.matrix(mat) || is.data.frame(mat))
  mat <- as.matrix(mat)

  missing_samples <- setdiff(colnames(mat), metadata[[sample_col]])
  if (length(missing_samples) > 0) {
    stop("Metadata is missing samples: ", paste(missing_samples, collapse = ", "))
  }

  pca <- stats::prcomp(t(mat), center = center_features, scale. = scale_features)
  pct <- (pca$sdev^2) / sum(pca$sdev^2) * 100

  pca_df <- tibble::as_tibble(pca$x, rownames = sample_col) |>
    dplyr::left_join(metadata, by = sample_col)

  x_col <- paste0("PC", pc_x)
  y_col <- paste0("PC", pc_y)

  p <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])
  )

  if (!is.null(color_col) && !is.null(shape_col)) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(color = .data[[color_col]], shape = .data[[shape_col]]),
      size = point_size
    )
  } else if (!is.null(color_col)) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(color = .data[[color_col]]),
      size = point_size
    )
  } else if (!is.null(shape_col)) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(shape = .data[[shape_col]]),
      size = point_size
    )
  } else {
    p <- p + ggplot2::geom_point(size = point_size)
  }

  if (!is.null(label_col)) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = .data[[label_col]]),
      size = label_size,
      max.overlaps = max_overlaps
    )
  }

  p <- p +
    ggplot2::labs(
      title = title,
      x = sprintf("%s (%.1f%%)", x_col, pct[pc_x]),
      y = sprintf("%s (%.1f%%)", y_col, pct[pc_y])
    )

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = width, height = height, dpi = dpi)
  }

  p
}

# ------------------------------------------------------------
# Sample-distance heatmap
# ------------------------------------------------------------
# Plot a sample-to-sample distance heatmap from a feature-by-sample matrix.
#
# Args:
#   mat: numeric matrix with features in rows and samples in columns.
#   metadata: data frame containing sample annotations.
#   sample_col: column in metadata matching colnames(mat).
#   annotation_cols: metadata columns to display as row/column annotation.
#   output_file: optional PNG path. If provided, the function opens/closes a
#     PNG device to match the original manuscript plotting behavior.
#
# Returns:
#   The pheatmap object, invisibly.
plot_sample_distance_heatmap <- function(mat,
                                         metadata,
                                         sample_col,
                                         annotation_cols = NULL,
                                         title = "Sample distances",
                                         distance_method = "euclidean",
                                         clustering_distance_rows = "euclidean",
                                         clustering_distance_cols = "euclidean",
                                         output_file = NULL,
                                         width = 1000,
                                         height = 900) {
  stopifnot(is.matrix(mat) || is.data.frame(mat))
  mat <- as.matrix(mat)

  missing_samples <- setdiff(colnames(mat), metadata[[sample_col]])
  if (length(missing_samples) > 0) {
    stop("Metadata is missing samples: ", paste(missing_samples, collapse = ", "))
  }

  metadata <- metadata |>
    dplyr::filter(.data[[sample_col]] %in% colnames(mat)) |>
    dplyr::slice(match(colnames(mat), .data[[sample_col]])) |>
    as.data.frame()
  rownames(metadata) <- NULL

  dmat <- stats::dist(t(mat), method = distance_method) |>
    as.matrix()

  ann <- NULL
  if (!is.null(annotation_cols)) {
    ann <- metadata |>
      dplyr::select(dplyr::all_of(c(sample_col, annotation_cols))) |>
      as.data.frame()
    rownames(ann) <- NULL
    ann <- tibble::column_to_rownames(ann, sample_col)
  }

  if (!is.null(output_file)) {
    grDevices::png(output_file, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  heatmap_obj <- pheatmap::pheatmap(
    dmat,
    annotation_col = ann,
    annotation_row = ann,
    clustering_distance_rows = clustering_distance_rows,
    clustering_distance_cols = clustering_distance_cols,
    main = title
  )

  invisible(heatmap_obj)
}

# ------------------------------------------------------------
# Volcano plot
# ------------------------------------------------------------
# Build a volcano plot from a differential-analysis result table.
#
# Args:
#   tbl: data frame containing fold-change and p-value columns.
#   lfc_col: log2 fold-change column.
#   p_col: p-value or adjusted p-value column.
#   lfc_cut, p_cut: significance thresholds.
#   title: plot title.
#   output_file: optional path where the plot should be saved.
#
# Returns:
#   A ggplot object, or NULL if there are fewer than min_points valid rows.
plot_volcano_qc <- function(tbl,
                            lfc_col,
                            p_col,
                            lfc_cut,
                            p_cut,
                            title,
                            min_points = 3,
                            point_size = 1.2,
                            point_alpha = 0.7,
                            base_size = 14,
                            colors = c(
                              "NS" = "grey70",
                              "Significant" = "red2",
                              "lfc_only" = "green4",
                              "p_only" = "royalblue3"
                            ),
                            output_file = NULL,
                            width = 8,
                            height = 4,
                            dpi = 200) {
  lfc_lab <- paste0("|log2FC| >= ", lfc_cut)
  p_lab <- paste0(p_col, " <= ", p_cut)
  both_lab <- "Significant"

  tbl_plot <- tbl |>
    dplyr::mutate(
      LFC = as.numeric(.data[[lfc_col]]),
      P_val = as.numeric(.data[[p_col]])
    ) |>
    dplyr::filter(is.finite(LFC), is.finite(P_val)) |>
    dplyr::mutate(
      negLog10P = -log10(P_val + 1e-300),
      group = dplyr::case_when(
        P_val <= p_cut & abs(LFC) >= lfc_cut ~ both_lab,
        P_val <= p_cut ~ p_lab,
        abs(LFC) >= lfc_cut ~ lfc_lab,
        TRUE ~ "NS"
      )
    )

  if (nrow(tbl_plot) < min_points) {
    return(NULL)
  }

  col_vals <- c("NS" = colors[["NS"]], "Significant" = colors[["Significant"]])
  col_vals[lfc_lab] <- colors[["lfc_only"]]
  col_vals[p_lab] <- colors[["p_only"]]

  p <- ggplot2::ggplot(tbl_plot, ggplot2::aes(x = LFC, y = negLog10P)) +
    ggplot2::geom_point(
      ggplot2::aes(color = group),
      alpha = point_alpha,
      size = point_size
    ) +
    ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 2) +
    ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = 2) +
    ggplot2::scale_color_manual(values = col_vals, name = NULL) +
    ggplot2::labs(
      title = title,
      x = "log2(Fold Change)",
      y = paste0("-log10(", p_col, ")")
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.position = "top"
    )

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = width, height = height, dpi = dpi)
  }

  p
}

# ------------------------------------------------------------
# MitoCarta/glycolysis volcano plot
# ------------------------------------------------------------
# Build volcano plots restricted to MitoCarta genes and glycolysis genes.
# The selected feature set is the union of:
# - MitoCarta symbols from the project MitoCarta file,
# - KEGG/Reactome glycolysis-related gene sets from msigdbr.
#
# The plotting function is omics-agnostic: callers provide the result files and
# the column names used for gene symbols, log2 fold-change, and adjusted p-value.

read_mitocarta_symbols <- function(mitocarta_file, symbol_col = "Symbol") {
  if (is.character(mitocarta_file) && length(mitocarta_file) == 1) {
    ext <- tolower(tools::file_ext(mitocarta_file))
    mitocarta_tbl <- switch(
      ext,
      "xlsx" = openxlsx::read.xlsx(mitocarta_file),
      "xls" = openxlsx::read.xlsx(mitocarta_file),
      "csv" = readr::read_csv(mitocarta_file, show_col_types = FALSE),
      "tsv" = readr::read_tsv(mitocarta_file, show_col_types = FALSE),
      stop("Unsupported MitoCarta file extension: ", ext)
    )
  } else {
    mitocarta_tbl <- mitocarta_file
  }

  if (!(symbol_col %in% colnames(mitocarta_tbl))) {
    stop("MitoCarta table is missing the symbol column: ", symbol_col)
  }

  unique(stats::na.omit(trimws(as.character(mitocarta_tbl[[symbol_col]]))))
}

get_msigdbr_glycolysis_symbols <- function(
    species = "Homo sapiens",
    collection = "C2",
    subcollections = c("CP:KEGG_LEGACY", "CP:REACTOME"),
    gene_set_pattern = "GLYCOLYSIS|GLUCONEOGENESIS") {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("The msigdbr package is required to build the glycolysis gene set.")
  }

  gene_sets <- lapply(subcollections, function(subcollection) {
    msigdbr::msigdbr(
      species = species,
      collection = collection,
      subcollection = subcollection
    )
  }) |>
    dplyr::bind_rows() |>
    dplyr::filter(grepl(gene_set_pattern, .data$gs_name, ignore.case = TRUE))

  if (nrow(gene_sets) == 0) {
    stop("No glycolysis gene sets were found in msigdbr with pattern: ",
         gene_set_pattern)
  }

  unique(stats::na.omit(trimws(as.character(gene_sets$gene_symbol))))
}

build_mitocarta_glycolysis_feature_set <- function(
    mitocarta_file,
    mitocarta_symbol_col = "Symbol",
    species = "Homo sapiens",
    glycolysis_subcollections = c("CP:KEGG_LEGACY", "CP:REACTOME"),
    glycolysis_gene_set_pattern = "GLYCOLYSIS|GLUCONEOGENESIS") {
  mitocarta_symbols <- read_mitocarta_symbols(
    mitocarta_file = mitocarta_file,
    symbol_col = mitocarta_symbol_col
  )

  glycolysis_symbols <- get_msigdbr_glycolysis_symbols(
    species = species,
    subcollections = glycolysis_subcollections,
    gene_set_pattern = glycolysis_gene_set_pattern
  )

  feature_set <- dplyr::bind_rows(
    tibble::tibble(Symbol = mitocarta_symbols, Source = "MitoCarta"),
    tibble::tibble(Symbol = glycolysis_symbols, Source = "Glycolysis")
  ) |>
    dplyr::distinct() |>
    dplyr::group_by(.data$Symbol) |>
    dplyr::summarise(
      Source = paste(sort(unique(.data$Source)), collapse = " + "),
      .groups = "drop"
    )

  feature_set
}

read_differential_result_table <- function(file) {
  ext <- tolower(tools::file_ext(file))
  switch(
    ext,
    "csv" = readr::read_csv(file, show_col_types = FALSE),
    "tsv" = readr::read_tsv(file, show_col_types = FALSE),
    "txt" = readr::read_tsv(file, show_col_types = FALSE),
    stop("Unsupported differential result file extension: ", ext)
  )
}

plot_mitocarta_glycolysis_volcano <- function(
    comparisons,
    mitocarta_file,
    out_dir,
    output_prefix,
    feature_col,
    lfc_col,
    padj_col,
    lfc_cut,
    padj_cut,
    title,
    subtitle = "MitoCarta and KEGG/Reactome glycolysis features",
    mitocarta_symbol_col = "Symbol",
    label_significant = TRUE,
    shade_significant_regions = TRUE,
    write_tables = TRUE,
    width = 9,
    height = 5,
    dpi = 300) {
  required_cols <- c("file", "label", "shape_group")
  missing_cols <- setdiff(required_cols, colnames(comparisons))
  if (length(missing_cols) > 0) {
    stop("Comparison table is missing columns: ",
         paste(missing_cols, collapse = ", "))
  }

  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("The ggrepel package is required for volcano labels.")
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  feature_set <- build_mitocarta_glycolysis_feature_set(
    mitocarta_file = mitocarta_file,
    mitocarta_symbol_col = mitocarta_symbol_col
  )

  plot_tbl <- lapply(seq_len(nrow(comparisons)), function(i) {
    info <- comparisons[i, , drop = FALSE]
    result_tbl <- read_differential_result_table(info$file)

    missing_result_cols <- setdiff(c(feature_col, lfc_col, padj_col), colnames(result_tbl))
    if (length(missing_result_cols) > 0) {
      stop("Result table is missing columns in ", info$file, ": ",
           paste(missing_result_cols, collapse = ", "))
    }

    result_tbl |>
      dplyr::transmute(
        Symbol = trimws(as.character(.data[[feature_col]])),
        logFC = as.numeric(.data[[lfc_col]]),
        padj = as.numeric(.data[[padj_col]]),
        Comparison = as.character(info$label),
        ShapeGroup = as.character(info$shape_group),
        File = basename(info$file)
      ) |>
      dplyr::filter(!is.na(.data$Symbol), .data$Symbol != "") |>
      dplyr::inner_join(feature_set, by = "Symbol")
  }) |>
    dplyr::bind_rows() |>
    dplyr::filter(is.finite(.data$logFC), is.finite(.data$padj)) |>
    dplyr::mutate(
      negLog10Padj = -log10(.data$padj + 1e-300),
      Direction = dplyr::case_when(
        .data$padj <= padj_cut & .data$logFC >= lfc_cut ~ "Up",
        .data$padj <= padj_cut & .data$logFC <= -lfc_cut ~ "Down",
        TRUE ~ "Not significant"
      ),
      Label = dplyr::if_else(
        label_significant & .data$Direction != "Not significant",
        .data$Symbol,
        NA_character_
      )
    )

  if (nrow(plot_tbl) == 0) {
    stop("No result rows overlap the MitoCarta/glycolysis feature set.")
  }

  if (isTRUE(write_tables)) {
    readr::write_tsv(
      feature_set,
      file.path(out_dir, paste0(output_prefix, "_feature_set.tsv"))
    )
    readr::write_tsv(
      plot_tbl,
      file.path(out_dir, paste0(output_prefix, "_plot_data.tsv"))
    )
  }

  shape_values <- c("CT" = 16, "CC" = 1)
  shape_labels <- c("CT" = "C/T", "CC" = "C/C")
  missing_shapes <- setdiff(unique(plot_tbl$ShapeGroup), names(shape_values))
  if (length(missing_shapes) > 0) {
    stop("Unsupported shape_group values: ", paste(missing_shapes, collapse = ", "),
         ". Use CT or CC.")
  }

  plot_tbl <- plot_tbl |>
    dplyr::mutate(
      Panel = factor(
        unname(shape_labels[.data$ShapeGroup]),
        levels = unname(shape_labels)
      )
    )

  p <- ggplot2::ggplot(
    plot_tbl,
    ggplot2::aes(x = .data$logFC, y = .data$negLog10Padj)
  )

  if (isTRUE(shade_significant_regions)) {
    p <- p +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = -lfc_cut,
        ymin = -log10(padj_cut), ymax = Inf,
        fill = "#d9edf7", alpha = 0.55
      ) +
      ggplot2::annotate(
        "rect",
        xmin = lfc_cut, xmax = Inf,
        ymin = -log10(padj_cut), ymax = Inf,
        fill = "#d9edf7", alpha = 0.55
      )
  }

  p <- p +
    ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.35) +
    ggplot2::geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", linewidth = 0.35) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$Direction),
      size = 2.5,
      alpha = 0.9
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = .data$Label),
      size = 3,
      max.overlaps = Inf,
      min.segment.length = 0,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$Panel), nrow = 1) +
    ggplot2::scale_color_manual(
      values = c("Up" = "#d73027", "Down" = "#2166ac", "Not significant" = "grey70"),
      breaks = c("Up", "Down", "Not significant"),
      name = NULL
    ) +
    ggplot2::labs(
      title = title,
      x = "log2 fold-change",
      y = "-log10(adjusted p-value)"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.position = "top",
      legend.box = "horizontal",
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey60"),
      strip.text = ggplot2::element_text(face = "bold")
    )

  png_file <- file.path(out_dir, paste0(output_prefix, ".png"))
  pdf_file <- file.path(out_dir, paste0(output_prefix, ".pdf"))
  ggplot2::ggsave(png_file, p, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(pdf_file, p, width = width, height = height)

  invisible(list(
    plot = p,
    plot_data = plot_tbl,
    feature_set = feature_set,
    png_file = png_file,
    pdf_file = pdf_file
  ))
}
