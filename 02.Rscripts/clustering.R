# ============================================================
# Multi-omics C/T vs C/C clustering helpers
# ============================================================
#
# This script contains reusable functions for the two-contrast multi-omics
# stimulation comparison used in the manuscript figures. The same functions are
# used for AC16 and iPSC-CM; dataset-specific files, labels, colors, and output
# folders are defined in config.R.

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

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# ------------------------------------------------------------
# Generic utilities
# ------------------------------------------------------------

safe_mkdir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

read_clustering_result_table <- function(file) {
  ext <- tolower(tools::file_ext(file))
  switch(
    ext,
    "csv" = readr::read_csv(file, show_col_types = FALSE),
    "tsv" = readr::read_tsv(file, show_col_types = FALSE),
    "txt" = readr::read_tsv(file, show_col_types = FALSE),
    stop("Unsupported result table extension for clustering: ", ext)
  )
}

first_existing_column <- function(tbl, candidates, file) {
  hit <- intersect(candidates, colnames(tbl))
  if (length(hit) == 0) {
    stop("None of these columns were found in ", file, ": ",
         paste(candidates, collapse = ", "))
  }
  hit[[1]]
}

scale_lfc_for_heatmap <- function(x, threshold, scale_max) {
  abs_x <- abs(x)
  y <- ifelse(
    abs_x <= threshold,
    abs_x / threshold,
    1 + (abs_x - threshold) / max(scale_max - threshold, 0.1)
  )
  pmax(pmin(sign(x) * y, 2), -2)
}


# ------------------------------------------------------------
# Input loading and cluster assignment
# ------------------------------------------------------------

read_multiomics_clustering_inputs <- function(analysis) {
  comparisons <- analysis$comparisons
  required_cols <- c("omic", "file", "contrast", "lfc_cut", "padj_cut")
  missing_cols <- setdiff(required_cols, colnames(comparisons))
  if (length(missing_cols) > 0) {
    stop("Clustering comparison table is missing columns: ",
         paste(missing_cols, collapse = ", "))
  }

  feature_candidates <- analysis$feature_candidates %||% c(
    "Genes", "Gene_name", "SYMBOL", "Metabolite", "protein_id",
    "Feature", "Name", "ID"
  )
  lfc_candidates <- analysis$lfc_candidates %||% c(
    "logFC", "log2FoldChange", "log2FC"
  )
  padj_candidates <- analysis$padj_candidates %||% c(
    "adj.P.Val", "padj", "FDR", "p.adj"
  )

  rows <- lapply(seq_len(nrow(comparisons)), function(i) {
    info <- comparisons[i, , drop = FALSE]
    if (!file.exists(info$file)) {
      stop("Clustering input file does not exist: ", info$file)
    }

    tbl <- read_clustering_result_table(info$file)
    if (nrow(tbl) == 0) return(NULL)

    feature_col <- first_existing_column(tbl, feature_candidates, info$file)
    lfc_col <- first_existing_column(tbl, lfc_candidates, info$file)
    padj_col <- first_existing_column(tbl, padj_candidates, info$file)

    tibble(
      Omic = as.character(info$omic),
      Feature = as.character(tbl[[feature_col]]),
      logFC = as.numeric(tbl[[lfc_col]]),
      adj.P.Val = as.numeric(tbl[[padj_col]]),
      Contrast = as.character(info$contrast),
      SourceFile = basename(info$file),
      lfc_cut = as.numeric(info$lfc_cut),
      padj_cut = as.numeric(info$padj_cut)
    ) |>
      filter(!is.na(.data$Feature), .data$Feature != "") |>
      mutate(
        FeatureID = paste(.data$Omic, .data$Feature, sep = " | "),
        signif = !is.na(.data$adj.P.Val) &
          !is.na(.data$logFC) &
          .data$adj.P.Val <= .data$padj_cut &
          abs(.data$logFC) >= .data$lfc_cut
      )
  })

  bind_rows(rows)
}

assign_two_contrast_clusters <- function(all_omics,
                                         ct_contrast = "CT",
                                         cc_contrast = "CC",
                                         cluster_labels,
                                         include_mixed = FALSE) {
  keep_ids <- all_omics |>
    filter(.data$signif) |>
    pull(.data$FeatureID) |>
    unique()

  wide_data <- all_omics |>
    filter(.data$FeatureID %in% keep_ids) |>
    select(FeatureID, Omic, Feature, Contrast, logFC, signif) |>
    tidyr::pivot_wider(
      names_from = Contrast,
      values_from = c(logFC, signif),
      values_fn = list(logFC = mean, signif = any),
      values_fill = list(signif = FALSE)
    ) |>
    mutate(across(starts_with("logFC"), ~replace_na(.x, 0.0)))

  ct_lfc <- paste0("logFC_", ct_contrast)
  cc_lfc <- paste0("logFC_", cc_contrast)
  ct_sig <- paste0("signif_", ct_contrast)
  cc_sig <- paste0("signif_", cc_contrast)

  required_cols <- c(ct_lfc, cc_lfc, ct_sig, cc_sig)
  missing_cols <- setdiff(required_cols, colnames(wide_data))
  if (length(missing_cols) > 0) {
    stop("Cluster assignment is missing contrast columns: ",
         paste(missing_cols, collapse = ", "))
  }

  annot_tbl <- wide_data |>
    mutate(
      Cluster = case_when(
        .data[[ct_sig]] & .data[[cc_sig]] &
          sign(.data[[ct_lfc]]) == sign(.data[[cc_lfc]]) &
          .data[[ct_lfc]] > 0 ~ "1",
        .data[[ct_sig]] & .data[[cc_sig]] &
          sign(.data[[ct_lfc]]) == sign(.data[[cc_lfc]]) &
          .data[[ct_lfc]] < 0 ~ "2",
        .data[[ct_sig]] & !.data[[cc_sig]] & .data[[ct_lfc]] > 0 ~ "3",
        .data[[ct_sig]] & !.data[[cc_sig]] & .data[[ct_lfc]] < 0 ~ "4",
        !.data[[ct_sig]] & .data[[cc_sig]] & .data[[cc_lfc]] > 0 ~ "5",
        !.data[[ct_sig]] & .data[[cc_sig]] & .data[[cc_lfc]] < 0 ~ "6",
        TRUE ~ "Mixed"
      ),
      ClusterLabel = dplyr::recode(.data$Cluster, !!!cluster_labels),
      PlotFeature = .data$Cluster != "Mixed" | isTRUE(include_mixed)
    )

  annot_tbl
}


# ------------------------------------------------------------
# Heatmap and legend outputs
# ------------------------------------------------------------

build_multiomics_cluster_heatmap <- function(annot_tbl,
                                             analysis,
                                             plot_settings,
                                             show_heatmap_legend = FALSE,
                                             show_annotation_legend = FALSE) {
  omic_order <- names(plot_settings$omic_colors)
  cluster_order <- names(plot_settings$cluster_colors)
  cluster_labels <- plot_settings$cluster_labels

  ht_list <- NULL
  for (omic_name in omic_order) {
    sub_df <- annot_tbl |>
      filter(.data$Omic == omic_name, .data$PlotFeature) |>
      mutate(
        Cluster = factor(.data$Cluster, levels = cluster_order),
        Omic = factor(.data$Omic, levels = omic_order)
      ) |>
      arrange(.data$Cluster)

    if (nrow(sub_df) == 0) next

    mat_raw <- as.matrix(sub_df[, c("logFC_CT", "logFC_CC")])
    colnames(mat_raw) <- analysis$column_labels

    threshold <- plot_settings$thresholds[[omic_name]]
    scale_max <- plot_settings$scale_max[[omic_name]]
    mat_scaled <- scale_lfc_for_heatmap(mat_raw, threshold, scale_max)

    col_fun <- circlize::colorRamp2(
      c(-2, -1, 0, 1, 2),
      plot_settings$heatmap_colors
    )

    row_ha <- ComplexHeatmap::rowAnnotation(
      OMIC = sub_df$Omic,
      Clusters = sub_df$Cluster,
      col = list(
        OMIC = plot_settings$omic_colors,
        Clusters = plot_settings$cluster_colors
      ),
      show_annotation_name = FALSE,
      show_legend = show_annotation_legend
    )

    ht <- ComplexHeatmap::Heatmap(
      mat_scaled,
      name = paste0(omic_name, " log2FC"),
      col = col_fun,
      column_labels = analysis$column_labels,
      row_split = sub_df$Cluster,
      row_title = NULL,
      row_title_rot = 0,
      left_annotation = row_ha,
      cluster_rows = TRUE,
      cluster_row_slices = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      show_heatmap_legend = show_heatmap_legend,
      column_title = paste0(omic_name, " (n=", nrow(sub_df), ")"),
      border = TRUE,
      height = grid::unit(min(plot_settings$max_panel_height_cm,
                              nrow(sub_df) * plot_settings$row_height_cm + 2), "cm"),
      heatmap_legend_param = list(
        at = c(-2, -1, 0, 1, 2),
        labels = c(
          paste0("-", scale_max),
          paste0("-", threshold),
          "0",
          as.character(threshold),
          as.character(scale_max)
        ),
        title = paste0(omic_name, "\nlog2FC")
      )
    )

    ht_list <- if (is.null(ht_list)) ht else ht_list %v% ht
  }

  if (is.null(ht_list)) {
    stop("No features are available for the multi-omics clustering heatmap.")
  }

  ht_list
}

make_multiomics_cluster_legends <- function(plot_settings) {
  cluster_legend <- ComplexHeatmap::Legend(
    labels = unname(plot_settings$cluster_legend_labels),
    legend_gp = grid::gpar(fill = plot_settings$cluster_colors),
    title = "Cluster legend",
    ncol = 1,
    labels_gp = grid::gpar(fontsize = 8),
    title_gp = grid::gpar(fontsize = 10, fontface = "bold")
  )

  omic_legend <- ComplexHeatmap::Legend(
    labels = names(plot_settings$omic_colors),
    legend_gp = grid::gpar(fill = plot_settings$omic_colors),
    title = "OMIC",
    nrow = 1,
    labels_gp = grid::gpar(fontsize = 8),
    title_gp = grid::gpar(fontsize = 10, fontface = "bold")
  )

  scale_legends <- lapply(names(plot_settings$omic_colors), function(omic_name) {
    threshold <- plot_settings$thresholds[[omic_name]]
    scale_max <- plot_settings$scale_max[[omic_name]]
    col_fun <- circlize::colorRamp2(
      c(-scale_max, -threshold, 0, threshold, scale_max),
      plot_settings$heatmap_colors
    )

    ComplexHeatmap::Legend(
      col_fun = col_fun,
      title = paste0(omic_name, "\nlog2FC"),
      at = c(-scale_max, -threshold, 0, threshold, scale_max),
      labels = c(
        paste0("-", scale_max),
        paste0("-", threshold),
        "0",
        as.character(threshold),
        as.character(scale_max)
      ),
      direction = "vertical",
      legend_height = grid::unit(2.6, "cm"),
      grid_width = grid::unit(0.4, "cm"),
      labels_gp = grid::gpar(fontsize = 7),
      title_gp = grid::gpar(fontsize = 8, fontface = "bold")
    )
  })

  ComplexHeatmap::packLegend(
    cluster_legend,
    omic_legend,
    ComplexHeatmap::packLegend(list = scale_legends, direction = "horizontal"),
    direction = "vertical",
    gap = grid::unit(4, "mm")
  )
}

save_multiomics_heatmap_outputs <- function(ht_list,
                                            legends,
                                            out_dir,
                                            output_prefix,
                                            heatmap_width = 8,
                                            heatmap_height = 14,
                                            legend_width = 9,
                                            legend_height = 4,
                                            dpi = 600) {
  safe_mkdir(out_dir)

  heatmap_pdf <- file.path(out_dir, paste0(output_prefix, "_Clustering_heatmap.pdf"))
  heatmap_png <- file.path(out_dir, paste0(output_prefix, "_Clustering_heatmap.png"))
  legend_pdf <- file.path(out_dir, paste0(output_prefix, "_Clustering_legends.pdf"))
  legend_png <- file.path(out_dir, paste0(output_prefix, "_Clustering_legends.png"))

  grDevices::pdf(heatmap_pdf, width = heatmap_width, height = heatmap_height)
  ComplexHeatmap::draw(ht_list, heatmap_legend_side = "right", annotation_legend_side = "right")
  grDevices::dev.off()

  grDevices::png(heatmap_png, width = heatmap_width, height = heatmap_height,
                 units = "in", res = dpi)
  ComplexHeatmap::draw(ht_list, heatmap_legend_side = "right", annotation_legend_side = "right")
  grDevices::dev.off()

  grDevices::pdf(legend_pdf, width = legend_width, height = legend_height)
  grid::grid.newpage()
  ComplexHeatmap::draw(legends, x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"))
  grDevices::dev.off()

  grDevices::png(legend_png, width = legend_width, height = legend_height,
                 units = "in", res = dpi)
  grid::grid.newpage()
  ComplexHeatmap::draw(legends, x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"))
  grDevices::dev.off()

  invisible(list(
    heatmap_pdf = heatmap_pdf,
    heatmap_png = heatmap_png,
    legend_pdf = legend_pdf,
    legend_png = legend_png
  ))
}


# ------------------------------------------------------------
# Composition and Reactome enrichment
# ------------------------------------------------------------

plot_multiomics_cluster_composition <- function(annot_tbl,
                                                out_dir,
                                                output_prefix,
                                                plot_settings,
                                                width = 8,
                                                height = 5,
                                                dpi = 600) {
  if (!requireNamespace("scatterpie", quietly = TRUE)) {
    message("scatterpie is not installed; skipping omic-composition plot.")
    return(NULL)
  }

  safe_mkdir(out_dir)
  omic_cols <- names(plot_settings$omic_colors)
  cluster_order <- names(plot_settings$cluster_colors)

  pie_tbl <- annot_tbl |>
    dplyr::filter(.data$PlotFeature) |>
    dplyr::count(.data$Cluster, .data$Omic, name = "n") |>
    tidyr::pivot_wider(names_from = Omic, values_from = n, values_fill = 0)

  for (omic_name in omic_cols) {
    if (!(omic_name %in% colnames(pie_tbl))) pie_tbl[[omic_name]] <- 0
  }

  pie_tbl <- pie_tbl |>
    mutate(
      total = rowSums(across(all_of(omic_cols))),
      x = as.numeric(factor(.data$Cluster, levels = cluster_order)),
      y = 1
    ) |>
    filter(.data$total > 0) |>
    mutate(
      r = 0.12 + (.data$total - min(.data$total)) /
        (max(.data$total) - min(.data$total) + 1e-9) * 0.25
    )

  readr::write_tsv(
    pie_tbl,
    file.path(out_dir, paste0(output_prefix, "_Omic_composition.tsv"))
  )

  p <- ggplot2::ggplot() +
    scatterpie::geom_scatterpie(
      data = pie_tbl,
      ggplot2::aes(x = .data$x, y = .data$y, r = .data$r),
      cols = omic_cols,
      color = "grey25",
      linewidth = 0.25
    ) +
    ggplot2::geom_text(
      data = pie_tbl,
      ggplot2::aes(x = .data$x, y = .data$y - .data$r - 0.08,
                   label = paste0("n=", .data$total)),
      size = 3,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = plot_settings$omic_colors, name = "OMIC") +
    ggplot2::scale_x_continuous(
      breaks = seq_along(cluster_order),
      labels = cluster_order,
      limits = c(0.5, length(cluster_order) + 0.5)
    ) +
    ggplot2::scale_y_continuous(breaks = NULL, limits = c(0.45, 1.35)) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Cluster", y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggplot2::ggsave(
    file.path(out_dir, paste0(output_prefix, "_Omic_composition.png")),
    p,
    width = width,
    height = height,
    dpi = dpi
  )
  ggplot2::ggsave(
    file.path(out_dir, paste0(output_prefix, "_Omic_composition.pdf")),
    p,
    width = width,
    height = height
  )

  invisible(p)
}

run_multiomics_cluster_reactome <- function(annot_tbl,
                                            out_dir,
                                            output_prefix,
                                            min_genes = 5,
                                            top_n = 5) {
  safe_mkdir(out_dir)

  ora_tbl <- annot_tbl |>
    filter(.data$PlotFeature, .data$Omic %in% c("Transcriptomic", "Proteomic")) |>
    mutate(Gene = stringr::str_trim(sub("^.*\\|\\s*", "", .data$Feature))) |>
    filter(!is.na(.data$Gene), .data$Gene != "")

  clusters <- sort(unique(ora_tbl$Cluster))
  enrich_list <- list()

  for (cluster_id in clusters) {
    genes <- ora_tbl |>
      filter(.data$Cluster == cluster_id) |>
      pull(.data$Gene) |>
      unique()

    if (length(genes) < min_genes) next
    if (!requireNamespace("gprofiler2", quietly = TRUE)) {
      message("gprofiler2 is not installed; skipping Reactome ORA.")
      break
    }

    gost_res <- tryCatch(
      gprofiler2::gost(
        query = genes,
        organism = "hsapiens",
        sources = "REAC",
        correction_method = "fdr",
        significant = TRUE,
        evcodes = TRUE
      ),
      error = function(e) {
        message("Reactome ORA failed for cluster ", cluster_id, ": ", e$message)
        NULL
      }
    )

    if (!is.null(gost_res$result) && nrow(gost_res$result) > 0) {
      enrich_list[[cluster_id]] <- as_tibble(gost_res$result) |>
        mutate(Cluster = cluster_id, query_size = length(genes)) |>
        select(Cluster, term_name, p_value, term_id,
               intersection_size, query_size, intersection)
    }
  }

  if (length(enrich_list) == 0) {
    message("No Reactome ORA terms were returned.")
    return(invisible(NULL))
  }

  final_tbl <- bind_rows(enrich_list) |>
    arrange(.data$Cluster, .data$p_value)

  readr::write_tsv(
    final_tbl,
    file.path(out_dir, paste0(output_prefix, "_Reactome_ORA_by_cluster.tsv"))
  )

  top_tbl <- final_tbl |>
    group_by(.data$Cluster) |>
    slice_min(order_by = .data$p_value, n = top_n, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      term_name_wrap = stringr::str_wrap(.data$term_name, width = 45),
      log10p = -log10(.data$p_value),
      Cluster = factor(.data$Cluster, levels = sort(unique(.data$Cluster)))
    )

  p <- ggplot2::ggplot(
    top_tbl,
    ggplot2::aes(x = .data$log10p, y = reorder(.data$term_name_wrap, .data$log10p),
                 fill = .data$Cluster)
  ) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~Cluster, scales = "free_y", ncol = 2) +
    ggplot2::labs(x = "-log10(FDR p-value)", y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "none")

  ggplot2::ggsave(
    file.path(out_dir, paste0(output_prefix, "_Reactome_ORA_top", top_n, ".png")),
    p,
    width = 11,
    height = 8,
    dpi = 300
  )
  ggplot2::ggsave(
    file.path(out_dir, paste0(output_prefix, "_Reactome_ORA_top", top_n, ".pdf")),
    p,
    width = 11,
    height = 8
  )

  invisible(final_tbl)
}


# ------------------------------------------------------------
# Top-level runner
# ------------------------------------------------------------

run_multiomics_ct_cc_clustering <- function(analysis,
                                            plot_settings,
                                            make_heatmap = TRUE,
                                            make_composition = TRUE,
                                            run_ora = TRUE) {
  safe_mkdir(analysis$out_dir)

  all_omics <- read_multiomics_clustering_inputs(analysis)
  annot_tbl <- assign_two_contrast_clusters(
    all_omics = all_omics,
    ct_contrast = "CT",
    cc_contrast = "CC",
    cluster_labels = plot_settings$cluster_labels,
    include_mixed = isTRUE(plot_settings$include_mixed)
  )

  readr::write_tsv(
    annot_tbl,
    file.path(analysis$out_dir, paste0(analysis$output_prefix, "_Cluster_assignment.tsv"))
  )

  ht_list <- NULL
  legends <- NULL
  if (isTRUE(make_heatmap)) {
    ht_list <- build_multiomics_cluster_heatmap(
      annot_tbl = annot_tbl,
      analysis = analysis,
      plot_settings = plot_settings,
      show_heatmap_legend = FALSE,
      show_annotation_legend = FALSE
    )
    legends <- make_multiomics_cluster_legends(plot_settings)
    save_multiomics_heatmap_outputs(
      ht_list = ht_list,
      legends = legends,
      out_dir = analysis$out_dir,
      output_prefix = analysis$output_prefix,
      heatmap_width = analysis$heatmap_width %||% 8,
      heatmap_height = analysis$heatmap_height %||% 14,
      legend_width = analysis$legend_width %||% 9,
      legend_height = analysis$legend_height %||% 4,
      dpi = analysis$dpi %||% 600
    )
  }

  composition_plot <- NULL
  if (isTRUE(make_composition)) {
    composition_plot <- plot_multiomics_cluster_composition(
      annot_tbl = annot_tbl,
      out_dir = analysis$out_dir,
      output_prefix = analysis$output_prefix,
      plot_settings = plot_settings,
      dpi = analysis$dpi %||% 600
    )
  }

  ora_tbl <- NULL
  if (isTRUE(run_ora)) {
    ora_tbl <- run_multiomics_cluster_reactome(
      annot_tbl = annot_tbl,
      out_dir = analysis$out_dir,
      output_prefix = analysis$output_prefix
    )
  }

  invisible(list(
    all_omics = all_omics,
    table = annot_tbl,
    heatmap = ht_list,
    legends = legends,
    composition_plot = composition_plot,
    ora = ora_tbl
  ))
}
