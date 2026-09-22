# ============================================================
# Transcriptomics helper functions
# ============================================================
#
# This script contains the transcriptomics-specific pieces of the analysis:
# sample metadata construction, count loading/filtering, DESeq2 normalization,
# gene annotation, and DEG table export.
#
# Generic plotting functions live in:
#   02.Rscripts/plots.R
#
# The plotting calls below intentionally keep the legacy output filenames and
# visual defaults used by the manuscript-producing code.

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

plots_file <- if (file.exists("02.Rscripts/plots.R")) {
  "02.Rscripts/plots.R"
} else if (file.exists("plots.R")) {
  "plots.R"
} else {
  NA_character_
}

if (!is.na(plots_file)) {
  source(plots_file)
}
rm(plots_file)

# ------------------------------------------------------------
# Sample metadata helpers
# ------------------------------------------------------------

get_sample_no <- function(x) {
  m <- sub("^([0-9]+)_.*", "\\1", x)
  ifelse(grepl("^[0-9]+$", m), as.integer(m), NA_integer_)
}

build_metadata <- function(counts_path) {
  hdr <- data.table::fread(counts_path, nrows = 0)
  all_cols <- colnames(hdr)
  nonsample <- c("Geneid","Chr","Start","End","Strand","Length")
  sample_cols <- setdiff(all_cols, nonsample)
  # Drop external/control samples from the expression matrix.
  sample_cols <- setdiff(sample_cols, grep("^Cneg|^Cpos", sample_cols, value = TRUE))

  sample_no <- vapply(sample_cols, get_sample_no, integer(1))

  # Sample-number groups defining dataset, line, and stimulation status.
  iPSC_6132_S   <- c(17,20,33,42)
  iPSC_6132_NS  <- c(18,19,26,48)
  iPSC_6135_S   <- c(25,32,37,47)
  iPSC_6135_NS  <- c(22,31,41,46)
  iPSC_6137_S   <- c(21,24,30,45)
  iPSC_6137_NS  <- c(29,36,40,44)
  iPSC_6135H1_S <- c(23,28,35,39)
  iPSC_6135H1_NS<- c(27,34,38,43)

  AC16_1184_NS <- c(1,3,5,11)
  AC16_1184_S  <- c(2,4,6,12)
  AC16_127_NS  <- c(7,9,13,15)
  AC16_127_S   <- c(8,10,14,16)

  flag <- function(x, set) x %in% set

  tibble(
    sample = sample_cols,
    sample_no = sample_no
  ) %>%
    mutate(
      dataset = case_when(
        flag(sample_no, c(iPSC_6132_S,iPSC_6132_NS,
                          iPSC_6135_S,iPSC_6135_NS,
                          iPSC_6137_S,iPSC_6137_NS,
                          iPSC_6135H1_S,iPSC_6135H1_NS)) ~ "iPSC",
        flag(sample_no, c(AC16_1184_NS,AC16_1184_S,AC16_127_NS,AC16_127_S)) ~ "AC16",
        TRUE ~ "UNKNOWN"
      ),
      line = case_when(
        flag(sample_no, c(iPSC_6132_S, iPSC_6132_NS)) ~ "6132",
        flag(sample_no, c(iPSC_6135_S, iPSC_6135_NS)) ~ "6135",
        flag(sample_no, c(iPSC_6137_S, iPSC_6137_NS)) ~ "6137",
        flag(sample_no, c(iPSC_6135H1_S, iPSC_6135H1_NS)) ~ "6135H1",
        flag(sample_no, c(AC16_1184_NS,AC16_1184_S)) ~ "AC16_1.184",
        flag(sample_no, c(AC16_127_NS,AC16_127_S)) ~ "AC16_1.27",
        TRUE ~ NA_character_
      ),
      stim = case_when(
        flag(sample_no, c(iPSC_6132_S, iPSC_6135_S, iPSC_6137_S, iPSC_6135H1_S,
                          AC16_1184_S, AC16_127_S)) ~ "S",
        flag(sample_no, c(iPSC_6132_NS, iPSC_6135_NS, iPSC_6137_NS, iPSC_6135H1_NS,
                          AC16_1184_NS, AC16_127_NS)) ~ "NS",
        TRUE ~ NA_character_
      ),
      celltype = case_when(
        dataset == "iPSC" ~ "iPSC-CM",
        dataset == "AC16" ~ "AC16",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(dataset, celltype, line, stim, sample_no)
}

# ------------------------------------------------------------
# Count loading and filtering helpers
# ------------------------------------------------------------

load_counts_subset <- function(counts_path, samples) {
  dt <- data.table::fread(counts_path)
  keep_cols <- c("Geneid","Length", samples)
  dt <- dt[, ..keep_cols]
  rn <- dt$Geneid
  mat <- as.matrix(dt[, setdiff(colnames(dt), c("Geneid","Length")), with=FALSE])
  rownames(mat) <- rn
  mat
}

filter_low_counts <- function(cts, min_count=10, min_samps=3) {
  keep <- rowSums(cts >= min_count) >= min_samps
  cts[keep, , drop=FALSE]
}

# ------------------------------------------------------------
# Normalization and QC
# ------------------------------------------------------------

normalize_and_qc <- function(counts_path, meta, dataset, out_dir="results") {
  meta_ds <- meta %>% filter(dataset == !!dataset)
  cts <- load_counts_subset(counts_path, meta_ds$sample)
  cts <- filter_low_counts(cts, 10, 3)

  coldata <- meta_ds %>% column_to_rownames("sample")
  dds <- DESeqDataSetFromMatrix(cts, coldata, design = ~ 1)
  dds <- estimateSizeFactors(dds)

  # Export normalized and variance-stabilized matrices.
  norm_dir <- file.path(out_dir, "normalized")
  qc_dir   <- file.path(out_dir, "qc")
  dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

  norm_counts <- counts(dds, normalized = TRUE)
  readr::write_tsv(as_tibble(norm_counts, rownames = "Geneid"),
                   file.path(norm_dir, paste0("normalized_counts_", dataset, ".tsv")))

  vst_mat <- assay(vst(dds, blind = TRUE))
  readr::write_tsv(as_tibble(vst_mat, rownames = "Geneid"),
                   file.path(norm_dir, paste0("vst_", dataset, ".tsv")))

  # PCA QC plot.
  p <- plot_pca_qc(
    mat = vst_mat,
    metadata = meta_ds,
    sample_col = "sample",
    color_col = "line",
    shape_col = "stim",
    label_col = "sample",
    title = paste0(dataset, " — PCA (VST)"),
    output_file = file.path(qc_dir, paste0("PCA_", dataset, ".png")),
    width = 7,
    height = 5,
    dpi = 200
  )
  print(p)

  # Sample-distance heatmap.
  plot_sample_distance_heatmap(
    mat = vst_mat,
    metadata = meta_ds,
    sample_col = "sample",
    annotation_cols = c("line", "stim"),
    title = paste0(dataset, " — sample distances"),
    output_file = file.path(qc_dir, paste0("sample_distance_", dataset, ".png")),
    width = 1000,
    height = 900
  )

  invisible(list(dds=dds, meta=meta_ds))
}

# ------------------------------------------------------------
# Gene annotation helpers
# ------------------------------------------------------------

add_gene_annotations <- function(tbl, species = c("human","rat","mouse")) {
  species <- match.arg(species)
  pkg <- switch(species,
    human = "org.Hs.eg.db",
    rat   = "org.Rn.eg.db",
    mouse = "org.Mm.eg.db"
  )
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing annotation package: ", pkg,
         ". Install it before running transcriptomic annotation.",
         call. = FALSE)
  }
  # Get the OrgDb object exported by the package, not only the namespace.
  db <- get(paste0(pkg), envir = asNamespace(pkg))   # e.g. org.Hs.eg.db::org.Hs.eg.db

  # Strip Ensembl version suffixes before annotation.
  ens_with_ver <- if ("Geneid" %in% names(tbl)) tbl$Geneid else rownames(tbl)
  ens <- sub("\\..*$", "", ens_with_ver)

  # Map to HGNC, RGD, or MGI symbols depending on species.
  sym <- AnnotationDbi::mapIds(db, keys = ens, keytype = "ENSEMBL",
                               column = "SYMBOL", multiVals = "first")

  # Keep the original Geneid column if present.
  if (!"Geneid" %in% names(tbl)) {
    tbl <- tibble::as_tibble(tbl, rownames = "Geneid")
  }
  tbl %>%
    dplyr::mutate(
      ENSEMBL_ID = ens,
      Gene_name  = unname(sym[ENSEMBL_ID])
    ) %>%
    dplyr::relocate(ENSEMBL_ID, Gene_name, .before = 1)
}


# ------------------------------------------------------------
# Public heart-tissue immune-axis expression exploration
# ------------------------------------------------------------
#
# These helpers generate public-data expression and gene-gene correlation
# heatmaps for the IFN-gamma/CXCR3 chemokine axis plus broad immune-cell
# markers. RNA-seq count tables are normalized with DESeq2 size factors.
# Agilent microarray Feature Extraction files are read from gProcessedSignal,
# log2-transformed, quantile-normalized with limma, and collapsed to gene
# symbols before plotting.

public_group_levels <- function() c("CCC", "Control", "DCM")

public_group_colors <- function() {
  c(
    "CCC" = "#B2182B",
    "Control" = "#2166AC",
    "DCM" = "#4DAF4A"
  )
}

infer_public_group <- function(sample) {
  dplyr::case_when(
    stringr::str_detect(sample, "^sevCCC|_CCC_") ~ "CCC",
    stringr::str_detect(sample, "^CTRL|_CTL_") ~ "Control",
    stringr::str_detect(sample, "^DCM|_DCM_") ~ "DCM",
    TRUE ~ "Other"
  )
}

public_row_center_scale <- function(mat) {
  scaled <- t(apply(mat, 1, function(x) {
    if (all(is.na(x))) return(rep(NA_real_, length(x)))

    x_sd <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(x_sd) || x_sd == 0) {
      return(rep(0, length(x)))
    }

    (x - mean(x, na.rm = TRUE)) / x_sd
  }))
  dimnames(scaled) <- dimnames(mat)
  scaled
}

select_public_ensembl_rows_by_symbol <- function(mat,
                                                 target_genes,
                                                 species = c("human", "rat", "mouse")) {
  species <- match.arg(species)
  pkg <- switch(species,
    human = "org.Hs.eg.db",
    rat   = "org.Rn.eg.db",
    mouse = "org.Mm.eg.db"
  )
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing annotation package: ", pkg,
         ". Install it before running public expression plots.",
         call. = FALSE)
  }

  db <- get(paste0(pkg), envir = asNamespace(pkg))
  ens <- sub("\\..*$", "", rownames(mat))
  sym <- suppressMessages(
    AnnotationDbi::mapIds(
      db,
      keys = ens,
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )
  )

  annot <- tibble::tibble(
    feature_id = rownames(mat),
    Gene_name = unname(sym[ens]),
    target_order = match(toupper(unname(sym[ens])), toupper(target_genes)),
    mean_expression = rowMeans(mat, na.rm = TRUE)
  ) |>
    dplyr::filter(!is.na(.data$target_order)) |>
    dplyr::arrange(.data$target_order, dplyr::desc(.data$mean_expression)) |>
    dplyr::distinct(.data$Gene_name, .keep_all = TRUE)

  target_mat <- matrix(
    NA_real_,
    nrow = length(target_genes),
    ncol = ncol(mat),
    dimnames = list(target_genes, colnames(mat))
  )

  if (nrow(annot) > 0) {
    target_mat[annot$Gene_name, ] <- mat[annot$feature_id, , drop = FALSE]
  }

  missing_genes <- setdiff(target_genes, annot$Gene_name)
  if (length(missing_genes) > 0) {
    message("Public expression table is missing target gene(s): ",
            paste(missing_genes, collapse = ", "))
  }

  target_mat
}

read_public_rnaseq_normalized_matrix <- function(input_path,
                                                 target_genes,
                                                 species = "human",
                                                 pseudo_count = 1) {
  if (!file.exists(input_path)) {
    stop("Public normalized RNA-seq file does not exist: ", input_path,
         call. = FALSE)
  }

  sample_cols <- strsplit(readLines(input_path, n = 1), "\t", fixed = TRUE)[[1]]
  sample_cols <- stringr::str_replace_all(sample_cols, '^"|"$', "")

  expr_tbl <- data.table::fread(
    input_path,
    skip = 1,
    header = FALSE,
    col.names = c("Geneid", sample_cols)
  )

  metadata <- tibble::tibble(
    sample = sample_cols,
    sample_label = sample_cols,
    group = infer_public_group(.data$sample),
    sample_no = readr::parse_number(.data$sample)
  ) |>
    dplyr::mutate(group = factor(.data$group, levels = c(public_group_levels(), "Other"))) |>
    dplyr::arrange(.data$group, .data$sample_no)

  expr_mat <- expr_tbl |>
    dplyr::select("Geneid", dplyr::all_of(metadata$sample)) |>
    tibble::column_to_rownames("Geneid") |>
    as.matrix()
  storage.mode(expr_mat) <- "double"
  expr_mat <- log2(expr_mat + pseudo_count)

  list(
    expression = select_public_ensembl_rows_by_symbol(expr_mat, target_genes, species),
    metadata = metadata
  )
}

read_public_rnaseq_count_matrix <- function(input_path,
                                            target_genes,
                                            species = "human",
                                            pseudo_count = 1) {
  if (!file.exists(input_path)) {
    stop("Public RNA-seq count file does not exist: ", input_path,
         call. = FALSE)
  }

  count_tbl <- readr::read_tsv(input_path, show_col_types = FALSE)
  if (!"Geneid" %in% colnames(count_tbl)) {
    stop("Public RNA-seq count table must contain a Geneid column: ",
         input_path, call. = FALSE)
  }

  sample_cols <- setdiff(colnames(count_tbl), "Geneid")
  metadata <- tibble::tibble(
    sample = sample_cols,
    sample_label = sample_cols,
    group = infer_public_group(.data$sample),
    sample_no = readr::parse_number(.data$sample)
  ) |>
    dplyr::mutate(group = factor(.data$group, levels = c(public_group_levels(), "Other"))) |>
    dplyr::arrange(.data$group, .data$sample_no)

  count_mat <- count_tbl |>
    dplyr::select("Geneid", dplyr::all_of(metadata$sample)) |>
    tibble::column_to_rownames("Geneid") |>
    as.matrix()
  storage.mode(count_mat) <- "double"

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = metadata |> tibble::column_to_rownames("sample"),
    design = ~ 1
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  expr_mat <- log2(DESeq2::counts(dds, normalized = TRUE) + pseudo_count)

  list(
    expression = select_public_ensembl_rows_by_symbol(expr_mat, target_genes, species),
    metadata = metadata
  )
}

parse_public_microarray_sample <- function(path) {
  stem <- basename(path) |>
    stringr::str_remove("\\.txt\\.gz$") |>
    stringr::str_remove("\\.txt$")
  parts <- strsplit(stem, "_", fixed = TRUE)[[1]]
  group_raw <- parts[3]

  tibble::tibble(
    sample = stem,
    sample_label = if (length(parts) >= 2) parts[2] else stem,
    group = dplyr::case_when(
      group_raw == "CCC" ~ "CCC",
      group_raw == "CTL" ~ "Control",
      group_raw == "DCM" ~ "DCM",
      TRUE ~ "Other"
    ),
    gsm = parts[1],
    raw_file = path
  )
}

read_public_agilent_feature_table <- function(path) {
  feature_tbl <- data.table::fread(
    path,
    skip = "FEATURES",
    select = c(
      "FEATURES", "ControlType", "GeneName", "SystematicName",
      "gProcessedSignal", "gIsFound"
    )
  )
  colnames(feature_tbl)[1] <- "RowType"

  feature_tbl |>
    dplyr::filter(
      .data$RowType == "DATA",
      .data$ControlType == 0,
      !is.na(.data$GeneName),
      .data$GeneName != "",
      is.finite(.data$gProcessedSignal)
    ) |>
    dplyr::mutate(
      Gene_name = stringr::str_trim(.data$GeneName),
      signal = pmax(as.numeric(.data$gProcessedSignal), 1)
    ) |>
    dplyr::group_by(.data$Gene_name) |>
    dplyr::summarise(signal = stats::median(.data$signal, na.rm = TRUE), .groups = "drop")
}

read_public_microarray_matrix <- function(input_dir,
                                          target_genes,
                                          groups_keep = public_group_levels()) {
  if (!dir.exists(input_dir)) {
    stop("Public microarray directory does not exist: ", input_dir,
         call. = FALSE)
  }

  files <- list.files(input_dir, pattern = "\\.txt(\\.gz)?$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No Agilent microarray text files found in: ", input_dir,
         call. = FALSE)
  }

  metadata <- purrr::map_dfr(files, parse_public_microarray_sample) |>
    dplyr::mutate(group = factor(.data$group, levels = c(public_group_levels(), "Other"))) |>
    dplyr::filter(.data$group %in% groups_keep) |>
    dplyr::arrange(.data$group, .data$sample_label) |>
    dplyr::mutate(sample_no = dplyr::row_number())

  if (nrow(metadata) == 0) {
    stop("No public microarray samples remained after group filtering.",
         call. = FALSE)
  }

  sample_vectors <- purrr::map(metadata$raw_file, read_public_agilent_feature_table)
  gene_universe <- Reduce(intersect, purrr::map(sample_vectors, ~ .x$Gene_name))
  if (length(gene_universe) == 0) {
    stop("No common microarray gene symbols across selected samples.",
         call. = FALSE)
  }

  signal_list <- purrr::map2(sample_vectors, metadata$sample, function(tbl, sample_id) {
    tbl |>
      dplyr::filter(.data$Gene_name %in% gene_universe) |>
      dplyr::arrange(.data$Gene_name) |>
      dplyr::pull(.data$signal)
  })
  signal_mat <- do.call(cbind, signal_list)

  rownames(signal_mat) <- sort(gene_universe)
  colnames(signal_mat) <- metadata$sample
  log_mat <- log2(signal_mat)
  norm_mat <- limma::normalizeBetweenArrays(log_mat, method = "quantile")

  target_mat <- matrix(
    NA_real_,
    nrow = length(target_genes),
    ncol = ncol(norm_mat),
    dimnames = list(target_genes, colnames(norm_mat))
  )
  present <- intersect(target_genes, rownames(norm_mat))
  if (length(present) > 0) {
    target_mat[present, ] <- norm_mat[present, , drop = FALSE]
  }

  missing_genes <- setdiff(target_genes, present)
  if (length(missing_genes) > 0) {
    message("Public microarray data is missing target gene(s): ",
            paste(missing_genes, collapse = ", "))
  }

  list(expression = target_mat, metadata = metadata)
}

filter_public_expression_groups <- function(expr_mat,
                                            metadata,
                                            groups_keep = public_group_levels()) {
  if (!"sample_no" %in% colnames(metadata)) {
    metadata$sample_no <- seq_len(nrow(metadata))
  }

  metadata <- metadata |>
    dplyr::filter(
      .data$sample %in% colnames(expr_mat),
      .data$group %in% groups_keep
    ) |>
    dplyr::mutate(group = factor(as.character(.data$group), levels = public_group_levels())) |>
    dplyr::arrange(.data$group, .data$sample_no, .data$sample_label)

  available_groups <- unique(as.character(metadata$group))
  missing_groups <- setdiff(groups_keep, available_groups)
  if (length(missing_groups) > 0) {
    message("No public samples available for group(s): ",
            paste(missing_groups, collapse = ", "))
  }

  list(
    expression = expr_mat[, metadata$sample, drop = FALSE],
    metadata = metadata
  )
}

plot_public_expression_heatmap <- function(expr_mat,
                                           metadata,
                                           out_dir,
                                           output_prefix,
                                           dataset_label,
                                           target_genes = rownames(expr_mat),
                                           width = 9.2,
                                           height = 6.2,
                                           dpi = 300) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expr_mat <- expr_mat[target_genes, metadata$sample, drop = FALSE]
  width <- max(width, 0.38 * ncol(expr_mat) + 2.8)
  scaled_mat <- public_row_center_scale(expr_mat)
  scaled_mat <- pmax(pmin(scaled_mat, 2), -2)

  heatmap_colors <- if (exists("CLUSTERING_PLOT_SETTINGS") &&
                        !is.null(CLUSTERING_PLOT_SETTINGS$heatmap_colors)) {
    CLUSTERING_PLOT_SETTINGS$heatmap_colors
  } else {
    c("#00265E", "#5E7F9D", "white", "#A85B61", "#67001E")
  }

  annot_colors <- public_group_colors()
  annot_colors <- annot_colors[intersect(names(annot_colors), levels(droplevels(metadata$group)))]

  top_annot <- ComplexHeatmap::HeatmapAnnotation(
    Phenotype = metadata$group,
    col = list(Phenotype = annot_colors),
    annotation_name_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    simple_anno_size = grid::unit(4, "mm")
  )

  output_file <- file.path(out_dir, paste0(output_prefix, "_expression_heatmap.png"))
  grDevices::png(output_file, width = width, height = height, units = "in", res = dpi)
  on.exit(grDevices::dev.off(), add = TRUE)

  ht <- ComplexHeatmap::Heatmap(
    scaled_mat,
    name = "Row z-score",
    col = circlize::colorRamp2(c(-2, -1, 0, 1, 2), heatmap_colors),
    na_col = "grey90",
    top_annotation = top_annot,
    column_split = metadata$group,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    column_labels = metadata$sample_label,
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Row z-score",
      at = c(-2, -1, 0, 1, 2)
    )
  )

  ComplexHeatmap::draw(
    ht,
    column_title = dataset_label,
    column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = grid::unit(c(10, 10, 10, 10), "mm")
  )

  invisible(list(heatmap = ht, file = output_file, matrix = scaled_mat))
}

plot_public_correlation_heatmaps <- function(expr_mat,
                                             metadata,
                                             out_dir,
                                             output_prefix,
                                             dataset_label,
                                             target_genes = rownames(expr_mat),
                                             width = 9.5,
                                             height = 5.4,
                                             dpi = 300) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expr_mat <- expr_mat[target_genes, metadata$sample, drop = FALSE]
  groups_present <- levels(droplevels(metadata$group))

  corr_df <- purrr::map_dfr(groups_present, function(group_name) {
    group_samples <- metadata$sample[metadata$group == group_name]
    group_mat <- expr_mat[, group_samples, drop = FALSE]
    corr_mat <- suppressWarnings(
      stats::cor(t(group_mat), use = "pairwise.complete.obs", method = "pearson")
    )
    diag(corr_mat) <- 1

    as.data.frame(as.table(corr_mat), stringsAsFactors = FALSE) |>
      tibble::as_tibble() |>
      dplyr::rename(gene_y = "Var1", gene_x = "Var2", correlation = "Freq") |>
      dplyr::mutate(group = group_name)
  })

  corr_df <- corr_df |>
    dplyr::mutate(
      gene_y = factor(.data$gene_y, levels = rev(target_genes)),
      gene_x = factor(.data$gene_x, levels = target_genes),
      group = factor(.data$group, levels = public_group_levels())
    )

  p <- ggplot2::ggplot(
    corr_df,
    ggplot2::aes(x = .data$gene_x, y = .data$gene_y, fill = .data$correlation)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradientn(
      colors = c("#2166AC", "white", "#B2182B"),
      limits = c(-1, 1),
      na.value = "grey90",
      name = "Pearson r"
    ) +
    ggplot2::facet_wrap(~group, nrow = 1) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0(dataset_label, " gene-gene correlations"),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      strip.text = ggplot2::element_text(face = "bold")
    )

  output_file <- file.path(out_dir, paste0(output_prefix, "_correlation_heatmaps.png"))
  ggplot2::ggsave(output_file, p, width = width, height = height, dpi = dpi)

  invisible(list(plot = p, file = output_file, plot_data = corr_df))
}

read_public_expression_analysis <- function(analysis) {
  input_type <- analysis$input_type
  if (identical(input_type, "rnaseq_normalized")) {
    read_public_rnaseq_normalized_matrix(
      input_path = analysis$input_path,
      target_genes = analysis$target_genes
    )
  } else if (identical(input_type, "rnaseq_counts")) {
    read_public_rnaseq_count_matrix(
      input_path = analysis$input_path,
      target_genes = analysis$target_genes
    )
  } else if (identical(input_type, "microarray_agilent")) {
    read_public_microarray_matrix(
      input_dir = analysis$input_path,
      target_genes = analysis$target_genes,
      groups_keep = analysis$groups_keep
    )
  } else {
    stop("Unsupported public expression input_type: ", input_type,
         call. = FALSE)
  }
}

run_public_expression_exploration <- function(analysis) {
  public_expr <- read_public_expression_analysis(analysis)
  filtered <- filter_public_expression_groups(
    expr_mat = public_expr$expression,
    metadata = public_expr$metadata,
    groups_keep = analysis$groups_keep
  )

  expression <- plot_public_expression_heatmap(
    expr_mat = filtered$expression,
    metadata = filtered$metadata,
    out_dir = analysis$out_dir,
    output_prefix = analysis$output_prefix,
    dataset_label = analysis$dataset_label,
    target_genes = analysis$target_genes
  )

  correlations <- plot_public_correlation_heatmaps(
    expr_mat = filtered$expression,
    metadata = filtered$metadata,
    out_dir = analysis$out_dir,
    output_prefix = analysis$output_prefix,
    dataset_label = analysis$dataset_label,
    target_genes = analysis$target_genes
  )

  readr::write_tsv(
    filtered$metadata,
    file.path(analysis$out_dir, paste0(analysis$output_prefix, "_sample_metadata.tsv"))
  )

  invisible(list(
    expression = expression,
    correlations = correlations,
    metadata = filtered$metadata
  ))
}

# Select top labels for optional volcano/MA-style annotations.
pick_top_labels <- function(tbl_plot, n = 10, y_col = c("padj","pvalue"),
                            by = c("p", "both_dir"), fc_col = "log2FoldChange") {
  y_col <- match.arg(y_col)
  by <- match.arg(by)

  # Fall back to stable IDs when gene symbols are missing.
  if (!"Gene_name" %in% names(tbl_plot)) tbl_plot$Gene_name <- tbl_plot$Geneid
  lbl <- dplyr::coalesce(tbl_plot$Gene_name, tbl_plot$ENSEMBL_ID, tbl_plot$Geneid)

  if (by == "p") {
    # Top n by significance.
    keep <- order(tbl_plot[[y_col]], decreasing = FALSE)
    unique(lbl[keep])[seq_len(min(n, length(lbl)))]
  } else {
    # Top n/2 up-regulated and n/2 down-regulated labels.
    sig <- tbl_plot |>
      dplyr::mutate(sig = (is.finite(.data[[y_col]]) & .data[[y_col]] < ifelse(y_col=="padj", 0.05, 0.001)))
    up  <- sig |>
      dplyr::filter(sig, .data[[fc_col]] > 0) |>
      dplyr::arrange(.data[[y_col]]) |>
      dplyr::slice_head(n = ceiling(n/2))
    dn  <- sig |>
      dplyr::filter(sig, .data[[fc_col]] < 0) |>
      dplyr::arrange(.data[[y_col]]) |>
      dplyr::slice_head(n = floor(n/2))
    unique(c(up$Gene_name %||% up$ENSEMBL_ID %||% up$Geneid,
             dn$Gene_name %||% dn$ENSEMBL_ID %||% dn$Geneid))
  }
}

# Small infix helper to coalesce vector-like values.
`%||%` <- function(a, b) ifelse(is.na(a) | a=="", b, a)


# ------------------------------------------------------------
# Differential expression analysis
# ------------------------------------------------------------

run_deg <- function(counts_path, meta, dataset, subset_expr, design, contrast, tag, out_dir="results") {
  lfc_cut = 1.5
  padj_cut = 0.05
  
  # Filter metadata to the requested dataset and comparison subset.
  meta_ds <- meta %>% 
    filter(dataset == !!dataset) %>% 
    filter(!!rlang::parse_expr(subset_expr))
  
  stopifnot(nrow(meta_ds) > 0)

  # Load and filter count matrix for the selected samples.
  cts <- load_counts_subset(counts_path, meta_ds$sample)
  cts <- filter_low_counts(cts, 10, 3)

  # Run DESeq2.
  coldata <- meta_ds %>% column_to_rownames("sample")
  dds <- DESeqDataSetFromMatrix(cts, coldata, design = as.formula(design))
  dds <- DESeq(dds)

  # Preserve the legacy shrinkage behavior used by the manuscript code:
  # shrink the last DESeq2 coefficient. The contrast argument is kept in the
  # function signature for compatibility with existing pipeline calls.
  res_names <- resultsNames(dds)
  coef_name <- res_names[length(res_names)] 
  shrinked  <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  
  tbl <- as_tibble(shrinked, rownames = "Geneid")
  tbl <- add_gene_annotations(tbl, species = "human")

  # Add a manuscript-style significance flag.
  tbl <- tbl %>%
    mutate(Significant = case_when(
      abs(log2FoldChange) >= lfc_cut & padj <= padj_cut ~ "YES",
      TRUE ~ "NO"
    ))

  # Select the statistical column used for display and volcano plotting.
  if (!"padj" %in% names(tbl))  tbl$padj  <- NA_real_
  if (!"pvalue" %in% names(tbl)) tbl$pvalue <- NA_real_
  y_col <- if (sum(is.finite(tbl$padj)) >= 10) "padj" else "pvalue"
  
  tbl <- tbl %>% arrange(.data[[y_col]], pvalue)

  # Export DEG result tables.
  out_dir_deg <- file.path(out_dir)
  dir.create(out_dir_deg, recursive = TRUE, showWarnings = FALSE)
  
  tbl$FC <- sign(tbl$log2FoldChange) * 2^abs(tbl$log2FoldChange)
  
  readr::write_tsv(tbl, file.path(out_dir_deg, paste0(tag, ".tsv")))
  writexl::write_xlsx(tbl, file.path(out_dir_deg, paste0(tag, ".xlsx")))

  # Volcano plot.
  plot_volcano_qc(
    tbl = tbl,
    lfc_col = "log2FoldChange",
    p_col = y_col,
    lfc_cut = lfc_cut,
    p_cut = padj_cut,
    title = tag,
    output_file = file.path(out_dir_deg, paste0(tag, "_volcano.png")),
    width = 8,
    height = 4,
    dpi = 200
  )

  message("DEG saved for: ", tag)
  invisible(tbl)
}
