# ============================================================
# Metabolomics DEM analysis helpers
# ============================================================
#
# This file contains reusable functions for metabolomics differential
# metabolite (DEM) table generation and QC plots. The implementation preserves
# the historical processing choices:
# - data are split by cell type before normalization,
# - standards are excluded,
# - intensities are log2-transformed with half-minimum positive pseudocount,
# - samples are median-centered,
# - outliers are removed after the initial QC plots and before limma.

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
# Generic utilities
# ------------------------------------------------------------

safe_mkdir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}


# ------------------------------------------------------------
# Metadata and input loading
# ------------------------------------------------------------

build_metabolomics_metadata <- function() {
  sample_id <- sprintf("Sample_%02d", 1:48)
  clones <- c("1184", "1.27", "6132", "6137", "6135", "6135-H1")

  pheno <- data.frame(
    SampleID = sample_id,
    CellType = c(rep("AC16", 2 * 8), rep("iPSCM", 4 * 8)),
    Clone = rep(clones, each = 8),
    Stimulation = rep(rep(c("Non_stimulated", "Stimulated"), each = 4), times = 6),
    Replicate = rep(1:4, times = 12),
    stringsAsFactors = FALSE
  )
  rownames(pheno) <- pheno$SampleID
  pheno
}

read_metabolomics_raw_data <- function(raw_data_file, qual_cut = 4) {
  raw_data <- data.table::fread(raw_data_file)
  raw_data[raw_data$Qual. < qual_cut, ]
}


# ------------------------------------------------------------
# Preprocessing
# ------------------------------------------------------------

prepare_metabolomics_matrix <- function(raw_data, metadata) {
  metab_data <- raw_data[Class != "Stds."]
  target_samples <- metadata$SampleID

  missing_samples <- setdiff(target_samples, colnames(metab_data))
  if (length(missing_samples) > 0) {
    stop("Metabolomics raw data is missing samples: ",
         paste(missing_samples, collapse = ", "))
  }

  mat <- as.matrix(metab_data[, ..target_samples])
  storage.mode(mat) <- "double"
  rownames(mat) <- metab_data$Name

  min_nonzero <- min(mat[mat > 0], na.rm = TRUE)
  log2_mat <- log2(mat + (min_nonzero / 2))
  sample_med <- matrixStats::colMedians(log2_mat, na.rm = TRUE)
  log2_norm <- sweep(log2_mat, 2, sample_med, "-")

  list(
    metabolite_table = metab_data,
    matrix = log2_norm
  )
}

subset_metabolomics_metadata <- function(metadata,
                                         celltype,
                                         target_clones = NULL) {
  ph <- metadata %>%
    filter(CellType == celltype)

  if (!is.null(target_clones)) {
    ph <- ph %>% filter(Clone %in% target_clones)
  }

  ph
}


# ------------------------------------------------------------
# QC plots
# ------------------------------------------------------------

run_metabolomics_qc_plots <- function(mat,
                                      metadata,
                                      out_dir,
                                      dataset_label,
                                      sample_col = "SampleID",
                                      color_col = "Clone",
                                      shape_col = "Stimulation",
                                      label_col = "SampleID") {
  safe_mkdir(out_dir)

  common_samples <- intersect(colnames(mat), metadata[[sample_col]])
  if (length(common_samples) == 0) {
    stop("No common samples between metabolomics matrix and metadata.")
  }

  mat <- mat[, common_samples, drop = FALSE]
  metadata <- metadata %>%
    filter(.data[[sample_col]] %in% common_samples) %>%
    slice(match(common_samples, .data[[sample_col]]))

  p <- plot_pca_qc(
    mat = mat,
    metadata = metadata,
    sample_col = sample_col,
    color_col = color_col,
    shape_col = shape_col,
    label_col = label_col,
    title = paste(dataset_label, "metabolomics QC PCA"),
    scale_features = TRUE,
    center_features = TRUE,
    point_size = 4,
    output_file = file.path(out_dir, "PCA_PC1_PC2.png"),
    width = 6,
    height = 6,
    dpi = 200
  )
  print(p)

  annotation <- metadata %>%
    select(all_of(c(sample_col, color_col, shape_col))) %>%
    as.data.frame() %>%
    `rownames<-`(NULL) %>%
    tibble::column_to_rownames(sample_col)

  grDevices::png(file.path(out_dir, "Correlation_Heatmap.png"), width = 800, height = 800)
  pheatmap::pheatmap(
    stats::cor(mat),
    annotation_col = annotation,
    annotation_row = annotation,
    main = paste(dataset_label, "metabolomics correlation")
  )
  grDevices::dev.off()

  plot_sample_distance_heatmap(
    mat = mat,
    metadata = metadata,
    sample_col = sample_col,
    annotation_cols = c(color_col, shape_col),
    title = paste(dataset_label, "metabolomics sample distances"),
    output_file = file.path(out_dir, "Sample_Distance_Heatmap.png"),
    width = 1000,
    height = 900
  )

  invisible(TRUE)
}

run_metabolomics_volcano_plots <- function(de_dir,
                                           out_dir,
                                           file_pattern = "\\.tsv$",
                                           lfc_cut = 0.58,
                                           padj_cut = 0.05) {
  if (!exists("plot_volcano_qc")) {
    stop("Generic volcano plotting helper is not available. Source 02.Rscripts/plots.R first.")
  }
  if (!dir.exists(de_dir)) {
    stop("DEM directory does not exist: ", de_dir)
  }

  safe_mkdir(out_dir)
  dem_files <- list.files(de_dir, pattern = file_pattern, full.names = TRUE)
  dem_files <- dem_files[!grepl("Summary", basename(dem_files), ignore.case = TRUE)]

  for (fpath in dem_files) {
    tag <- tools::file_path_sans_ext(basename(fpath))
    dem_tbl <- readr::read_tsv(fpath, show_col_types = FALSE)

    p <- plot_volcano_qc(
      tbl = dem_tbl,
      lfc_col = "logFC",
      p_col = "adj.P.Val",
      lfc_cut = lfc_cut,
      p_cut = padj_cut,
      title = tag,
      output_file = file.path(out_dir, paste0(tag, "_volcano.png")),
      width = 7,
      height = 5,
      dpi = 200
    )

    if (!is.null(p)) print(p)
  }

  invisible(TRUE)
}


# ------------------------------------------------------------
# DEM analysis
# ------------------------------------------------------------

format_metabolomics_contrast_label <- function(x) {
  gsub("\\.", "_", x)
}

orient_metabolomics_clone_pair <- function(pair) {
  preferred_orders <- list(
    c("1.27", "1184"),
    c("6135", "6135-H1"),
    c("6137", "6135-H1"),
    c("6135", "6137")
  )

  for (preferred in preferred_orders) {
    if (setequal(pair, preferred)) {
      return(list(case = preferred[1], control = preferred[2]))
    }
  }

  list(case = pair[2], control = pair[1])
}

build_metabolomics_contrasts <- function(metadata, design) {
  clones <- unique(metadata$Clone)
  stims <- unique(metadata$Stimulation)
  contrast_list <- list()
  contrast_info <- list()

  add_contrast <- function(name,
                           expression,
                           case,
                           control,
                           comparison_type,
                           stimulation_context = NA_character_) {
    contrast_list[[name]] <<- expression
    contrast_info[[name]] <<- data.frame(
      contrast = name,
      comparison_type = comparison_type,
      stimulation_context = stimulation_context,
      comparison_case = case,
      comparison_control = control,
      positive_logFC_higher_in = case,
      negative_logFC_higher_in = control,
      stringsAsFactors = FALSE
    )
  }

  for (clone_id in clones) {
    stim_group <- make.names(paste0(clone_id, "_Stimulated"))
    ns_group <- make.names(paste0(clone_id, "_Non_stimulated"))

    if (all(c(stim_group, ns_group) %in% colnames(design))) {
      contrast_name <- paste0(
        "Stim_vs_NS_",
        format_metabolomics_contrast_label(clone_id)
      )
      add_contrast(
        name = contrast_name,
        expression = paste(stim_group, "-", ns_group),
        case = paste(clone_id, "Stimulated"),
        control = paste(clone_id, "Non_stimulated"),
        comparison_type = "Within-clone stimulation effect",
        stimulation_context = "Stimulated_vs_Non_stimulated"
      )
    }
  }

  if (length(clones) > 1) {
    clone_pairs <- utils::combn(clones, 2, simplify = FALSE)
    for (stim_state in stims) {
      stim_label <- ifelse(stim_state == "Non_stimulated", "NS", "Stim")
      for (pair in clone_pairs) {
        ordered_pair <- orient_metabolomics_clone_pair(pair)
        group_case <- make.names(paste(ordered_pair$case, stim_state, sep = "_"))
        group_control <- make.names(paste(ordered_pair$control, stim_state, sep = "_"))

        if (all(c(group_case, group_control) %in% colnames(design))) {
          tag <- paste0(
            format_metabolomics_contrast_label(ordered_pair$case),
            "_vs_",
            format_metabolomics_contrast_label(ordered_pair$control),
            "_",
            stim_label
          )
          add_contrast(
            name = tag,
            expression = paste(group_case, "-", group_control),
            case = paste(ordered_pair$case, stim_state),
            control = paste(ordered_pair$control, stim_state),
            comparison_type = "Between-clone phenotype effect",
            stimulation_context = stim_state
          )
        }
      }
    }
  }

  attr(contrast_list, "contrast_info") <- dplyr::bind_rows(contrast_info)
  contrast_list
}

run_metabolomics_dem_by_celltype <- function(celltype,
                                             raw_data_file,
                                             out_dir,
                                             target_clones = NULL,
                                             output_exclude_pattern = NULL,
                                             outliers = character(),
                                             qual_cut = 4,
                                             lfc_cut = 0.58,
                                             padj_cut = 0.05,
                                             make_qc = TRUE,
                                             make_volcano = TRUE) {
  message("Running metabolomics DEM analysis for: ", celltype)
  safe_mkdir(out_dir)

  qc_dir <- file.path(out_dir, "QC_plots")
  de_dir <- file.path(out_dir, "DE_tables")
  volcano_dir <- file.path(out_dir, "volcano")
  safe_mkdir(qc_dir)
  safe_mkdir(de_dir)

  metadata_full <- build_metabolomics_metadata()
  metadata <- subset_metabolomics_metadata(
    metadata = metadata_full,
    celltype = celltype,
    target_clones = target_clones
  )

  raw_data <- read_metabolomics_raw_data(raw_data_file, qual_cut = qual_cut)
  prepared <- prepare_metabolomics_matrix(raw_data, metadata)
  log2_norm <- prepared$matrix
  metab_data <- prepared$metabolite_table

  if (isTRUE(make_qc)) {
    run_metabolomics_qc_plots(
      mat = log2_norm,
      metadata = metadata,
      out_dir = qc_dir,
      dataset_label = celltype
    )
  }

  if (length(outliers) > 0) {
    message("Removing metabolomics outliers: ", paste(outliers, collapse = ", "))
    log2_norm <- log2_norm[, !(colnames(log2_norm) %in% outliers), drop = FALSE]
    metadata <- metadata %>% filter(!SampleID %in% outliers)
  }

  metadata <- metadata %>%
    filter(SampleID %in% colnames(log2_norm)) %>%
    slice(match(colnames(log2_norm), SampleID)) %>%
    mutate(Group = make.names(paste(Clone, Stimulation, sep = "_")))

  design <- model.matrix(~ 0 + Group, data = metadata)
  colnames(design) <- levels(factor(metadata$Group))
  rownames(design) <- metadata$SampleID

  contrast_list <- build_metabolomics_contrasts(metadata, design)
  contrast_info <- attr(contrast_list, "contrast_info")
  if (length(contrast_list) == 0) {
    stop("No valid metabolomics contrasts could be built for: ", celltype)
  }

  contrast_matrix <- limma::makeContrasts(contrasts = unlist(contrast_list), levels = design)
  fit <- limma::lmFit(log2_norm, design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix))

  export_contrasts <- colnames(contrast_matrix)
  if (!is.null(output_exclude_pattern)) {
    export_contrasts <- export_contrasts[
      !grepl(output_exclude_pattern, export_contrasts, ignore.case = FALSE)
    ]
  }
  contrast_info <- contrast_info %>%
    dplyr::filter(.data$contrast %in% export_contrasts)

  metab_annot <- as.data.frame(
    metab_data[, .(Class, Name, Formula, `m/z`, Adduct, Polarity, Peak, `Qual.`)]
  )
  rownames(metab_annot) <- metab_annot$Name

  summary_df <- data.frame()
  for (contrast_name in export_contrasts) {
    info <- contrast_info %>%
      dplyr::filter(.data$contrast == contrast_name) %>%
      dplyr::slice_head(n = 1)

    tt <- limma::topTable(fit2, coef = contrast_name, number = Inf) %>%
      mutate(
        logFC = as.numeric(.data$logFC),
        FC = dplyr::case_when(
          .data$logFC > 0 ~ 2^.data$logFC,
          .data$logFC < 0 ~ -2^abs(.data$logFC),
          TRUE ~ 1
        ),
        significant = ifelse(abs(.data$logFC) >= lfc_cut & adj.P.Val <= padj_cut, TRUE, FALSE)
      )

    final_table <- cbind(metab_annot[rownames(tt), ], tt) %>%
      dplyr::mutate(
        contrast = contrast_name,
        comparison_type = info$comparison_type,
        stimulation_context = info$stimulation_context,
        comparison_case = info$comparison_case,
        comparison_control = info$comparison_control,
        positive_logFC_higher_in = info$positive_logFC_higher_in,
        negative_logFC_higher_in = info$negative_logFC_higher_in,
        .before = "logFC"
      ) %>%
      dplyr::relocate("FC", .after = "logFC")

    openxlsx::write.xlsx(final_table, file.path(de_dir, paste0(contrast_name, ".xlsx")))
    write.table(final_table, file.path(de_dir, paste0(contrast_name, ".tsv")),
                sep = "\t", row.names = FALSE)
    write.table(final_table, file.path(de_dir, paste0(contrast_name, ".txt")),
                sep = "\t", row.names = FALSE)

    summary_df <- rbind(summary_df, data.frame(
      file_name = contrast_name,
      comparison_type = info$comparison_type,
      stimulation_context = info$stimulation_context,
      comparison_case = info$comparison_case,
      comparison_control = info$comparison_control,
      positive_logFC_higher_in = info$positive_logFC_higher_in,
      negative_logFC_higher_in = info$negative_logFC_higher_in,
      nb_DEM_total = sum(tt$significant),
      nb_DEM_up = sum(tt$significant & tt$logFC > 0),
      nb_DEM_down = sum(tt$significant & tt$logFC < 0)
    ))
  }

  openxlsx::write.xlsx(summary_df, file.path(de_dir, "Summary_Table.xlsx"))
  write.table(summary_df, file.path(de_dir, "Summary_Table.tsv"),
              sep = "\t", row.names = FALSE)
  write.table(summary_df, file.path(de_dir, "Summary_Table.txt"),
              sep = "\t", row.names = FALSE)

  if (isTRUE(make_volcano)) {
    run_metabolomics_volcano_plots(
      de_dir = de_dir,
      out_dir = volcano_dir,
      lfc_cut = lfc_cut,
      padj_cut = padj_cut
    )
  }

  invisible(list(
    summary = summary_df,
    fit = fit2,
    normalized_matrix = log2_norm,
    metadata = metadata,
    de_dir = de_dir,
    qc_dir = qc_dir
  ))
}
