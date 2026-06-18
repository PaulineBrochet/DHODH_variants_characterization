# ============================================================
# Proteomics DEP analysis helpers
# ============================================================
#
# This file contains reusable functions for proteomics differential
# enrichment/abundance protein (DEP) table generation. Plotting, enrichment,
# and manuscript figure code are intentionally kept outside this file.
#
# Historical analysis choices are preserved for reproducibility:
# - AC16 is a single-batch analysis and uses BPCA-imputed normalized data for
#   limma, as in the old manuscript code.
# - iPSC combines several batches, uses BPCA/batch correction only for QC
#   matrices, and runs limma on the merged scale-normalized non-imputed matrix
#   with Batch as a covariate.

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

write_matrix_csv <- function(mat, file) {
  safe_mkdir(dirname(file))
  readr::write_csv(
    as.data.frame(mat) %>% tibble::rownames_to_column("Genes"),
    file
  )
}

safe_write_tsv <- function(df, file) {
  safe_mkdir(dirname(file))
  readr::write_tsv(df, file)
}

safe_write_xlsx <- function(df, file) {
  safe_mkdir(dirname(file))
  writexl::write_xlsx(df, file)
}

write_text <- function(lines, file) {
  safe_mkdir(dirname(file))
  writeLines(lines, con = file)
}


# ------------------------------------------------------------
# Raw report cleaning
# ------------------------------------------------------------

get_sample_cols <- function(report) {
  cols <- grep("\\.mzML$", colnames(report), value = TRUE)
  if (length(cols) == 0) {
    stop("No .mzML sample columns found in the proteomics report.")
  }
  cols
}

clean_ips_report_genes <- function(report) {
  report[["Genes"]] <- sub(";.*", "", as.character(report[["Genes"]]))
  if (!("Protein.Names" %in% colnames(report))) {
    report$Protein.Names <- NA_character_
  }

  report %>%
    mutate(
      Protein.Names = as.character(Protein.Names),
      Genes = as.character(Genes),
      Genes = if_else(
        is.na(Genes) | trimws(Genes) == "",
        if_else(is.na(Protein.Names), NA_character_, sub("_.*", "", Protein.Names)),
        Genes
      )
    ) %>%
    filter(!is.na(Genes) & Genes != "")
}

read_ips_raw_report <- function(file, trim_tail = 0) {
  report <- read.delim(file, sep = "\t", header = TRUE, check.names = FALSE)
  if (trim_tail > 0) {
    if (nrow(report) <= trim_tail) {
      stop("Report has fewer rows than trim_tail for: ", file)
    }
    report <- report[seq_len(nrow(report) - trim_tail), , drop = FALSE]
  }
  clean_ips_report_genes(report)
}

build_gene_matrix_from_report <- function(report, sample_cols) {
  stopifnot(all(sample_cols %in% colnames(report)))

  df <- report %>%
    select(Genes, all_of(sample_cols)) %>%
    mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
    mutate(across(all_of(sample_cols), ~ ifelse(is.na(.x) | .x <= 0, NA_real_, .x)))

  # Multiple protein groups can map to one gene symbol. The iPSC historical
  # code collapsed them by median intensity before log transformation.
  df2 <- df %>%
    group_by(Genes) %>%
    summarise(across(all_of(sample_cols), ~ median(.x, na.rm = TRUE)), .groups = "drop")

  mat <- as.matrix(df2 %>% tibble::column_to_rownames("Genes"))
  storage.mode(mat) <- "double"
  mat[!is.finite(mat)] <- NA_real_
  mat
}

clean_ac16_report <- function(report) {
  genes_crap <- c(
    "cRAP-KRT15", "cRAP-SRPP", "cRAP-ALB", "cRAP-AHSG", "cRAP-FETUB",
    "OPG085", "ALB.cRAP.ALB", "L1", "LYZ.cRAP.LYZ", "CUSTOS"
  )

  report <- report[
    !(report$Genes %in% genes_crap) &
      !(report$Genes == "") &
      !grepl("^X\\.", report$Genes),
    ,
    drop = FALSE
  ]

  # Manual duplicate fixes from the historical AC16 proteomics analysis.
  idx_tmpo <- which(report$Genes == "TMPO")
  if (length(idx_tmpo) == 2) report$Genes[idx_tmpo] <- c("TMPO2A", "TMPO2B")

  idx_polr2m <- which(report$Genes == "POLR2M")
  if (length(idx_polr2m) == 2) report$Genes[idx_polr2m] <- c("POLR2M1A", "POLR2M1AD")

  report
}


# ------------------------------------------------------------
# Matrix filtering and imputation
# ------------------------------------------------------------

filter_present_per_group <- function(mat, groups, min_rep = 1, min_groups = length(groups)) {
  missing_groups <- setdiff(unlist(groups, use.names = FALSE), colnames(mat))
  if (length(missing_groups) > 0) {
    stop("Some group samples are missing from the matrix: ",
         paste(missing_groups, collapse = ", "))
  }

  presence <- sapply(groups, function(group_samples) {
    rowSums(!is.na(mat[, group_samples, drop = FALSE])) >= min_rep
  })
  if (is.null(dim(presence))) presence <- matrix(presence, ncol = 1)

  n_groups <- rowSums(presence)
  keep <- n_groups >= min_groups

  list(
    matrix = mat[keep, , drop = FALSE],
    presence = presence[keep, , drop = FALSE],
    n_groups = n_groups[keep]
  )
}

filter_by_group_presence <- function(mat, meta, group_col = "Phenotype", min_n = 2) {
  stopifnot(all(meta$sample_uid %in% colnames(mat)))
  mat <- mat[, meta$sample_uid, drop = FALSE]

  groups <- unique(meta[[group_col]])
  counts <- sapply(groups, function(group_name) {
    cols <- meta$sample_uid[meta[[group_col]] == group_name]
    rowSums(!is.na(mat[, cols, drop = FALSE]))
  })
  if (is.null(dim(counts))) counts <- matrix(counts, ncol = 1)

  keep <- apply(counts, 1, function(x) all(x >= min_n))
  mat[keep, , drop = FALSE]
}

median_impute_matrix <- function(mat) {
  x <- as.matrix(mat)
  storage.mode(x) <- "double"
  x[!is.finite(x)] <- NA_real_

  for (i in seq_len(nrow(x))) {
    row_median <- median(x[i, ], na.rm = TRUE)
    if (!is.finite(row_median)) row_median <- 0
    x[i, is.na(x[i, ])] <- row_median
  }

  x
}

bpca_impute_safe <- function(mat, nPcs = 5, seed = 123) {
  x <- as.matrix(mat)
  storage.mode(x) <- "double"
  x[!is.finite(x)] <- NA_real_

  keep1 <- rowSums(!is.na(x)) > 0
  x <- x[keep1, , drop = FALSE]
  row_var <- apply(x, 1, function(v) var(v, na.rm = TRUE))
  keep2 <- is.finite(row_var) & row_var > 0
  x <- x[keep2, , drop = FALSE]

  max_npcs <- max(1, min(ncol(x) - 1, nrow(x) - 1))
  nPcs_use <- min(nPcs, max_npcs)

  set.seed(seed)
  obj <- pcaMethods::pca(x, method = "bpca", nPcs = nPcs_use, scale = "none")
  pcaMethods::completeObs(obj)
}


# ------------------------------------------------------------
# iPSC metadata helpers
# ------------------------------------------------------------

make_ips_metadata_240 <- function(n = 30) {
  idx <- seq_len(n)
  genotype <- dplyr::case_when(
    idx %in% 1:10 ~ "6135_H1",
    idx %in% 11:20 ~ "6137",
    idx %in% 21:30 ~ "6135_CCC",
    TRUE ~ NA_character_
  )
  stim <- dplyr::case_when(
    idx %in% 1:5 ~ "NS",
    idx %in% 6:10 ~ "S",
    idx %in% 11:15 ~ "NS",
    idx %in% 16:20 ~ "S",
    idx %in% 21:25 ~ "NS",
    idx %in% 26:30 ~ "S",
    TRUE ~ NA_character_
  )

  tibble(
    sample_index = idx,
    sample_uid = paste0("240_", idx),
    Sample_ID = as.character(idx),
    Replicate = rep(paste0("R", 1:5), times = 6),
    Genotype = genotype,
    Stimulation = stim,
    Phenotype = paste0(genotype, "_", stim),
    Batch = "240_Nov2025",
    Date = "Nov 2025",
    Spectrometer = "Exploris240"
  )
}

make_ips_metadata_480 <- function(n = 40, exclude_lines = "6132") {
  idx <- seq_len(n)
  batch <- ifelse(idx <= 20, "480_Oct2025", "480_Dec2025")
  date <- ifelse(idx <= 20, "Oct 2025", "Dec 2025")

  genotype <- dplyr::case_when(
    idx %in% 1:10 ~ "6137",
    idx %in% 11:20 ~ "6135_CCC",
    idx %in% 21:30 ~ "6135_H1",
    idx %in% 31:40 ~ "6132",
    TRUE ~ NA_character_
  )
  stim <- dplyr::case_when(
    idx %in% 1:5 ~ "NS",
    idx %in% 6:10 ~ "S",
    idx %in% 11:15 ~ "NS",
    idx %in% 16:20 ~ "S",
    idx %in% 21:25 ~ "NS",
    idx %in% 26:30 ~ "S",
    idx %in% 31:35 ~ "S",
    idx %in% 36:40 ~ "NS",
    TRUE ~ NA_character_
  )

  tibble(
    sample_index = idx,
    sample_uid = paste0("480_", idx),
    Sample_ID = as.character(idx),
    Replicate = rep(paste0("R", 1:5), times = 8),
    Genotype = genotype,
    Stimulation = stim,
    Phenotype = paste0(genotype, "_", stim),
    Batch = batch,
    Date = date,
    Spectrometer = "Exploris480"
  ) %>%
    filter(!Genotype %in% exclude_lines)
}


# ------------------------------------------------------------
# Shared limma DEP exporter
# ------------------------------------------------------------

run_proteomics_limma_dep <- function(mat,
                                     design,
                                     contrast_defs,
                                     annotation_tbl = NULL,
                                     out_dir,
                                     tables_dir = file.path(out_dir, "tables"),
                                     lfc_cut = 0.58,
                                     padj_cut = 0.05,
                                     sort_by = "P",
                                     lmfit_method = "ls",
                                     pre_contrast_ebayes = FALSE,
                                     post_contrast_trend = TRUE,
                                     fc_mode = c("signed_ratio", "signed_power"),
                                     write_summary = TRUE,
                                     summary_control_case = NULL,
                                     file_prefix = "",
                                     csv_output = FALSE,
                                     xlsx_suffix = ".xlsx",
                                     significant_col = "Significant") {
  fc_mode <- match.arg(fc_mode)
  safe_mkdir(out_dir)
  safe_mkdir(tables_dir)

  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"

  if (!all(colnames(mat) %in% rownames(design))) {
    rownames(design) <- colnames(mat)
  }
  design <- design[colnames(mat), , drop = FALSE]

  fit <- limma::lmFit(mat, design, method = lmfit_method)
  if (isTRUE(pre_contrast_ebayes)) {
    fit <- limma::eBayes(fit, trend = TRUE)
  }

  contrast_matrix <- limma::makeContrasts(contrasts = contrast_defs, levels = design)
  colnames(contrast_matrix) <- names(contrast_defs)

  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit3 <- limma::eBayes(fit2, trend = post_contrast_trend)

  results <- list()
  summary_rows <- list()

  for (tag in names(contrast_defs)) {
    message("Running proteomics DEP contrast: ", tag)

    res <- limma::topTable(fit3, coef = tag, number = Inf, sort.by = sort_by) %>%
      tibble::rownames_to_column("Genes")

    if (!is.null(annotation_tbl)) {
      res <- res %>% left_join(annotation_tbl, by = "Genes")
      front_cols <- intersect(c("Genes", "Protein.Group", "First.Protein.Description"), names(res))
      res <- res %>% select(all_of(front_cols), everything())
    }

    res <- res %>%
      mutate(
        logFC = as.numeric(logFC),
        P.Value = as.numeric(P.Value),
        adj.P.Val = as.numeric(adj.P.Val),
        FC = dplyr::case_when(
          fc_mode == "signed_power" ~ sign(logFC) * 2^abs(logFC),
          logFC >= 0 ~ 2^logFC,
          TRUE ~ -1 / (2^logFC)
        )
      )

    res[[significant_col]] <- ifelse(
      abs(res$logFC) >= lfc_cut & res$adj.P.Val <= padj_cut,
      "Yes",
      "No"
    )

    out_base <- file.path(tables_dir, paste0(file_prefix, tag))
    if (isTRUE(csv_output)) {
      write.csv(res, file = paste0(out_base, ".csv"), row.names = FALSE)
    } else {
      safe_write_tsv(res, paste0(out_base, ".tsv"))
    }
    safe_write_xlsx(res, paste0(out_base, xlsx_suffix))

    dep_tbl <- res %>%
      filter(
        is.finite(adj.P.Val),
        adj.P.Val <= padj_cut,
        is.finite(logFC),
        abs(logFC) >= lfc_cut
      )

    cc <- if (is.null(summary_control_case)) {
      c(Control = NA_character_, Case = NA_character_)
    } else {
      summary_control_case(tag)
    }

    summary_rows[[tag]] <- tibble(
      Comparison = tag,
      Control = unname(cc["Control"]),
      Case = unname(cc["Case"]),
      Nb_DEP = nrow(dep_tbl),
      Nb_DEP_UP = sum(dep_tbl$logFC > 0, na.rm = TRUE),
      Nb_DEP_DOWN = sum(dep_tbl$logFC < 0, na.rm = TRUE)
    )

    results[[tag]] <- res
  }

  summary_tbl <- bind_rows(summary_rows)
  if (isTRUE(write_summary)) {
    safe_write_tsv(summary_tbl, file.path(out_dir, "Summary_DEP.tsv"))
    safe_write_xlsx(summary_tbl, file.path(out_dir, "Summary_DEP.xlsx"))
  }

  invisible(list(
    fit = fit3,
    contrast_matrix = contrast_matrix,
    results = results,
    summary = summary_tbl
  ))
}


# ------------------------------------------------------------
# AC16 DEP wrapper
# ------------------------------------------------------------

run_ac16_proteomics_dep <- function(report_file,
                                    out_dir,
                                    sample_groups,
                                    contrast_defs,
                                    normalized_matrix_file = "normalized_matrix_log2.csv",
                                    min_rep_per_group = 1,
                                    min_groups_required = length(sample_groups),
                                    nPcs_impute = 5,
                                    seed = 123,
                                    lfc_cut = 0.58,
                                    padj_cut = 0.05) {
  safe_mkdir(out_dir)

  message("Loading AC16 proteomics report: ", report_file)
  report <- read.delim(report_file, sep = "\t", header = TRUE)
  report <- clean_ac16_report(report)

  abundance_cols <- grep("julia_.*\\.mzML", colnames(report), value = TRUE)
  if (length(abundance_cols) == 0) {
    stop("No AC16 abundance columns matching 'julia_*.mzML' were found.")
  }

  matrix_abund <- dplyr::select(report, dplyr::all_of(abundance_cols))
  matrix_abund <- apply(matrix_abund, 2, function(x) as.numeric(gsub(",", ".", x)))
  rownames(matrix_abund) <- report$Genes

  matrix_abund[matrix_abund <= 0] <- NA_real_
  matrix_abund_log2 <- log2(matrix_abund)

  colnames(matrix_abund_log2) <- colnames(matrix_abund_log2) %>%
    basename() %>%
    stringr::str_extract("julia_\\d+[A-Z]?") %>%
    stringr::str_replace("julia_", "")

  configured_samples <- unlist(sample_groups, use.names = FALSE)
  missing_samples <- setdiff(configured_samples, colnames(matrix_abund_log2))
  if (length(missing_samples) > 0) {
    stop("Some configured AC16 samples are missing from the raw matrix: ",
         paste(missing_samples, collapse = ", "))
  }
  # Preserve raw report column order after excluding the high-dose samples.
  # BPCA imputation can be order-sensitive, so this reproduces the historical
  # no-50.50 AC16 normalized matrix and downstream DEP tables.
  selected_samples <- colnames(matrix_abund_log2)[colnames(matrix_abund_log2) %in% configured_samples]
  matrix_abund_log2 <- matrix_abund_log2[, selected_samples, drop = FALSE]

  message("AC16 missing values before filtering: ",
          round(sum(is.na(matrix_abund_log2)) / length(matrix_abund_log2) * 100, 2), "%")

  filtered <- filter_present_per_group(
    matrix_abund_log2,
    sample_groups,
    min_rep = min_rep_per_group,
    min_groups = min_groups_required
  )
  filtered_proteins <- filtered$matrix

  message("AC16 proteins retained after group filter: ", nrow(filtered_proteins))
  message("AC16 missing values after filtering: ",
          round(sum(is.na(filtered_proteins)) / length(filtered_proteins) * 100, 2), "%")

  matrix_norm_global <- limma::normalizeBetweenArrays(filtered_proteins, method = "scale")
  write_matrix_csv(matrix_norm_global, file.path(out_dir, normalized_matrix_file))

  # The historical AC16 script ran a BPCA diagnostic before the final BPCA
  # object. Keeping that call avoids changing RNG-sensitive behavior.
  set.seed(seed)
  invisible(pcaMethods::pca(matrix_norm_global, method = "bpca", nPcs = 10, scale = "none"))
  imputed_obj <- pcaMethods::pca(matrix_norm_global, method = "bpca", nPcs = nPcs_impute, scale = "none")
  imputed_bpca <- pcaMethods::completeObs(imputed_obj)

  sample_to_condition <- unlist(
    lapply(names(sample_groups), function(condition) {
      stats::setNames(rep(condition, length(sample_groups[[condition]])), sample_groups[[condition]])
    })
  )
  condition <- unname(sample_to_condition[colnames(imputed_bpca)])
  if (any(is.na(condition))) {
    stop("Some AC16 samples could not be mapped to a condition: ",
         paste(colnames(imputed_bpca)[is.na(condition)], collapse = ", "))
  }

  design <- model.matrix(~ 0 + factor(condition))
  colnames(design) <- make.names(levels(factor(condition)))
  rownames(design) <- colnames(imputed_bpca)

  annotation_tbl <- report %>%
    select(any_of(c("Genes", "Protein.Group", "Protein.Names", "First.Protein.Description")))

  run_proteomics_limma_dep(
    mat = imputed_bpca,
    design = design,
    contrast_defs = contrast_defs,
    annotation_tbl = annotation_tbl,
    out_dir = out_dir,
    tables_dir = out_dir,
    lfc_cut = lfc_cut,
    padj_cut = padj_cut,
    sort_by = "B",
    lmfit_method = "robust",
    pre_contrast_ebayes = TRUE,
    post_contrast_trend = FALSE,
    fc_mode = "signed_power",
    write_summary = TRUE,
    summary_control_case = infer_ac16_control_case,
    file_prefix = "limma_",
    csv_output = TRUE,
    xlsx_suffix = "_excel.xlsx",
    significant_col = "significant"
  )
}

infer_ac16_control_case <- function(tag) {
  if (grepl("^[0-9]+_.*_vs_NS$", tag)) {
    parts <- strsplit(tag, "_")[[1]]
    clone <- parts[1]
    dose <- sub(paste0("^", clone, "_"), "", sub("_vs_NS$", "", tag))
    return(c(Control = paste0(clone, "_NS"), Case = paste0(clone, "_", dose)))
  }
  if (grepl("^127_.*_vs_1184_.*$", tag)) {
    parts <- strsplit(tag, "_vs_")[[1]]
    return(c(Control = parts[2], Case = parts[1]))
  }
  c(Control = NA_character_, Case = NA_character_)
}


# ------------------------------------------------------------
# iPSC preprocessing and DEP wrappers
# ------------------------------------------------------------

run_ips_proteomics_preprocessing <- function(report240_file,
                                             report480_file,
                                             out_base,
                                             trim_tail_240 = 3,
                                             trim_tail_480 = 7,
                                             min_n_per_group = 1,
                                             nPcs_impute = 5,
                                             exclude_lines = "6132") {
  safe_mkdir(out_base)
  mat_dir <- file.path(out_base, "matrices")
  meta_dir <- file.path(out_base, "metadata")
  safe_mkdir(mat_dir)
  safe_mkdir(meta_dir)

  message("Reading iPSC proteomics raw reports.")
  rep240 <- read_ips_raw_report(report240_file, trim_tail = trim_tail_240)
  rep480 <- read_ips_raw_report(report480_file, trim_tail = trim_tail_480)

  mz240 <- get_sample_cols(rep240)
  mz480 <- get_sample_cols(rep480)

  if (length(mz240) < 30) stop("Expected at least 30 mzML columns in 240 report.")
  if (length(mz480) < 40) stop("Expected at least 40 mzML columns in 480 report.")

  mz240 <- mz240[1:30]
  mz480 <- mz480[1:40]

  meta240 <- make_ips_metadata_240(30) %>% mutate(mzml_col = mz240)
  meta480 <- make_ips_metadata_480(40, exclude_lines = exclude_lines) %>%
    mutate(mzml_col = mz480[sample_index])

  meta_oct <- meta480 %>% filter(Batch == "480_Oct2025")
  meta_dec <- meta480 %>% filter(Batch == "480_Dec2025")

  write_csv(meta240, file.path(meta_dir, "metadata_240_Nov2025.csv"))
  write_csv(meta_oct, file.path(meta_dir, "metadata_480_Oct2025.csv"))
  write_csv(meta_dec, file.path(meta_dir, "metadata_480_Dec2025.csv"))
  write_csv(bind_rows(meta240, meta480), file.path(meta_dir, "metadata_all_from_raw.csv"))

  m240_lin <- build_gene_matrix_from_report(rep240, sample_cols = meta240$mzml_col)
  m480_lin <- build_gene_matrix_from_report(rep480, sample_cols = meta480$mzml_col)

  colnames(m240_lin) <- meta240$sample_uid
  colnames(m480_lin) <- meta480$sample_uid

  mOct_lin <- m480_lin[, meta_oct$sample_uid, drop = FALSE]
  mDec_lin <- m480_lin[, meta_dec$sample_uid, drop = FALSE]

  m240 <- log2(m240_lin + 1)
  mOct <- log2(mOct_lin + 1)
  mDec <- log2(mDec_lin + 1)
  m240[!is.finite(m240)] <- NA_real_
  mOct[!is.finite(mOct)] <- NA_real_
  mDec[!is.finite(mDec)] <- NA_real_

  shared_all3 <- Reduce(intersect, list(rownames(m240), rownames(mOct), rownames(mDec)))
  message("iPSC shared genes across all three batches: ", length(shared_all3))

  m240_s <- m240[shared_all3, , drop = FALSE]
  mOct_s <- mOct[shared_all3, , drop = FALSE]
  mDec_s <- mDec[shared_all3, , drop = FALSE]

  m240_k <- filter_by_group_presence(m240_s, meta240, "Phenotype", min_n_per_group)
  mOct_k <- filter_by_group_presence(mOct_s, meta_oct, "Phenotype", min_n_per_group)
  mDec_k <- filter_by_group_presence(mDec_s, meta_dec, "Phenotype", min_n_per_group)

  keep <- Reduce(intersect, list(rownames(m240_k), rownames(mOct_k), rownames(mDec_k)))
  message("iPSC proteins retained after batch/phenotype filter: ", length(keep))

  merged_raw <- cbind(
    m240_s[keep, , drop = FALSE],
    mOct_s[keep, , drop = FALSE],
    mDec_s[keep, , drop = FALSE]
  )

  meta_all <- bind_rows(meta240, meta_oct, meta_dec) %>%
    mutate(
      Batch = factor(Batch, levels = c("240_Nov2025", "480_Oct2025", "480_Dec2025")),
      Phenotype = factor(Phenotype)
    ) %>%
    arrange(Batch, Phenotype, as.integer(Sample_ID))

  merged_raw <- merged_raw[, meta_all$sample_uid, drop = FALSE]

  write_matrix_csv(merged_raw, file.path(mat_dir, "matrix_log2_merged_raw_sharedFiltered.csv"))
  write_csv(meta_all, file.path(meta_dir, "metadata_merged_3batches.csv"))

  merged_norm <- limma::normalizeBetweenArrays(merged_raw, method = "scale")
  write_matrix_csv(merged_norm, file.path(mat_dir, "matrix_log2_merged_scaleNormalized.csv"))

  # These imputed and batch-corrected matrices are retained for QC continuity,
  # but they are not used by run_ips_proteomics_dep().
  merged_imp <- bpca_impute_safe(merged_norm, nPcs = nPcs_impute, seed = 123)
  merged_imp <- merged_imp[, meta_all$sample_uid, drop = FALSE]
  write_matrix_csv(merged_imp, file.path(mat_dir, "matrix_log2_merged_norm_BPCAimputed_QConly.csv"))

  design <- model.matrix(~ 0 + Phenotype, data = meta_all)
  merged_bc <- limma::removeBatchEffect(merged_imp, batch = meta_all$Batch, design = design)
  write_matrix_csv(merged_bc, file.path(mat_dir, "matrix_log2_merged_norm_imputed_batchCorrected_QConly.csv"))

  invisible(list(
    metadata = meta_all,
    merged_raw = merged_raw,
    merged_norm = merged_norm,
    merged_imputed_QConly = merged_imp,
    merged_batchCorrected_QConly = merged_bc,
    matrix_file = file.path(mat_dir, "matrix_log2_merged_scaleNormalized.csv"),
    metadata_file = file.path(meta_dir, "metadata_merged_3batches.csv")
  ))
}

run_ips_proteomics_dep <- function(matrix_file,
                                   metadata_file,
                                   report_files,
                                   out_dir,
                                   condition_levels,
                                   contrast_defs,
                                   exclude_lines = "6132",
                                   lfc_cut = 0.58,
                                   padj_cut = 0.05) {
  tables_dir <- file.path(out_dir, "tables")
  qc_dir <- file.path(out_dir, "QC")
  safe_mkdir(tables_dir)
  safe_mkdir(qc_dir)

  message("Loading iPSC proteomics normalized matrix: ", matrix_file)
  mat_df <- readr::read_csv(matrix_file, show_col_types = FALSE)
  if (!("Genes" %in% colnames(mat_df))) {
    stop("Expected a 'Genes' column in: ", matrix_file)
  }

  mat <- mat_df %>%
    tibble::column_to_rownames("Genes") %>%
    as.matrix()
  storage.mode(mat) <- "double"

  sample_info <- readr::read_csv(metadata_file, show_col_types = FALSE) %>%
    mutate(
      sample_uid = as.character(sample_uid),
      Batch = as.character(Batch),
      Genotype = as.character(Genotype),
      Stimulation = as.character(Stimulation),
      Phenotype = paste0(Genotype, "_", Stimulation)
    ) %>%
    filter(!Genotype %in% exclude_lines)

  common_samples <- intersect(colnames(mat), sample_info$sample_uid)
  if (length(common_samples) == 0) {
    stop("No common samples between the iPSC matrix and metadata after filtering.")
  }

  mat <- mat[, common_samples, drop = FALSE]
  sample_info <- sample_info %>% slice(match(common_samples, sample_uid))
  sample_info$Batch <- factor(sample_info$Batch)

  unknown_levels <- setdiff(unique(sample_info$Phenotype), condition_levels)
  if (length(unknown_levels) > 0) {
    stop("Unexpected iPSC proteomics condition labels: ",
         paste(unknown_levels, collapse = ", "))
  }

  condition_factor <- factor(sample_info$Phenotype, levels = condition_levels)
  design <- model.matrix(~ 0 + condition_factor + Batch, data = sample_info)
  colnames(design) <- make.names(colnames(design))
  rownames(design) <- sample_info$sample_uid

  qr_rank <- qr(design)$rank
  if (qr_rank < ncol(design)) {
    write_text(
      c(
        "WARNING: design matrix is not full rank.",
        paste0("rank = ", qr_rank, " < ncol = ", ncol(design)),
        "This can happen if some conditions are absent in some batches."
      ),
      file.path(qc_dir, "QC_design_rank_warning.txt")
    )
  }

  annotation_tbl <- load_ips_protein_annotations(report_files)

  dep <- run_proteomics_limma_dep(
    mat = mat,
    design = design,
    contrast_defs = contrast_defs,
    annotation_tbl = annotation_tbl,
    out_dir = out_dir,
    tables_dir = tables_dir,
    lfc_cut = lfc_cut,
    padj_cut = padj_cut,
    sort_by = "P",
    lmfit_method = "ls",
    pre_contrast_ebayes = TRUE,
    post_contrast_trend = TRUE,
    fc_mode = "signed_ratio",
    write_summary = TRUE,
    summary_control_case = infer_ips_control_case,
    file_prefix = "",
    csv_output = FALSE,
    xlsx_suffix = ".xlsx",
    significant_col = "Significant"
  )

  qc_design <- tibble(
    sample_uid = sample_info$sample_uid,
    Batch = sample_info$Batch,
    Phenotype_raw = sample_info$Phenotype,
    condition_used = as.character(condition_factor)
  )
  safe_write_tsv(qc_design, file.path(qc_dir, "DEP_sample_annotations.tsv"))

  invisible(dep)
}

load_ips_protein_annotations <- function(report_files) {
  reports <- lapply(report_files, function(file) {
    read.delim(file, sep = "\t", header = TRUE, check.names = FALSE) %>%
      clean_ips_report_genes()
  })

  meta_cols <- unique(unlist(lapply(reports, colnames)))
  meta_cols <- intersect(
    c("Genes", "Protein.Group", "Protein.Names", "First.Protein.Description"),
    meta_cols
  )

  bind_rows(lapply(reports, function(report) select(report, any_of(meta_cols)))) %>%
    distinct(Genes, .keep_all = TRUE)
}

infer_ips_control_case <- function(tag) {
  if (grepl("^[0-9A-Za-z]+_S_vs_NS$", tag)) {
    clone <- sub("_S_vs_NS$", "", tag)
    return(c(Control = paste0(clone, "_NS"), Case = paste0(clone, "_S")))
  }
  if (grepl("^NS_[0-9A-Za-z]+_vs_[0-9A-Za-z]+$", tag)) {
    parts <- strsplit(sub("^NS_", "", tag), "_vs_")[[1]]
    return(c(Control = paste0(parts[2], "_NS"), Case = paste0(parts[1], "_NS")))
  }
  if (grepl("^S_[0-9A-Za-z]+_vs_[0-9A-Za-z]+$", tag)) {
    parts <- strsplit(sub("^S_", "", tag), "_vs_")[[1]]
    return(c(Control = paste0(parts[2], "_S"), Case = paste0(parts[1], "_S")))
  }
  c(Control = NA_character_, Case = NA_character_)
}


# ------------------------------------------------------------
# Proteomics QC and volcano plot wrappers
# ------------------------------------------------------------

read_proteomics_matrix <- function(matrix_file) {
  mat_df <- readr::read_csv(matrix_file, show_col_types = FALSE)
  if (!("Genes" %in% colnames(mat_df))) {
    stop("Expected a 'Genes' column in: ", matrix_file)
  }

  mat <- mat_df %>%
    tibble::column_to_rownames("Genes") %>%
    as.matrix()
  storage.mode(mat) <- "double"
  mat
}

build_ac16_proteomics_metadata <- function(sample_groups) {
  bind_rows(lapply(names(sample_groups), function(condition_name) {
    parts <- strsplit(condition_name, "_", fixed = TRUE)[[1]]
    tibble(
      sample = sample_groups[[condition_name]],
      clone = parts[1],
      condition = parts[2],
      phenotype = condition_name
    )
  }))
}

read_ips_proteomics_metadata <- function(metadata_file, exclude_lines = "6132") {
  readr::read_csv(metadata_file, show_col_types = FALSE) %>%
    mutate(
      sample_uid = as.character(sample_uid),
      Batch = as.character(Batch),
      Genotype = as.character(Genotype),
      Stimulation = as.character(Stimulation),
      Phenotype = paste0(Genotype, "_", Stimulation)
    ) %>%
    filter(!Genotype %in% exclude_lines)
}

run_proteomics_matrix_qc_plots <- function(matrix_file,
                                           metadata,
                                           sample_col,
                                           color_col,
                                           shape_col = NULL,
                                           label_col = NULL,
                                           out_dir,
                                           dataset_label,
                                           annotation_cols = NULL,
                                           pca_scale_features = TRUE) {
  if (!exists("plot_pca_qc") || !exists("plot_sample_distance_heatmap")) {
    stop("Generic plotting helpers are not available. Source 02.Rscripts/plots.R first.")
  }

  safe_mkdir(out_dir)
  mat <- read_proteomics_matrix(matrix_file)

  common_samples <- intersect(colnames(mat), metadata[[sample_col]])
  if (length(common_samples) == 0) {
    stop("No common samples between proteomics matrix and metadata.")
  }

  mat <- mat[, common_samples, drop = FALSE]
  metadata <- metadata %>%
    filter(.data[[sample_col]] %in% common_samples) %>%
    slice(match(common_samples, .data[[sample_col]]))

  mat_complete <- median_impute_matrix(mat)

  p <- plot_pca_qc(
    mat = mat_complete,
    metadata = metadata,
    sample_col = sample_col,
    color_col = color_col,
    shape_col = shape_col,
    label_col = label_col,
    title = paste0(dataset_label, " - proteomics PCA"),
    scale_features = pca_scale_features,
    center_features = TRUE,
    output_file = file.path(out_dir, paste0("PCA_", dataset_label, "_proteomics.png")),
    width = 7,
    height = 5,
    dpi = 200
  )
  print(p)

  plot_sample_distance_heatmap(
    mat = mat_complete,
    metadata = metadata,
    sample_col = sample_col,
    annotation_cols = annotation_cols,
    title = paste0(dataset_label, " - proteomics sample distances"),
    output_file = file.path(out_dir, paste0("sample_distance_", dataset_label, "_proteomics.png")),
    width = 1000,
    height = 900
  )

  invisible(list(matrix = mat, matrix_complete = mat_complete, metadata = metadata, pca = p))
}

read_dep_table <- function(file) {
  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    readr::read_csv(file, show_col_types = FALSE)
  } else {
    readr::read_tsv(file, show_col_types = FALSE)
  }
}

run_proteomics_volcano_plots <- function(dep_dir,
                                         out_dir,
                                         file_pattern = "\\.(tsv|csv)$",
                                         lfc_cut = 0.58,
                                         padj_cut = 0.05,
                                         lfc_col = "logFC",
                                         padj_col = "adj.P.Val",
                                         skip_pattern = "(Summary|summary)") {
  if (!exists("plot_volcano_qc")) {
    stop("Generic volcano plotting helper is not available. Source 02.Rscripts/plots.R first.")
  }
  if (!dir.exists(dep_dir)) {
    stop("DEP directory does not exist: ", dep_dir)
  }

  safe_mkdir(out_dir)
  dep_files <- list.files(dep_dir, pattern = file_pattern, full.names = TRUE)
  dep_files <- dep_files[!grepl(skip_pattern, basename(dep_files))]

  for (fpath in dep_files) {
    tag <- tools::file_path_sans_ext(basename(fpath))
    dep_tbl <- read_dep_table(fpath)

    p <- plot_volcano_qc(
      tbl = dep_tbl,
      lfc_col = lfc_col,
      p_col = padj_col,
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
