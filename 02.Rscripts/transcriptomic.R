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
