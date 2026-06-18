# ============================================================
# Enrichment helper functions
# ============================================================
#
# This script contains reusable enrichment functions for omics result tables.
# The top-level functions are intentionally generic:
# - run_universal_gsea(): Reactome GSEA from a ranked fold-change list.
# - run_universal_gprofiler(): Reactome over-representation via g:Profiler.
# - run_enrichment_for_deg_directory(): batch runner over DEG/DEP TSV files.
#
# Project-specific directories and thresholds should be defined in config.R.

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
# Reactome gene-set dictionaries
# ------------------------------------------------------------
# Pre-load Reactome gene sets for Homo sapiens. GSEA uses TERM2GENE tables,
# either with HGNC symbols or Ensembl gene IDs depending on the input table.
reactome_hgnc <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>%
  dplyr::select(term = gs_name, gene = gene_symbol) %>% distinct()

reactome_ensembl <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>%
  dplyr::select(term = gs_name, gene = ensembl_gene) %>% distinct()

# ------------------------------------------------------------
# Reactome GSEA
# ------------------------------------------------------------
# Run Reactome GSEA from a result table containing gene identifiers and a
# fold-change column. The full GSEA table is saved, and a dotplot is generated
# for significant positive/negative pathways when available.
run_universal_gsea <- function(df, tag, out_dir, 
                               gene_col = "SYMBOL", 
                               fc_col = "log2FoldChange",
                               gene_type = "HGNC") {
  
  message("--- Running Reactome GSEA: ", tag, " [Type: ", gene_type, "]")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Select the identifier dictionary matching the input DEG/DEP table.
  selected_t2g <- if(toupper(gene_type) == "HGNC") reactome_hgnc else reactome_ensembl
  
  # Build a ranked list. Duplicate identifiers are collapsed by mean FC.
  ranked_tbl <- df %>%
    filter(!is.na(.data[[gene_col]]), !is.na(.data[[fc_col]])) %>%
    group_by(.data[[gene_col]]) %>%
    summarise(fc = mean(.data[[fc_col]]), .groups = "drop") %>%
    arrange(desc(fc))
  
  geneList <- ranked_tbl$fc
  names(geneList) <- ranked_tbl[[gene_col]]
  
  # Run GSEA and keep all pathways in the output table.
  gsea_res <- tryCatch({
    clusterProfiler::GSEA(geneList = geneList, 
                          TERM2GENE = selected_t2g,
                          pvalueCutoff = 1,
                          verbose = FALSE)
  }, error = function(e) { 
    message("  ! GSEA Error: ", e$message)
    return(NULL) 
  })
  
  if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) {
    message("  ! No GSEA results found.")
    return(NULL)
  }
  
  # Save the complete result table.
  write_tsv(as.data.frame(gsea_res), file.path(out_dir, paste0(tag, "_GSEA_Reactome.tsv")))
  writexl::write_xlsx(as.data.frame(gsea_res), file.path(out_dir, paste0(tag, "_GSEA_Reactome.xlsx")))

  # Plot the top significant up/down pathways when available.
  plot_df <- as.data.frame(gsea_res) %>%
    filter(p.adjust <= 0.05) %>%
    group_by(sign(NES)) %>%
    slice_max(order_by = abs(NES), n = 10) %>%
    ungroup() %>%
    mutate(Description = str_to_title(str_replace_all(str_remove(Description, "REACTOME_"), "_", " ")),
           Description = reorder(Description, NES))

  if(nrow(plot_df) > 0) {
    p <- ggplot(plot_df, aes(x = NES, y = Description, size = -log10(p.adjust), color = NES)) +
      geom_point() +
      scale_color_gradient2(low = "blue3", mid = "white", high = "red3", midpoint = 0) +
      labs(title = paste("Reactome GSEA:", tag), y = NULL, x = "Normalized Enrichment Score (NES)") +
      theme_bw() +
      theme(axis.text.y = element_text(size = 9))
    
    ggsave(file.path(out_dir, paste0(tag, "_GSEA_dotplot.png")), p, width = 11, height = 7)
  }
  
  return(gsea_res)
}


# ------------------------------------------------------------
# g:Profiler Reactome enrichment
# ------------------------------------------------------------
# Run Reactome over-representation analysis separately for up- and
# down-regulated genes using g:Profiler. Results and a signed barplot are saved
# when significant terms are returned.
run_universal_gprofiler <- function(df, tag, out_dir, 
                                    gene_col = "SYMBOL", 
                                    fc_col = "log2FoldChange", 
                                    padj_col = "padj",
                                    lfc_cut = 1.5, 
                                    padj_cut = 0.05,
                                    gene_type = "HGNC") {
  
  message("--- Running g:Profiler: ", tag, " [Type: ", gene_type, "]")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Select significant up- and down-regulated genes.
  up_genes <- df %>% filter(.data[[padj_col]] <= padj_cut, .data[[fc_col]] >= lfc_cut) %>% pull(.data[[gene_col]])
  dn_genes <- df %>% filter(.data[[padj_col]] <= padj_cut, .data[[fc_col]] <= -lfc_cut) %>% pull(.data[[gene_col]])
  
  # Query helper. g:Profiler is skipped for very small gene lists.
  query_gp <- function(genes, direction) {
    if(length(genes) < 3) return(NULL)
    res <- gprofiler2::gost(query = genes, 
                            organism = "hsapiens", 
                            sources = "REAC", 
                            evcodes = TRUE)$result
    if(!is.null(res)) res$direction <- direction
    return(res)
  }
  
  # Run up/down Reactome queries.
  res_up <- query_gp(up_genes, "Up")
  res_dn <- query_gp(dn_genes, "Down")
  all_res <- bind_rows(res_up, res_dn)
  
  if(!is.null(all_res) && nrow(all_res) > 0) {
    # Save result tables.
    write_tsv(all_res, file.path(out_dir, paste0(tag, "_gProfiler_REAC.tsv")))
    writexl::write_xlsx(all_res, file.path(out_dir, paste0(tag, "_gProfiler_REAC.xlsx")))

    # Build a signed summary barplot.
    plot_df <- all_res %>%
      group_by(direction) %>%
      slice_min(order_by = p_value, n = 10) %>%
      ungroup() %>%
      mutate(logP = -log10(p_value),
             val = ifelse(direction == "Up", logP, -logP),
             term_name = reorder(term_name, val))
    
    p <- ggplot(plot_df, aes(x = val, y = term_name, fill = direction)) +
      geom_col() +
      scale_fill_manual(values = c("Up" = "red2", "Down" = "blue2")) +
      geom_vline(xintercept = 0, color = "black") +
      labs(title = paste("g:Profiler Reactome:", tag),
           x = "-log10(p-value) [Negative for Down-regulated]", y = NULL) +
      theme_bw()
    
    ggsave(file.path(out_dir, paste0(tag, "_gProfiler_barplot.png")), p, width = 10, height = 6)
  } else {
    message("  ! No significant g:Profiler enrichment found.")
  }
  
  return(all_res)
}


# ------------------------------------------------------------
# Batch enrichment runner for DEG/DEP directories
# ------------------------------------------------------------
# Iterate over result tables in a directory and run both enrichment
# methods. This is the generic replacement for repeated notebook loops.
read_enrichment_input_table <- function(file) {
  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    readr::read_csv(file, show_col_types = FALSE)
  } else {
    readr::read_tsv(file, show_col_types = FALSE)
  }
}

run_enrichment_for_deg_directory <- function(deg_dir,
                                             out_dir,
                                             gene_col = "Gene_name",
                                             fc_col = "log2FoldChange",
                                             padj_col = "padj",
                                             lfc_cut = 1.5,
                                             padj_cut = 0.05,
                                             gene_type = "HGNC",
                                             file_pattern = "\\.tsv$",
                                             skip_pattern = "(summary|resume)",
                                             run_gsea = TRUE,
                                             run_gprofiler = TRUE) {
  if (!dir.exists(deg_dir)) {
    stop("DEG directory does not exist: ", deg_dir)
  }

  result_files <- list.files(deg_dir, pattern = file_pattern, full.names = TRUE)
  if (!is.null(skip_pattern)) {
    result_files <- result_files[!grepl(skip_pattern, basename(result_files), ignore.case = TRUE)]
  }

  if (length(result_files) == 0) {
    message("No DEG/DEP result files found in: ", deg_dir)
    return(invisible(NULL))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (fpath in result_files) {
    comparison_tag <- tools::file_path_sans_ext(basename(fpath))
    message("\nProcessing enrichment for: ", comparison_tag)

    deg_data <- read_enrichment_input_table(fpath)
    if (nrow(deg_data) == 0) {
      message("  ! Skip: file is empty")
      next
    }

    if (isTRUE(run_gsea)) {
      run_universal_gsea(
        df = deg_data,
        tag = comparison_tag,
        out_dir = file.path(out_dir, "GSEA"),
        gene_col = gene_col,
        fc_col = fc_col,
        gene_type = gene_type
      )
    }

    if (isTRUE(run_gprofiler)) {
      run_universal_gprofiler(
        df = deg_data,
        tag = comparison_tag,
        out_dir = file.path(out_dir, "gProfiler"),
        gene_col = gene_col,
        fc_col = fc_col,
        padj_col = padj_col,
        lfc_cut = lfc_cut,
        padj_cut = padj_cut,
        gene_type = gene_type
      )
    }
  }

  message("\nFinished enrichment. Results written to: ", out_dir)
  invisible(TRUE)
}
# ------------------------------------------------------------
# Parse transcriptomics DEG filenames
# ------------------------------------------------------------
# Convert DEG table filenames into lightweight metadata used by MitoCarta
# plotting. The parser supports both current iPSC and AC16 transcriptomics
# naming conventions.
parse_transcriptomics_deg_filename <- function(path, dataset = NA_character_) {
  fn <- basename(path)
  if (!str_detect(fn, "\\.tsv$")) return(NULL)

  core <- fn %>% str_remove("\\.tsv$")

  # Skip summary-like tables if they are present as TSV files.
  if (str_detect(core, regex("DEGs_resume|summary", ignore_case = TRUE))) return(NULL)

  # Within-line stimulation effect: <LINE>_S_vs_NS.
  if (str_detect(core, "^[0-9A-Za-z_]+_S_vs_NS$")) {
    lhs <- str_remove(core, "_S_vs_NS$")
    return(tibble(
      file = path,
      cellline = dataset,
      contrast = core,
      context = "within",
      lhs = lhs,
      rhs = NA_character_,
      group_key = lhs,
      label_pretty = paste0(lhs, " (S vs NS)")
    ))
  }

  # Between-line iPSC effect: <LINE_A><NS/S>_vs_<LINE_B><NS/S>.
  if (str_detect(core, "^[0-9A-Za-z-]+(NS|S)_vs_[0-9A-Za-z-]+(NS|S)$")) {
    parts <- str_match(core, "^([0-9A-Za-z-]+?)(NS|S)_vs_([0-9A-Za-z-]+?)(NS|S)$")
    lhs <- parts[, 2]
    ctx <- parts[, 3]
    rhs <- parts[, 4]

    return(tibble(
      file = path,
      cellline = dataset,
      contrast = core,
      context = ctx,
      lhs = lhs,
      rhs = rhs,
      group_key = lhs,
      label_pretty = paste0(ctx, ": ", lhs, " vs ", rhs)
    ))
  }

  # Between-line AC16 effect: AC16_<LINE_A><NS/S>_vs_<LINE_B><NS/S>.
  if (str_detect(core, "^AC16_[0-9.]+(NS|S)_vs_[0-9.]+(NS|S)$")) {
    parts <- str_match(core, "^AC16_([0-9.]+)(NS|S)_vs_([0-9.]+)(NS|S)$")
    lhs <- paste0("AC16_", parts[, 2])
    ctx <- parts[, 3]
    rhs <- paste0("AC16_", parts[, 4])

    return(tibble(
      file = path,
      cellline = dataset,
      contrast = core,
      context = ctx,
      lhs = lhs,
      rhs = rhs,
      group_key = lhs,
      label_pretty = paste0(ctx, ": ", lhs, " vs ", rhs)
    ))
  }

  # Between-line proteomics-style effect: (NS|S)_<A>_vs_<B>.
  if (str_detect(core, "^(NS|S)_[0-9A-Za-z_]+_vs_[0-9A-Za-z_]+$")) {
    ctx <- str_match(core, "^(NS|S)_")[,2]
    pair <- str_remove(core, "^(NS|S)_")
    parts <- str_split(pair, "_vs_", simplify = TRUE)
    lhs <- parts[,1]
    rhs <- parts[,2]

    return(tibble(
      file = path,
      cellline = dataset,
      contrast = core,
      context = ctx,
      lhs = lhs,
      rhs = rhs,
      group_key = lhs,
      label_pretty = paste0(ctx, ": ", lhs, " vs ", rhs)
    ))
  }

  # Fallback for uncommon filenames.
  tibble(
    file = path,
    cellline = dataset,
    contrast = core,
    context = "unknown",
    lhs = NA_character_,
    rhs = NA_character_,
    group_key = "mixed",
    label_pretty = core
  )
}

# Backward-compatible alias for older callers.
parse_ips_rna_deg_filename <- function(path) {
  parse_transcriptomics_deg_filename(path, dataset = "IPS")
}


expand_mitocarta_hierarchy <- function(mitocarta_df, levels_keep = 2) {
  stopifnot(levels_keep %in% c(1, 2, 3))

  df <- mitocarta_df %>%
    dplyr::transmute(
      Symbol = Symbol,
      pathways_raw = `MitoCarta3.0_MitoPathways`
    ) %>%
    dplyr::filter(!is.na(Symbol), !is.na(pathways_raw), pathways_raw != "") %>%
    tidyr::separate_rows(pathways_raw, sep = "\\s*\\|\\s*") %>%
    dplyr::mutate(pathways_raw = stringr::str_squish(pathways_raw))

  # Split into L1-L4, but make sure columns always exist
  df <- df %>%
    tidyr::separate(
      pathways_raw,
      into = c("L1", "L2", "L3", "L4"),
      sep = "\\s*>\\s*",
      fill = "right",
      extra = "merge",
      remove = TRUE
    ) %>%
    # Guarantee columns exist even if tidyr behaves weirdly on some inputs
    dplyr::mutate(
      L1 = dplyr::coalesce(L1, NA_character_),
      L2 = dplyr::coalesce(L2, NA_character_),
      L3 = dplyr::coalesce(L3, NA_character_)
    ) %>%
    # Clean values + drop junk ("0", "NA", "", "N/A")
    dplyr::mutate(
      dplyr::across(
        c(L1, L2, L3),
        ~ {
          x <- stringr::str_squish(.x)
          x <- dplyr::na_if(x, "")
          x <- dplyr::na_if(x, "0")
          x <- dplyr::na_if(x, "NA")
          x <- dplyr::na_if(x, "N/A")
          x
        }
      )
    ) %>%
    dplyr::filter(!is.na(L1))

  # Build nodes
  out <- if (levels_keep == 1) {
    df %>%
      dplyr::transmute(Symbol, level = 1L, node = L1)

  } else if (levels_keep == 2) {
    dplyr::bind_rows(
      df %>% dplyr::transmute(Symbol, level = 1L, node = L1),
      df %>%
        dplyr::filter(!is.na(L2)) %>%
        dplyr::transmute(Symbol, level = 2L, node = paste(L1, L2, sep = " > "))
    )

  } else {
    dplyr::bind_rows(
      df %>% dplyr::transmute(Symbol, level = 1L, node = L1),
      df %>%
        dplyr::filter(!is.na(L2)) %>%
        dplyr::transmute(Symbol, level = 2L, node = paste(L1, L2, sep = " > ")),
      df %>%
        dplyr::filter(!is.na(L2), !is.na(L3)) %>%
        dplyr::transmute(Symbol, level = 3L, node = paste(L1, L2, L3, sep = " > "))
    )
  }

  out %>%
    dplyr::filter(!is.na(node), node != "") %>%
    dplyr::distinct(Symbol, level, node)
}
# ------------------------------------------------------------
# IPS: parse DEP output filenames (TSV tables)
# Handles:
#   6132_S_vs_NS
#   NS_6135CCC_vs_6135H1
#   S_6135CCC_vs_6132
# ------------------------------------------------------------
parse_ips_prot_limma_filename <- function(path) {
  fn <- basename(path)
  if (!str_detect(fn, "\\.tsv$")) return(NULL)

  core <- fn %>% str_remove("\\.tsv$")

  # 1) within clone: <CLONE>_S_vs_NS
  if (str_detect(core, "^[0-9A-Za-z_]+_S_vs_NS$")) {
    lhs <- str_remove(core, "_S_vs_NS$")
    return(tibble(
      file = path,
      cellline = "IPS",
      contrast = core,
      context = "within",
      lhs = lhs,
      rhs = NA_character_,
      group_key = lhs,                                  # for coloring
      label_pretty = paste0(lhs, " (S vs NS)")
    ))
  }

  # 2) between clones in NS or S context: (NS|S)_<A>_vs_<B>
  if (str_detect(core, "^(NS|S)_[0-9A-Za-z_]+_vs_[0-9A-Za-z_]+$")) {
    ctx <- str_match(core, "^(NS|S)_")[,2]
    pair <- str_remove(core, "^(NS|S)_")
    parts <- str_split(pair, "_vs_", simplify = TRUE)
    lhs <- parts[,1]
    rhs <- parts[,2]

    return(tibble(
      file = path,
      cellline = "IPS",
      contrast = core,
      context = ctx,                                    # "NS" or "S"
      lhs = lhs,
      rhs = rhs,
      group_key = lhs,                                  # color by LHS clone
      label_pretty = paste0(ctx, ": ", lhs, " vs ", rhs)
    ))
  }

  # fallback
  tibble(
    file = path,
    cellline = "IPS",
    contrast = core,
    context = "unknown",
    lhs = NA_character_,
    rhs = NA_character_,
    group_key = "mixed",
    label_pretty = core
  )
}

# ------------------------------------------------------------
# Read limma-like DEP TSV robustly
# ------------------------------------------------------------
read_limma_csv_ips <- function(path) {
  df <- readr::read_tsv(path, show_col_types = FALSE)
  cn <- names(df)

  gene_col <- dplyr::case_when(
    "Genes"     %in% cn ~ "Genes",
    "Gene_name" %in% cn ~ "Gene_name",
    "SYMBOL"    %in% cn ~ "SYMBOL",
    "Symbol"    %in% cn ~ "Symbol",
    TRUE ~ NA_character_
  )
  if (is.na(gene_col)) stop("No gene column found in: ", path)

  lfc_col <- dplyr::case_when(
    "logFC"          %in% cn ~ "logFC",
    "log2FC"         %in% cn ~ "log2FC",
    "log2FoldChange" %in% cn ~ "log2FoldChange",
    TRUE ~ NA_character_
  )
  if (is.na(lfc_col)) stop("No logFC/log2FC column found in: ", path)

  padj_col <- dplyr::case_when(
    "adj.P.Val" %in% cn ~ "adj.P.Val",
    "padj"      %in% cn ~ "padj",
    "FDR"       %in% cn ~ "FDR",
    TRUE ~ NA_character_
  )
  if (is.na(padj_col)) stop("No adj.P.Val/padj/FDR column found in: ", path)

  df %>%
    transmute(
      Gene_name = .data[[gene_col]],
      log2FC    = as.numeric(.data[[lfc_col]]),
      padj      = as.numeric(.data[[padj_col]])
    ) %>%
    filter(!is.na(Gene_name), !is.na(log2FC), !is.na(padj))
}

# Expand MitoCarta pathways into hierarchical nodes (level 1 + 2 + optional 3)
# MitoCarta3.0_MitoPathways format: "A > B > C | A > D | ..."
expand_mitocarta_hierarchy <- function(mitocarta_df, levels_keep = 2) {
  stopifnot(levels_keep %in% c(1, 2, 3))

  df <- mitocarta_df %>%
    dplyr::transmute(
      Symbol = Symbol,
      pathways_raw = `MitoCarta3.0_MitoPathways`
    ) %>%
    dplyr::filter(!is.na(Symbol), !is.na(pathways_raw), pathways_raw != "") %>%
    tidyr::separate_rows(pathways_raw, sep = "\\s*\\|\\s*") %>%
    dplyr::mutate(pathways_raw = stringr::str_squish(pathways_raw))

  # Split into L1-L4, but make sure columns always exist
  df <- df %>%
    tidyr::separate(
      pathways_raw,
      into = c("L1", "L2", "L3", "L4"),
      sep = "\\s*>\\s*",
      fill = "right",
      extra = "merge",
      remove = TRUE
    ) %>%
    # Guarantee columns exist even if tidyr behaves weirdly on some inputs
    dplyr::mutate(
      L1 = dplyr::coalesce(L1, NA_character_),
      L2 = dplyr::coalesce(L2, NA_character_),
      L3 = dplyr::coalesce(L3, NA_character_)
    ) %>%
    # Clean values + drop junk ("0", "NA", "", "N/A")
    dplyr::mutate(
      dplyr::across(
        c(L1, L2, L3),
        ~ {
          x <- stringr::str_squish(.x)
          x <- dplyr::na_if(x, "")
          x <- dplyr::na_if(x, "0")
          x <- dplyr::na_if(x, "NA")
          x <- dplyr::na_if(x, "N/A")
          x
        }
      )
    ) %>%
    dplyr::filter(!is.na(L1))

  # Build nodes
  out <- if (levels_keep == 1) {
    df %>%
      dplyr::transmute(Symbol, level = 1L, node = L1)

  } else if (levels_keep == 2) {
    dplyr::bind_rows(
      df %>% dplyr::transmute(Symbol, level = 1L, node = L1),
      df %>%
        dplyr::filter(!is.na(L2)) %>%
        dplyr::transmute(Symbol, level = 2L, node = paste(L1, L2, sep = " > "))
    )

  } else {
    dplyr::bind_rows(
      df %>% dplyr::transmute(Symbol, level = 1L, node = L1),
      df %>%
        dplyr::filter(!is.na(L2)) %>%
        dplyr::transmute(Symbol, level = 2L, node = paste(L1, L2, sep = " > ")),
      df %>%
        dplyr::filter(!is.na(L2), !is.na(L3)) %>%
        dplyr::transmute(Symbol, level = 3L, node = paste(L1, L2, L3, sep = " > "))
    )
  }

  out %>%
    dplyr::filter(!is.na(node), node != "") %>%
    dplyr::distinct(Symbol, level, node)
}

# Build numbered labels like:
# "1 Metabolism"
# "1.1 Amino acid metabolism"
# We do numbering based on alphabetical order of categories within each parent (stable).
make_numbered_labels <- function(nodes_df) {

  df <- nodes_df %>%
    mutate(
      parent = if_else(level == 1, NA_character_,
                       str_remove(node, "\\s*>\\s*[^>]+$")),
      leaf = if_else(level == 1, node,
                     str_remove(node, "^.*\\s*>\\s*"))
    ) %>%
    mutate(
      parent = if_else(is.na(parent), NA_character_, str_squish(parent)),
      leaf   = str_squish(leaf)
    )

  # Level 1 numbering (alphabetical)
  lvl1 <- df %>%
    filter(level == 1) %>%
    distinct(leaf) %>%
    arrange(leaf) %>%
    mutate(idx1 = row_number())

  # Attach idx1 to everyone based on their TOP parent (for level2 that's 'parent')
  df2 <- df %>%
    mutate(parent_l1 = if_else(level == 1, leaf, parent)) %>%
    left_join(lvl1, by = c("parent_l1" = "leaf"))

  # Level 2 numbering within each parent_l1 (alphabetical)
  df2 <- df2 %>%
    group_by(parent_l1) %>%
    mutate(idx2 = if_else(level == 2, dense_rank(leaf), NA_integer_)) %>%
    ungroup()

  df2 %>%
    mutate(
      label_num = case_when(
        level == 1 ~ sprintf("%d", idx1),
        level == 2 ~ paste0(idx1, ".", idx2),
        TRUE ~ NA_character_
      ),
      # sort key to keep hierarchical ordering: 01, 01.01, 01.02, 02, 02.01 ...
      sort_key = case_when(
        level == 1 ~ sprintf("%02d", idx1),
        level == 2 ~ paste0(sprintf("%02d", idx1), ".", sprintf("%02d", idx2)),
        TRUE ~ NA_character_
      ),
      is_main = (level == 1)
    ) %>%
    select(level, node, parent, leaf, idx1, idx2, label_num, sort_key, is_main)
}

# Build counts up/down per pathway node for one DEG table
build_pathway_counts <- function(deg_df, mit_nodes, lfc_cut = 0.5, padj_cut = 0.05) {
  deg_sig <- deg_df %>%
    filter(padj <= padj_cut, abs(log2FC) >= lfc_cut) %>%
    mutate(Direction = if_else(log2FC > 0, "Upregulated", "Downregulated")) %>%
    select(Gene_name, Direction)

  j <- deg_sig %>%
    # inner_join(mit_nodes, by = c("Gene_name" = "Symbol"))
      inner_join(mit_nodes, by = c("Gene_name" = "Genes"))


  if (nrow(j) == 0) return(tibble())

  wide <- j %>%
    distinct(Gene_name, Direction, level, node) %>%
    dplyr::count(level, node, Direction, name = "n") %>%
    tidyr::pivot_wider(names_from = Direction, values_from = n, values_fill = 0)

  # Ensure both columns exist even if one direction is absent
  if (!"Upregulated" %in% names(wide)) wide$Upregulated <- 0L
  if (!"Downregulated" %in% names(wide)) wide$Downregulated <- 0L

  wide %>%
    rename(up = Upregulated, down = Downregulated) %>%
    mutate(
      up = as.integer(up),
      down = as.integer(down)
    )
}

# Plot signed up/down counts for one comparison.
plot_mitocarta_bar <- function(pathway_counts_tbl, clone_label, clone_color) {
  # pathway_counts_tbl must contain Pathway_label, up, and down columns.

  data_long <- pathway_counts_tbl %>%
    mutate(
      up = as.integer(up),
      down = as.integer(down)
    ) %>%
    pivot_longer(cols = c(up, down), names_to = "dir", values_to = "count") %>%
    mutate(
      Direction = if_else(dir == "up", "Upregulated", "Downregulated"),
      signed = if_else(Direction == "Downregulated", -count, count)
    )

  max_abs <- max(abs(data_long$signed), na.rm = TRUE)
  max_abs <- max_abs + 1

  ggplot(data_long, aes(x = signed, y = Pathway_label, fill = Direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, colour = "grey40", linetype = "dashed") +
    scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
    scale_x_continuous(
      limits = c(-max_abs, max_abs),
      breaks = pretty(c(-max_abs, max_abs)),
      labels = function(x) abs(x)
    ) +
    labs(
      x = "Number of genes (down / up)",
      y = NULL,
      fill = NULL
      #title = paste0(clone_label, " NS vs S")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", color = clone_color),
      axis.text.y = ggtext::element_markdown(size = 9),
      legend.position = "top",
      panel.grid.major.y = element_blank()
    )
}


# ------------------------------------------------------------
# Main runner (IPS)
# ------------------------------------------------------------
run_mitocarta_plots_proteomics_IPS <- function(
    deg_dir,
    mitocarta_file,
    out_dir = deg_dir,
    contrast_keep = NULL,
    lfc_cut = 0.58,
    padj_cut = 0.05,
    levels_keep = 2,
    group_colors = c(
      "6132"    = "#1b9e77",
      "6135CCC" = "#d95f02",
      "6135_CCC"= "#d95f02",
      "6135H1"  = "#7570b3",
      "6135_H1" = "#7570b3",
      "6137"    = "#e7298a",
      "mixed"   = "grey30"
    )
) {
  plot_dir <- file.path(out_dir, "mitocarta_barplots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  # Load mapping: Ensure the mapping table uses "Genes" as the key
  mit_nodes <- expand_mitocarta_hierarchy(mitocarta_file, levels_keep = levels_keep)
  # If expand_mitocarta_hierarchy uses 'Symbol', rename it here to match
  if("Symbol" %in% colnames(mit_nodes)) mit_nodes <- rename(mit_nodes, Genes = Symbol)

  dep_files <- list.files(deg_dir, pattern = "\\.tsv$", full.names = TRUE)
  dep_files <- dep_files[!grepl("Summary_DEP", basename(dep_files), ignore.case = TRUE)]

  meta <- bind_rows(lapply(dep_files, parse_ips_prot_limma_filename))

  if (!is.null(contrast_keep)) {
    meta <- meta %>% filter(str_detect(contrast, fixed(contrast_keep)))
  }
  
  if (nrow(meta) == 0) stop("No .tsv found matching contrast_keep in: ", deg_dir)

  for (i in seq_len(nrow(meta))) {
    info <- meta[i, ]
    
    # 1. Load data and force "Genes" to be a column (not rownames)
    deg_df <- read_limma_csv_ips(info$file)
    colnames(deg_df)[1] <- 'Genes'
    # if (!"Genes" %in% colnames(deg_df)) {
    #   deg_df <- deg_df %>% tibble::rownames_to_column("Genes")
    # }

    # 2. Create the Detailed Gene Table
    # One row per Gene per Category
    gene_detail <- deg_df %>%
      inner_join(mit_nodes, by = "Genes") %>%
      mutate(
        # Symmetric Fold Change (Standard Bioinfo)
        FC = if_else(log2FC >= 0, 2^log2FC, -1/(2^log2FC)),
        # Fixed Significance logic
        Significant = if_else(abs(log2FC) >= lfc_cut & padj <= padj_cut, "Yes", "No")
      ) %>%
      select(
        Gene = Genes, 
        log2FC = log2FC, 
        FC, 
        padj = padj, 
        Mitocarta_category = node, 
        Significant
      ) %>%
      arrange(Mitocarta_category, padj)

    # Save Detailed Report
    gene_detail <- gene_detail[gene_detail$Significant == "Yes",]
    detail_path <- file.path(plot_dir, paste0("gene_detail_mitocarta_IPS_PROT_", info$contrast))
    readr::write_tsv(gene_detail, paste0(detail_path, ".tsv"))
    writexl::write_xlsx(gene_detail, paste0(detail_path, ".xlsx"))

    # 3. Barplot Logic
    colnames(deg_df)[1] <- 'Gene_name'

    counts <- build_pathway_counts(deg_df, mit_nodes, lfc_cut = lfc_cut, padj_cut = padj_cut)
    
    if (nrow(counts) == 0) {
      message("No MitoCarta genes passing thresholds for: ", info$contrast)
      next
    }

    labeled <- make_numbered_labels(counts %>% distinct(level, node))

    counts2 <- counts %>%
      left_join(labeled %>% select(level, node, leaf, label_num, sort_key, is_main),
                by = c("level", "node")) %>%
      filter(!is.na(label_num)) %>%
      mutate(
        label = paste(label_num, leaf),
        label_display = if_else(level == 1, label, paste0("&nbsp;&nbsp;&nbsp;&nbsp;", label)),
        label_display = if_else(is_main, paste0("<b>", label_display, "</b>"), label_display),
        Pathway_label = forcats::fct_rev(forcats::fct_inorder(label_display))
      )

    # Plotting
    key <- as.character(info$group_key)
    group_color <- group_colors[key]
    if (length(group_color) == 0 || is.na(group_color)) group_color <- group_colors[["mixed"]]

    p <- plot_mitocarta_bar(
      pathway_counts_tbl = counts2 %>% select(Pathway_label, up, down),
      clone_label = as.character(info$label_pretty),
      clone_color = unname(group_color)
    )

    ggsave(file.path(plot_dir, paste0("barplot_mitocarta_IPS_PROT_", info$contrast, ".png")),
           p, width = 7, height = 5, dpi = 300)

    message("Finished processing: ", info$contrast)
  }

  invisible(TRUE)
}


# ------------------------------------------------------------
# Read DESeq2-like transcriptomics tables
# ------------------------------------------------------------
# Standardize a DEG table to the columns needed by MitoCarta plotting:
# Gene_name, log2FC, and padj.
read_deseq2_tsv <- function(path) {
  df <- readr::read_tsv(path, show_col_types = FALSE)
  cn <- names(df)

  gene_col <- dplyr::case_when(
    "Gene_name" %in% cn ~ "Gene_name",
    "Genes"     %in% cn ~ "Genes",
    "gene"      %in% cn ~ "gene",
    "Gene"      %in% cn ~ "Gene",
    "symbol"    %in% cn ~ "symbol",
    "SYMBOL"    %in% cn ~ "SYMBOL",
    "Symbol"    %in% cn ~ "Symbol",
    TRUE ~ NA_character_
  )
  if (is.na(gene_col)) stop("No gene column found in: ", path)

  lfc_col <- dplyr::case_when(
    "log2FoldChange" %in% cn ~ "log2FoldChange",
    "logFC"          %in% cn ~ "logFC",
    "log2FC"         %in% cn ~ "log2FC",
    TRUE ~ NA_character_
  )
  if (is.na(lfc_col)) stop("No log2FoldChange/logFC column found in: ", path)

  padj_col <- dplyr::case_when(
    "padj"      %in% cn ~ "padj",
    "adj.P.Val" %in% cn ~ "adj.P.Val",
    "FDR"       %in% cn ~ "FDR",
    TRUE ~ NA_character_
  )
  if (is.na(padj_col)) stop("No padj/adj.P.Val/FDR column found in: ", path)

  df %>%
    transmute(
      Gene_name = .data[[gene_col]],
      log2FC    = as.numeric(.data[[lfc_col]]),
      padj      = as.numeric(.data[[padj_col]])
    ) %>%
    filter(!is.na(Gene_name), !is.na(log2FC), !is.na(padj))
}

# Backward-compatible alias for older callers.
read_deseq2_tsv_ips <- read_deseq2_tsv

# ------------------------------------------------------------
# Transcriptomics MitoCarta runner
# ------------------------------------------------------------
# For each DEG table in a directory, this function:
# - reads DESeq2-like transcriptomics results,
# - intersects significant genes with MitoCarta pathway annotations,
# - writes a detailed gene/category table,
# - writes pathway up/down count summaries,
# - saves a signed pathway barplot.
#
# mitocarta_file may be either:
# - a data frame already loaded in R, or
# - a path to a CSV file.
run_mitocarta_plots_transcriptomics <- function(
  deg_dir,
  mitocarta_file,
  out_dir = deg_dir,
  dataset = NA_character_,
  output_prefix = "RNA",
  contrast_keep = NULL,
  lfc_cut = 1.5,
  padj_cut = 0.05,
  levels_keep = 2,
  group_colors = c(
    "AC16_1184" = "#1b9e77",
    "AC16_127" = "#d95f02",
    "6132"    = "#1b9e77",
    "6135"    = "#d95f02",
    "6135CCC" = "#d95f02",
    "6135_CCC"= "#d95f02",
    "6135H1"  = "#7570b3",
    "6135_H1" = "#7570b3",
    "6137"    = "#e7298a",
    "mixed"   = "grey30"
  )
) {
  plot_dir <- out_dir
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  # Load MitoCarta annotations and standardize the join key.
  if (is.character(mitocarta_file) && length(mitocarta_file) == 1) {
    mitocarta_file <- openxlsx::read.xlsx(mitocarta_file)

    bad <- mitocarta_file$MitoCarta3.0_MitoPathways == "MIM" &
      !is.na(mitocarta_file$X7)

    mitocarta_file$Description[bad] <- paste(
      mitocarta_file$Description[bad],
      mitocarta_file$MitoCarta3.0_SubMitoLocalization[bad],
      sep = ", "
    )

    mitocarta_file$MitoCarta3.0_SubMitoLocalization[bad] <- mitocarta_file$MitoCarta3.0_MitoPathways[bad]
    mitocarta_file$MitoCarta3.0_MitoPathways[bad] <- mitocarta_file$X7[bad]
  }
  mit_nodes <- expand_mitocarta_hierarchy(mitocarta_file, levels_keep = levels_keep)
  if ("Symbol" %in% colnames(mit_nodes)) mit_nodes <- rename(mit_nodes, Genes = Symbol)

  deg_files <- list.files(deg_dir, pattern = "\\.tsv$", full.names = TRUE)
  deg_files <- deg_files[!grepl("DEGs_resume|Summary", basename(deg_files), ignore.case = TRUE)]

  meta <- bind_rows(lapply(deg_files, parse_transcriptomics_deg_filename, dataset = dataset)) %>%
    filter(!is.na(file))

  if (!is.null(contrast_keep)) {
    meta <- meta %>% filter(str_detect(contrast, fixed(contrast_keep)))
  }
  if (nrow(meta) == 0) stop("No DEG .tsv found matching contrast_keep in: ", deg_dir)

  for (i in seq_len(nrow(meta))) {
    info <- meta[i, ]

    # Load DEG data as Gene_name/log2FC/padj.
    deg_df <- read_deseq2_tsv(info$file)

    # Create one row per significant gene and MitoCarta category.
    gene_detail <- deg_df %>%
      inner_join(mit_nodes, by = c("Gene_name" = "Genes")) %>%
      mutate(
        FC = if_else(log2FC >= 0, 2^log2FC, -1/(2^log2FC)),
        Significant = if_else(abs(log2FC) >= lfc_cut & padj <= padj_cut, "Yes", "No")
      ) %>%
      select(
        Gene = Gene_name, 
        log2FC, 
        FC, 
        padj, 
        Mitocarta_category = node, 
        Significant
      ) %>%
      filter(Significant == "Yes") %>%
      arrange(Mitocarta_category, padj)

    # Save detailed gene/category tables.
    detail_path <- file.path(plot_dir, paste0("gene_detail_mitocarta_", output_prefix, "_", info$contrast))
    readr::write_tsv(gene_detail, paste0(detail_path, ".tsv"))
    writexl::write_xlsx(gene_detail, paste0(detail_path, ".xlsx"))

    # Count significant up/down genes per MitoCarta pathway node.
    counts <- build_pathway_counts(deg_df, mit_nodes, lfc_cut = lfc_cut, padj_cut = padj_cut)
    
    if (nrow(counts) == 0) {
      message("No MitoCarta genes passing thresholds for: ", info$contrast)
      next
    }

    labeled <- make_numbered_labels(counts %>% distinct(level, node))

    counts2 <- counts %>%
      left_join(labeled %>% select(level, node, leaf, label_num, sort_key, is_main),
                by = c("level", "node")) %>%
      filter(!is.na(label_num), !is.na(leaf)) %>%
      mutate(
        label = paste(label_num, leaf),
        label_display = if_else(level == 1, label, paste0("&nbsp;&nbsp;&nbsp;&nbsp;", label)),
        label_display = if_else(is_main, paste0("<b>", label_display, "</b>"), label_display)
      ) %>%
      arrange(sort_key) %>%
      mutate(Pathway_label = forcats::fct_rev(forcats::fct_inorder(label_display)))

    # Plot signed up/down pathway counts.
    key <- as.character(info$group_key)
    group_color <- group_colors[key]
    if (length(group_color) == 0 || is.na(group_color)) group_color <- group_colors[["mixed"]]
    group_color <- unname(group_color)

    p <- plot_mitocarta_bar(
      pathway_counts_tbl = counts2 %>% select(Pathway_label, up, down),
      clone_label = as.character(info$label_pretty),
      clone_color = group_color
    )

    out_png <- file.path(plot_dir, paste0("barplot_mitocarta_", output_prefix, "_", info$contrast, ".png"))
    ggsave(out_png, p, width = 7, height = 5, dpi = 300)

    # Save the plotted pathway-count summary.
    out_tsv_counts <- file.path(plot_dir, paste0("counts_mitocarta_", output_prefix, "_", info$contrast, ".tsv"))
    readr::write_tsv(counts2 %>% dplyr::select(Pathway = label, up, down), out_tsv_counts)

    message("Finished MitoCarta transcriptomics plot: ", info$contrast)
  }

  invisible(TRUE)
}

# Backward-compatible alias for older iPSC transcriptomics callers.
run_mitocarta_plots_transcriptomics_IPS <- function(...) {
  run_mitocarta_plots_transcriptomics(..., dataset = "IPS", output_prefix = "IPS_RNA")
}
