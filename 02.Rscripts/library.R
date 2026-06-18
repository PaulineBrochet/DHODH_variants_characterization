# ============================================================
# Project package loader
# ============================================================
#
# This is the only script in the cleaned workflow that should call library().
# Other scripts define analysis functions and assume this loader has already
# been sourced. Keeping package loading centralized makes the notebook and
# GitHub workflow easier to audit and reproduce.

PROJECT_PACKAGE_GROUPS <- list(
  core_io = c(
    "tidyverse", "data.table", "readr", "tibble", "dplyr", "tidyr",
    "stringr", "forcats", "purrr"
  ),
  plotting = c(
    "ggplot2", "ggrepel", "pheatmap", "ComplexHeatmap", "circlize",
    "EnhancedVolcano", "scatterpie", "grid"
  ),
  statistics = c(
    "DESeq2", "limma", "vsn", "matrixStats", "pcaMethods"
  ),
  annotation_enrichment = c(
    "AnnotationDbi", "org.Hs.eg.db", "org.Rn.eg.db", "org.Mm.eg.db",
    "msigdbr", "clusterProfiler", "gprofiler2"
  ),
  export = c(
    "openxlsx", "writexl"
  )
)

required_project_packages <- function(groups = names(PROJECT_PACKAGE_GROUPS)) {
  unique(unlist(PROJECT_PACKAGE_GROUPS[groups], use.names = FALSE))
}

check_project_packages <- function(groups = names(PROJECT_PACKAGE_GROUPS)) {
  packages <- required_project_packages(groups)
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0) {
    stop(
      "Missing required R package(s): ",
      paste(missing, collapse = ", "),
      "\nInstall the missing CRAN/Bioconductor packages before running the pipeline.",
      call. = FALSE
    )
  }

  invisible(packages)
}

load_project_libraries <- function(groups = names(PROJECT_PACKAGE_GROUPS),
                                   quiet = TRUE) {
  packages <- check_project_packages(groups)
  cache_key <- paste(sort(packages), collapse = "|")

  if (identical(getOption("dhodh.loaded_packages"), cache_key)) {
    return(invisible(packages))
  }

  loader <- function(pkg) {
    library(pkg, character.only = TRUE)
  }

  if (isTRUE(quiet)) {
    suppressPackageStartupMessages(invisible(lapply(packages, loader)))
  } else {
    invisible(lapply(packages, loader))
  }

  options(dhodh.loaded_packages = cache_key)
  invisible(packages)
}
