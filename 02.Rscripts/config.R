# ============================================================
# Project configuration
# ============================================================
#
# Centralized paths and analysis settings used by the notebook callers.
# Analysis functions should stay generic; project-specific choices such as
# input files, output folders, excluded lines, and requested comparisons live
# here.

# ------------------------------------------------------------
# Raw data paths
# ------------------------------------------------------------

# Transcriptomic

PATH.TRANSCRIPTOMOC_RAW_DATA <- '01.Data/01.Transcriptomic/Raw_counts.tab'
PATH.TRANSCRIPTOMIC_RAW_DATA <- PATH.TRANSCRIPTOMOC_RAW_DATA
PATH.TRANSCRIPTOMIC_SAMPLE_METADATA <- '01.Data/01.Transcriptomic/samples.tsv'
PATH.PUBLIC_RNASEQ_COUNTS <- "01.Data/04.Public_data/01.RNAseq/GSE191081_Count_table.txt"
PATH.PUBLIC_RNASEQ_NORMALIZED <- "01.Data/04.Public_data/01.RNAseq/RNA_seq_norm_data.txt"
PATH.PUBLIC_MICROARRAY_DIR <- "01.Data/04.Public_data/02.Microarray"
PATH.MITOCARTA_DB <- "01.Data/mitocarta.xlsx"

# Metabolomics raw table
PATH.METABOLOMIC_RAW_DATA <- "01.Data/03.Metabolomic/Raw_norm_peak_area.csv"

# Proteomics raw reports
PATH.PROTEOMIC_AC16_REPORT <- "01.Data/02.Proteomic/report.pg_matrix.tsv"
PATH.PROTEOMIC_IPS_REPORT_240 <- "01.Data/02.Proteomic/report.pg_matrix_240.tsv"
PATH.PROTEOMIC_IPS_REPORT_480 <- "01.Data/02.Proteomic/report.pg_matrix_480.tsv"


# ------------------------------------------------------------
# Results paths
# ------------------------------------------------------------

# AC16
PATH.AC16_results <- '03.Results/01.AC16/'

# IPS
PATH.IPS_results <- '03.Results/02.IPS/'

# Public data
PATH.PUBLIC_results <- "03.Results/03.Public_data/"


# ------------------------------------------------------------
# Proteomics: shared thresholds
# ------------------------------------------------------------

PROTEOMIC_LFC_CUT <- 0.58
PROTEOMIC_PADJ_CUT <- 0.05


# ------------------------------------------------------------
# Proteomics: AC16 DEP settings
# ------------------------------------------------------------

# AC16 is a single-batch proteomics dataset. The 50.50 condition is excluded
# before filtering and normalization because it is not part of the retained
# manuscript DEP analysis.
PROTEOMIC_AC16_SAMPLE_GROUPS <- list(
  `1184_NS` = c("01B", "10", "16"),
  `1184_10.5` = c("02B", "11", "17"),
  `127_NS` = c("04B", "07B", "13"),
  `127_10.5` = c("05B", "08", "14")
)

PROTEOMIC_AC16_CONTRASTS <- c(
  # Within-clone stimulation effects
  `1184_10.5_vs_NS` = "X1184_10.5 - X1184_NS",
  `127_10.5_vs_NS` = "X127_10.5 - X127_NS",

  # Between-clone effects, expressed as 127 minus 1184
  `127_NS_vs_1184_NS` = "X127_NS - X1184_NS",
  `127_10.5_vs_1184_10.5` = "X127_10.5 - X1184_10.5"
)

PROTEOMIC_AC16_ANALYSIS <- list(
  report_file = PATH.PROTEOMIC_AC16_REPORT,
  out_dir = file.path(PATH.AC16_results, "02.Proteomic", "DEP"),
  sample_groups = PROTEOMIC_AC16_SAMPLE_GROUPS,
  contrast_defs = PROTEOMIC_AC16_CONTRASTS,
  normalized_matrix_file = "normalized_matrix_log2_no5050.csv",
  min_rep_per_group = 1,
  min_groups_required = length(PROTEOMIC_AC16_SAMPLE_GROUPS),
  nPcs_impute = 5,
  lfc_cut = PROTEOMIC_LFC_CUT,
  padj_cut = PROTEOMIC_PADJ_CUT
)


# ------------------------------------------------------------
# Proteomics: iPSC DEP settings
# ------------------------------------------------------------

PROTEOMIC_IPS_EXCLUDED_LINES <- "6132"

PROTEOMIC_IPS_CONDITION_LEVELS <- c(
  "6135_CCC_NS", "6135_CCC_S",
  "6135_H1_NS", "6135_H1_S",
  "6137_NS", "6137_S"
)

PROTEOMIC_IPS_CONTRASTS <- c(
  # Within-line stimulation effect
  `6135CCC_S_vs_NS` = "condition_factor6135_CCC_S - condition_factor6135_CCC_NS",
  `6135H1_S_vs_NS` = "condition_factor6135_H1_S - condition_factor6135_H1_NS",
  `6137_S_vs_NS` = "condition_factor6137_S - condition_factor6137_NS",

  # Baseline line effects
  `NS_6135CCC_vs_6135H1` = "condition_factor6135_CCC_NS - condition_factor6135_H1_NS",
  `NS_6135H1_vs_6137` = "condition_factor6135_H1_NS - condition_factor6137_NS",
  `NS_6135CCC_vs_6137` = "condition_factor6135_CCC_NS - condition_factor6137_NS",

  # Stimulated line effects
  `S_6135CCC_vs_6135H1` = "condition_factor6135_CCC_S - condition_factor6135_H1_S",
  `S_6135H1_vs_6137` = "condition_factor6135_H1_S - condition_factor6137_S",
  `S_6135CCC_vs_6137` = "condition_factor6135_CCC_S - condition_factor6137_S"
)

PROTEOMIC_IPS_PREPROCESSING <- list(
  report240_file = PATH.PROTEOMIC_IPS_REPORT_240,
  report480_file = PATH.PROTEOMIC_IPS_REPORT_480,
  out_base = file.path(PATH.IPS_results, "02.Proteomic", "QC_batch_effect_correction"),
  trim_tail_240 = 3,
  trim_tail_480 = 7,
  min_n_per_group = 1,
  nPcs_impute = 5,
  exclude_lines = PROTEOMIC_IPS_EXCLUDED_LINES
)

PROTEOMIC_IPS_ANALYSIS <- list(
  matrix_file = file.path(
    PROTEOMIC_IPS_PREPROCESSING$out_base,
    "matrices",
    "matrix_log2_merged_scaleNormalized.csv"
  ),
  metadata_file = file.path(
    PROTEOMIC_IPS_PREPROCESSING$out_base,
    "metadata",
    "metadata_merged_3batches.csv"
  ),
  report_files = c(PATH.PROTEOMIC_IPS_REPORT_240, PATH.PROTEOMIC_IPS_REPORT_480),
  out_dir = file.path(PATH.IPS_results, "02.Proteomic", "DEP"),
  condition_levels = PROTEOMIC_IPS_CONDITION_LEVELS,
  contrast_defs = PROTEOMIC_IPS_CONTRASTS,
  exclude_lines = PROTEOMIC_IPS_EXCLUDED_LINES,
  lfc_cut = PROTEOMIC_LFC_CUT,
  padj_cut = PROTEOMIC_PADJ_CUT
)


# ------------------------------------------------------------
# Proteomics: plot settings
# ------------------------------------------------------------

PROTEOMIC_AC16_PLOT_ANALYSIS <- list(
  matrix_file = file.path(
    PROTEOMIC_AC16_ANALYSIS$out_dir,
    PROTEOMIC_AC16_ANALYSIS$normalized_matrix_file
  ),
  dep_dir = PROTEOMIC_AC16_ANALYSIS$out_dir,
  qc_out_dir = file.path(PATH.AC16_results, "02.Proteomic", "QC"),
  volcano_out_dir = file.path(PATH.AC16_results, "02.Proteomic", "DEP", "volcano"),
  dataset_label = "AC16",
  file_pattern = "^limma_.*\\.csv$",
  lfc_cut = PROTEOMIC_LFC_CUT,
  padj_cut = PROTEOMIC_PADJ_CUT
)

PROTEOMIC_IPS_PLOT_ANALYSIS <- list(
  matrix_file = PROTEOMIC_IPS_ANALYSIS$matrix_file,
  metadata_file = PROTEOMIC_IPS_ANALYSIS$metadata_file,
  dep_dir = file.path(PROTEOMIC_IPS_ANALYSIS$out_dir, "tables"),
  qc_out_dir = file.path(PATH.IPS_results, "02.Proteomic", "QC"),
  volcano_out_dir = file.path(PROTEOMIC_IPS_ANALYSIS$out_dir, "volcano"),
  dataset_label = "IPS",
  file_pattern = "\\.tsv$",
  exclude_lines = PROTEOMIC_IPS_EXCLUDED_LINES,
  lfc_cut = PROTEOMIC_LFC_CUT,
  padj_cut = PROTEOMIC_PADJ_CUT
)


# ------------------------------------------------------------
# Proteomics: MitoCarta/glycolysis volcano settings
# ------------------------------------------------------------

PROTEOMIC_AC16_GLYCO_MITO_VOLCANO <- list(
  comparisons = data.frame(
    file = c(
      file.path(PROTEOMIC_AC16_ANALYSIS$out_dir, "limma_127_10.5_vs_NS.csv"),
      file.path(PROTEOMIC_AC16_ANALYSIS$out_dir, "limma_1184_10.5_vs_NS.csv")
    ),
    label = c("AC16 1.27 10.5 vs NS", "AC16 1.184 10.5 vs NS"),
    shape_group = c("CT", "CC"),
    stringsAsFactors = FALSE
  ),
  out_dir = file.path(PATH.AC16_results, "02.Proteomic", "glycolysis_mitocarta_volcano"),
  output_prefix = "AC16_proteomics_glycolysis_mitocarta_volcano",
  title = "AC16 proteomics",
  feature_col = "Genes",
  lfc_col = "logFC",
  padj_col = "adj.P.Val",
  lfc_cut = PROTEOMIC_LFC_CUT,
  padj_cut = PROTEOMIC_PADJ_CUT
)

PROTEOMIC_IPS_GLYCO_MITO_VOLCANO <- list(
  comparisons = data.frame(
    file = c(
      file.path(PROTEOMIC_IPS_ANALYSIS$out_dir, "tables", "6135CCC_S_vs_NS.tsv"),
      file.path(PROTEOMIC_IPS_ANALYSIS$out_dir, "tables", "6135H1_S_vs_NS.tsv")
    ),
    label = c("iPSC-CM 6135 S vs NS", "iPSC-CM 6135-H1 S vs NS"),
    shape_group = c("CT", "CC"),
    stringsAsFactors = FALSE
  ),
  out_dir = file.path(PATH.IPS_results, "02.Proteomic", "glycolysis_mitocarta_volcano"),
  output_prefix = "IPS_proteomics_glycolysis_mitocarta_volcano",
  title = "iPSC-CM proteomics",
  feature_col = "Genes",
  lfc_col = "logFC",
  padj_col = "adj.P.Val",
  lfc_cut = PROTEOMIC_LFC_CUT,
  padj_cut = PROTEOMIC_PADJ_CUT
)


# ------------------------------------------------------------
# Proteomics: enrichment settings
# ------------------------------------------------------------

PROTEOMIC_ENRICHMENT_ANALYSES <- data.frame(
  dataset = c("AC16", "IPS"),
  dep_dir = c(
    PROTEOMIC_AC16_ANALYSIS$out_dir,
    file.path(PROTEOMIC_IPS_ANALYSIS$out_dir, "tables")
  ),
  out_dir = c(
    file.path(PATH.AC16_results, "02.Proteomic", "Enrichment"),
    file.path(PATH.IPS_results, "02.Proteomic", "Enrichment")
  ),
  file_pattern = c("^limma_.*\\.csv$", "\\.tsv$"),
  gene_col = c("Genes", "Genes"),
  fc_col = c("logFC", "logFC"),
  padj_col = c("adj.P.Val", "adj.P.Val"),
  lfc_cut = c(PROTEOMIC_LFC_CUT, PROTEOMIC_LFC_CUT),
  padj_cut = c(PROTEOMIC_PADJ_CUT, PROTEOMIC_PADJ_CUT),
  gene_type = c("HGNC", "HGNC"),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# Metabolomics: DEM settings
# ------------------------------------------------------------

METABOLOMIC_LFC_CUT <- 0.58
METABOLOMIC_PADJ_CUT <- 0.05
METABOLOMIC_QUAL_CUT <- 4

METABOLOMIC_AC16_ANALYSIS <- list(
  celltype = "AC16",
  raw_data_file = PATH.METABOLOMIC_RAW_DATA,
  out_dir = file.path(PATH.AC16_results, "03.Metabolomic"),
  target_clones = c("1184", "1.27"),
  outliers = character(),
  qual_cut = METABOLOMIC_QUAL_CUT,
  lfc_cut = METABOLOMIC_LFC_CUT,
  padj_cut = METABOLOMIC_PADJ_CUT
)

METABOLOMIC_IPS_ANALYSIS <- list(
  celltype = "iPSCM",
  raw_data_file = PATH.METABOLOMIC_RAW_DATA,
  out_dir = file.path(PATH.IPS_results, "03.Metabolomic"),
  target_clones = c("6132", "6137", "6135", "6135-H1"),
  output_exclude_pattern = "6132",
  outliers = c("Sample_32", "Sample_41"),
  qual_cut = METABOLOMIC_QUAL_CUT,
  lfc_cut = METABOLOMIC_LFC_CUT,
  padj_cut = METABOLOMIC_PADJ_CUT
)


# ------------------------------------------------------------
# Transcriptomics: AC16 DEG comparisons
# ------------------------------------------------------------

TRANSCRIPTOMIC_AC16_DATASET <- "AC16"

# Each row defines one AC16 DESeq2 caller. The columns are passed directly to
# run_deg() in the notebook caller, following the same structure as the iPSC
# comparison table below.
TRANSCRIPTOMIC_AC16_DEG_COMPARISONS <- data.frame(
  section = c(
    rep("Within-line stimulation effect", 2),
    rep("Between-line phenotype effect", 2)
  ),
  tag = c(
    "AC16_1184_S_vs_NS",
    "AC16_127_S_vs_NS",
    "AC16_127NS_vs_1184NS",
    "AC16_127S_vs_1184S"
  ),
  subset_expr = c(
    "line == 'AC16_1.184'",
    "line == 'AC16_1.27'",
    "stim == 'NS' & line %in% c('AC16_1.27', 'AC16_1.184')",
    "stim == 'S' & line %in% c('AC16_1.27', 'AC16_1.184')"
  ),
  design = c(
    rep("~ stim", 2),
    rep("~ line", 2)
  ),
  contrast_factor = c(
    rep("stim", 2),
    rep("line", 2)
  ),
  contrast_case = c(
    "S",
    "S",
    "AC16_1.27",
    "AC16_1.27"
  ),
  contrast_control = c(
    "NS",
    "NS",
    "AC16_1.184",
    "AC16_1.184"
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# Transcriptomics: iPSC DEG comparisons
# ------------------------------------------------------------

TRANSCRIPTOMIC_IPS_DATASET <- "iPSC"
TRANSCRIPTOMIC_IPS_EXCLUDED_LINES <- "6132"

# Each row defines one DESeq2 caller. The columns are passed directly to
# run_deg() in the notebook caller. The contrast columns are kept explicit so
# the final contrast vector remains visible while avoiding repeated calls.
TRANSCRIPTOMIC_IPS_DEG_COMPARISONS <- data.frame(
  section = c(
    rep("Within-line stimulation effect", 3),
    rep("Between-line phenotype effect", 6)
  ),
  tag = c(
    "6135_S_vs_NS",
    "6137_S_vs_NS",
    "6135H1_S_vs_NS",
    "6135NS_vs_6137NS",
    "6135S_vs_6137S",
    "6135NS_vs_6135H1NS",
    "6135S_vs_6135H1S",
    "6137NS_vs_6135H1NS",
    "6137S_vs_6135H1S"
  ),
  subset_expr = c(
    "line == '6135'",
    "line == '6137'",
    "line == '6135H1'",
    "stim == 'NS' & line %in% c('6135', '6137')",
    "stim == 'S' & line %in% c('6135', '6137')",
    "stim == 'NS' & line %in% c('6135', '6135H1')",
    "stim == 'S' & line %in% c('6135', '6135H1')",
    "stim == 'NS' & line %in% c('6137', '6135H1')",
    "stim == 'S' & line %in% c('6137', '6135H1')"
  ),
  design = c(
    rep("~ stim", 3),
    rep("~ line", 6)
  ),
  contrast_factor = c(
    rep("stim", 3),
    rep("line", 6)
  ),
  contrast_case = c(
    "S",
    "S",
    "S",
    "6135",
    "6135",
    "6135",
    "6135",
    "6137",
    "6137"
  ),
  contrast_control = c(
    "NS",
    "NS",
    "NS",
    "6137",
    "6137",
    "6135H1",
    "6135H1",
    "6135H1",
    "6135H1"
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# Transcriptomics: enrichment settings
# ------------------------------------------------------------

# Each row defines one transcriptomics enrichment caller. The DEG directories
# are produced by the AC16/iPSC transcriptomics chunks above.
TRANSCRIPTOMIC_ENRICHMENT_ANALYSES <- data.frame(
  dataset = c(
    TRANSCRIPTOMIC_AC16_DATASET,
    TRANSCRIPTOMIC_IPS_DATASET
  ),
  deg_dir = c(
    file.path(PATH.AC16_results, "01.Transcriptomic", "deg"),
    file.path(PATH.IPS_results, "01.Transcriptomic", "deg")
  ),
  out_dir = c(
    file.path(PATH.AC16_results, "01.Transcriptomic", "02.Enrichment"),
    file.path(PATH.IPS_results, "01.Transcriptomic", "02.Enrichment")
  ),
  gene_col = c("Gene_name", "Gene_name"),
  fc_col = c("log2FoldChange", "log2FoldChange"),
  padj_col = c("padj", "padj"),
  lfc_cut = c(1.5, 1.5),
  padj_cut = c(0.05, 0.05),
  gene_type = c("HGNC", "HGNC"),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# Transcriptomics: MitoCarta settings
# ------------------------------------------------------------

TRANSCRIPTOMIC_MITOCARTA_ANALYSES <- data.frame(
  dataset = c(
    TRANSCRIPTOMIC_AC16_DATASET,
    TRANSCRIPTOMIC_IPS_DATASET
  ),
  deg_dir = c(
    file.path(PATH.AC16_results, "01.Transcriptomic", "deg"),
    file.path(PATH.IPS_results, "01.Transcriptomic", "deg")
  ),
  out_dir = c(
    file.path(PATH.AC16_results, "01.Transcriptomic", "mitocarta_barplots"),
    file.path(PATH.IPS_results, "01.Transcriptomic", "mitocarta_barplots")
  ),
  mitocarta_path = c(PATH.MITOCARTA_DB, PATH.MITOCARTA_DB),
  output_prefix = c("AC16_RNA", "IPS_RNA"),
  lfc_cut = c(1.5, 1.5),
  padj_cut = c(0.05, 0.05),
  levels_keep = c(2, 2),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# Transcriptomics: MitoCarta/glycolysis volcano settings
# ------------------------------------------------------------

TRANSCRIPTOMIC_AC16_GLYCO_MITO_VOLCANO <- list(
  comparisons = data.frame(
    file = c(
      file.path(PATH.AC16_results, "01.Transcriptomic", "deg", "AC16_127_S_vs_NS.tsv"),
      file.path(PATH.AC16_results, "01.Transcriptomic", "deg", "AC16_1184_S_vs_NS.tsv")
    ),
    label = c("AC16 1.27 S vs NS", "AC16 1.184 S vs NS"),
    shape_group = c("CT", "CC"),
    stringsAsFactors = FALSE
  ),
  out_dir = file.path(PATH.AC16_results, "01.Transcriptomic", "glycolysis_mitocarta_volcano"),
  output_prefix = "AC16_RNA_glycolysis_mitocarta_volcano",
  title = "AC16 transcriptomics",
  feature_col = "Gene_name",
  lfc_col = "log2FoldChange",
  padj_col = "padj",
  lfc_cut = 1.5,
  padj_cut = 0.05
)

TRANSCRIPTOMIC_IPS_GLYCO_MITO_VOLCANO <- list(
  comparisons = data.frame(
    file = c(
      file.path(PATH.IPS_results, "01.Transcriptomic", "deg", "6135_S_vs_NS.tsv"),
      file.path(PATH.IPS_results, "01.Transcriptomic", "deg", "6135H1_S_vs_NS.tsv")
    ),
    label = c("iPSC-CM 6135 S vs NS", "iPSC-CM 6135-H1 S vs NS"),
    shape_group = c("CT", "CC"),
    stringsAsFactors = FALSE
  ),
  out_dir = file.path(PATH.IPS_results, "01.Transcriptomic", "glycolysis_mitocarta_volcano"),
  output_prefix = "IPS_RNA_glycolysis_mitocarta_volcano",
  title = "iPSC-CM transcriptomics",
  feature_col = "Gene_name",
  lfc_col = "log2FoldChange",
  padj_col = "padj",
  lfc_cut = 1.5,
  padj_cut = 0.05
)


# ------------------------------------------------------------
# Public data: heart-tissue immune-axis expression settings
# ------------------------------------------------------------

PUBLIC_IMMUNE_AXIS_GENES <- c(
  "IFNG",
  "CXCR3",
  "CXCL9",
  "CXCL10",
  "CXCL11",
  "TBX21",
  "CD3D",
  "CD3E",
  "CD3G"
)

PUBLIC_RNASEQ_CCC_CTRL_ANALYSIS <- list(
  input_type = "rnaseq_normalized",
  input_path = PATH.PUBLIC_RNASEQ_NORMALIZED,
  out_dir = file.path(PATH.PUBLIC_results, "01.RNAseq_CCC_CTRL"),
  output_prefix = "public_rnaseq_CCC_CTRL_immune_axis",
  dataset_label = "Public RNA-seq heart tissue CCC vs Control",
  target_genes = PUBLIC_IMMUNE_AXIS_GENES,
  groups_keep = c("CCC", "Control")
)

PUBLIC_RNASEQ_ALL_PHENOTYPES_ANALYSIS <- list(
  input_type = "rnaseq_counts",
  input_path = PATH.PUBLIC_RNASEQ_COUNTS,
  out_dir = file.path(PATH.PUBLIC_results, "02.RNAseq_all_phenotypes"),
  output_prefix = "public_rnaseq_all_phenotypes_immune_axis",
  dataset_label = "Public RNA-seq heart tissue all phenotypes",
  target_genes = PUBLIC_IMMUNE_AXIS_GENES,
  groups_keep = c("CCC", "Control", "DCM")
)

PUBLIC_MICROARRAY_CCC_CTRL_ANALYSIS <- list(
  input_type = "microarray_agilent",
  input_path = PATH.PUBLIC_MICROARRAY_DIR,
  out_dir = file.path(PATH.PUBLIC_results, "03.Microarray_CCC_CTRL"),
  output_prefix = "public_microarray_CCC_CTRL_immune_axis",
  dataset_label = "Public microarray heart tissue CCC vs Control",
  target_genes = PUBLIC_IMMUNE_AXIS_GENES,
  groups_keep = c("CCC", "Control")
)

PUBLIC_MICROARRAY_ALL_PHENOTYPES_ANALYSIS <- list(
  input_type = "microarray_agilent",
  input_path = PATH.PUBLIC_MICROARRAY_DIR,
  out_dir = file.path(PATH.PUBLIC_results, "04.Microarray_all_phenotypes"),
  output_prefix = "public_microarray_all_phenotypes_immune_axis",
  dataset_label = "Public microarray heart tissue all available phenotypes",
  target_genes = PUBLIC_IMMUNE_AXIS_GENES,
  groups_keep = c("CCC", "Control", "DCM")
)


# ------------------------------------------------------------
# Multi-omics C/T vs C/C clustering settings
# ------------------------------------------------------------

# These colors and scale limits are shared by AC16 and iPSC-CM so the two
# heatmaps can use a single OMIC, cluster, and log2FC scale legend.
CLUSTERING_OMIC_COLORS <- c(
  "Transcriptomic" = "#1b9e77",
  "Proteomic" = "#d95f02",
  "Metabolomic" = "#7570b3"
)

CLUSTERING_CLUSTER_COLORS <- c(
  "1" = "#d73027",
  "2" = "#4575b4",
  "3" = "#fc8d59",
  "4" = "#91bfdb",
  "5" = "#b2182b",
  "6" = "#2166ac"
)

CLUSTERING_CLUSTER_LABELS <- c(
  "1" = "1 Conserved up",
  "2" = "2 Conserved down",
  "3" = "3 C/T-specific up",
  "4" = "4 C/T-specific down",
  "5" = "5 C/C-specific up",
  "6" = "6 C/C-specific down",
  "Mixed" = "Opposite / mixed"
)

CLUSTERING_CLUSTER_LEGEND_LABELS <- c(
  "1" = "1 Features upregulated in both stimulated cells",
  "2" = "2 Features downregulated in both stimulated cells",
  "3" = "3 Features upregulated only in stimulated C/T cells",
  "4" = "4 Features downregulated only in stimulated C/T cells",
  "5" = "5 Features upregulated only in stimulated C/C cells",
  "6" = "6 Features downregulated only in stimulated C/C cells"
)

CLUSTERING_PLOT_SETTINGS <- list(
  omic_colors = CLUSTERING_OMIC_COLORS,
  cluster_colors = CLUSTERING_CLUSTER_COLORS,
  cluster_labels = CLUSTERING_CLUSTER_LABELS,
  cluster_legend_labels = CLUSTERING_CLUSTER_LEGEND_LABELS,
  thresholds = list(
    Transcriptomic = 1.5,
    Proteomic = PROTEOMIC_LFC_CUT,
    Metabolomic = METABOLOMIC_LFC_CUT
  ),
  scale_max = list(
    Transcriptomic = 17,
    Proteomic = 6,
    Metabolomic = 13
  ),
  heatmap_colors = c("#00265E", "#5E7F9D", "white", "#A85B61", "#67001E"),
  include_mixed = FALSE,
  row_height_cm = 0.018,
  max_panel_height_cm = 9
)

MULTIOMICS_CLUSTERING_AC16 <- list(
  dataset = "AC16",
  out_dir = file.path(PATH.AC16_results, "04.Clustering_CT_CC"),
  output_prefix = "AC16_CT_CC",
  column_labels = c("AC16\n1.27 C/T\nS/NS", "AC16\n1.184 C/C\nS/NS"),
  heatmap_width = 7,
  heatmap_height = 11,
  legend_width = 9,
  legend_height = 4,
  dpi = 600,
  comparisons = data.frame(
    omic = c(
      "Transcriptomic", "Transcriptomic",
      "Proteomic", "Proteomic",
      "Metabolomic", "Metabolomic"
    ),
    file = c(
      file.path(PATH.AC16_results, "01.Transcriptomic", "deg", "AC16_127_S_vs_NS.tsv"),
      file.path(PATH.AC16_results, "01.Transcriptomic", "deg", "AC16_1184_S_vs_NS.tsv"),
      file.path(PATH.AC16_results, "02.Proteomic", "DEP", "limma_127_10.5_vs_NS.csv"),
      file.path(PATH.AC16_results, "02.Proteomic", "DEP", "limma_1184_10.5_vs_NS.csv"),
      file.path(PATH.AC16_results, "03.Metabolomic", "DE_tables", "Stim_vs_NS_1_27.tsv"),
      file.path(PATH.AC16_results, "03.Metabolomic", "DE_tables", "Stim_vs_NS_1184.tsv")
    ),
    contrast = c("CT", "CC", "CT", "CC", "CT", "CC"),
    lfc_cut = c(1.5, 1.5, PROTEOMIC_LFC_CUT, PROTEOMIC_LFC_CUT, METABOLOMIC_LFC_CUT, METABOLOMIC_LFC_CUT),
    padj_cut = c(0.05, 0.05, PROTEOMIC_PADJ_CUT, PROTEOMIC_PADJ_CUT, METABOLOMIC_PADJ_CUT, METABOLOMIC_PADJ_CUT),
    stringsAsFactors = FALSE
  )
)

MULTIOMICS_CLUSTERING_IPS <- list(
  dataset = "iPSC-CM",
  out_dir = file.path(PATH.IPS_results, "04.Clustering_CT_CC"),
  output_prefix = "IPS_CT_CC",
  column_labels = c("iPSC-CM\n6135 C/T\nS/NS", "iPSC-CM\n6135-H1 C/C\nS/NS"),
  heatmap_width = 7,
  heatmap_height = 13,
  legend_width = 9,
  legend_height = 4,
  dpi = 600,
  comparisons = data.frame(
    omic = c(
      "Transcriptomic", "Transcriptomic",
      "Proteomic", "Proteomic",
      "Metabolomic", "Metabolomic"
    ),
    file = c(
      file.path(PATH.IPS_results, "01.Transcriptomic", "deg", "6135_S_vs_NS.tsv"),
      file.path(PATH.IPS_results, "01.Transcriptomic", "deg", "6135H1_S_vs_NS.tsv"),
      file.path(PATH.IPS_results, "02.Proteomic", "DEP", "tables", "6135CCC_S_vs_NS.tsv"),
      file.path(PATH.IPS_results, "02.Proteomic", "DEP", "tables", "6135H1_S_vs_NS.tsv"),
      file.path(PATH.IPS_results, "03.Metabolomic", "DE_tables", "Stim_vs_NS_6135.tsv"),
      file.path(PATH.IPS_results, "03.Metabolomic", "DE_tables", "Stim_vs_NS_6135-H1.tsv")
    ),
    contrast = c("CT", "CC", "CT", "CC", "CT", "CC"),
    lfc_cut = c(1.5, 1.5, PROTEOMIC_LFC_CUT, PROTEOMIC_LFC_CUT, METABOLOMIC_LFC_CUT, METABOLOMIC_LFC_CUT),
    padj_cut = c(0.05, 0.05, PROTEOMIC_PADJ_CUT, PROTEOMIC_PADJ_CUT, METABOLOMIC_PADJ_CUT, METABOLOMIC_PADJ_CUT),
    stringsAsFactors = FALSE
  )
)
