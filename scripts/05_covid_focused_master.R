## ========== 05_covid_focused_master.R ==========
## Driver script: sets configuration, loads Seurat object, and sources modules in order.

## Choose which utils to use:
##  (A) Your project utils
##  source("/home/isadeghi/projects/covid_kids/scripts/00_utils.R")
##  (B) Local minimal utils (included here)
source("00_utils_local.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(stringdist)
  library(igraph); library(ggraph); library(ComplexHeatmap); library(FNN)
  library(Seurat); library(tibble); library(scales); library(ggrepel)
  library(patchwork); library(vegan); library(readr); library(circlize)
  library(ggridges); library(effsize); library(tidytext)
})

## ---- Global config (edit) ----
CONF_LEVEL <- "low"   # choose: "low" | "medium" | "high"

## ---- Load merged Seurat ----
## If you used project utils, out_dir/fig_dir/tab_dir come from there.
## Otherwise, 00_utils_local.R created them.
seu <- readRDS(file.path(out_dir, "merged_seurat_final.rds"))
stopifnot(inherits(seu, "Seurat"))

## ---- Source modules ----
source("05a_define_covid_hit.R")
source("05b_build_topN_and_tables.R")
source("05c_epitope_annotation_and_mix.R")
source("05d_DE_medium_by_lineage.R")
source("05e_tcr_beta_similarity_network.R")
source("05f_effect_sizes_OR.R")
source("05g_cd8_DE_by_cluster.R")
source("05h_clonal_expansion_violin.R")
source("05i_vj_skew_enrichment.R")
source("05j_trbvj_delta_heatmap.R")
source("05k_topCT_dotplot_and_umap.R")
source("05l_spectratype_ridges_and_ecdf.R")
source("05m_topCT_marker_heatmap.R")
source("05n_ct_epitope_barplot.R")

message("COVID-focused pipeline complete for CONF_LEVEL = ", CONF_LEVEL)
