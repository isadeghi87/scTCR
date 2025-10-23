## ======== 08_autoaggression_signature.R (updated) ========
## Expression of human auto-aggression signature in COVID-hit vs Other
## -------------------------------------------------------------------
## Inputs:
##  - object: out/merged_seurat_final.rds
##  - Excel gene list: tabs/Top70 genes - human auto-aggression gene signature.xlsx
## Outputs:
##  figs/
##    autoagg_module_score_violin_CD4.png
##    autoagg_module_score_violin_CD8.png
##    autoagg_signature_umap_diverging.png
##    autoagg_sentinels_umap.png
##  tabs/
##    autoagg_module_score_stats.tsv
##    autoagg_marker_list_used.tsv
## -------------------------------------------------------------------

setwd("/home/isadeghi/projects/covid_kids")
source(file.path("scripts", "00_utils.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(readxl); library(patchwork)
})

seu <- readRDS(file.path(out_dir, "merged_seurat_final.rds"))

## Use SCT for scoring/violin; switch to 'integrated' only for UMAP drawing
DefaultAssay(seu) <- if ("SCT" %in% Assays(seu)) "SCT" else DefaultAssay(seu)

## COVID-hit status (LOW confidence = best_per)
best_per <- readr::read_tsv(file.path(tab_dir, "covid_best_per_cell_CD4_CD8.tsv"), show_col_types = FALSE)
cov_barcodes <- intersect(best_per$barcode, colnames(seu))
seu$covid_hit <- ifelse(colnames(seu) %in% cov_barcodes, "COVID-hit", "Other")

## Read Excel list (first sheet, first column)
auto_xlsx <- file.path(tab_dir, "Top70 genes - human auto-aggression gene signature.xlsx")
stopifnot(file.exists(auto_xlsx))
auto_tbl   <- readxl::read_xlsx(auto_xlsx, sheet = 1)
auto_genes <- as.character(auto_tbl[[1]])
auto_genes <- sort(unique(intersect(auto_genes, rownames(seu))))
if (length(auto_genes) < 5) warning("Auto-aggression list small after intersection: ", length(auto_genes))
write_tsv(tibble(gene=auto_genes), file.path(tab_dir, "autoagg_marker_list_used.tsv"))

## Score (module)
seu <- AddModuleScore(seu, features = list(auto_genes), name = "score_AutoAgg")
score_col <- grep("^score_AutoAgg", colnames(seu@meta.data), value=TRUE)[1]
seu$AutoAggZ <- scale(seu@meta.data[[score_col]])[,1]   # centered at 0

## Violin by lineage
seu_cd4 <- subset(seu, subset = celltype == "CD4")
seu_cd8 <- subset(seu, subset = celltype == "CD8")

plot_vln <- function(obj, col, title){
  ggplot(obj@meta.data, aes(x=covid_hit, y=.data[[col]], fill=covid_hit)) +
    geom_violin(scale="width", trim=TRUE, alpha=.85, color=NA) +
    geom_boxplot(width=.12, outlier.size=.4, fill="white") +
    scale_fill_manual(values=c("Other"="grey70","COVID-hit"="tomato")) +
    labs(x=NULL, y="Auto-aggression module score", title=title) +
    theme_bw(base_size=11) + theme(legend.position="none")
}

p_cd4 <- plot_vln(seu_cd4, score_col, "Auto-aggression signature (CD4)")
p_cd8 <- plot_vln(seu_cd8, score_col, "Auto-aggression signature (CD8)")
ggsave(plot = p_cd4, file.path(fig_dir, "autoagg_module_score_violin_CD4.png"), width = 5.2, height = 3.8,dpi =300)
ggsave(plot = p_cd8, file.path(fig_dir, "autoagg_module_score_violin_CD8.png"), width = 5.2, height = 3.8,dpi = 300)


## Stats
w_cd4 <- wilcox.test(seu_cd4@meta.data[[score_col]] ~ seu_cd4$covid_hit)
w_cd8 <- wilcox.test(seu_cd8@meta.data[[score_col]] ~ seu_cd8$covid_hit)
bind_rows(
  tibble(test="Wilcoxon", lineage="CD4", p_value=w_cd4$p.value, n=ncol(seu_cd4)),
  tibble(test="Wilcoxon", lineage="CD8", p_value=w_cd8$p.value, n=ncol(seu_cd8))
) |> write_tsv(file.path(tab_dir, "autoagg_module_score_stats.tsv"))

## UMAP of the z-scored signature (diverging; centered at 0)
if ("umap" %in% names(seu@reductions)) {
  DefaultAssay(seu) <- "integrated"
  p_umap_sig <- FeaturePlot(seu, features="AutoAggZ", pt.size=0.3) +
    scale_color_gradient2(low="#6A51A3", mid="white", high="#F16913", midpoint=0) +
    labs(title="Auto-aggression signature (z-score)")
  ggsave(plot = p_umap_sig, file.path(fig_dir, "autoagg_signature_umap_diverging.png"), 
  width = 6.5, height = 5.2,dpi = 300)

  # Optional: UMAP of sentinels (not z-scored, standard FeaturePlot colors)
  sentinels <- intersect(c("CX3CR1","GZMB","PRF1","NKG7","TBX21","IFNG","ITGAL","TIGIT"), rownames(seu))
  if (length(sentinels) > 0) {
    p_sentinels <- FeaturePlot(seu, features=sentinels, combine=TRUE) &
      theme(plot.title = element_text(size=9))
    ggsave(plot = p_sentinels, file.path(fig_dir, "autoagg_sentinels_umap.png"),
    width =  8, height = 6,dpi = 300)
  }
  DefaultAssay(seu) <- "SCT"
}

message("08_autoaggression_signature.R done.")
