## ========== 05k_topCT_dotplot_and_umap.R ==========
## Purpose: DotPlot panel for canonical markers across top clonotypes + UMAPs for covid/epitopes.

cov_cells_keep <- colnames(seu)[seu$covid_hit=="COVID-hit"]
assign_vec <- setNames(rep("OtherCOVID", ncol(seu)), colnames(seu))
for (i in seq_len(nrow(topN))) {
  bc <- colnames(seu)[seu$CTaa == topN$clonotype[i] & colnames(seu) %in% cov_cells_keep]
  if (length(bc)) assign_vec[bc] <- sprintf("CT:%02d", i)
}
Idents(seu) <- factor(assign_vec, levels=c(sprintf("CT:%02d", seq_len(nrow(topN))), "OtherCOVID"))

panel <- intersect(c(
  "IL2RA","HLA-DRA","PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","TNFRSF9","ICOS","CXCL13",
  "GZMK","GZMB","PRF1","NKG7","GNLY","CCL5","IFNG",
  "CCR7","SELL","IL7R","TCF7","LEF1","BCL2","CXCR5","BCL6","FOXP3","IKZF2"
), rownames(seu))

DefaultAssay(seu) <- if ("RNA" %in% Assays(seu)) "RNA" else DefaultAssay(seu)
if (DefaultAssay(seu)=="RNA" && !"data" %in% Layers(seu[["RNA"]])) seu <- NormalizeData(seu, verbose=FALSE)

p_dot <- DotPlot(seu, features=panel, cols=c("lightgrey","firebrick")) +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  ggtitle(sprintf("Canonical markers in top COVID clonotypes (%s)", CONF_LEVEL))
sp(p_dot, sprintf("covid_%s_topCT_dotplot.png", CONF_LEVEL), 12, 6)

# UMAPs
if ("umap" %in% names(seu@reductions)) {
  sp(DimPlot(seu, group.by="covid_hit") + ggtitle(sprintf("COVID-hit (%s) vs Other", CONF_LEVEL)),
     sprintf("covid_%s_umap_hit.png", CONF_LEVEL), 6, 5)

  cov_ep <- readr::read_tsv(file.path(tab_dir, "covid_best_per_cell_CD4_CD8.tsv"), show_col_types=FALSE) %>% dplyr::select(barcode, epitope) %>% dplyr::distinct()
  seu$covid_epitope <- NA_character_
  idx <- match(cov_ep$barcode, colnames(seu)); seu$covid_epitope[idx[!is.na(idx)]] <- cov_ep$epitope
  top_epi <- seu@meta.data %>% dplyr::filter(covid_hit=="COVID-hit", !is.na(covid_epitope)) %>%
    dplyr::count(covid_epitope, sort=TRUE) %>% dplyr::slice_head(n=5) %>% dplyr::pull(covid_epitope)
  seu$epi_top5 <- ifelse(!is.na(seu$covid_epitope) & seu$covid_epitope %in% top_epi, seu$covid_epitope, "other")
  sp(DimPlot(seu, group.by="epi_top5") + ggtitle("Top COVID epitopes (UMAP)"),
     sprintf("covid_%s_umap_top_epitopes.png", CONF_LEVEL), 6, 5)
}
