## ========== 05d_DE_medium_by_lineage.R ==========
## Purpose: differential expression for "medium" tier within lineage & dominant cluster.

Idents(seu) <- seu$celltype
de_medium_by_lineage <- function(lineage, min.cells = 10) {
  bcs_med <- colnames(seu)[seu$celltype == lineage & seu$covid_like == "medium"]
  if (length(bcs_med) < min.cells) return(NULL)
  dom_cluster <- names(sort(table(seu$seurat_clusters[bcs_med]), decreasing=TRUE))[1]
  bcs_ctrl <- colnames(seu)[seu$celltype == lineage & seu$seurat_clusters == dom_cluster & seu$covid_like == "none"]
  if (length(bcs_ctrl) < min.cells) return(NULL)
  seu_sub <- subset(seu, cells = c(bcs_med, bcs_ctrl))
  seu_sub$covid_tier <- ifelse(colnames(seu_sub) %in% bcs_med, "medium", "none")
  Idents(seu_sub) <- seu_sub$covid_tier
  markers <- tryCatch(FindMarkers(seu_sub, "medium","none", test.use="wilcox", logfc.threshold=.25, min.pct=.1), error=function(e) NULL)
  if (is.null(markers) || nrow(markers)==0) return(NULL)
  markers <- markers |> tibble::rownames_to_column("gene") |>
    dplyr::mutate(lineage=lineage, cluster=dom_cluster, p_adj=p.adjust(p_val, "BH")) |>
    dplyr::arrange(p_adj, dplyr::desc(avg_log2FC))
  out_file <- file.path(tab_dir, paste0("DE_medium_vs_none_", lineage, "_cluster", dom_cluster, ".tsv"))
  readr::write_tsv(markers, out_file)
  p <- ggplot2::ggplot(markers, ggplot2::aes(avg_log2FC, -log10(p_adj))) + ggplot2::geom_point(alpha=.6) +
    ggrepel::geom_text_repel(data=head(markers, 15), ggplot2::aes(label=gene), max.overlaps=50, size=3) +
    ggplot2::labs(title=paste0("Medium vs None (", lineage, ", cluster ", dom_cluster, ")")) + ggplot2::theme_bw()
  sp(p, paste0("volcano_medium_", lineage, "_cluster", dom_cluster, ".png"), 6.5, 5)
  markers
}
res_cd4 <- de_medium_by_lineage("CD4")
res_cd8 <- de_medium_by_lineage("CD8")
