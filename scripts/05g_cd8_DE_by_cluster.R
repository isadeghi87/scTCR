## ========== 05g_cd8_DE_by_cluster.R ==========
## Purpose: DE for CD8 medium vs none within each cluster with enough cells.

Idents(seu) <- seu$celltype
bcs_med_cd8 <- colnames(seu)[seu$celltype=="CD8" & seu$covid_like=="medium"]
if (length(bcs_med_cd8) > 0) {
  dom_by_cluster <- sort(table(seu$seurat_clusters[bcs_med_cd8]), decreasing=TRUE)
  for (cl in names(dom_by_cluster)) {
    b_med  <- colnames(seu)[seu$celltype=="CD8" & seu$covid_like=="medium" & seu$seurat_clusters==cl]
    b_ctrl <- colnames(seu)[seu$celltype=="CD8" & seu$covid_like=="none"   & seu$seurat_clusters==cl]
    if (length(b_med) < 10 || length(b_ctrl) < 20) next
    seu_sub <- subset(seu, cells=c(b_med, b_ctrl))
    seu_sub$group <- ifelse(colnames(seu_sub) %in% b_med, "medium", "none")
    Idents(seu_sub) <- seu_sub$group
    markers <- FindMarkers(seu_sub, "medium", "none", test.use="wilcox", logfc.threshold=.25, min.pct=.1)
    if (!nrow(markers)) next
    markers <- markers |> tibble::rownames_to_column("gene") |> dplyr::mutate(cluster=cl, p_adj=p.adjust(p_val, "BH")) |>
      dplyr::arrange(p_adj, dplyr::desc(avg_log2FC))
    outf <- file.path(tab_dir, paste0("DE_medium_vs_none_CD8_cluster", cl, ".tsv"))
    readr::write_tsv(markers, outf)
    p <- ggplot2::ggplot(markers, ggplot2::aes(avg_log2FC, -log10(p_adj))) + ggplot2::geom_point(alpha=.6) +
      ggplot2::labs(title=paste0("CD8 medium vs none (cluster ", cl, ")")) + ggplot2::theme_bw()
    sp(p, paste0("volcano_CD8_medium_cluster", cl, ".png"), 6.5, 5)
  }
}
