## ========== 05m_topCT_marker_heatmap.R ==========
## Purpose: heatmap of canonical marker averages per identity (CTs + OtherCOVID).

keep_ids <- names(which(table(Idents(seu)) >= 5))
seu_sub  <- subset(seu, idents = keep_ids)

panel <- intersect(c(
  "IL2RA","HLA-DRA","PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","TNFRSF9","ICOS","CXCL13",
  "GZMK","GZMB","PRF1","NKG7","GNLY","CCL5","IFNG",
  "CCR7","SELL","IL7R","TCF7","LEF1","BCL2","CXCR5","BCL6","FOXP3","IKZF2"
), rownames(seu_sub))

DefaultAssay(seu_sub) <- if ("RNA" %in% Assays(seu_sub)) "RNA" else DefaultAssay(seu_sub)
if (DefaultAssay(seu_sub)=="RNA" && !"data" %in% Layers(seu_sub[["RNA"]])) seu_sub <- NormalizeData(seu_sub, verbose=FALSE)

avg <- AverageExpression(seu_sub, features = panel, assays = DefaultAssay(seu_sub), slot = "data")[[1]]
avg_z <- t(scale(t(avg)))
col_fun <- circlize::colorRamp2(c(-1.5, 0, 1.5), c("#2166AC","white","#B2182B"))
ht <- ComplexHeatmap::Heatmap(avg_z, name="z", col=col_fun, cluster_rows=TRUE, cluster_columns=TRUE,
                              column_title="Identities (CTs vs OtherCOVID)")
png(file.path(fig_dir, "topCT_marker_heatmap.png"), 1400, 900, res=150); draw(ht); dev.off()
