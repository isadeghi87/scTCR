## ========== 06_liver_trm_genes.R ==========
## Liver-homing / tissue-resident memory (TRM) signatures in COVID-linked T cells
## Uses SCT assay for expression; integrated only for UMAP visuals if needed.

source("/home/isadeghi/projects/covid_kids/scripts/00_utils.R")
setwd("/home/isadeghi/projects/covid_kids/")
# ---- Setup
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
  library(scales); library(ggrepel); library(readr);library(ComplexHeatmap)
})

seu <- readRDS(file.path(out_dir, "merged_seurat_final.rds"))


# ---------- 1) choose assay & cells
seu_cd8 <- subset(seu, subset = celltype == "CD8")
DefaultAssay(seu_cd8) <- "SCT"
seu_cd8$covid_hit <- factor(seu_cd8$covid_hit, levels = c("COVID-hit","Other"))

# ---------- 2) your gene sets
gene_sets <- list(
  Chemokine_Receptors = c("CXCR3","CXCR6","CCR5","CCR6","CX3CR1"),
  Integrins_Adhesion  = c("ITGAL","ITGB2","ITGA4","ITGB1","ITGAE","ITGA1","CD69"),
  TRM_TFs             = c("PRDM1","ZNF683","RUNX3","TBX21","EOMES","BHLHE40"),
  Survival_Metabolic  = c("CD44","S1PR1","FABP4","FABP5"),
  Checkpoint          = c("TIGIT"),
  Effector_Cytokines  = c("GZMB","IFNG","TNF","IL17A","IL17F")
)

present <- rownames(seu_cd8[["SCT"]])
gene_sets <- lapply(gene_sets, function(v) intersect(v, present))
panel <- unique(unlist(gene_sets))

# ---------- 3) matrices: mean expression & % expressing
avg <- AverageExpression(seu_cd8, features = panel, assays = "SCT",
                         group.by = "covid_hit", slot = "data")$SCT
avg <- avg[, c("COVID-hit","Other"), drop = FALSE]                 # enforce order

# Row z-score and cap to keep colors comparable
avg_z <- t(scale(t(avg)))
avg_z[!is.finite(avg_z)] <- 0
avg_z <- pmin(pmax(avg_z, -2), 2)                                  # cap at ±2

mat <- GetAssayData(seu_cd8, assay = "SCT", slot = "data")[panel, ]
pct <- sapply(colnames(avg), function(g){
  cells <- colnames(seu_cd8)[seu_cd8$covid_hit == g]
  if (!length(cells)) return(rep(NA_real_, nrow(mat)))
  Matrix::rowMeans(mat[, cells, drop = FALSE] > 0) * 100
})
colnames(pct) <- colnames(avg)
pct <- pct[rownames(avg_z), , drop = FALSE]

# ---------- 4) row annotations: category stripes & row splitting
stacked <- utils::stack(gene_sets) # values=gene, ind=category
cat_by_gene <- setNames(stacked$ind, stacked$values)[rownames(avg_z)]

pal_cat <- c(
  Chemokine_Receptors="#1f77b4", Integrins_Adhesion="#2ca02c", TRM_TFs="#9467bd",
  Survival_Metabolic="#8c564b",  Checkpoint="#e377c2",         Effector_Cytokines="#ff7f0e"
)

ha_row <- rowAnnotation(
  Category = cat_by_gene,
  col = list(Category = pal_cat),
  gp = gpar(col = NA), width = unit(4, "mm")
)

# ---------- 5) column annotation (group colors)
pal_grp <- c("COVID-hit"="#E4572E", "Other"="#4C78A8")
ha_col <- HeatmapAnnotation(
  Group = factor(colnames(avg_z), levels = c("COVID-hit","Other")),
  col = list(Group = pal_grp)
)

# ---------- 6) colors
col_expr <- colorRamp2(c(-1, 0, 1), c("blue","white","red"))   # purple–white–orange
col_pct  <- colorRamp2(c(0, 50, 100), c("white", "#9ecae1", "#08519c"))

# ---------- 7) build heatmaps (split rows by category, cluster within split)
## --- 1) Build the expression heatmap (NOT drawn yet) ---
ht_expr <- Heatmap(
  avg_z,                              # z-scored expression (genes x groups)
  name = "Expr (z)",
  col  = col_expr,                    # purple → white → orange
  cluster_rows    = TRUE,
  cluster_columns = FALSE,
  row_split = factor(cat_by_gene, levels = names(gene_sets)),  # your category split
  top_annotation  = ha_col,           # column annotation (COVID-hit vs Other)
  left_annotation = ha_row,           # row-side category bars
  row_names_gp    = gpar(fontsize = 9),
  # column_title    = "CD8 T cells",
  border = TRUE,
  show_row_names = T,             # Hide row names (category names)
  row_title       = NULL,             # Hide the row title (category names)
)

## --- 2) Draw just to compute the final clustering order ---
ht_drawn <- draw(ht_expr)

## Extract the final row order (respecting row_split)
ord_list <- row_order(ht_drawn)                   # a list, one integer vector per split
ord_all  <- unlist(ord_list, use.names = FALSE)   # flatten
genes_ordered <- rownames(avg_z)[ord_all]         # names are safer than integer positions

## --- 3) Reorder the % positive matrix to match and build the 2nd heatmap ---
pct_reordered <- pct[genes_ordered, , drop = FALSE]

ht_pct <- Heatmap(
  pct_reordered,
  name = "Percent positive",
  col  = col_pct,                  # white → blue
  width = unit(12, "mm"),
  cluster_rows    = FALSE,         # order already enforced
  cluster_columns = FALSE,
  show_row_names  = T,
  border = TRUE
)

## --- 4) Now combine the *original* ht_expr with ht_pct and draw ONCE ---
ht_list <- ht_expr + ht_pct        # combine the not-drawn ht_expr

png(file.path(fig_dir, "TRM_liver_markers_CD8_heatmap_nice.png"), width = 1800, height = 1200, res = 150)
draw(ht_list,
     heatmap_legend_side     = "right",
     annotation_legend_side  = "right")

dev.off()

DefaultAssay(seu_cd8) <- "SCT"

p_dot <- DotPlot(
  seu_cd8, features = panel, assay = "SCT", group.by = "covid_hit"
) +
  ## 3-color mapping centered at 0 (scaled exp in DotPlot is z-like)
  scale_color_gradient2(
    low  = "#6A00A8",   # purple
    mid  = "white",
    high = "#F9844A",   # orange
    midpoint = 0
  ) +
  labs(title = "Liver-homing/TRM markers in CD8 T cells",
       color = "Avg exp (scaled)", size = "Percent expressed") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.minor = element_blank())

p_dot

ggsave(plot = p_dot, file.path(fig_dir, "TRM_liver_markers_CD8_dotplot.png"), width = 11, height = 7,dpi = 300)


p_dot <- DotPlot(
  seu_cd8, features = panel, assay = "SCT", group.by = "covid_hit"
) +
  ## 3-color mapping centered at 0 (scaled exp in DotPlot is z-like)
  scale_color_gradient2(
    low  = "#6A00A8",   # purple
    mid  = "white",
    high = "#F9844A",   # orange
    midpoint = 0
  ) +
  labs(title = "Liver-homing/TRM markers in CD8 T cells",
       color = "Avg exp (scaled)", size = "Percent expressed") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.minor = element_blank())

p_dot

ggsave(plot = p_dot, file.path(fig_dir, "TRM_liver_markers_CD8_dotplot.png"), width = 11, height = 7,dpi = 300)


# Stats for the violin
df_cd8 <- FetchData(seu_cd8, vars = c("score_LiverTRM1","covid_hit"))
wil_p  <- wilcox.test(score_LiverTRM1 ~ covid_hit, data = df_cd8)$p.value

write_tsv(tibble::tibble(test="CD8 Wilcoxon", feature="score_LiverTRM1", p_value=wil_p),
          file.path(tab_dir, "TRM_module_score_CD8_stats.tsv"))

###

for (genes in names(gene_sets)) {
  gene_panel = gene_sets[[genes]]
  
  # violin (left)
  p_vln <- VlnPlot(seu_cd8, features = gene_panel, 
                   assay = "SCT", 
                   group.by = "covid_hit",
                   pt.size = 0) +
    labs(title = genes, x = NULL, y = "Module score (SCT)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
  vplot_name = paste0(genes,'_vlnPlot.png')
  ggsave(plot = p_vln,filename = file.path(fig_dir,vplot_name),
         width = 8,height = 7,dpi =300)
  
  # UMAP (right): FeaturePlot can read meta columns via FetchData
  p_umap <- FeaturePlot(seu_cd8, features = gene_panel,split.by = 'covid_hit', reduction = "umap") +
    # scale_color_gradientn(colors = expr_cols, name = "Score") +
    theme_bw(base_size = 11) +
    theme(legend.position = "right") +
    ggtitle(genes)
  
  umap_name = paste0(genes,'_umapPlot.png')
  ggsave(plot = p_umap,filename = file.path(fig_dir,umap_name),
         width = 8,height = 9,dpi =300)
  
  
  # allplots <- p_vln + p_umap + plot_layout(widths = c(1, 1.2))
  
  # print(allplots)
}
