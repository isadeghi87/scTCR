## ========== 08_naive_Tcells.R ==========
## Naïve T-cell module scoring, UMAPs (purple–white–orange), and group-level stats

## ---- Project setup ----
## (edit to your local path if needed)
setwd("/home/isadeghi/projects/covid_kids")

source("scripts/00_utils.R")  # defines out_dir, fig_dir, tab_dir, sp(), etc.

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
  library(readr);  library(patchwork); library(stringr)
})

## ---- Load Seurat ----
seu <- readRDS(file.path(out_dir, "merged_seurat_final.rds"))

## ---- Define groups from covid_like ----
stopifnot("covid_like" %in% colnames(seu@meta.data))
seu$covid_hit <- ifelse(seu$covid_like %in% c("low","medium"), "COVID-hit", "Other")

## ---- Choose assay for expression (SCT) ----
DefaultAssay(seu) <- if ("SCT" %in% Assays(seu)) "SCT" else DefaultAssay(seu)

## ---- Naïve panel (RNA-level) from Biocompare page ----
## Positive markers (naïve-high)
naive_panel_pos <- c("CCR7","SELL","IL7R","TCF7","LEF1","FOXP1","CD27","CD28")
## Optional negatives (activation/memory-high) – used only for QC/plots
naive_panel_neg <- c("CD44","HLA-DRA","IL2RA")

## Harmonize to available genes (case-insensitive)
gene_lookup <- rownames(seu)
fix_genes <- function(v) {
  m <- match(toupper(v), toupper(gene_lookup))
  v[!is.na(m)] <- gene_lookup[m[!is.na(m)]]
  unique(v[!is.na(m)])
}
pos_use <- fix_genes(naive_panel_pos)
neg_use <- fix_genes(naive_panel_neg)

if (length(pos_use) < 3) {
  stop("Too few naïve markers found in object: ", paste(pos_use, collapse=", "))
}

## ---- Module score (positive-only; standard practice) ----
## Store as NaiveScore1
seu <- AddModuleScore(seu, features = list(pos_use), name = "NaiveScore", assay = DefaultAssay(seu))
col_score <- grep("^NaiveScore", colnames(seu@meta.data), value = TRUE)[1]
colnames(seu@meta.data)[colnames(seu@meta.data)==col_score] <- "NaiveScore1"

## Z-scale within object for visualization (center across all cells)
seu$NaiveScoreZ <- scale(seu$NaiveScore1)[,1]

## ---- Cluster-level summary (which clusters look naïve?) ----
naive_by_cluster <- seu@meta.data |>
  dplyr::mutate(seurat_clusters = as.character(seurat_clusters)) |>
  dplyr::group_by(seurat_clusters) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_NaiveScoreZ = mean(NaiveScoreZ, na.rm = TRUE)
  ) |>
  dplyr::arrange(dplyr::desc(mean_NaiveScoreZ))

readr::write_tsv(naive_by_cluster, file.path(tab_dir, "naive_score_by_cluster.tsv"))

## Simple barplot (positive values = more naïve-like)
p_bar <- ggplot(naive_by_cluster, aes(x = reorder(seurat_clusters, mean_NaiveScoreZ),
                                      y = mean_NaiveScoreZ,
                                      fill = mean_NaiveScoreZ > 0)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE"="#F97316","FALSE"="grey70"), guide = "none") +
  labs(x = "Seurat cluster", y = "Mean NaïveScore (Z)", title = "Naïve score by cluster") +
  theme_bw(base_size = 12)
ggsave(plot = p_bar, filename = file.path(fig_dir, "naive_score_by_cluster.png"), width = 9, height = 6)

## ---- Violin: naïve score by lineage and group ----
Idents(seu) <- seu$celltype
for (lin in intersect(c("CD4","CD8"), unique(seu$celltype))) {
  cells_lin <- colnames(seu)[seu$celltype==lin]
  if (length(cells_lin) < 50) next
  p <- VlnPlot(subset(seu, cells = cells_lin),
               features = "NaiveScoreZ", group.by = "group",
               pt.size = 0.2, cols = c("tomato","grey60")) +
    labs(title = paste("Naïve signature (", lin, ")", sep=""),
         y = "Naïve module score (Z)", x = NULL) +
    theme_bw(base_size = 12)
  print(p)
  ggsave(plot = p, filename = file.path(fig_dir, paste0("naive_signature_", lin, ".png")), 8, 5)
  
  ## Wilcoxon test (COVID-hit vs Other)
  df <- FetchData(subset(seu, cells = cells_lin), vars = c("NaiveScoreZ","group"))
  wt <- tryCatch(wilcox.test(NaiveScoreZ ~ group, data = df), error = function(e) NULL)
  if (!is.null(wt)) {
    tibble::tibble(test="Wilcoxon", lineage=lin, p_value=wt$p.value, n=nrow(df)) |>
      readr::write_tsv(file.path(tab_dir, paste0("naive_signature_stats_", lin, ".tsv")))
  }
}

## ---- UMAP FeaturePlots with purple → white → orange ----
## What the scale bar shows: SCT expression (log-normalized residuals).
umap_feats <- unique(c(pos_use, neg_use))
if ("umap" %in% names(seu@reductions)) {
  DefaultAssay(seu) <- "SCT"  # use SCT values for color
  gp <- FeaturePlot(seu, features = umap_feats, combine = FALSE, order = TRUE,
                    min.cutoff = "q01", max.cutoff = "q99")
  gp <- lapply(gp, function(g) {
    g + scale_colour_gradientn(
      colours = c("#5e3c99", "white", "#e66101")   # purple → white → orange
    ) +
      labs(colour = "Expr (SCT)") +
      theme(plot.title = element_text(face="bold"))
  })
  ## Save in pages of 6
  n <- length(gp); k <- 6
  for (i in seq(1, n, by = k)) {
    j <- min(i+k-1, n)
    pg <- patchwork::wrap_plots(gp[i:j], ncol = 2)
    print(pg)
    ggsave(plot = pg, filename = file.path(fig_dir, sprintf("naive_markers_umap_%02d.png", (i-1)%/%k + 1)), width = 8, height = 6.5)
  }
}

pos_vln = VlnPlot(seu,features = naive_panel_pos)
ggsave(plot = pos_vln, filename = file.path(fig_dir, paste0("naive_signature_positive_vlnPlot.png")), width = 8, height = 7)

cov_pos = VlnPlot(seu,features = naive_panel_neg, group.by = 'covid_hit')
ggsave(plot = cov_pos, filename = file.path(fig_dir, paste0("naive_signature_positive_covid_hit_vlnPlot.png")), width = 8, height = 5)

## ---- Per-lineage percent naïve (simple threshold on Z>0) ----
seu$naive_high <- seu$NaiveScoreZ > 0
tab_lin <- seu@meta.data |>
  dplyr::count(celltype, group, naive_high, name="n") |>
  dplyr::group_by(celltype, group) |>
  dplyr::mutate(pct = 100*n/sum(n)) |>
  dplyr::ungroup()

readr::write_tsv(tab_lin, file.path(tab_dir, "naive_high_by_lineage_group.tsv"))

seu$naive_high <- seu$NaiveScoreZ > 0

## -------- Overall 2x2 Fisher --------
tab_overall <- table(seu$naive_high, seu$covid_hit)
ft_overall  <- fisher.test(tab_overall)  # enrichment of naive_high in COVID-hit vs Other
OR_overall  <- unname(ft_overall$estimate)

print(tab_overall)
cat(sprintf("Overall OR = %.3f, p = %.3g\n", OR_overall, ft_overall$p.value))


## -------- Plot % naïve-high by lineage (95% Wilson CI), grouped by covid_hit --------
library(dplyr); library(binom); library(ggplot2)
p_df = as.data.frame(tab_overall)
colnames(p_df) = c('naïve','covid_hit','Freq')

p <- ggplot(p_df, aes(x = covid_hit, y = Freq, fill = naïve)) +
  geom_col(position = position_dodge(width=.7), width=.6, color="black") +
  labs(x=NULL, y="Number of cells",
       title="Naïve-high by lineage)",
       subtitle=sprintf("Overall OR≈%.2f (Fisher p=%.2g)", OR_overall, ft_overall$p.value)) +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "naive_high_by_lineage_covidHIT.png"), p, width=7, height=5, dpi=300)
