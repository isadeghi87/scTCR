## ========== 05_covid_focused.R ==========
source("/home/isadeghi/projects/covid_kids/scripts/00_utils.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(stringdist)
  library(igraph); library(ggraph); library(ComplexHeatmap); library(FNN)
  library(Seurat); library(tibble); library(scales)
  library(ggrepel);library(patchwork)
  library(vegan)
})

## ---- Load merged Seurat and COVID hit tables ----
seu <- readRDS(file.path(out_dir, "merged_seurat_final.rds"))

## Confidence level → which cells count as “COVID-hit”
## choose one of: "high" (paired same epitope) / "medium" (paired any) / "low" (best single-chain)
CONF_LEVEL <- "medium"

# Build cov_set from your saved hit tables
best_per    <- readr::read_tsv(file.path(tab_dir, "covid_best_per_cell_CD4_CD8.tsv"), show_col_types = FALSE)
paired_any  <- readr::read_tsv(file.path(tab_dir, "covid_paired_any_CD4_CD8.tsv"), show_col_types = FALSE)
paired_same <- readr::read_tsv(file.path(tab_dir, "covid_paired_same_epitope_CD4_CD8.tsv"), show_col_types = FALSE)

cov_cells_low    <- best_per$barcode
cov_cells_medium <- paired_any$barcode[paired_any$paired_any]
cov_cells_high   <- paired_same$barcode[paired_same$paired_same_epitope]

## What you have now (LOW confidence: best_per = single-chain best hit)
table(seu$covid_hit)
addmargins(table(seu$covid_hit, seu$celltype))

## If you want MEDIUM / HIGH instead, flip the source of barcodes:
# medium = paired_any$paired_any; high = paired_same$paired_same_epitope
CONF_LEVEL <- "low"  # "low" | "medium" | "high"
cov_barcodes <- switch(
  CONF_LEVEL,
  low    = intersect(best_per$barcode, colnames(seu)),
  medium = intersect(paired_any$barcode[paired_any$paired_any], colnames(seu)),
  high   = intersect(paired_same$barcode[paired_same$paired_same_epitope], colnames(seu))
)
seu$covid_hit <- ifelse(colnames(seu) %in% cov_barcodes, "COVID-hit", "Other")


## --- Build topN (top COVID clonotypes) robustly ---

# Make sure we know which barcodes are COVID-hit
stopifnot("CTaa" %in% colnames(seu@meta.data))
cov_barcodes <- intersect(best_per$barcode, colnames(seu))
seu$covid_hit <- ifelse(colnames(seu) %in% cov_barcodes, "COVID-hit", "Other")

# 2) Choose a clonotype label to count (prefer CTaa; fall back if needed)
clono_col <- dplyr::case_when(
  "CTaa"     %in% colnames(seu@meta.data) ~ "CTaa",
  "clone_id" %in% colnames(seu@meta.data) ~ "clone_id",
  TRUE                                     ~ "TRB_cdr3"   # last-resort fallback
)

# 3) Collect COVID-hit cells with a non-empty clonotype label
cov_meta <- seu@meta.data |>
  tibble::rownames_to_column("barcode") |>
  dplyr::filter(barcode %in% cov_barcodes,
                !is.na(.data[[clono_col]]),
                .data[[clono_col]] != "")

if (nrow(cov_meta) == 0) {
  warning("No COVID-hit cells with a non-empty clonotype label — cannot build topN.")
  topN <- tibble::tibble(clonotype = character(), n_cells = integer())
} else {
  # 4) Count and pick the top 20 clonotypes
  topN <- cov_meta |>
    dplyr::count(.data[[clono_col]], name = "n_cells") |>
    dplyr::arrange(dplyr::desc(n_cells)) |>
    dplyr::slice_head(n = 20) |>
    dplyr::rename(clonotype = 1)
  
  # 5) Save for downstream steps that expect topN
  readr::write_tsv(topN, file.path(tab_dir, "covid_top20_clonotypes.tsv"))
}

# (Optional) quick peek
print(head(topN, 10))

# quick tables
addmargins(table(seu$covid_like, seu$celltype))
round(100*prop.table(table(seu$covid_like, seu$celltype), 2), 1)

addmargins(table(seu$covid_like, seu$seurat_clusters))
round(100*prop.table(table(seu$covid_like, seu$seurat_clusters), 2), 1)

# enrichment by lineage and cluster (Fisher tests)
enrich <- function(fac){
  tab <- table(seu$covid_like != "none", fac)
  data.frame(
    level = colnames(tab),
    p = apply(tab, 2, \(col) fisher.test(matrix(c(col[2], sum(col)-col[2], col[1], sum(col)-col[1]), 2))$p.value)
  ) |> dplyr::mutate(q = p.adjust(p, "BH")) |> dplyr::arrange(q)
}
enr_lineage <- enrich(seu$celltype)
enr_cluster <- enrich(seu$seurat_clusters)
readr::write_tsv(enr_lineage, file.path(tab_dir, "enrichment_lineage_covid_like.tsv"))
readr::write_tsv(enr_cluster, file.path(tab_dir, "enrichment_cluster_covid_like.tsv"))

#===================

# load the hits we saved in 03_ script
hits <- readr::read_tsv(file.path(tab_dir, "covid_hits_CD4_CD8.tsv"), show_col_types = FALSE)

# pick 1 epitope per barcode with a simple priority: TRB exact > TRA exact
prio <- dplyr::mutate(hits, priority = dplyr::case_when(
  match_chain=="TRB" & match_type=="exact" ~ 3L,
  match_chain=="TRA" & match_type=="exact" ~ 2L,
  TRUE ~ 0L
))
best_epi <- prio |>
  dplyr::group_by(barcode) |>
  dplyr::slice_max(priority, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(barcode, epitope, match_chain, match_type)

seu$covid_epitope <- NA_character_
m <- match(best_epi$barcode, colnames(seu))
seu$covid_epitope[m[!is.na(m)]] <- best_epi$epitope

# epitope mix in low vs medium, split by lineage
epi_mix <- seu@meta.data |>
  tibble::rownames_to_column("barcode") |>
  dplyr::filter(covid_like %in% c("low","medium"), !is.na(covid_epitope)) |>
  dplyr::count(celltype, covid_like, covid_epitope, name="n") |>
  dplyr::group_by(celltype, covid_like) |>
  dplyr::mutate(pct = 100*n/sum(n)) |> dplyr::ungroup()

readr::write_tsv(epi_mix, file.path(tab_dir, "covid_epitope_mix_by_lineage_tier.tsv"))

#==========plots 
ggplot(epi_mix, aes(x=reorder(covid_epitope, pct), y=pct, fill=covid_like)) +
  geom_col(position="dodge") + coord_flip() + facet_wrap(~celltype, scales="free_y") +
  labs(x="Epitope", y="% of COVID-like", title="COVID epitopes by tier & lineage") +
  theme_bw(base_size=10) -> p_epi
sp(p_epi, "covid_epitope_mix.png", 9, 6)


## Marker genes for medium (paired exact) — controlled comparison

## We’ll do DE within lineage (CD4 or CD8) and dominant cluster of the medium cells to reduce state confounding.
Idents(seu) <- seu$celltype

de_medium_by_lineage <- function(lineage, min.cells = 10) {
  bcs_med <- colnames(seu)[seu$celltype == lineage & seu$covid_like == "medium"]
  if (length(bcs_med) < min.cells) return(NULL)
  dom_cluster <- names(sort(table(seu$seurat_clusters[bcs_med]), decreasing=TRUE))[1]
  
  # comparison set: same lineage & cluster, covid_like == "none"
  bcs_ctrl <- colnames(seu)[seu$celltype == lineage &
                              seu$seurat_clusters == dom_cluster &
                              seu$covid_like == "none"]
  if (length(bcs_ctrl) < min.cells) return(NULL)
  
  seu_sub <- subset(seu, cells = c(bcs_med, bcs_ctrl))
  seu_sub$covid_tier <- ifelse(colnames(seu_sub) %in% bcs_med, "medium", "none")
  Idents(seu_sub) <- seu_sub$covid_tier
  
  markers <- tryCatch({
    FindMarkers(seu_sub, ident.1 = "medium", ident.2 = "none",
                test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1)
  }, error = function(e) NULL)
  
  if (is.null(markers) || nrow(markers)==0) return(NULL)
  
  markers <- markers |> tibble::rownames_to_column("gene") |>
    dplyr::mutate(lineage = lineage, cluster = dom_cluster,
                  p_adj = p.adjust(p_val, method="BH")) |>
    dplyr::arrange(p_adj, dplyr::desc(avg_log2FC))
  
  out_file <- file.path(tab_dir, paste0("DE_medium_vs_none_", lineage, "_cluster", dom_cluster, ".tsv"))
  readr::write_tsv(markers, out_file)
  
  # quick volcano
  top_n <- 15
  
  # protect from p_adj==0
  markers$padj_plot <- pmax(markers$p_adj, .Machine$double.xmin)
  
  markers <- markers |>
    dplyr::arrange(p_adj) |>
    dplyr::mutate(
      is_top = dplyr::row_number() <= top_n,
      dir    = dplyr::case_when(
        is_top & avg_log2FC >  0 ~ "top_up",
        is_top & avg_log2FC <= 0 ~ "top_down",
        TRUE                      ~ "other"
      ),
      label  = ifelse(is_top, gene, NA_character_)
    )
  
  ggplot(markers, aes(x = avg_log2FC, y = -log10(padj_plot))) +
    geom_point(aes(color = dir), alpha = .75, size = 1.8) +
    ggrepel::geom_text_repel(
      data = subset(markers, is_top),
      aes(label = label),
      max.overlaps = 50, size = 3, seed = 1
    ) +
    scale_color_manual(
      values = c(other = "grey70", top_up = "#D55E00", top_down = "#0072B2"),
      guide = "none"
    ) +
    labs(
      x = "avg_log2FC",
      y = "-log10(FDR)",
      title = sprintf("Medium vs None (%s, cluster %s)", lineage, dom_cluster)
    ) +
    theme_bw() -> p_volc
  sp(p_volc, paste0("volcano_medium_", lineage, "_cluster", dom_cluster, ".png"), 6.5, 5)
  
  markers
}

res_cd4 <- de_medium_by_lineage("CD4")
res_cd8 <- de_medium_by_lineage("CD8")

### TCR community among COVID-like (β CDR3 similarity graph)

suppressPackageStartupMessages({ library(igraph); library(stringdist) })

cov_betas <- seu@meta.data |>
  tibble::rownames_to_column("barcode") |>
  dplyr::filter(covid_like %in% c("low","medium"), !is.na(TRB_cdr3)) |>
  dplyr::distinct(TRB_cdr3, .keep_all=TRUE) |>
  dplyr::pull(TRB_cdr3)

# connect if same length and Levenshtein <= 1
edges <- lapply(seq_along(cov_betas), function(i){
  q <- cov_betas[i]
  d <- stringdist(q, cov_betas, method="lv")
  idx <- which(d == 1 & nchar(cov_betas)==nchar(q))
  if (!length(idx)) return(NULL)
  data.frame(from=q, to=cov_betas[idx])
}) |> dplyr::bind_rows()

if (!is.null(edges) && nrow(edges)>0) {
  g <- graph_from_data_frame(edges, directed=FALSE)
  comps <- igraph::components(g)$membership
  net_tab <- data.frame(TRB_cdr3 = names(comps), community = comps)
  readr::write_tsv(net_tab, file.path(tab_dir, "tcr_beta_lev1_communities.tsv"))
}

### Report effect sizes (odds ratios + CIs)
# lineage ORs
or_lineage <- lapply(split(seu@meta.data, seu$celltype), function(df){
  tab <- table(df$covid_like != "none")
  # success = covid-like
  m <- matrix(c(tab["TRUE"], tab["FALSE"], 
                sum(seu$covid_like!="none")-tab["TRUE"], 
                sum(seu$covid_like=="none")-tab["FALSE"]), nrow=2, byrow=TRUE)
  ft <- fisher.test(m)
  data.frame(level = unique(df$celltype), OR = unname(ft$estimate),
             CI_lo = ft$conf.int[1], CI_hi = ft$conf.int[2], p = ft$p.value)
}) |> dplyr::bind_rows()
readr::write_tsv(or_lineage, file.path(tab_dir, "covid_like_lineage_OR.tsv"))

# cluster ORs (top 10 by OR)
or_cluster <- lapply(split(seu@meta.data, seu$seurat_clusters), function(df){
  tab <- table(df$covid_like != "none")
  m <- matrix(c(tab["TRUE"], tab["FALSE"],
                sum(seu$covid_like!="none")-tab["TRUE"],
                sum(seu$covid_like=="none")-tab["FALSE"]), 2, byrow=TRUE)
  ft <- fisher.test(m)
  data.frame(cluster = unique(df$seurat_clusters), OR = unname(ft$estimate),
             CI_lo = ft$conf.int[1], CI_hi = ft$conf.int[2], p = ft$p.value)
}) |> dplyr::bind_rows() |> dplyr::arrange(dplyr::desc(OR))
readr::write_tsv(or_cluster, file.path(tab_dir, "covid_like_cluster_OR.tsv"))

##=========== CD8-focused DE for paired exact (your “medium”)
# Compare medium CD8 to CD8, same cluster, covid_like==none (controls). If a cluster has <10 medium cells, it will skip it.
## Tip: expect activation/cytotoxic signals (e.g., GZMB/PRF1/IFNG/CCL5), depending on which clusters your medium cells sit in.

Idents(seu) <- seu$celltype
bcs_med_cd8 <- colnames(seu)[seu$celltype=="CD8" & seu$covid_like=="medium"]
if (length(bcs_med_cd8) > 0) {
  # run per dominant cluster among medium CD8
  dom_by_cluster <- sort(table(seu$seurat_clusters[bcs_med_cd8]), decreasing=TRUE)
  for (cl in names(dom_by_cluster)) {
    b_med <- colnames(seu)[seu$celltype=="CD8" & seu$covid_like=="medium" & seu$seurat_clusters==cl]
    if (length(b_med) < 10) next
    b_ctrl <- colnames(seu)[seu$celltype=="CD8" & seu$covid_like=="none" & seu$seurat_clusters==cl]
    if (length(b_ctrl) < 20) next
    
    seu_sub <- subset(seu, cells = c(b_med, b_ctrl))
    seu_sub$group <- ifelse(colnames(seu_sub) %in% b_med, "medium", "none")
    Idents(seu_sub) <- seu_sub$group
    
    markers <- FindMarkers(seu_sub, ident.1="medium", ident.2="none",
                           test.use="wilcox", logfc.threshold=0.25, min.pct=0.1)
    if (nprow(markers)==0) next
    print(cl)
    markers <- markers |> tibble::rownames_to_column("gene") |>
      mutate(cluster=cl, p_adj=p.adjust(p_val, "BH")) |>
      arrange(p_adj, desc(avg_log2FC))
    outf <- file.path(tab_dir, paste0("DE_medium_vs_none_CD8_cluster", cl, ".tsv"))
    readr::write_tsv(markers, outf)
    
    # quick volcano
    p <- ggplot(markers, aes(avg_log2FC, -log10(p_adj))) + geom_point(alpha=.6) +
      labs(title=paste0("CD8 medium vs none (cluster ", cl, ")")) + theme_bw()
    sp(p, paste0("volcano_CD8_medium_cluster", cl, ".png"), 6.5, 5)
  }
}


#========== Are COVID-linked cells more expanded?
##============clonal expression of covid like cells
# compare clone_size distributions
df_exp <- seu@meta.data |>
  tibble::rownames_to_col
umn("barcode") |>
  dplyr::mutate(clone_size = as.integer(ifelse(is.na(clone_id), 0, ave(clone_id, clone_id, FUN=length)))) |>
  dplyr::mutate(is_covid = factor(covid_like %in% c("low","medium"), levels=c(FALSE,TRUE), labels=c("none","covid-like")))

wilcox.test(clone_size ~ is_covid, data = df_exp)

p_clone <- ggplot(df_exp, aes(x = is_covid, y = clone_size, fill = is_covid)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7, color = NA) +
  geom_boxplot(width = .15, outlier.size = 0.6, fatten = 1.2,
               fill = "white", color = "black") +
  scale_y_continuous(trans = "log1p") +
  scale_fill_brewer(palette = "Set2") +  # or: scale_fill_viridis_d()
  labs(x = NULL, y = "Clone size (log1p)",
       title = "Clonal expansion in COVID-like vs others") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank()
  )

sp(p_clone, "clone_size_covid_like_colored.png", 6, 4.5)


##=====V/J skew in COVID-linked cells (TRB)
## ===== V/J skew among COVID-like cells (TRB) =====
library(dplyr); library(ggplot2); library(readr)

# 1) Harmonize barcodes for joining: strip trailing "_1"/"_2" from Seurat barcodes
seumeta <- seu@meta.data %>%
  tibble::rownames_to_column("barcode_merged") %>%
  mutate(base_bc = sub("_\\d+$", "", barcode_merged))

# 2) Make a deduplicated V/J table from 'pair' (comes from CD4+CD8 concatenation)
pair_CD4 <- readr::read_tsv(file.path(tab_dir, "pair_CD4.tsv"), show_col_types = FALSE)
pair_CD8 <- readr::read_tsv(file.path(tab_dir, "pair_CD8.tsv"), show_col_types = FALSE)
pair <- dplyr::bind_rows(pair_CD4, pair_CD8)

pair_vj <- pair %>%
  transmute(base_bc = barcode, TRB_v, TRB_j, TRA_v, TRA_j) %>%
  distinct(base_bc, .keep_all = TRUE)

# 3) Join and define COVID group
vj_df <- seumeta %>%
  left_join(pair_vj, by = "base_bc") %>%
  mutate(group = ifelse(covid_like %in% c("low","medium"), "covid", "other"),
         lineage = celltype) %>%
  filter(!is.na(TRB_v), !is.na(TRB_j))

# 4) Enrichment helpers (per-feature 2x2 Fisher + OR + BH q)
enrich_feature <- function(vec, grp, fname){
  tab <- table(vec, grp)
  # keep only features observed in >= 3 cells overall (optional guard)
  keep <- rowSums(tab) >= 3
  tab  <- tab[keep, , drop = FALSE]
  if (nrow(tab) == 0) return(tibble::tibble())
  ress <- lapply(seq_len(nrow(tab)), function(i){
    r <- tab[i, ]
    # 2x2 against the rest
    mat <- matrix(c(r["covid"], sum(r) - r["covid"],
                    sum(tab[,"covid"]) - r["covid"],
                    sum(tab[,"other"]) - (sum(r) - r["covid"])), nrow=2, byrow=TRUE)
    ft <- suppressWarnings(fisher.test(mat))
    data.frame(
      feature = rownames(tab)[i],
      covid = as.integer(r["covid"]),
      other = as.integer(r["other"]),
      OR = unname(ft$estimate),
      p  = ft$p.value,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  ress$q <- p.adjust(ress$p, method = "BH")
  ress <- arrange(ress, q)
  write_tsv(ress, file.path(tab_dir, paste0(fname, "_enrichment.tsv")))
  ress
}

# 5) Run enrichment (all cells together); optionally also per lineage
enr_trbv <- enrich_feature(vj_df$TRB_v, vj_df$group, "TRBV")
enr_trbj <- enrich_feature(vj_df$TRB_j, vj_df$group, "TRBJ")

# 6) Plots: top enriched TRBV/TRBJ (by q, capped at 15)
plot_top <- function(df, feat, tit, n = 15){
  if (nrow(df) == 0) return(NULL)
  d <- df |>
    dplyr::slice_min(q, n = min(n, nrow(df))) |>
    dplyr::mutate(sig = factor(ifelse(q < 0.05, "FDR < 0.05", "Not significant"),
                               levels = c("FDR < 0.05", "Not significant")))
  
  ggplot(d, aes(x = reorder(feature, OR), y = OR, fill = sig)) +
    geom_col() +
    coord_flip() +
    geom_hline(yintercept = 1, linetype = 2) +
    scale_fill_manual(
      name = "Legend",
      values = c("FDR < 0.05" = "tomato", "Not significant" = "grey70")
    ) +
    labs(x = feat, y = "Odds ratio (covid vs other)", title = tit) +
    theme_bw(base_size = 11)
}


p_trbv <- plot_top(enr_trbv, "TRBV", "TRBV enrichment in COVID-like")
p_trbj <- plot_top(enr_trbj, "TRBJ", "TRBJ enrichment in COVID-like")
if (!is.null(p_trbv)) sp(p_trbv, "TRBV_enrichment_covid.png", 6.5, 5)
if (!is.null(p_trbj)) sp(p_trbj, "TRBJ_enrichment_covid.png", 6.5, 5)

# 7) (Optional) V–J combo skew: normalized heatmaps for covid vs other
# combo_tab <- vj_df %>%
#   count(group, TRB_v, TRB_j, name = "n") %>%
#   group_by(group) %>% mutate(frac = n / sum(n)) %>% ungroup()
# 
# combo_mat <- function(g){
#   m <- combo_tab %>% filter(group == g) %>%
#     select(TRB_v, TRB_j, frac) %>%
#     tidyr::pivot_wider(names_from = TRB_j, values_from = frac, values_fill = 0)
#   mat <- as.matrix(data.frame(m[,-1], row.names = m[[1]]))
#   mat
# }
# 
# mat_covid <- combo_mat("covid")
# mat_other <- combo_mat("other")
# 
# # Δ fraction heatmap (covid - other), clipped for display
# if (nrow(mat_covid) && ncol(mat_covid)) {
#   common_rows <- intersect(rownames(mat_covid), rownames(mat_other))
#   common_cols <- intersect(colnames(mat_covid), colnames(mat_other))
#   M <- mat_covid[common_rows, common_cols, drop=FALSE] - mat_other[common_rows, common_cols, drop=FALSE]
#   M[abs(M) < 1e-5] <- 0
#   # base R image plot (no ComplexHeatmap dependency)
#   png(file.path(fig_dir, "TRBVJ_delta_fraction_covid_vs_other.png"), width=1400, height=1000, res=160)
#   op <- par(mar=c(8,10,2,2))
#   image(t(M[nrow(M):1, , drop=FALSE]), axes=FALSE)
#   axis(1, at = seq(0,1,length.out=ncol(M)), labels = colnames(M), las=2, cex.axis=.7)
#   axis(2, at = seq(0,1,length.out=nrow(M)), labels = rev(rownames(M)), las=2, cex.axis=.7)
#   title("TRBV–TRBJ Δ fraction (covid - other)")
#   par(op)
#   dev.off()
# }

library(dplyr); library(tidyr)
library(ComplexHeatmap); library(circlize)

# Rebuild count tables to filter by abundance
ct <- vj_df %>% count(group, TRB_v, TRB_j, name="n")
wide_counts <- ct %>% group_by(TRB_v, TRB_j) %>% summarise(n_tot = sum(n), .groups="drop")
keep_pairs <- wide_counts %>% filter(n_tot >= 10) %>% select(TRB_v, TRB_j)

cov <- ct %>% filter(group=="covid") %>%
  right_join(keep_pairs, by=c("TRB_v","TRB_j")) %>%
  replace_na(list(n=0)) %>%
  group_by(TRB_v) %>% mutate(frac = n/sum(n)) %>% ungroup() %>% rename(frac = "frac_covid")

oth <- ct %>% filter(group=="other") %>%
  right_join(keep_pairs, by=c("TRB_v","TRB_j")) %>%
  replace_na(list(n=0)) %>%
  group_by(TRB_v) %>% mutate(frac = n/sum(n)) %>% ungroup()%>% rename(frac = "frac_other")

## 1) Build Δ (covid - other) *only* and pivot to V x J
Delta <- cov %>%
  inner_join(oth, by = c("TRB_v","TRB_j")) %>%
  transmute(TRB_v, TRB_j, delta = frac_covid - frac_other) %>%
  pivot_wider(names_from = TRB_j, values_from = delta, values_fill = 0)

mat <- as.matrix(Delta[,-1, drop = FALSE])
rownames(mat) <- Delta$TRB_v

## 2) If some TRBV rows repeat, collapse duplicates (mean). This also makes rownames unique.
if (any(duplicated(rownames(mat)))) {
  grp <- factor(rownames(mat))
  mat <- rowsum(mat, grp) / as.vector(table(grp))
}

## 3) OPTIONAL: row-wise z-score with guard for sd=0 (avoids NaN)
row_scale <- function(M){
  mu <- rowMeans(M, na.rm = TRUE)
  sdv <- apply(M, 1, sd, na.rm = TRUE)
  Z <- sweep(M, 1, mu, "-")
  Z <- sweep(Z, 1, ifelse(sdv > 0, sdv, 1), "/")   # divide by 1 if sd==0
  Z[!is.finite(Z)] <- 0                            # replace any Inf/NaN
  Z
}
mat_z <- row_scale(mat)   # or skip this line and use `mat` directly

## 4) Color range robust to NA/NaN
rng <- range(mat_z[is.finite(mat_z)], na.rm = TRUE)
rng
# e.g., use ComplexHeatmap
# col_fun <- circlize::colorRamp2(c(rng[1], 0, rng[2]), c("#2166AC","white","#B2182B"))
col_fun <- colorRamp2(c(min(mat_z), 0, max(mat_z)), c("#043915", "#FFF8E8", "#FF3F7F"))

ht <- Heatmap(
  mat_z,
  name = "Δ frac",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = TRUE, show_column_dend = TRUE,
  row_title = "TRBV", 
  column_title = "TRBJ",
  column_title_side = 'bottom'
)
png(file.path(fig_dir, "TRBVJ_delta_fraction_covid_vs_other.png"), width=1400, height=1100, res=160)
draw(ht)
dev.off()


## ---- 6) Canonical marker panel across top clonotypes (DotPlot) ----
# make sure cov_set exists (e.g., from your CONF_LEVEL choice)
cov_set <- rownames(seu)[seu$covid_like %in% c("low","medium")]
cov_cells_keep <- colnames(seu)[seu$covid_hit == "COVID-hit"]
assign_vec <- rep("OtherCOVID", ncol(seu))
names(assign_vec) <- colnames(seu)

for (i in seq_len(nrow(topN))) {
  bc <- colnames(seu)[seu$CTaa == topN$clonotype[i] & colnames(seu) %in% cov_cells_keep]
  if (length(bc)) assign_vec[bc] <- paste0("CT:", sprintf("%02d", i))
}

Idents(seu) <- factor(assign_vec,
                      levels = c(paste0("CT:", sprintf("%02d", seq_len(nrow(topN)))), "OtherCOVID"))

# 2) use an expression assay for color (integrated/SCT are centered)
DefaultAssay(seu) <- if ("RNA" %in% Assays(seu)) "RNA" else DefaultAssay(seu)

# 3) if needed, ensure RNA has log-normalized data
if (DefaultAssay(seu) == "RNA" && !"data" %in% Layers(seu[["RNA"]])) {
  seu <- NormalizeData(seu, verbose = FALSE)
}

# 4) DotPlot (same marker panel you already defined)
p_dot <- DotPlot(seu, features = panel, assay = DefaultAssay(seu), cols = c("lightgrey","firebrick")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle(sprintf("Canonical markers in top COVID clonotypes (%s)", CONF_LEVEL))
print(p_dot)

ggsave(filename = "./out_combined_tcr_covid/figs/covid_low_topCT_dotplot.pdf",plot = p_dot,
       width = 12,height = 12)

## ---- 7) UMAP: highlight COVID-hit & top epitopes ----
if ("umap" %in% names(seu@reductions)) {
  sp(
    DimPlot(seu, group.by = "covid_hit") +
      ggtitle(sprintf("COVID-hit (%s) vs Other", CONF_LEVEL)),
    sprintf("covid_%s_umap_hit.png", CONF_LEVEL), 6, 5
  )
  
  # label top 5 epitopes across COVID-hit cells
  cov_ep <- best_per %>% dplyr::select(barcode, epitope) %>% dplyr::distinct()
  seu$covid_epitope <- NA_character_
  idx <- match(cov_ep$barcode, colnames(seu))
  seu$covid_epitope[idx[!is.na(idx)]] <- cov_ep$epitope
  
  top_epi <- seu@meta.data %>%
    dplyr::filter(covid_hit == "COVID-hit", !is.na(covid_epitope)) %>%
    dplyr::count(covid_epitope, sort = TRUE) %>%
    dplyr::slice_head(n = 5) %>% dplyr::pull(covid_epitope)
  
  seu$epi_top5 <- ifelse(!is.na(seu$covid_epitope) & seu$covid_epitope %in% top_epi,
                         seu$covid_epitope, "other")
  
  sp(
    DimPlot(seu, group.by = "epi_top5", label = FALSE) +
      ggtitle("Top COVID epitopes (UMAP)"),
    sprintf("covid_%s_umap_top_epitopes.png", CONF_LEVEL), 6, 5
  )
}


# 8) Also save the joined table for reproducibility
write_tsv(vj_df %>% select(barcode_merged, base_bc, group, lineage, TRB_v, TRB_j),
          file.path(tab_dir, "vj_joined_table.tsv"))

# Attach labels: "CT:01  (n=14)  CAANTGGF..._CASSLGH..."
# make sure identities are set to CT:01..CT:20 and OtherCOVID first
Idents(seu) <- factor(assign_vec,
                      levels = c(paste0("CT:", sprintf("%02d", 1:nrow(topN))), "OtherCOVID"))

# lab_map: named vector like c("CT:01"="CT:01 (n=14)  CAANTGGF..._CASSLGH…", ...)
# (what you printed is perfect)

# rename identities safely
seu <- RenameIdents(seu, lab_map)
min_cells <- 5

keep_ids <- names(which(table(Idents(seu)) >= min_cells))
seu_sub   <- subset(seu, idents = keep_ids)

panel <- intersect(c(
  "IL2RA","HLA-DRA","PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","TNFRSF9","ICOS","CXCL13",
  "GZMK","GZMB","PRF1","NKG7","GNLY","CCL5","IFNG",
  "CCR7","SELL","IL7R","TCF7","LEF1","BCL2","CXCR5","BCL6","FOXP3","IKZF2"), rownames(seu_sub))

p_dot <- DotPlot(seu_sub, features = panel, cols = c("lightgrey","firebrick")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle(sprintf("Canonical markers in top COVID clonotypes (%s)", CONF_LEVEL))
sp(p_dot, sprintf("covid_%s_topCT_dotplot.png", CONF_LEVEL), 12, 6)


### spectratype ridges
suppressPackageStartupMessages(library(ggridges))

spec_df <- seu@meta.data |>
  tibble::rownames_to_column("barcode") |>
  mutate(
    group = paste0(celltype, "_", ifelse(covid_like %in% c("low","medium"), "COVID", "other")),
    len_trb = nchar(TRB_cdr3)
  ) |>
  filter(!is.na(len_trb))



p_ridge <- ggplot(spec_df, aes(x = len_trb, y = group, fill = group)) +
  geom_density_ridges(alpha = .6, scale = 1.2, rel_min_height = .01) +
  labs(x = "TRB CDR3 length (aa)", y = NULL, title = "Spectratype of TRB by lineage × COVID status") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
print(p_ridge)
sp(p_ridge, "spectratype_trb_lineage_covid.png", 8, 6)

### test
cells <- seu@meta.data |>
  tibble::rownames_to_column("barcode") |>
  dplyr::select(barcode, celltype, covid_hit, TRA_cdr3, TRB_cdr3)

spec <- rbind(
  cells |> filter(!is.na(TRA_cdr3)) |> transmute(chain="TRA", group=paste(celltype, ifelse(covid_hit=="COVID-hit","COVID","other"), sep="_"), len=nchar(TRA_cdr3)),
  cells |> filter(!is.na(TRB_cdr3)) |> transmute(chain="TRB", group=paste(celltype, ifelse(covid_hit=="COVID-hit","COVID","other"), sep="_"), len=nchar(TRB_cdr3))
)


## “Empirical cumulative distribution functions (ECDFs) of 
# TRB CDR3 amino-acid length for COVID-linked (‘COVID’) and non-linked (‘other’) clonotypes

library(effsize)

x8  <- spec$len[spec$group=="CD8_COVID"]
y8  <- spec$len[spec$group=="CD8_other"]
x4  <- spec$len[spec$group=="CD4_COVID"]
y4  <- spec$len[spec$group=="CD4_other"]

## 1) Permutation KS (handles ties/discrete)
perm_ks <- function(x, y, B=5000, seed=1){
  set.seed(seed)
  Dobs <- suppressWarnings(ks.test(x, y)$statistic)
  z <- c(x,y); n <- length(x)
  Dperm <- replicate(B, {
    i <- sample(seq_along(z), n)
    suppressWarnings(ks.test(z[i], z[-i])$statistic)
  })
  p <- mean(Dperm >= Dobs)
  c(D = unname(Dobs), p_perm = p)
}
perm_ks(x8,y8); perm_ks(x4,y4)

## 2) Location-based tests (Welch t / Wilcoxon) on length
t.test(x8, y8); wilcox.test(x8, y8, exact=FALSE)
t.test(x4, y4); wilcox.test(x4, y4, exact=FALSE)

## 3) Effect size (Cliff’s delta; robust for ordinal/discrete)
cliff.delta(x8, y8)
cliff.delta(x4, y4)

## 4) Simple summaries for context
summ <- \(a,b) data.frame(nA=length(a), nB=length(b),
                          medA=median(a), medB=median(b),
                          IQR_A=IQR(a), IQR_B=IQR(b))
summ(x8,y8); summ(x4,y4)

### plot
plot_ecdf <- function(x, y, labA, labB){
  dx <- data.frame(len=x, g=labA); dy <- data.frame(len=y, g=labB)
  ggplot(rbind(dx,dy), aes(len, colour=g)) +
    stat_ecdf(geom="step") +
    scale_x_continuous(breaks=seq(8,22,2)) +
    labs(x="TRB CDR3 length (aa)", y="ECDF", colour=NULL) +
    theme_bw()
}
p1 <-plot_ecdf(x8,y8,"CD8_COVID","CD8_other")
p2 <- plot_ecdf(x4,y4,"CD4_COVID","CD4_other")

#compute p-values
p_wilcox_cd8 <- wilcox.test(x8, y8, exact = FALSE)$p.value
p_wilcox_cd4 <- wilcox.test(x4, y4, exact = FALSE)$p.value

# pretty labels
lab8 <- sprintf("Wilcoxon p = %.2g", p_wilcox_cd8)
lab4 <- sprintf("Wilcoxon p = %.2g", p_wilcox_cd4)

# add as subtitle (cleanest)
p1a <- p1 + labs(title = "CD8 vs Other", subtitle = lab8)
p2a <- p2 + labs(title = "CD4 vs Other", subtitle = lab4)

p_stats <- p1a+p2a
sp(p_stats, "ecdf_trb_len_by_lineage_covid_status.png", 8, 6)

### Marker genes in each top clonotype vs other COVID-hit (fast visual)
# Average expression per identity (CTs + OtherCOVID) for the panel
avg <- AverageExpression(seu_sub, features = panel, assays = DefaultAssay(seu_sub), slot = "data")[[1]]
# z-score by gene
avg_z <- t(scale(t(avg)))

# visualize
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#2166AC","white","#B2182B"))
ht <- Heatmap(avg_z, name="z", col=col_fun, cluster_rows=TRUE, cluster_columns=TRUE,
              column_title="Identities (CTs vs OtherCOVID)")

print(ht)
png(file.path(fig_dir, "topCT_marker_heatmap.png"), 1400, 900, res=150); draw(ht); dev.off()


## Clonotype frequency vs epitope barplot (COVID only):
library(tidytext)
cov_only <- subset(seu, subset = covid_hit == "COVID-hit")
cov_only@meta.data %>%
  count(CTaa, covid_epitope, sort=TRUE) %>% slice_head(n=10) %>%
  ggplot(aes(reorder_within(CTaa, n, covid_epitope), n, fill=covid_epitope)) +
  geom_col(show.legend=FALSE) + coord_flip() +
  scale_x_reordered() + labs(x="Clonotype", y="# cells", title="Top CTs per epitope (COVID-hit)")


##=========== chordiagram
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

# returns a data.frame with columns: base_bc, TRB_v, TRB_j, lineage, covid_group
build_trb_table <- function(seu, tab_dir){
  # prefer your existing pair tables
  pair_cd4 <- readr::read_tsv(file.path(tab_dir, "pair_CD4.tsv"), show_col_types = FALSE)
  pair_cd8 <- readr::read_tsv(file.path(tab_dir, "pair_CD8.tsv"), show_col_types = FALSE)
  pair     <- bind_rows(pair_cd4, pair_cd8) %>%
    transmute(base_bc = barcode, TRB_v, TRB_j,
              lineage = ifelse(barcode %in% pair_cd4$barcode, "CD4", "CD8"))
  
  # join COVID status from Seurat meta (allowing for *_1/_2 suffixes)
  mv <- seu@meta.data %>%
    tibble::rownames_to_column("barcode_merged") %>%
    mutate(base_bc = sub("_\\d+$", "", barcode_merged)) %>%
    select(base_bc, covid_like)
  
  pair %>%
    left_join(mv, by = "base_bc") %>%
    mutate(covid_group = ifelse(covid_like %in% c("low","medium"), "covid", "other")) %>%
    filter(!is.na(TRB_v), !is.na(TRB_j))
}

suppressPackageStartupMessages(library(circlize))

# Make a chord diagram from TRB V-J usage
make_vj_chord <- function(df,
                          title        = "TRB V–J usage chord diagram",
                          file_out     = "vj_chord_TRB.png",
                          min_pair_count = 10,
                          normalize_rows = FALSE,
                          from_contigs = FALSE,
                          V_col = NULL, J_col = NULL,         # only used when from_contigs=TRUE
                          seed = 1) {
  
  set.seed(seed)
  
  if (from_contigs) {
    stopifnot(!is.null(V_col), !is.null(J_col))
    tab <- df %>%
      transmute(TRB_v = .data[[V_col]], TRB_j = .data[[J_col]]) %>%
      filter(!is.na(TRB_v), !is.na(TRB_j)) %>%
      count(TRB_v, TRB_j, name = "n")
  } else {
    tab <- df %>% count(TRB_v, TRB_j, name = "n")
  }
  
  # prune to keep the plot readable
  tab <- tab %>% filter(n >= min_pair_count)
  if (!nrow(tab)) {
    message("[chord] No pairs above min_pair_count = ", min_pair_count, " → skipping: ", file_out)
    return(invisible(NULL))
  }
  
  # wide matrix (rows = V, cols = J)
  mat <- tab %>%
    tidyr::pivot_wider(names_from = TRB_j, values_from = n, values_fill = 0) %>%
    as.data.frame()
  rownames(mat) <- mat$TRB_v
  mat$TRB_v <- NULL
  m <- as.matrix(mat)
  
  if (normalize_rows) {
    rs <- rowSums(m); rs[rs == 0] <- 1
    m <- sweep(m, 1, rs, "/")
  }
  
  ## ---- Make sector names unique & order them nicely
  V_names <- rownames(m)
  J_names <- colnames(m)
  V_lab   <- paste0("V:", V_names)
  J_lab   <- paste0("J:", J_names)
  dimnames(m) <- list(V_lab, J_lab)
  
  # order by total weights for a cleaner layout
  rord <- order(rowSums(m), decreasing = TRUE)
  cord <- order(colSums(m), decreasing = TRUE)
  m <- m[rord, cord, drop = FALSE]
  V_lab <- rownames(m); J_lab <- colnames(m)
  
  ## ---- Palettes
  # V side: many distinct hues (wraps automatically)
  palV <- hue_pal()(max(8, length(V_lab)))
  names(palV) <- V_lab
  # J side: another palette (salmon/orange tone family)
  palJ <- hue_pal(h = c(10, 30))(max(6, length(J_lab)))
  names(palJ) <- J_lab
  grid_col <- c(palV, palJ)
  
  ## ---- Optional: color each link by its V gene to help track V-specific preferences
  # Build long edge list from m
  edges <- as.data.frame(as.table(m), stringsAsFactors = FALSE)
  colnames(edges) <- c("from", "to", "value")
  edges <- edges %>% filter(value > 0)
  
  # link color = color of the V sector (from)
  link_col <- palV[match(edges$from, names(palV))]
  # use some transparency so overlaps are visible
  link_col <- alpha(link_col, 0.6)
  
  ## ---- nice spacing: larger gap between V block and J block
  gap_after <- c(rep(2, length(V_lab) - 1), 10, rep(2, length(J_lab) - 1), 10)
  
  png(file_out, width = 1400, height = 1200, res = 180)
  circos.clear()
  circos.par(gap.after = gap_after, start.degree = 90, track.margin = c(0, 0.01))
  par(mar = c(1, 1, 3, 1))
  
  # We plot from edge list so we can pass per-link colors
  chordDiagram(
    x = edges,
    grid.col = grid_col,
    col = link_col,
    transparency = 0.15,           # lower transparency = more vivid
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.08)
  )
  
  # Label sectors with small text that follows the circle
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(
      get.cell.meta.data("xcenter"),
      get.cell.meta.data("ylim")[1],
      sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.6
    )
  }, bg.border = NA)
  
  title(title, cex.main = 1.25)
  dev.off()
  
  message("[chord] saved: ", file_out)
  invisible(m)
}


## run
# build the base table
trb_df <- build_trb_table(seu, tab_dir)

# ---- Overall (all lineages combined)
make_vj_chord(
  df       = trb_df,
  title    = "TRB V–J usage chord diagram (all T cells)",
  file_out = file.path(fig_dir, "vj_chord_TRB_all.png"),
  min_pair_count = 15
)

# ---- By lineage
make_vj_chord(
  df       = filter(trb_df, lineage == "CD8"),
  title    = "TRB V–J usage chord diagram (CD8)",
  file_out = file.path(fig_dir, "vj_chord_TRB_CD8.png"),
  min_pair_count = 10
)
make_vj_chord(
  df       = filter(trb_df, lineage == "CD4"),
  title    = "TRB V–J usage chord diagram (CD4)",
  file_out = file.path(fig_dir, "vj_chord_TRB_CD4.png"),
  min_pair_count = 10
)

# ---- COVID vs other (CD8 only, to match your focus)
make_vj_chord(
  df       = filter(trb_df, lineage == "CD8", covid_group == "covid"),
  title    = "TRB V–J usage chord diagram (CD8, COVID-like)",
  file_out = file.path(fig_dir, "vj_chord_TRB_CD8_covid.png"),
  min_pair_count = 6
)
make_vj_chord(
  df       = filter(trb_df, lineage == "CD8", covid_group == "other"),
  title    = "TRB V–J usage chord diagram (CD8, other)",
  file_out = file.path(fig_dir, "vj_chord_TRB_CD8_other.png"),
  min_pair_count = 6
)
