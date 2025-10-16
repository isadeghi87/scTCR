## ========== 04_scirpy_like_analyses.R ==========
source("/home/isadeghi/projects/covid_kids/scripts/00_utils.R")

suppressPackageStartupMessages({
  library(igraph)
  library(ggraph)
  library(Biostrings)
  library(ggseqlogo)
  library(ggpubr)
  library(reshape2)
  library(ggalluvial)
  library(vegan)      # Morisita-Horn / Jaccard
})

seu <- readRDS(file.path(out_dir, "merged_seurat_final.rds"))

## Helper: safe NA-to-empty
na0 <- function(x) ifelse(is.na(x), "", x)

## Build per-cell table we’ll re-use
cells <- seu@meta.data |>
  tibble::rownames_to_column("barcode") |>
  dplyr::transmute(
    barcode, celltype, cluster = seurat_clusters,
    TRA = toupper(TRA_cdr3), TRB = toupper(TRB_cdr3),
    CTaa = na0(CTaa), clone_id = na0(clone_id), clone_size = as.integer(clone_size)
  )

## --------------------------------------------------------------------------
## A) Clonal expansion categories (like scirpy.tl.clonal_expansion)
## --------------------------------------------------------------------------
cut_expansion <- function(n) {
  if (is.na(n) || n <= 1) return("Singleton")
  if (n <= 5) return("Small (2-5)")
  if (n <= 20) return("Medium (6-20)")
  return("Large (21+)")
}
cells$expansion <- vapply(cells$clone_size, cut_expansion, character(1))

pA1 <- cells |>
  dplyr::count(celltype, expansion) |>
  dplyr::group_by(celltype) |>
  dplyr::mutate(pct = 100 * n / sum(n)) |>
  ggplot(aes(expansion, pct, fill = celltype)) +
  geom_col(position = "dodge") +
  coord_flip() + ylab("% of cells") + xlab("Clonal expansion") +
  ggtitle("Clonal expansion by lineage (CD4 vs CD8)") +
  theme_bw(base_size = 11)
sp(pA1, "scirpy_like_expansion_by_lineage.png", 7, 4)

pA2 <- cells |>
  dplyr::count(cluster = factor(cluster), expansion) |>
  dplyr::group_by(cluster) |>
  dplyr::mutate(pct = 100 * n / sum(n)) |>
  ggplot(aes(expansion, pct, fill = expansion)) +
  geom_col() + coord_flip() +
  facet_wrap(~ cluster, ncol = 4, scales = "free_y") +
  ylab("% of cells") + xlab("Clonal expansion") +
  ggtitle("Clonal expansion by Seurat cluster") +
  theme_bw(base_size = 10) + theme(legend.position = "none")
sp(pA2, "scirpy_like_expansion_by_cluster.png", 9, 6)

## --------------------------------------------------------------------------
## B) Repertoire overlap metrics (like scirpy.pl.repertoire_overlap)
##    Morisita-Horn & Jaccard on clonotype presence across groups
## --------------------------------------------------------------------------
# presence/absence of CTaa per group
ct_by_group <- cells |>
  dplyr::filter(CTaa != "") |>
  dplyr::count(group = celltype, CTaa) |>
  tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0)
mat_ct <- as.matrix(ct_by_group[,-1]); rownames(mat_ct) <- ct_by_group$CTaa

if (ncol(mat_ct) >= 2) {
  # Morisita-Horn distance -> similarity
  mh <- as.matrix(vegdist(t(mat_ct), method = "horn"))
  mh_sim <- 1 - mh
  png(file.path(fig_dir, "scirpy_like_overlap_morisita_horn.png"), 800, 700, res = 140)
  Heatmap(mh_sim, name = "Morisita-Horn\nsimilarity",
          cluster_rows = FALSE, cluster_columns = FALSE)
  dev.off()
  
  # Jaccard on presence/absence
  bin_ct <- (mat_ct > 0) * 1
  jac <- as.matrix(vegdist(t(bin_ct), method = "jaccard"))
  jac_sim <- 1 - jac
  png(file.path(fig_dir, "scirpy_like_overlap_jaccard.png"), 800, 700, res = 140)
  Heatmap(jac_sim, name = "Jaccard\nsimilarity",
          cluster_rows = FALSE, cluster_columns = FALSE)
  dev.off()
}

## --------------------------------------------------------------------------
## C) Spectratype (CDR3 length distributions) per chain & lineage
## --------------------------------------------------------------------------
len_df <- rbind(
  cells |> dplyr::filter(!is.na(TRA)) |> dplyr::transmute(lineage = celltype, chain = "TRA", len = nchar(TRA)),
  cells |> dplyr::filter(!is.na(TRB)) |> dplyr::transmute(lineage = celltype, chain = "TRB", len = nchar(TRB))
)
pC <- ggplot(len_df, aes(x = len, y = ..density.., fill = lineage)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.5) +
  facet_wrap(~ chain, nrow = 1) +
  xlab("CDR3 length (aa)") + ylab("Density") +
  ggtitle("CDR3 spectratype by chain & lineage") +
  theme_bw(base_size = 11)
sp(pC, "scirpy_like_spectratype.png", 8, 4)

## --------------------------------------------------------------------------
## D) Clonotype similarity network (like scirpy.pl.tcr_network)
##    Nodes = clonotypes (TRB), edges if Levenshtein dist ≤ 1 (same length)
## --------------------------------------------------------------------------
library(stringdist)

# Build clonotype table using TRB as key (most public datasets are TRB-centric)
clones_trb <- cells |>
  dplyr::filter(TRB != "") |>
  dplyr::group_by(TRB) |>
  dplyr::summarise(
    size = n(),
    n_cd4 = sum(celltype == "CD4"),
    n_cd8 = sum(celltype == "CD8"),
    clusters = paste(sort(unique(cluster)), collapse = ","),
    .groups = "drop"
  )

trb_vec <- clones_trb$TRB
if (length(trb_vec) > 1) {
  # Compute pairwise distances efficiently by length bucket
  edges <- list()
  by_len <- split(seq_along(trb_vec), nchar(trb_vec))
  for (ix in by_len) {
    if (length(ix) < 2) next
    seqs <- trb_vec[ix]
    D <- stringdistmatrix(seqs, seqs, method = "lv")
    diag(D) <- 99
    pairs <- which(D <= 1, arr.ind = TRUE)
    if (nrow(pairs)) {
      edges[[length(edges) + 1]] <- data.frame(
        from = clones_trb$TRB[ix[pairs[,1]]],
        to   = clones_trb$TRB[ix[pairs[,2]]],
        dist = D[pairs],
        stringsAsFactors = FALSE
      )
    }
  }
  edges <- if (length(edges)) dplyr::bind_rows(edges) else data.frame()
  edges <- unique(edges[edges$from < edges$to, , drop = FALSE])
  
  g <- graph_from_data_frame(d = edges,
                             vertices = clones_trb |> dplyr::mutate(
                               frac_cd4 = n_cd4 / pmax(size, 1),
                               frac_cd8 = n_cd8 / pmax(size, 1)
                             ),
                             directed = FALSE)
  
  V(g)$size_viz <- scales::rescale(V(g)$size, to = c(3, 16))
  V(g)$color <- ifelse(V(g)$n_cd8 > V(g)$n_cd4, "#1f77b4", "#2ca02c")  # CD8-ish vs CD4-ish
  
  pD <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(alpha = 1 - dist), show.legend = FALSE) +
    geom_node_point(aes(size = size_viz, color = color)) +
    scale_color_identity() +
    ggtitle("TRB clonotype similarity network (Levenshtein ≤ 1)") +
    theme_void(base_size = 11)
  sp(pD, "scirpy_like_trb_network.png", 8, 6)
  
  # Export components (potential “public-like” neighborhoods)
  comps <- components(g)
  comp_tab <- data.frame(TRB = names(comps$membership),
                         comp_id = comps$membership,
                         stringsAsFactors = FALSE) |>
    dplyr::left_join(clones_trb, by = c("TRB"))
  readr::write_tsv(comp_tab, file.path(tab_dir, "trb_network_components.tsv"))
}

## --------------------------------------------------------------------------
## E) Public / shared clonotypes across lineages (like scirpy.tl.get_shared_clonotypes)
## --------------------------------------------------------------------------
shared_trb <- cells |>
  dplyr::filter(TRB != "") |>
  dplyr::distinct(TRB, celltype) |>
  dplyr::count(TRB) |>
  dplyr::filter(n > 1) |>
  dplyr::left_join(
    cells |>
      dplyr::group_by(TRB) |>
      dplyr::summarise(n_cells = n(), cd4 = sum(celltype == "CD4"), cd8 = sum(celltype == "CD8"), .groups = "drop"),
    by = "TRB"
  ) |>
  dplyr::arrange(desc(n_cells))
readr::write_tsv(shared_trb, file.path(tab_dir, "public_trb_cd4_cd8.tsv"))

## --------------------------------------------------------------------------
## F) TRA–TRB pairing alluvial (like scirpy.pl.umap with pair categories)
## --------------------------------------------------------------------------
pair_cat <- cells |>
  dplyr::mutate(pairing = dplyr::case_when(
    TRA != "" & TRB != "" ~ "paired",
    TRA != "" & TRB == "" ~ "alpha-only",
    TRA == "" & TRB != "" ~ "beta-only",
    TRUE ~ "none"
  ))
pF <- ggplot(pair_cat,
             aes(y = ..count.., axis1 = celltype, axis2 = pairing)) +
  geom_alluvium(aes(fill = pairing), alpha = 0.8) +
  geom_stratum(width = 1/8) + geom_text(stat = "stratum", infer.label = TRUE, size = 3) +
  scale_x_discrete(limits = c("Lineage", "Pairing"), expand = c(.1, .1)) +
  ylab("# cells") + ggtitle("TRA–TRB pairing by lineage") +
  theme_bw(base_size = 11)
sp(pF, "scirpy_like_pairing_alluvial.png", 7, 4)

## --------------------------------------------------------------------------
## G) Top clones localization across clusters (like scirpy.pl.clonal_expansion with pies)
## --------------------------------------------------------------------------
topN <- cells |>
  dplyr::filter(CTaa != "") |>
  dplyr::count(CTaa, name = "n_cells") |>
  dplyr::arrange(desc(n_cells)) |>
  dplyr::slice_head(n = 20)

loc <- cells |>
  dplyr::semi_join(topN, by = "CTaa") |>
  dplyr::count(CTaa, cluster) |>
  dplyr::group_by(CTaa) |>
  dplyr::mutate(pct = 100 * n / sum(n))

pG <- loc |>
  ggplot(aes(x = reorder(CTaa, -pct), y = pct, fill = cluster)) +
  geom_col() + coord_flip() +
  labs(x = "Top clonotypes", y = "% cells per cluster",
       title = "Top-20 clonotypes: cluster composition") +
  theme_bw(base_size = 9)
sp(pG, "scirpy_like_top_clones_cluster_mix.png", 8, 7)

## --------------------------------------------------------------------------
## H) Motif logo for dominant TRB CDR3s (like scirpy’s motif views)
## --------------------------------------------------------------------------
# Take TRB of top expanded clones and draw a sequence logo (aligned by length mode)
top_trb <- cells |>
  dplyr::filter(TRB != "") |>
  dplyr::count(TRB, name = "n") |>
  dplyr::arrange(desc(n)) |>
  dplyr::slice_head(n = 200)

# use the most common length to avoid gapping in logo
mode_len <- as.integer(names(sort(table(nchar(top_trb$TRB)), decreasing = TRUE))[1])
seqs <- top_trb$TRB[nchar(top_trb$TRB) == mode_len]
if (length(seqs) >= 10) {
  pH <- ggseqlogo(seqs, method = "prob") + ggtitle(sprintf("TRB CDR3 motif (top, len=%d)", mode_len))
  sp(pH, "scirpy_like_trb_logo.png", 7, 3)
}

## --------------------------------------------------------------------------
## I) UMAP “clonal neighborhood” (fraction of neighbors sharing clonotype)
## --------------------------------------------------------------------------
# Build a neighbor graph from the integrated reduction and compute per-cell fraction
if ("integrated" %in% Assays(seu) && "umap" %in% names(seu@reductions)) {
  emb <- Embeddings(seu, "umap")
  # kNN in UMAP space
  k <- 15
  nn_idx <- FNN::get.knn(emb, k = k)$nn.index  # needs FNN
  frac_same <- numeric(nrow(emb))
  ctvec <- seu$CTaa
  for (i in seq_len(nrow(emb))) {
    nns <- nn_idx[i,]
    frac_same[i] <- mean(ctvec[nns] == ctvec[i], na.rm = TRUE)
  }
  seu$clonal_neighborhood <- frac_same
  pI <- FeaturePlot(seu, features = "clonal_neighborhood", reduction = "umap") +
    ggtitle("Clonal neighborhood (UMAP kNN fraction sharing CTaa)")
  sp(pI, "scirpy_like_clonal_neighborhood.png", 6, 5)
}

message("Done: see figs/*scirpy_like_*.png and tables/* for new outputs.")
