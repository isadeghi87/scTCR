## ========== 02_tcr_integration_diversity.R ==========
source("/home/isadeghi/projects/covid_kids/scripts/00_utils.R")
seu <- readRDS(file.path(out_dir, "merged_seurat.rds"))
ct  <- readRDS(file.path(out_dir, "contigs_raw.rds"))

## Build clono lists for scRepertoire
to_scR <- function(df, sample){
  d <- df; names(d) <- make.names(names(d))
  out <- list(d); names(out) <- sample
  scRepertoire::combineTCR(out, samples = sample, ID = paste0(sample,"_patient"))
}
cl4 <- to_scR(c4, "CD4")
cl8 <- to_scR(c8, "CD8")
clono <- c(cl4, cl8)

## Strip any prefixes in combineTCR barcodes (keep after last underscore)
for(i in seq_along(clono)){
  if ("barcode" %in% names(clono[[i]])){
    clono[[i]]$barcode <- sub("^.*_(?=[^_]+$)", "", clono[[i]]$barcode, perl=TRUE)
  }
}

## Pick clonotype column
clone_col <- if ("CTaa" %in% names(clono[[1]])) "CTaa" else if ("CTnt" %in% names(clono[[1]])) "CTnt" else "CTgene"

## Combine into Seurat
## 1) Collect the suffixes used in the merged object (e.g. "_1", "_2")
suffixes <- unique(sub(".*(-\\d+)(_\\d+)$", "\\2", colnames(seu)))
stopifnot(all(grepl("^_\\d+$", suffixes)))  # sanity

## 2) For each clono list element, choose the suffix that maximizes matches
pick_suffix <- function(bcs, all_cells, cand_suffixes) {
  # if they already match as-is (e.g., you re-ran merge with add.cell.ids), keep as-is
  if (mean(bcs %in% all_cells) > 0.5) return("")  # already good enough
  scores <- sapply(cand_suffixes, function(sfx) mean(paste0(bcs, sfx) %in% all_cells))
  cand_suffixes[ which.max(scores) ]
}

for (nm in names(clono)) {
  if (!"barcode" %in% names(clono[[nm]])) next
  sfx <- pick_suffix(clono[[nm]]$barcode, colnames(seu), suffixes)
  if (nzchar(sfx)) clono[[nm]]$barcode <- paste0(clono[[nm]]$barcode, sfx)
  cat(sprintf("Mapped %s with suffix '%s' (match rate: %.1f%%)\n",
              nm, sfx,
              100*mean(clono[[nm]]$barcode %in% colnames(seu))))
}

## 3) Re-run combineExpression (handle groupBy vs group.by depending on version)
seu <- combineExpression( clono,                      # contig list
                          seu,                        # Seurat object
                          cloneCall = "aa",
                          chain = "both",
                          group.by = 'sample',
                          filterNA  = TRUE)


## 4) Quick confirmations
print(intersect(c("CTaa","CTnt","CTgene","cloneType"), colnames(seu@meta.data)))
cat("Mapped (non-NA CTaa): ", sum(!is.na(seu$CTaa)), " / ", ncol(seu), "\n", sep = "")

## Diversity (downsample by celltype for fairness)
if ("clonalDiversity" %in% ls(getNamespace("scRepertoire"))) {
  div_ct <- scRepertoire::clonalDiversity(seu, cloneCall=clone_col, group="celltype", n.boots=200)
  # save figure or table depending on version
  if (inherits(div_ct, "gg")) {
    sp(div_ct + ggtitle(paste0("Clonal diversity by celltype (", clone_col, ")")), "diversity_by_celltype.png", 6, 5)
  } else {
    readr::write_tsv(as.data.frame(div_ct), file.path(tab_dir,"diversity_by_celltype.tsv"))
  }
}

# Homeostasis (expansion categories) on merged object
p_homeo <- scRepertoire::clonalHomeostasis(seu, cloneCall = "aa")
sp(p_homeo, "clonal_homeostasis_merged.png")

# TRB V/J usage (whole sample)
## ---- Manual TRB V/J usage & V–J heatmap (version-agnostic) ----
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(ComplexHeatmap)})

contigs_all <- dplyr::bind_rows(CD4 = ct$CD4, CD8 = ct$CD8, .id = "lineage")

trb <- contigs_all %>%
  dplyr::filter(productive %in% c(TRUE, "True")) %>%
  dplyr::mutate(locus = ifelse(!is.na(chain), chain, locus)) %>%
  dplyr::filter(grepl("TRB|TCRB|Beta", locus, ignore.case = TRUE)) %>%
  dplyr::transmute(lineage, V = v_gene, J = j_gene) %>%
  dplyr::filter(!is.na(V), !is.na(J))

## 1) TRB V usage per lineage
v_counts <- trb %>% count(lineage, V, name = "n") %>% group_by(lineage) %>%
  mutate(pct = 100*n/sum(n)) %>% ungroup()

p_v <- ggplot(v_counts, aes(x = reorder(V, pct), y = pct, fill = lineage)) +
  geom_col(position = "dodge") + coord_flip() +
  labs(x = "TRB V gene", y = "% of TRB reads", title = "TRB V usage by lineage") +
  theme_bw(base_size = 11)
sp(p_v, "vj_usage_TRB_V_manual.png", w = 8, h = 6)

## 2) TRB J usage per lineage
j_counts <- trb %>% count(lineage, J, name = "n") %>% group_by(lineage) %>%
  mutate(pct = 100*n/sum(n)) %>% ungroup()

p_j <- ggplot(j_counts, aes(x = reorder(J, pct), y = pct, fill = lineage)) +
  geom_col(position = "dodge") + coord_flip() +
  labs(x = "TRB J gene", y = "% of TRB reads", title = "TRB J usage by lineage") +
  theme_bw(base_size = 11)
sp(p_j, "vj_usage_TRB_J_manual.png", w = 8, h = 6)

## 3) TRB V–J heatmap (pooled or per lineage)
vj_mat <- trb %>%
  count(V, J, name = "n") %>%
  tidyr::pivot_wider(names_from = J, values_from = n, values_fill = 0) %>%
  as.data.frame()
row.names(vj_mat) <- vj_mat$V; vj_mat$V <- NULL
m <- as.matrix(vj_mat)

png(file.path(fig_dir, "vj_heatmap_TRB_manual.png"), width = 1200, height = 900, res = 150)
Heatmap(m, name = "count", column_title = "TRB V–J usage (all lineages)",
        cluster_rows = TRUE, cluster_columns = TRUE, show_row_dend = FALSE, show_column_dend = FALSE)
dev.off()

## (Optional) split heatmaps by lineage and save both
for (lin in unique(trb$lineage)) {
  vj_l <- trb %>% filter(lineage == lin) %>%
    count(V, J, name = "n") %>%
    tidyr::pivot_wider(names_from = J, values_from = n, values_fill = 0) %>%
    as.data.frame()
  if (nrow(vj_l) > 0) {
    row.names(vj_l) <- vj_l$V; vj_l$V <- NULL
    png(file.path(fig_dir, paste0("vj_heatmap_TRB_", lin, "_manual.png")), width = 1200, height = 900, res = 150)
    Heatmap(as.matrix(vj_l), name = "count", column_title = paste0("TRB V–J (", lin, ")"),
            cluster_rows = TRUE, cluster_columns = TRUE, show_row_dend = FALSE, show_column_dend = FALSE)
    dev.off()
  }
}

# Save updated Seurat
saveRDS(seu, file.path(out_dir, "merged_seurat_with_tcr.rds"))
