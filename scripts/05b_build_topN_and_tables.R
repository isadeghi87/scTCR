## ========== 05b_build_topN_and_tables.R ==========
## Purpose: build the 'topN' clonotypes among COVID-hit cells and write quick enrichment tables.

stopifnot("CTaa" %in% colnames(seu@meta.data))

# pick clonotype label (prefer CTaa)
clono_col <- dplyr::case_when(
  "CTaa" %in% colnames(seu@meta.data) ~ "CTaa",
  "clone_id" %in% colnames(seu@meta.data) ~ "clone_id",
  TRUE ~ "TRB_cdr3"
)

cov_meta <- seu@meta.data |>
  tibble::rownames_to_column("barcode") |>
  dplyr::filter(barcode %in% colnames(seu)[seu$covid_hit=="COVID-hit"],
                !is.na(.data[[clono_col]]), .data[[clono_col]]!="")

topN <- if (nrow(cov_meta)) {
  cov_meta |>
    dplyr::count(.data[[clono_col]], name="n_cells") |>
    dplyr::arrange(dplyr::desc(n_cells)) |>
    dplyr::slice_head(n=20) |>
    dplyr::rename(clonotype = 1)
} else tibble::tibble(clonotype=character(), n_cells=integer())

readr::write_tsv(topN, file.path(tab_dir, "covid_top20_clonotypes.tsv"))
message("topN saved: ", nrow(topN), " clonotypes")

# Optional enrichment tables if covid_like present
if ("covid_like" %in% colnames(seu@meta.data)) {
  enrich <- function(fac){
    tab <- table(seu$covid_like != "none", fac)
    data.frame(
      level = colnames(tab),
      p = apply(tab, 2, \(col) fisher.test(matrix(c(col[2], sum(col)-col[2], col[1], sum(col)-col[1]), 2))$p.value)
    ) |> dplyr::mutate(q = p.adjust(p, "BH")) |> dplyr::arrange(q)
  }
  readr::write_tsv(enrich(seu$celltype),        file.path(tab_dir, "enrichment_lineage_covid_like.tsv"))
  readr::write_tsv(enrich(seu$seurat_clusters), file.path(tab_dir, "enrichment_cluster_covid_like.tsv"))
}
