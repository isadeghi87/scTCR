## ========== 05e_tcr_beta_similarity_network.R ==========
## Purpose: build TRB CDR3 similarity graph (Levenshtein ≤ 1) among COVID-like cells
cov_betas <- seu@meta.data |> tibble::rownames_to_column("barcode") |>
  dplyr::filter(covid_like %in% c("low","medium"), !is.na(TRB_cdr3)) |>
  dplyr::distinct(TRB_cdr3, .keep_all=TRUE) |> dplyr::pull(TRB_cdr3)

edges <- lapply(split(seq_along(cov_betas), nchar(cov_betas)), function(ix){
  if (length(ix)<2) return(NULL)
  D <- stringdist::stringdistmatrix(cov_betas[ix], cov_betas[ix], method="lv"); diag(D) <- 99
  prs <- which(D<=1, arr.ind=TRUE)
  if (!nrow(prs)) return(NULL)
  tibble::tibble(from=cov_betas[ix[prs[,1]]], to=cov_betas[ix[prs[,2]]], dist=D[prs])
}) |> dplyr::bind_rows()

if (!is.null(edges) && nrow(edges)>0) {
  g <- igraph::graph_from_data_frame(edges, directed=FALSE)
  comps <- igraph::components(g)$membership
  net_tab <- data.frame(TRB_cdr3 = names(comps), community = comps)
  readr::write_tsv(net_tab, file.path(tab_dir, "tcr_beta_lev1_communities.tsv"))
}
