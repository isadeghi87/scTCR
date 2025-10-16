## ========== 05i_vj_skew_enrichment.R ==========
## Purpose: TRBV/TRBJ enrichment among COVID-like (Fisher-OR) + barplots.
pair_CD4 <- readr::read_tsv(file.path(tab_dir, "pair_CD4.tsv"), show_col_types=FALSE)
pair_CD8 <- readr::read_tsv(file.path(tab_dir, "pair_CD8.tsv"), show_col_types=FALSE)
pair <- dplyr::bind_rows(pair_CD4, pair_CD8)

seumeta <- seu@meta.data |> tibble::rownames_to_column("barcode_merged") |> dplyr::mutate(base_bc=sub("_\\d+$","", barcode_merged))
pair_vj <- pair |> dplyr::transmute(base_bc = barcode, TRB_v, TRB_j, TRA_v, TRA_j) |> dplyr::distinct(base_bc, .keep_all=TRUE)

vj_df <- seumeta |>
  dplyr::left_join(pair_vj, by="base_bc") |>
  dplyr::mutate(group = ifelse(covid_like %in% c("low","medium"), "covid", "other"), lineage=celltype) |>
  dplyr::filter(!is.na(TRB_v), !is.na(TRB_j))

enrich_feature <- function(vec, grp, fname){
  tab <- table(vec, grp); keep <- rowSums(tab) >= 3; tab <- tab[keep,, drop=FALSE]
  if (!nrow(tab)) return(tibble::tibble())
  ress <- lapply(seq_len(nrow(tab)), function(i){
    r <- tab[i,]
    mat <- matrix(c(r["covid"], sum(r)-r["covid"], sum(tab[,"covid"])-r["covid"], sum(tab[,"other"])-(sum(r)-r["covid"])), 2, byrow=TRUE)
    ft <- suppressWarnings(fisher.test(mat))
    data.frame(feature=rownames(tab)[i], covid=as.integer(r["covid"]), other=as.integer(r["other"]), OR=unname(ft$estimate), p=ft$p.value)
  }) |> dplyr::bind_rows() |> dplyr::mutate(q=p.adjust(p,"BH")) |> dplyr::arrange(q)
  readr::write_tsv(ress, file.path(tab_dir, paste0(fname, "_enrichment.tsv"))); ress
}
enr_trbv <- enrich_feature(vj_df$TRB_v, vj_df$group, "TRBV")
enr_trbj <- enrich_feature(vj_df$TRB_j, vj_df$group, "TRBJ")

plot_top <- function(df, feat, tit, n=15){
  if (!nrow(df)) return(NULL)
  d <- df |> dplyr::slice_min(q, n=min(n, nrow(df))) |>
    dplyr::mutate(sig = factor(ifelse(q < 0.05, "FDR < 0.05", "Not significant"), levels=c("FDR < 0.05","Not significant")))
  ggplot2::ggplot(d, ggplot2::aes(x=reorder(feature, OR), y=OR, fill=sig)) +
    ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::geom_hline(yintercept=1, linetype=2) +
    ggplot2::scale_fill_manual(name="Legend", values=c("FDR < 0.05"="tomato","Not significant"="grey70")) +
    ggplot2::labs(x=feat, y="Odds ratio (covid vs other)", title=tit) + ggplot2::theme_bw(11)
}
p_trbv <- plot_top(enr_trbv, "TRBV", "TRBV enrichment in COVID-like")
p_trbj <- plot_top(enr_trbj, "TRBJ", "TRBJ enrichment in COVID-like")
if (!is.null(p_trbv)) sp(p_trbv, "TRBV_enrichment_covid.png", 6.5, 5)
if (!is.null(p_trbj)) sp(p_trbj, "TRBJ_enrichment_covid.png", 6.5, 5)
