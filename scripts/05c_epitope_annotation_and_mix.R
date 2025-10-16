## ========== 05c_epitope_annotation_and_mix.R ==========
## Purpose: assign best epitope per barcode and plot lineage composition among COVID-hit.

hits <- readr::read_tsv(file.path(tab_dir, "covid_hits_CD4_CD8.tsv"), show_col_types = FALSE)
prio <- dplyr::mutate(hits, priority = dplyr::case_when(
  match_chain=="TRB" & match_type=="exact" ~ 3L,
  match_chain=="TRA" & match_type=="exact" ~ 2L,
  TRUE ~ 0L
))
best_epi <- prio |> dplyr::group_by(barcode) |> dplyr::slice_max(priority, with_ties = FALSE) |>
  dplyr::ungroup() |> dplyr::select(barcode, epitope)

seu$covid_epitope <- NA_character_
m <- match(best_epi$barcode, colnames(seu)); seu$covid_epitope[m[!is.na(m)]] <- best_epi$epitope

epi_mix <- seu@meta.data |>
  tibble::rownames_to_column("barcode") |>
  dplyr::filter(covid_hit=="COVID-hit", !is.na(covid_epitope)) |>
  dplyr::count(celltype, covid_epitope, name="n") |>
  dplyr::group_by(celltype) |>
  dplyr::mutate(pct = 100*n/sum(n)) |> dplyr::ungroup()

readr::write_tsv(epi_mix, file.path(tab_dir, "covid_epitope_mix_by_lineage.tsv"))

p_epi <- ggplot2::ggplot(epi_mix, ggplot2::aes(x=reorder(covid_epitope, pct), y=pct, fill=celltype)) +
  ggplot2::geom_col(position="dodge") + ggplot2::coord_flip() +
  ggplot2::labs(x="Epitope", y="% of COVID-hit", title="COVID epitopes by lineage") +
  ggplot2::theme_bw(base_size=10)
sp(p_epi, "covid_epitope_mix.png", 9, 6)
