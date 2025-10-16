## ========== 05n_ct_epitope_barplot.R ==========
## Purpose: show the top clonotypes per epitope among COVID-hit cells only.

cov_only <- subset(seu, subset = covid_hit == "COVID-hit")
p <- cov_only@meta.data %>%
  dplyr::count(CTaa, covid_epitope, sort=TRUE) %>% dplyr::slice_head(n=10) %>%
  ggplot2::ggplot(ggplot2::aes(tidytext::reorder_within(CTaa, n, covid_epitope), n, fill=covid_epitope)) +
  ggplot2::geom_col(show.legend=FALSE) + ggplot2::coord_flip() +
  tidytext::scale_x_reordered() + ggplot2::labs(x="Clonotype", y="# cells", title="Top CTs per epitope (COVID-hit)") +
  ggplot2::theme_bw()
sp(p, "top_cts_per_epitope_covid.png", 7, 6)
