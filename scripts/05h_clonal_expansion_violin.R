## ========== 05h_clonal_expansion_violin.R ==========
## Purpose: compare clone size distributions between COVID-like and others (log1p scale).

df_exp <- seu@meta.data |>
  tibble::rownames_to_column("barcode") |>
  dplyr::mutate(clone_size = as.integer(ifelse(is.na(clone_id), 0, ave(clone_id, clone_id, FUN=length)))) |>
  dplyr::mutate(is_covid = factor(covid_like %in% c("low","medium"), levels=c(FALSE,TRUE), labels=c("none","covid-like")))
suppressWarnings(print(wilcox.test(clone_size ~ is_covid, data = df_exp)))

p_clone <- ggplot2::ggplot(df_exp, ggplot2::aes(x=is_covid, y=clone_size, fill=is_covid)) +
  ggplot2::geom_violin(scale="width", trim=TRUE, alpha=.7, color=NA) +
  ggplot2::geom_boxplot(width=.15, outlier.size=.6, fatten=1.2, fill="white", color="black") +
  ggplot2::scale_y_continuous(trans="log1p") + ggplot2::scale_fill_brewer(palette="Set2") +
  ggplot2::labs(x=NULL, y="Clone size (log1p)", title="Clonal expansion in COVID-like vs others") +
  ggplot2::theme_bw(base_size=12) + ggplot2::theme(legend.position="none", panel.grid.minor=element_blank())
sp(p_clone, "clone_size_covid_like_colored.png", 6, 4.5)
