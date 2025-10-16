## ========== 05l_spectratype_ridges_and_ecdf.R ==========
## Purpose: (1) Spectratype ridges; (2) ECDFs of TRB CDR3 length with Wilcoxon p-values.

spec_df <- seu@meta.data |>
  tibble::rownames_to_column("barcode") |>
  dplyr::mutate(group = paste0(celltype, "_", ifelse(covid_hit=="COVID-hit","COVID","other")),
                len_trb = nchar(TRB_cdr3)) |>
  dplyr::filter(!is.na(len_trb))

p_ridge <- ggplot2::ggplot(spec_df, ggplot2::aes(x=len_trb, y=group, fill=group)) +
  ggridges::geom_density_ridges(alpha=.6, scale=1.2, rel_min_height=.01) +
  ggplot2::labs(x="TRB CDR3 length (aa)", y=NULL, title="Spectratype of TRB by lineage × COVID status") +
  ggplot2::theme_bw(base_size=11) + ggplot2::theme(legend.position="none")
sp(p_ridge, "spectratype_trb_lineage_covid.png", 8, 6)

cells <- seu@meta.data |> tibble::rownames_to_column("barcode") |> dplyr::select(barcode, celltype, covid_hit, TRA_cdr3, TRB_cdr3)
spec <- rbind(
  cells |> dplyr::filter(!is.na(TRA_cdr3)) |> dplyr::transmute(chain="TRA", group=paste(celltype, ifelse(covid_hit=="COVID-hit","COVID","other"), sep="_"), len=nchar(TRA_cdr3)),
  cells |> dplyr::filter(!is.na(TRB_cdr3)) |> dplyr::transmute(chain="TRB", group=paste(celltype, ifelse(covid_hit=="COVID-hit","COVID","other"), sep="_"), len=nchar(TRB_cdr3))
)
x8 <- spec$len[spec$group=="CD8_COVID"];  y8 <- spec$len[spec$group=="CD8_other"]
x4 <- spec$len[spec$group=="CD4_COVID"];  y4 <- spec$len[spec$group=="CD4_other"]

plot_ecdf <- function(x,y,labA,labB){
  dx <- data.frame(len=x, g=labA); dy <- data.frame(len=y, g=labB)
  ggplot2::ggplot(rbind(dx,dy), ggplot2::aes(len, colour=g)) +
    ggplot2::stat_ecdf(geom="step") +
    ggplot2::scale_x_continuous(breaks=seq(8,22,2)) +
    ggplot2::labs(x="TRB CDR3 length (aa)", y="ECDF", colour=NULL) +
    ggplot2::theme_bw()
}
lab8 <- sprintf("Wilcoxon p = %.2g", suppressWarnings(wilcox.test(x8, y8, exact=FALSE)$p.value))
lab4 <- sprintf("Wilcoxon p = %.2g", suppressWarnings(wilcox.test(x4, y4, exact=FALSE)$p.value))
p_stats <- plot_ecdf(x8,y8,"CD8_COVID","CD8_other")+ggplot2::labs(title="CD8 vs Other", subtitle=lab8) +
           plot_ecdf(x4,y4,"CD4_COVID","CD4_other")+ggplot2::labs(title="CD4 vs Other", subtitle=lab4)
sp(p_stats, "ecdf_trb_len_by_lineage_covid_status.png", 8, 6)
