## ========== 05a_define_covid_hit.R ==========
## Purpose: compute 'covid_hit' tag on cells based on desired confidence level.
## Inputs: best_per, paired_any, paired_same tables saved by previous scripts.
## Output: adds `seu$covid_hit` as "COVID-hit" vs "Other".

best_per    <- readr::read_tsv(file.path(tab_dir, "covid_best_per_cell_CD4_CD8.tsv"), show_col_types = FALSE)
paired_any  <- readr::read_tsv(file.path(tab_dir, "covid_paired_any_CD4_CD8.tsv"), show_col_types = FALSE)
paired_same <- readr::read_tsv(file.path(tab_dir, "covid_paired_same_epitope_CD4_CD8.tsv"), show_col_types = FALSE)

cov_barcodes <- switch(
  CONF_LEVEL,
  low    = intersect(best_per$barcode, colnames(seu)),
  medium = intersect(paired_any$barcode[paired_any$paired_any], colnames(seu)),
  high   = intersect(paired_same$barcode[paired_same$paired_same_epitope], colnames(seu))
)

seu$covid_hit <- ifelse(colnames(seu) %in% cov_barcodes, "COVID-hit", "Other")
message("covid_hit set. nCOVID = ", sum(seu$covid_hit=="COVID-hit"))
