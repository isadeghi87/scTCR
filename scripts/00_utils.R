## ========== 00_utils.R ==========
suppressWarnings(try(setwd('/home/isadeghi/projects/covid_kids/'), silent = TRUE))

## I/O
gex_cd4  <- "./data/CD4/sample_filtered_feature_bc_matrix.h5"
vdj_cd4  <- "./data/CD4/filtered_contig_annotations.csv"
gex_cd8  <- "./data/CD8/sample_filtered_feature_bc_matrix.h5"
vdj_cd8  <- "./data/CD8/filtered_contig_annotations.csv"
vdj_ref  <- "./reference/vdjdb_cdr3_covid_ref.tsv"   # your portal export

script_dir<- "./scripts"
out_dir  <- "out_combined_tcr_covid"
fig_dir  <- file.path(out_dir, "figs")
tab_dir  <- file.path(out_dir, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

## Plot saver
sp <- function(p, file, w=7, h=5) ggsave(filename=file.path(fig_dir, file), plot=p, width=w, height=h, dpi=300)

`%||%` <- function(x, y) if (!is.null(x)) x else y
