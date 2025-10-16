## ========== 03_covid_mapping_reporting.R ==========
source("/home/isadeghi/projects/covid_kids/scripts/00_utils.R")
seu <- readRDS(file.path(out_dir, "merged_seurat_with_tcr.rds"))
pair_CD4 <- readr::read_tsv(file.path(tab_dir,"pair_CD4.tsv"), show_col_types=FALSE)
pair_CD8 <- readr::read_tsv(file.path(tab_dir,"pair_CD8.tsv"), show_col_types=FALSE)
pair <- dplyr::bind_rows(pair_CD4, pair_CD8)

## --- Read VDJdb export and harmonize ---
stopifnot(file.exists(vdj_ref))
vdjdb <- readr::read_tsv(vdj_ref, show_col_types = FALSE)
names(vdjdb) <- tolower(names(vdjdb))
vdj0 <- vdjdb %>%
  dplyr::transmute(
    gene    = toupper(gene),
    cdr3    = toupper(cdr3),
    epitope = as.character(epitope),
    species = tolower(species),
    mhc_cls = `mhc class`,
    mhc_a   = `mhc a`,
    mhc_b   = `mhc b`,
    method  = method,
    score   = suppressWarnings(as.numeric(score))
  ) %>%
  dplyr::filter(grepl("homo", species), nzchar(cdr3))

## Split refs: CD4→MHCII; CD8→MHCI
split_refs <- function(vdj, mhc_class){
  v <- vdj %>% dplyr::filter(toupper(mhc_cls)==mhc_class)
  list(
    TRB = v %>% dplyr::filter(grepl("TRB|BETA", gene)) %>% dplyr::distinct(cdr3, .keep_all=TRUE),
    TRA = v %>% dplyr::filter(grepl("TRA|ALPHA", gene)) %>% dplyr::distinct(cdr3, .keep_all=TRUE)
  )
}
ref_cd4 <- split_refs(vdj0, "MHCII")
ref_cd8 <- split_refs(vdj0, "MHCI")

## Build clone table from merged Seurat (+ V/J if available)
vj <- try(pair %>% dplyr::select(barcode, TRA_v, TRA_j, TRB_v, TRB_j), silent=TRUE)
if (inherits(vj,"try-error")) vj <- NULL

clones_df <- seu@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::select(barcode, orig.ident, celltype, seurat_clusters, CTaa, clone_id, TRA_cdr3, TRB_cdr3) %>%
  dplyr::mutate(TRA_cdr3 = toupper(TRA_cdr3), TRB_cdr3 = toupper(TRB_cdr3)) %>%
  { if (!is.null(vj)) dplyr::left_join(., vj, by="barcode") else . }

## Matcher (exact + optional fuzzy), with safe joins
match_chain <- function(df_clones, ref_tbl, chain=c("TRB","TRA"), max_dist=1){
  chain <- match.arg(chain)
  cdr3_col <- if (chain=="TRB") "TRB_cdr3" else "TRA_cdr3"
  exact <- df_clones %>%
    dplyr::filter(!is.na(.data[[cdr3_col]])) %>%
    dplyr::inner_join(ref_tbl %>% dplyr::select(cdr3, epitope, mhc_cls, mhc_a, mhc_b),
                      by = setNames("cdr3", cdr3_col)) %>%
    dplyr::mutate(match_chain = chain, match_type="exact", cdr3_ref = .data[[cdr3_col]], dist=0L)
  
  fuzzy <- NULL
  if (max_dist > 0 && nrow(ref_tbl) > 0) {
    qs <- df_clones %>% dplyr::filter(!is.na(.data[[cdr3_col]])) %>%
      dplyr::distinct(.data[[cdr3_col]]) %>% dplyr::pull()
    fuzz_core <- lapply(qs, function(q){
      d <- stringdist::stringdist(q, ref_tbl$cdr3, method="lv")
      idx <- which(d <= max_dist & nchar(ref_tbl$cdr3) == nchar(q))
      if (!length(idx)) return(NULL)
      tibble::tibble(!!cdr3_col := q,
                     cdr3_ref = ref_tbl$cdr3[idx],
                     dist     = d[idx],
                     epitope  = ref_tbl$epitope[idx],
                     mhc_cls  = ref_tbl$mhc_cls[idx],
                     mhc_a    = ref_tbl$mhc_a[idx],
                     mhc_b    = ref_tbl$mhc_b[idx],
                     match_chain = chain,
                     match_type  = paste0("fuzzy", max_dist))
    }) %>% dplyr::bind_rows()
    if (!is.null(fuzz_core) && nrow(fuzz_core)) {
      map <- df_clones %>% dplyr::select(barcode, !!rlang::sym(cdr3_col)) %>% dplyr::distinct()
      fuzzy <- dplyr::inner_join(map, fuzz_core, by = cdr3_col,
                                 relationship = if ("relationship" %in% names(formals(dplyr::inner_join))) "many-to-many" else NULL)
    }
  }
  dplyr::bind_rows(exact, fuzzy)
}

## Lineage-aware matching rules:
##   CD4 (MHCII): TRA exact, TRB exact±fuzzy1
##   CD8 (MHCI):  TRA exact, TRB exact±fuzzy1
run_lineage <- function(df, refs, lineage){
  df_lin <- df %>% dplyr::filter(celltype == lineage)
  h_trb  <- match_chain(df_lin, refs$TRB, chain="TRB", max_dist=1)
  h_tra  <- match_chain(df_lin, refs$TRA, chain="TRA", max_dist=0)
  dplyr::bind_rows(h_trb, h_tra) %>% dplyr::mutate(lineage = lineage) %>% dplyr::distinct()
}

hits_cd4 <- run_lineage(clones_df, ref_cd4, "CD4")
hits_cd8 <- run_lineage(clones_df, ref_cd8, "CD8")
hits_all <- dplyr::bind_rows(hits_cd4, hits_cd8) %>% dplyr::distinct()

## Confidence ladder & paired flags
hits_all <- hits_all %>%
  dplyr::mutate(priority = dplyr::case_when(
    match_chain=="TRB" & match_type=="exact" ~ 3L,
    match_chain=="TRA" & match_type=="exact" ~ 2L,
    match_chain=="TRB" & grepl("^fuzzy", match_type) ~ 1L,
    TRUE ~ 0L
  ))
paired_any <- hits_all %>% dplyr::group_by(barcode) %>%
  dplyr::summarise(paired_any = all(c("TRB","TRA") %in% match_chain), .groups="drop")
paired_same <- hits_all %>% dplyr::group_by(barcode) %>%
  dplyr::summarise(paired_same_epitope = {
    e_trb <- unique(epitope[match_chain=="TRB"]); e_tra <- unique(epitope[match_chain=="TRA"])
    length(intersect(e_trb, e_tra)) > 0
  }, .groups="drop")
best_per_cell <- hits_all %>% dplyr::group_by(barcode) %>% dplyr::slice_max(priority, with_ties=FALSE) %>% dplyr::ungroup()

## Annotate merged Seurat
seu$covid_like <- "none"
seu$covid_like[colnames(seu) %in% best_per_cell$barcode] <- "low"
seu$covid_like[colnames(seu) %in% paired_any$barcode[paired_any$paired_any]] <- "medium"
seu$covid_like[colnames(seu) %in% paired_same$barcode[paired_same$paired_same_epitope]] <- "high"

## Save outputs
readr::write_tsv(hits_all,        file.path(tab_dir, "covid_hits_CD4_CD8.tsv"))
readr::write_tsv(best_per_cell,   file.path(tab_dir, "covid_best_per_cell_CD4_CD8.tsv"))
readr::write_tsv(paired_any,      file.path(tab_dir, "covid_paired_any_CD4_CD8.tsv"))
readr::write_tsv(paired_same,     file.path(tab_dir, "covid_paired_same_epitope_CD4_CD8.tsv"))

## Plots
sp(DimPlot(seu, group.by="covid_like") + ggtitle("COVID-like confidence (CD4+CD8)"),
   "umap_covid_like_confidence_CD4_CD8.png", 6,5)
sp(hits_all %>% dplyr::count(lineage, match_chain, match_type) %>%
     ggplot(aes(interaction(lineage,match_chain), n, fill=match_type)) +
     geom_col() + xlab("Lineage × Chain") + ylab("# matches") +
     ggtitle("SARS-CoV-2 matches by lineage & chain"),
   "covid_hits_by_lineage_chain.png", 7,4)

## Minimal report index
report <- c(
  "# Combined CD4+CD8 TCR/COVID analysis",
  "",
  "## Figures",
  "- CD4 UMAP: figs/umap_CD4.png",
  "- CD8 UMAP: figs/umap_CD8.png",
  "- Merged UMAP (celltype): figs/umap_merged_celltype.png",
  "- Homeostasis: figs/clonal_homeostasis_merged.png",
  "- V usage TRB: figs/vj_usage_TRB_V_merged.png",
  "- J usage TRB: figs/vj_usage_TRB_J_merged.png",
  "- COVID confidence (CD4+CD8): figs/umap_covid_like_confidence_CD4_CD8.png",
  "- COVID lineage×chain: figs/covid_hits_by_lineage_chain.png",
  "",
  "## Tables",
  "- CD4 pairs: tables/pair_CD4.tsv",
  "- CD8 pairs: tables/pair_CD8.tsv",
  "- COVID hits (all): tables/covid_hits_CD4_CD8.tsv",
  "- Best per cell: tables/covid_best_per_cell_CD4_CD8.tsv",
  "- Paired any: tables/covid_paired_any_CD4_CD8.tsv",
  "- Paired same epitope: tables/covid_paired_same_epitope_CD4_CD8.tsv"
)
writeLines(report, file.path(out_dir, "REPORT.md"))

## Save final object
saveRDS(seu, file.path(out_dir, "merged_seurat_final.rds"))
message("Done. See ", out_dir, "/REPORT.md and figs/tables.")
