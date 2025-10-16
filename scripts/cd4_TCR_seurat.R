## ====== 0) Setup ======
## wd
suppressWarnings(try(setwd('/home/isadeghi/projects/covid_kids/'), silent = TRUE))

## Paths to your uploaded CD4 files (adjust if needed)
vdj_csv   <- "./data/CD4/filtered_contig_annotations.csv"
gex_h5    <- "./data/CD4/sample_filtered_feature_bc_matrix.h5"

## Output folders
out_dir   <- "out_cd4_tcr_covid"
fig_dir   <- file.path(out_dir, "figs")
tab_dir   <- file.path(out_dir, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

## Packages
req_pkgs <- c(
  "Seurat","SeuratObject","dplyr","tidyr","ggplot2","patchwork","readr","stringr",
  "Matrix","data.table","scRepertoire","immunarch","stringdist","ComplexHeatmap","circlize"
)
for(p in req_pkgs){
  if(!requireNamespace(p, quietly = TRUE)){
    message("Installing ", p, " ...")
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
library(Seurat); library(SeuratObject); library(dplyr); library(tidyr); library(ggplot2)
library(patchwork); library(readr); library(stringr); library(Matrix); library(data.table)
library(scRepertoire); library(immunarch); library(stringdist); library(ComplexHeatmap); library(circlize)

theme_set(theme_bw(base_size = 12))

## Helper to save plots
sp <- function(p, file, w=7, h=5){ ggsave(filename=file.path(fig_dir, file), plot=p, width=w, height=h, dpi=150) }

## ====== 1) Load data (GEX + VDJ) ======
# 10x filtered feature-barcode matrix (H5)
seu <- tryCatch({
  Seurat::Read10X_h5(gex_h5)
}, error=function(e){
  stop("Could not read H5: ", e$message, "\nMake sure it's a 10x filtered_feature_bc_matrix .h5")
}) |> CreateSeuratObject(project = "CD4", min.cells = 3, min.features = 200)

# Basic QC
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
VlnPlot(seu, c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3) |> sp("qc_violin.png")

# Reasonable CD4 T filtering (keep lenient; adjust to your dataset)
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

# Normalize & reduce
seu <- SCTransform(seu, vst.flavor = "v2", verbose = FALSE) |>
  RunPCA(verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30)
seu <- FindNeighbors(seu, dims = 1:30) |> FindClusters(resolution = 0.6)

p1 <- DimPlot(seu, group.by="seurat_clusters", label=TRUE) + ggtitle("CD4 clustering")
sp(p1, "umap_clusters.png")

## ====== 2) Annotate CD4 T-cell states (quick markers) ======
# Score canonical CD4 states (naive, TCM, TEM, Tfh, Treg, activated/exhausted)
gene_sets <- list(
  CD4_naive = c("CCR7","LEF1","TCF7","IL7R","LST1"),
  CD4_TCM   = c("CCR7","IL7R","MAL","SELL"),
  CD4_TEM   = c("CXCR3","GZMK","IFNG","CCL5"),
  Tfh       = c("CXCR5","BCL6","PDCD1","ICOS"),
  Treg      = c("FOXP3","IL2RA","CTLA4","IKZF2"),
  Act_Exh   = c("PDCD1","HAVCR2","LAG3","TIGIT","GZMB","PRF1")
)
# Add module scores
for(nm in names(gene_sets)){
  genes <- intersect(gene_sets[[nm]], rownames(seu))
  if(length(genes)>=2) seu <- AddModuleScore(seu, list(genes), name = paste0(nm,"_score"), search = FALSE)
}
# Pick best label per cell
score_cols <- grep("_score1$", colnames(seu@meta.data), value=TRUE)
if(length(score_cols)>0){
  ms <- seu@meta.data[,score_cols,drop=FALSE]
  labs <- colnames(ms)[max.col(ms, ties.method = "first")]
  labs <- gsub("_score1$","",labs)
  seu$cd4_state <- labs
} else {
  seu$cd4_state <- "CD4_T"
}
p2 <- DimPlot(seu, group.by="cd4_state", label=FALSE) + ggtitle("CD4 state (module score heuristic)")
sp(p2, "umap_cd4_state.png")

## ====== 3) Read and prepare V(D)J contigs ======
contigs <- fread(vdj_csv) |> as.data.frame()
# Expect columns like: barcode, chain, cdr3, cdr3_nt, v_gene, j_gene, productive, is_cell, umis, reads, locus, etc.
# Keep only high-confidence, productive
contigs <- contigs |>
  filter(is_cell == "true" | is_cell == TRUE | is.na(is_cell)) |>
  filter(productive == "True" | productive == TRUE)

# Resolve multiple chains: keep top-UMI per locus
contigs <- contigs |>
  mutate(locus = ifelse(!is.na(chain), chain, locus)) |>
  group_by(barcode, locus) |>
  slice_max(order_by = as.numeric(umis %||% reads %||% 0), n = 1, with_ties = FALSE) |>
  ungroup()

# Make per-barcode paired clonotype keys (TRA + TRB CDR3aa if present; fallback to TRB only)
tra <- contigs |> filter(grepl("TRA|TCRA|Alpha", locus, ignore.case = TRUE)) |>
  select(barcode, TRA_cdr3 = cdr3, TRA_v=v_gene, TRA_j=j_gene)
trb <- contigs |> filter(grepl("TRB|TCRB|Beta",  locus, ignore.case = TRUE)) |>
  select(barcode, TRB_cdr3 = cdr3, TRB_v=v_gene, TRB_j=j_gene)

pair  <- full_join(tra, trb, by="barcode")
pair <- pair |>
  mutate(clone_def = case_when(
    !is.na(TRA_cdr3) & !is.na(TRB_cdr3) ~ paste0("TRA:",TRA_cdr3,"|TRB:",TRB_cdr3),
    is.na(TRA_cdr3) & !is.na(TRB_cdr3)  ~ paste0("TRB:",TRB_cdr3),
    TRUE ~ NA_character_
  ))
pair$clone_id <- as.character(factor(pair$clone_def))
## Normalize VDJ barcodes to include "-1"
if (!any(grepl("-\\d+$", pair$barcode))) {
  pair$barcode <- paste0(pair$barcode, "-1")
}

## Attach to Seurat (simple exact match)
hits <- intersect(colnames(seu), pair$barcode)
idx  <- match(hits, pair$barcode)
seu$TRA_cdr3 <- NA_character_
seu$TRB_cdr3 <- NA_character_
seu$clone_id <- NA_character_

row.names(pair) <- pair$barcode
## Assign by position (coerce to character just in case)
seu$TRA_cdr3[hits] <- as.character(pair$TRA_cdr3[idx])
seu$TRB_cdr3[hits] <- as.character(pair$TRB_cdr3[idx])
seu$clone_id[hits] <- as.character(pair$clone_id[idx])

## Derive clone_size and quick check
seu$clone_size <- ifelse(is.na(seu$clone_id), 0, ave(seu$clone_id, seu$clone_id, FUN=length))
cat("Matched cells:", length(hits), "\n")
head(seu@meta.data[, c("TRA_cdr3","TRB_cdr3","clone_id","clone_size")])

# Quick TCR QC summaries
qc_tcr <- list(
  n_cells            = ncol(seu),
  n_with_TRB         = sum(!is.na(seu$TRB_cdr3)),
  n_with_TRA_TRB     = sum(!is.na(seu$TRB_cdr3) & !is.na(seu$TRA_cdr3)),
  n_clonotypes       = n_distinct(na.omit(seu$clone_id))
)

write_tsv(as.data.frame(qc_tcr), file.path(tab_dir, "tcr_qc_summary.tsv"))

## Visualize clonal expansion on UMAP
seu$clone_size <- ifelse(is.na(seu$clone_id), 0, ave(seu$clone_id, seu$clone_id, FUN=length))
p3 <- DimPlot(seu, group.by="clone_size") + ggtitle("Clone size on UMAP")
sp(p3, "umap_clone_size.png")

## ====== 4) scRepertoire integration (abundance, diversity, V/J usage) ======
# Build a scRepertoire-like object from contigs
# scRepertoire expects a list per sample; we have one sample ("CD4")
contigs_cd4 <- contigs
names(contigs_cd4) <- make.names(names(contigs_cd4))
contigs_list <- list(CD4 = contigs_cd4)

contigs_sr <- contigs %>%
  mutate(
    chain   = if ("chain"   %in% names(.)) chain   else locus,
    cdr3    = if ("cdr3"    %in% names(.)) cdr3    else NA_character_,
    cdr3_nt = if ("cdr3_nt" %in% names(.)) cdr3_nt else NA_character_,
    v_gene  = if ("v_gene"  %in% names(.)) v_gene  else NA_character_,
    j_gene  = if ("j_gene"  %in% names(.)) j_gene  else NA_character_
  ) %>%
  select(any_of(c("barcode","chain","cdr3","cdr3_nt","v_gene","j_gene")))


# 1) Strip any sample prefix: keep everything AFTER the last underscore
for (i in seq_along(clono)) {
  if ("barcode" %in% names(clono[[i]])) {
    clono[[i]]$barcode <- sub("^.*_(?=[^_]+$)", "", clono[[i]]$barcode, perl = TRUE)
    # sanity: show first few after stripping
    if (i == 1) print(head(clono[[i]]$barcode))
  }
}

# 2) Choose the actual clonotype column present (your print shows these exist)
clone_col <- if ("CTaa" %in% names(clono[[1]])) "CTaa" else
  if ("CTnt" %in% names(clono[[1]])) "CTnt" else
    if ("CTgene" %in% names(clono[[1]])) "CTgene" else
      stop("No clonotype column (CTaa/CTnt/CTgene) found in clono[[1]].")

# 3) Ensure grouping column exists on Seurat
if (!"orig.ident" %in% colnames(seu@meta.data)) seu$orig.ident <- "CD4"

# 4) Detect correct group arg for your installed version
ce_formals <- names(formals(scRepertoire::combineExpression))
group_arg  <- if ("groupBy" %in% ce_formals) "groupBy" else "group.by"

# 5) Call combineExpression with REAL clonotype column and normalized barcodes
ce_args <- list(
  clono,                      # contig list
  seu,                        # Seurat object
  cloneCall = clone_col,
  filterNA  = TRUE
)
ce_args[[group_arg]] <- "orig.ident"

seu <- combineExpression( clono,                      # contig list
                          seu,                        # Seurat object
                          cloneCall = "CTaa",
                          chain = "both",
                          filterNA  = TRUE)

# 6) Quick sanity checks
print(intersect(c("CTaa","CTnt","CTgene","cloneType"), colnames(seu@meta.data)))
cat("Mapped cells (non-NA ", clone_col, "): ",
    sum(!is.na(seu@meta.data[[clone_col]])), " / ", ncol(seu), "\n", sep = "")

## --- Sanity checks ---
# Did scRepertoire add CTaa?
intersect(c("CTaa","CTnt","CTgene","cloneType"), colnames(seu@meta.data))

# How many cells have TRB and paired TRA+TRB from our own attachment?
cat("With TRB:", sum(!is.na(seu$TRB_cdr3)),
    " | Paired:", sum(!is.na(seu$TRB_cdr3) & !is.na(seu$TRA_cdr3)),
    " | Unique clonotypes:", dplyr::n_distinct(na.omit(seu$clone_id)), "\n")


# --- Build a tidy clono table and compare TRA/TRB with Seurat attachment ---
cl1 <- clono[[1]][, c("barcode","TCR1","cdr3_aa1","TCR2","cdr3_aa2","CTaa")]
# strip the "CD4_CD4_patient_" (anything up to last "_")
cl1$barcode <- sub("^.*_(?=[^_]+$)", "", cl1$barcode, perl = TRUE)

# determine which cdr3 belongs to TRA vs TRB
is_alpha1 <- grepl("^TRAV|^TRA", cl1$TCR1)
is_beta1  <- grepl("^TRBV|^TRB", cl1$TCR1)
is_alpha2 <- grepl("^TRAV|^TRA", cl1$TCR2)
is_beta2  <- grepl("^TRBV|^TRB", cl1$TCR2)

cl1$TRA_from_clono <- ifelse(is_alpha1, cl1$cdr3_aa1,
                             ifelse(is_alpha2, cl1$cdr3_aa2, NA))
cl1$TRB_from_clono <- ifelse(is_beta1,  cl1$cdr3_aa1,
                             ifelse(is_beta2,  cl1$cdr3_aa2,  NA))

meta <- seu@meta.data |> tibble::rownames_to_column("barcode")

# TRB agreement
chk_trb <- dplyr::left_join(meta[, c("barcode","TRB_cdr3")], 
                            cl1[, c("barcode","TRB_from_clono")], by="barcode")
sub_trb <- subset(chk_trb, !is.na(TRB_cdr3) | !is.na(TRB_from_clono))
agree_trb <- sum(sub_trb$TRB_cdr3 == sub_trb$TRB_from_clono, na.rm = TRUE)
cat("TRB CDR3 agreement: ", agree_trb, "/", nrow(sub_trb), " (",
    round(100*agree_trb/nrow(sub_trb), 2), "%)\n", sep = "")
if (agree_trb < nrow(sub_trb)) print(head(sub_trb[sub_trb$TRB_cdr3 != sub_trb$TRB_from_clono, ], 10))

# TRA agreement
chk_tra <- dplyr::left_join(meta[, c("barcode","TRA_cdr3")], 
                            cl1[, c("barcode","TRA_from_clono")], by="barcode")
sub_tra <- subset(chk_tra, !is.na(TRA_cdr3) | !is.na(TRA_from_clono))
agree_tra <- sum(sub_tra$TRA_cdr3 == sub_tra$TRA_from_clono, na.rm = TRUE)
cat("TRA CDR3 agreement: ", agree_tra, "/", nrow(sub_tra), " (",
    round(100*agree_tra/nrow(sub_tra), 2), "%)\n", sep = "")
if (agree_tra < nrow(sub_tra)) print(head(sub_tra[sub_tra$TRA_cdr3 != sub_tra$TRA_from_clono, ], 10))

# Each clone_id should map to a single CTaa (sanity vs scRepertoire)
tmp <- meta |> dplyr::filter(!is.na(clone_id)) |> dplyr::count(clone_id, CTaa)
ambig <- tmp |> dplyr::count(clone_id) |> dplyr::filter(n > 1)
cat("Clones with >1 CTaa label: ", nrow(ambig), "\n", sep = "")
if (nrow(ambig) > 0) print(head(dplyr::semi_join(tmp, ambig, by="clone_id"), 20))

## --- Diversity & plots (exports to out_cd4_tcr_covid) ---
set.seed(123)
clone_call <- "CTaa"  # you confirmed this exists
grp <- "seurat_clusters"

p_div <- scRepertoire::clonalDiversity(
  seu,
  cloneCall = clone_call,
  group     = grp,
  skip.boots = TRUE  # no downsampling since you're within one sample
)
# If your version returns a ggplot, just save it:
if (inherits(p_div, "gg")) {
  sp(p_div + ggtitle("Shannon diversity by cluster (CTaa)"),
     "diversity_shannon_by_cluster.png", w=6, h=8)
} else {
  # Some versions return a data.frame; make the plot yourself:
  div_tab <- as.data.frame(p_div)
  readr::write_tsv(div_tab, file.path(tab_dir, "diversity_shannon_by_cluster.tsv"))
  gp <- ggplot(div_tab, aes(x=Group, y=Shannon, color=Group)) +
    geom_point(size=2) + theme_bw() +
    ylab("Shannon Index Score") + ggtitle("Shannon diversity by cluster (CTaa)")
  sp(gp, "diversity_shannon_by_cluster.png", w=6, h=8)
}

p_homeo <- scRepertoire::clonalHomeostasis(seu, cloneCall = "aa")
sp(p_homeo, "clonal_homeostasis_cd4.png")

p_v_trb <- try(scRepertoire::vizGenes(seu, gene="V", chain="TRB", plot=FALSE), silent=TRUE)
if(!inherits(p_v_trb,"try-error")) sp(p_v_trb, "vj_usage_TRB_V.png", w=8, h=6)

p_j_trb <- try(scRepertoire::vizGenes(seu, gene="J", chain="TRB", plot=FALSE), silent=TRUE)
if(!inherits(p_j_trb,"try-error")) sp(p_j_trb, "vj_usage_TRB_J.png", w=8, h=6)


# CDR3 length distribution (TRB AA)
cdr3_len <- seu@meta.data |> filter(!is.na(TRB_cdr3)) |>
  transmute(len = nchar(TRB_cdr3))
p_len <- ggplot(cdr3_len, aes(len)) + geom_histogram(binwidth = 1) +
  xlab("TRB CDR3 length (AA)") + ylab("Cells") + ggtitle("CDR3β spectratype")
sp(p_len, "cdr3_len_TRB.png")

## ====== 5) Map clonotypes to SARS-CoV-2 (COVID-19) specificity ======
# Strategy:
#  A) Exact AA match of TRB CDR3 to VDJdb SARS-CoV-2 entries
#  B) Fuzzy match (1–2 AA distance) as supportive evidence (not definitive)
#
# immunarch ships a 'vdjdb' dataset; if not found, it will skip gracefully.
vdj_tsv = './reference/vdjdb_cdr3_covid_ref.tsv'
vdjdb <- read_tsv(vdj_tsv)

# Harmonize names → lower, and standardize key fields
names(vdjdb) <- tolower(names(vdjdb))
names(vdj) <- tolower(names(vdj))
vdj <- vdj %>%
  transmute(
    gene = toupper(gene),
    cdr3 = toupper(cdr3),
    v_gene =  v,
    j_gene =  j,
    species = tolower(species),
    epitope = as.character( epitope),
    epitope_gene = `epitope gene`,
    epitope_species = `epitope species`,
    mhc_a = `mhc a`,
    mhc_b = `mhc b`,
    mhc_cls = `mhc class`,
    reference =reference,
    method = method
  ) %>%
  filter(grepl("homo", species, ignore.case = TRUE), nzchar(cdr3))


# Expect: "gene","cdr3","v","j","species","mhc a","mhc b","mhc class","epitope","epitope gene","epitope species",...
if (!all(c("gene","cdr3","species","epitope") %in% names(vdjdb))) {
  stop("Unexpected VDJdb columns: ", paste(names(vdjdb), collapse=", "))
}

## ==== Keep both TRB & TRA when mapping to VDJdb (exact + fuzzy) ====
## ---- Tightened matching with confidence ladder ----
library(dplyr); library(tidyr); library(stringdist); library(tibble); library(ggplot2)

# 0) Clean ref (your vdjdb already read)
ref0 <- vdjdb
names(ref0) <- tolower(names(ref0))
ref0 <- ref0 %>%
  transmute(
    gene  = toupper(gene),
    cdr3  = toupper(cdr3),
    epitope = as.character(epitope),
    species = tolower(species),
    mhc_a = `mhc a`, mhc_b = `mhc b`, mhc_cls = `mhc class`,
    method = method, score = suppressWarnings(as.numeric(score))
  ) %>%
  filter(grepl("homo", species), nzchar(cdr3))

## ---- Tight, CD4-aware SARS-CoV-2 mapping (MHCII focus) ----
library(dplyr); library(tidyr); library(stringdist); library(tibble); library(ggplot2)

# knobs
MAX_DIST_TRB <- 0   # TRB exact-only (raise to 1 if you want a small fuzzy allowance)
MAX_DIST_TRA <- 0   # TRA exact-only
KEEP_SINGLECELL <- TRUE
RESTRICT_MHC <- "MHCII"    # for CD4; set NULL to allow both classes
REQUIRE_VJ_CONCORDANCE <- TRUE
PATIENT_HLA <- NULL        # e.g., c("HLA-DRB1*15:01"); keep NULL if unknown

# 0) Clean/reference table from your vdjdb data.frame
if (!is.null(RESTRICT_MHC)) {
  ref0 <- ref0 %>% filter(toupper(mhc_cls) == toupper(RESTRICT_MHC))
}

if (KEEP_SINGLECELL) {
  ref0 <- ref0 %>% filter(grepl("single", method, ignore.case = TRUE))
}
if (!is.null(PATIENT_HLA)) {
  any_hla <- function(x) any(na.omit(x) %in% PATIENT_HLA)
  ref0 <- ref0 %>% rowwise() %>% filter(any_hla(c(mhc_a, mhc_b))) %>% ungroup()
}

ref_trb <- ref0 %>% filter(grepl("TRB|BETA", gene)) %>% distinct(cdr3, .keep_all = TRUE)
ref_tra <- ref0 %>% filter(grepl("TRA|ALPHA", gene)) %>% distinct(cdr3, .keep_all = TRUE)

# 1) Clones from Seurat + V/J (if you still have them in 'pair')
# If you kept 'pair' with TRA_v/TRA_j/TRB_v/TRB_j earlier:
vj <- try(pair %>% select(barcode, TRA_v, TRA_j, TRB_v, TRB_j), silent = TRUE)
if (!inherits(vj, "try-error")) {
  row.names(vj) <- vj$barcode
} else {
  vj <- NULL
}

clones_df <- seu@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  select(barcode, cd4_state, seurat_clusters, CTaa, clone_id, TRA_cdr3, TRB_cdr3) %>%
  mutate(TRA_cdr3 = toupper(TRA_cdr3), TRB_cdr3 = toupper(TRB_cdr3)) %>%
  { if (!is.null(vj)) left_join(., vj, by = "barcode") else . }

# 2) Matching (exact-only per knobs)
match_chain_strict <- function(df_clones, ref_tbl, chain = c("TRB","TRA")){
  chain <- match.arg(chain)
  cdr3_col <- if (chain=="TRB") "TRB_cdr3" else "TRA_cdr3"
  v_col    <- if (chain=="TRB") "TRB_v"    else "TRA_v"
  j_col    <- if (chain=="TRB") "TRB_j"    else "TRA_j"
  
  res <- df_clones %>%
    filter(!is.na(.data[[cdr3_col]])) %>%
    inner_join(ref_tbl %>% select(cdr3, epitope, mhc_cls, mhc_a, mhc_b, v_ref, j_ref),
               by = setNames("cdr3", cdr3_col)) %>%
    mutate(match_chain = chain, match_type = "exact", cdr3_ref = .data[[cdr3_col]], dist = 0L)
  
  if (REQUIRE_VJ_CONCORDANCE && !is.null(vj)) {
    res <- res %>%
      filter(is.na(v_ref) | is.na(.data[[v_col]]) | grepl(substr(v_ref, 1, 5), .data[[v_col]])) %>%
      filter(is.na(j_ref) | is.na(.data[[j_col]]) | grepl(substr(j_ref, 1, 5), .data[[j_col]]))
  }
  res
}

hits_trb <- match_chain_strict(clones_df, ref_trb, chain = "TRB")
hits_tra <- match_chain_strict(clones_df, ref_tra, chain = "TRA")

# (Optional) allow TRB fuzzy1 on top of exact-only baseline
if (MAX_DIST_TRB > 0 && nrow(ref_trb) > 0) {
  uniq_trb <- clones_df %>% filter(!is.na(TRB_cdr3)) %>% distinct(TRB_cdr3) %>% pull()
  fuzzy_list <- lapply(uniq_trb, function(q){
    d <- stringdist::stringdist(q, ref_trb$cdr3, method = "lv")
    idx <- which(d <= MAX_DIST_TRB & nchar(ref_trb$cdr3) == nchar(q))
    if (!length(idx)) return(NULL)
    tibble(
      barcode = clones_df$barcode[clones_df$TRB_cdr3 == q],
      TRB_cdr3 = q,
      cdr3_ref = ref_trb$cdr3[idx],
      dist = d[idx],
      epitope = ref_trb$epitope[idx],
      mhc_cls = ref_trb$mhc_cls[idx],
      mhc_a   = ref_trb$mhc_a[idx],
      mhc_b   = ref_trb$mhc_b[idx],
      match_chain = "TRB", match_type = paste0("fuzzy", MAX_DIST_TRB)
    )
  }) %>% bind_rows()
  if (!is.null(fuzzy_list) && nrow(fuzzy_list)) {
    hits_trb <- bind_rows(hits_trb, fuzzy_list) %>% distinct()
  }
}

hits <- bind_rows(hits_trb, hits_tra) %>% distinct()

# 3) Paired evidence & confidence ladder
paired_any <- hits %>% group_by(barcode) %>% summarise(paired_any = all(c("TRB","TRA") %in% match_chain), .groups="drop")
paired_same_epi <- hits %>% group_by(barcode) %>%
  summarise(paired_same_epitope = {
    e_trb <- unique(epitope[match_chain=="TRB"])
    e_tra <- unique(epitope[match_chain=="TRA"])
    length(intersect(e_trb, e_tra)) > 0
  }, .groups="drop")

hits <- hits %>%
  mutate(priority = case_when(
    match_chain=="TRB" & match_type=="exact" ~ 3L,
    match_chain=="TRA" & match_type=="exact" ~ 2L,
    match_chain=="TRB" & grepl("^fuzzy", match_type) ~ 1L,
    TRUE ~ 0L
  ))

best_per_cell <- hits %>% group_by(barcode) %>% slice_max(priority, with_ties = FALSE) %>% ungroup()

# 4) Annotate Seurat
seu$covid_like <- "none"
seu$covid_like[colnames(seu) %in% best_per_cell$barcode] <- "low"
seu$covid_like[colnames(seu) %in% paired_any$barcode[paired_any$paired_any]] <- "medium"
seu$covid_like[colnames(seu) %in% paired_same_epi$barcode[paired_same_epi$paired_same_epitope]] <- "high"

# 5) Updated summary
sum_tab2 <- list(
  n_hits_trb_exact = sum(hits$match_chain=="TRB" & hits$match_type=="exact"),
  n_hits_trb_fuzzy = sum(hits$match_chain=="TRB" & grepl("^fuzzy", hits$match_type)),
  n_hits_tra_exact = sum(hits$match_chain=="TRA" & hits$match_type=="exact"),
  cells_any        = n_distinct(hits$barcode),
  cells_paired_any = sum(paired_any$paired_any),
  cells_paired_same_epitope = sum(paired_same_epi$paired_same_epitope)
)
print(sum_tab2)

match_bp = hits %>% count(match_chain, match_type) %>%
  ggplot(aes(match_chain, n, fill=match_type)) + geom_col() +
  ylab("# matches") + ggtitle("SARS-CoV-2 matches (MHCII, strict)")

# Plots/tables
ggsave(file.path(fig_dir, "covid_hits_by_chain_cd4_mhcii.png"),match_bp,
       width=6, height=4, dpi=150)

ggsave(file.path(fig_dir, "umap_covid_like_confidence_cd4.png"),
       DimPlot(seu, group.by="covid_like") + ggtitle("COVID-like confidence (CD4, MHCII)"),
       width=6, height=5, dpi=150)

readr::write_tsv(hits,             file.path(tab_dir, "covid_hits_cd4_MHCII_strict.tsv"))
readr::write_tsv(best_per_cell,    file.path(tab_dir, "covid_best_call_per_cell_cd4_MHCII.tsv"))
readr::write_tsv(paired_any,       file.path(tab_dir, "covid_paired_any_cd4_MHCII.tsv"))
readr::write_tsv(paired_same_epi,  file.path(tab_dir, "covid_paired_same_epitope_cd4_MHCII.tsv"))

# helpers
`%NA%` <- function(x, y) ifelse(is.na(x) | x == "", y, x)
`%||%` <- function(x, y) if (!is.null(x)) x else y

# split references
ref_trb <- vdj %>% filter(grepl("TRB|BETA", gene, ignore.case = TRUE)) %>% distinct(cdr3, .keep_all = TRUE)
ref_tra <- vdj %>% filter(grepl("TRA|ALPHA", gene, ignore.case = TRUE)) %>% distinct(cdr3, .keep_all = TRUE)

# 1) Build your clonotype table (with both α/β)
clones_df <- seu@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::select(barcode, cd4_state, seurat_clusters, CTaa, clone_id, TRA_cdr3, TRB_cdr3) %>%
  dplyr::mutate(TRA_cdr3 = toupper(TRA_cdr3), TRB_cdr3 = toupper(TRB_cdr3))

# 2) A function to get exact + fuzzy matches for a chain
## Fixed matcher: exact + fuzzy (≤1 AA), with correct joins
match_chain <- function(df_clones, ref_tbl, chain = c("TRB","TRA"), max_dist = 1){
  chain <- match.arg(chain)
  cdr3_col <- if (chain == "TRB") "TRB_cdr3" else "TRA_cdr3"
  
  # ---- Exact matches (join by CDR3 only)
  exact <- df_clones %>%
    dplyr::filter(!is.na(.data[[cdr3_col]])) %>%
    dplyr::inner_join(
      ref_tbl %>% dplyr::select(cdr3, epitope, epitope_gene, epitope_species, mhc_a, mhc_b, mhc_cls, reference, method),
      by = setNames("cdr3", cdr3_col)
    ) %>%
    dplyr::mutate(match_chain = chain, match_type = "exact")
  
  # ---- Fuzzy matches (≤ max_dist, same length), then attach barcode by CDR3
  uniq_q <- df_clones %>% dplyr::filter(!is.na(.data[[cdr3_col]])) %>%
    dplyr::distinct(.data[[cdr3_col]]) %>% dplyr::pull()
  fuzzy <- NULL
  if (length(uniq_q) > 0 && nrow(ref_tbl) > 0) {
    fuzzy_list <- lapply(uniq_q, function(q){
      d <- stringdist::stringdist(q, ref_tbl$cdr3, method = "lv")
      idx <- which(d <= max_dist & nchar(ref_tbl$cdr3) == nchar(q))
      if (length(idx) == 0) return(NULL)
      tibble::tibble(!!cdr3_col := q,
                     cdr3 = ref_tbl$cdr3[idx],
                     dist = d[idx],
                     epitope = ref_tbl$epitope[idx],
                     epitope_gene = ref_tbl$epitope_gene[idx],
                     epitope_species = ref_tbl$epitope_species[idx],
                     mhc_a = ref_tbl$mhc_a[idx],
                     mhc_b = ref_tbl$mhc_b[idx],
                     mhc_cls = ref_tbl$mhc_cls[idx],
                     reference = ref_tbl$reference[idx],
                     method = ref_tbl$method[idx],
                     match_chain = chain,
                     match_type = paste0("fuzzy", max_dist))
    })
    fuzzy_core <- dplyr::bind_rows(fuzzy_list)
    if (!is.null(fuzzy_core) && nrow(fuzzy_core) > 0) {
      # map CDR3 -> barcodes (attach barcode AFTER fuzzy match)
      map <- df_clones %>% dplyr::select(barcode, !!rlang::sym(cdr3_col))
      fuzzy <- dplyr::inner_join(map, fuzzy_core, by = cdr3_col)
    }
  }
  
  dplyr::bind_rows(exact, fuzzy)
}


# 3) Run α and β matching
# refs already made as ref_trb / ref_tra
hits_trb <- match_chain(clones_df, ref_trb, chain = "TRB", max_dist = 1)
hits_tra <- match_chain(clones_df, ref_tra, chain = "TRA", max_dist = 1)


# 4) Combine & summarize
covid_hits_both <- bind_rows(hits_trb, hits_tra) %>% distinct()

# 5) “Either-chain” flag and “paired-confirmed” flags
cov_any_by_cell <- covid_hits_both %>% group_by(barcode) %>% summarise(any_hit = any(!is.na(match_chain)))
cov_paired_same_epitope <- covid_hits_both %>%
  group_by(barcode) %>%
  summarise(paired_same_epi = any(match_chain == "TRB") & any(match_chain == "TRA") &
              any(outer(epitope[match_chain=="TRB"], epitope[match_chain=="TRA"], FUN = "=="), na.rm=TRUE))
cov_paired_any_epitope <- covid_hits_both %>%
  group_by(barcode) %>%
  summarise(paired_any = any(match_chain == "TRB") & any(match_chain == "TRA"))

cov_flags <- cov_any_by_cell %>%
  left_join(cov_paired_same_epitope, by = "barcode") %>%
  left_join(cov_paired_any_epitope,   by = "barcode") %>%
  mutate(across(everything(), ~replace_na(., FALSE)))

# 6) Write outputs
if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(covid_hits_both, file.path(tab_dir, "covid_hits_TRA_TRB_cd4.tsv"))


## ====== 6) Export clonotype tables and top clones ======
# Per-clone abundance
cl_tab <- clones_df |>
  group_by(clone_id, TRB_cdr3, TRA_cdr3) |>
  summarise(n_cells = n(), .groups = "drop") |>
  arrange(desc(n_cells))
write_tsv(cl_tab, file.path(tab_dir, "clonotypes_cd4.tsv"))

topN <- head(cl_tab, 20)
p_top <- ggplot(topN, aes(x=reorder(TRB_cdr3, n_cells), y=n_cells)) +
  geom_col() + coord_flip() + xlab("TRB CDR3 (top20)") + ylab("# cells") +
  ggtitle("Top CD4 clonotypes")
sp(p_top, "top20_clones.png")

## ====== 7) Nice-to-have: V-J chord diagram for TRB (sample-level) ======
vj_mat <- contigs |>
  filter(grepl("TRB|Beta", locus, ignore.case = TRUE)) |>
  mutate(V = v_gene, J = j_gene) |>
  filter(!is.na(V), !is.na(J)) |>
  count(V, J) |>
  tidyr::pivot_wider(names_from = J, values_from = n, values_fill = 0) |>
  as.data.frame()

if(nrow(vj_mat)>0){
  rownames(vj_mat) <- vj_mat$V; vj_mat$V <- NULL
  m <- as.matrix(vj_mat)
  png(file.path(fig_dir, "vj_chord_TRB.png"), width=1400, height=1200, res=200)
  circos.clear()
  chordDiagram(m, transparency = 0.6)
  title("TRB V-J usage chord diagram (CD4)")
  dev.off()
}

## ====== 8) Save Seurat object with all metadata ======
saveRDS(seu, file.path(out_dir, "cd4_tcr_covid_seurat.rds"))

## ====== 9) (Optional) Minimal HTML report ======
# Write a tiny report index so you can click through outputs
report_lines <- c(
  "# CD4 TCR/COVID analysis report",
  "",
  "## Figures",
  "- UMAP clusters: figs/umap_clusters.png",
  "- CD4 state: figs/umap_cd4_state.png",
  "- Clone size UMAP: figs/umap_clone_size.png",
  "- Clonal homeostasis: figs/clonal_homeostasis_cd4.png",
  "- TRB CDR3 length: figs/cdr3_len_TRB.png",
  "- Top20 clones: figs/top20_clones.png",
  if(file.exists(file.path(fig_dir,"vj_usage_TRB_V.png"))) "- TRB V usage: figs/vj_usage_TRB_V.png" else NULL,
  if(file.exists(file.path(fig_dir,"vj_usage_TRB_J.png"))) "- TRB J usage: figs/vj_usage_TRB_J.png" else NULL,
  if(file.exists(file.path(fig_dir,"covid_hits_bar.png"))) "- COVID hits bar: figs/covid_hits_bar.png" else NULL,
  if(file.exists(file.path(fig_dir,"umap_covid_like.png"))) "- COVID-like UMAP: figs/umap_covid_like.png" else NULL,
  if(file.exists(file.path(fig_dir,"vj_chord_TRB.png"))) "- VJ chord (TRB): figs/vj_chord_TRB.png" else NULL,
  "",
  "## Tables",
  "- TCR QC: tables/tcr_qc_summary.tsv",
  "- Diversity: tables/diversity_cd4.tsv",
  "- Clonotypes: tables/clonotypes_cd4.tsv",
  if(file.exists(file.path(tab_dir,"covid_trb_hits_cd4.tsv"))) "- COVID matches: tables/covid_trb_hits_cd4.tsv" else NULL,
  if(file.exists(file.path(tab_dir,"covid_like_by_state.tsv"))) "- COVID-like × state: tables/covid_like_by_state.tsv" else NULL,
  if(file.exists(file.path(tab_dir,"covid_like_state_enrichment.tsv"))) "- State enrichment (Fisher): tables/covid_like_state_enrichment.tsv" else NULL
)
writeLines(report_lines[!sapply(report_lines, is.null)], con = file.path(out_dir, "REPORT.md"))
message("Done. See ", out_dir, "/REPORT.md and figs/tables subfolders.")
