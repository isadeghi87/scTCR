## ========== 01_load_merge.R ==========
source("/home/isadeghi/projects/covid_kids/scripts/00_utils.R")
suppressPackageStartupMessages(library(future))
plan("sequential")                               # avoid parallel futures for small/medium data
options(future.globals.maxSize = 8 * 1024^3)  

## --- Load CD4 ---
mat4 <- Read10X_h5(gex_cd4)
seu4 <- CreateSeuratObject(mat4, project="CD4", min.cells=3, min.features=200)
seu4[["percent.mt"]] <- PercentageFeatureSet(seu4, "^MT-")
seu4 <- subset(seu4, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
seu4 <- SCTransform(seu4, vst.flavor="v2", verbose=FALSE) |> RunPCA() |> RunUMAP(dims=1:30) |>
  FindNeighbors(dims=1:30) |> FindClusters(resolution=0.6)
seu4$orig.ident <- "CD4"; seu4$celltype <- "CD4"
sp(DimPlot(seu4, group.by="seurat_clusters", label=TRUE)+ggtitle("CD4 clusters"), "umap_CD4.png")

## --- Load CD8 ---
mat8 <- Read10X_h5(gex_cd8)
seu8 <- CreateSeuratObject(mat8, project="CD8", min.cells=3, min.features=200)
seu8[["percent.mt"]] <- PercentageFeatureSet(seu8, "^MT-")
seu8 <- subset(seu8, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
seu8 <- SCTransform(seu8, vst.flavor="v2", verbose=FALSE) |> RunPCA() |> RunUMAP(dims=1:30) |>
  FindNeighbors(dims=1:30) |> FindClusters(resolution=0.6)
seu8$orig.ident <- "CD8"; seu8$celltype <- "CD8"
sp(DimPlot(seu8, group.by="seurat_clusters", label=TRUE)+ggtitle("CD8 clusters"), "umap_CD8.png")

## --- Read VDJ contigs, keep productive & top-UMI per locus ---
prep_contigs <- function(vdj_csv){
  df <- data.table::fread(vdj_csv) |> as.data.frame()
  df <- df |>
    dplyr::filter(is_cell %in% c(TRUE,"true",NA), productive %in% c(TRUE,"True")) |>
    dplyr::mutate(locus = ifelse(!is.na(chain), chain, locus))
  df <- df |>
    dplyr::group_by(barcode, locus) |>
    dplyr::slice_max(order_by = as.numeric(umis %||% reads %||% 0), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  df
}
c4 <- prep_contigs(vdj_cd4)
c8 <- prep_contigs(vdj_cd8)

## --- Build TRA/TRB tables & attach to Seurat meta ---
pair_from_contigs <- function(contigs){
  tra <- contigs |> dplyr::filter(grepl("TRA|TCRA|Alpha", locus, ignore.case=TRUE)) |>
    dplyr::select(barcode, TRA_cdr3=cdr3, TRA_v=v_gene, TRA_j=j_gene)
  trb <- contigs |> dplyr::filter(grepl("TRB|TCRB|Beta",  locus, ignore.case=TRUE)) |>
    dplyr::select(barcode, TRB_cdr3=cdr3, TRB_v=v_gene, TRB_j=j_gene)
  pair <- dplyr::full_join(tra, trb, by="barcode")
  if (!any(grepl("-\\d+$", pair$barcode))) pair$barcode <- paste0(pair$barcode, "-1")
  pair$clone_def <- dplyr::case_when(
    !is.na(pair$TRA_cdr3) & !is.na(pair$TRB_cdr3) ~ paste0("TRA:",pair$TRA_cdr3,"|TRB:",pair$TRB_cdr3),
    is.na(pair$TRA_cdr3) & !is.na(pair$TRB_cdr3)  ~ paste0("TRB:",pair$TRB_cdr3),
    TRUE ~ NA_character_
  )
  pair$clone_id <- as.character(factor(pair$clone_def))
  rownames(pair) <- pair$barcode
  pair
}
p4 <- pair_from_contigs(c4); p8 <- pair_from_contigs(c8)

attach_pair <- function(seu, pair){
  hits <- intersect(colnames(seu), pair$barcode)
  idx  <- match(hits, pair$barcode)
  seu$TRA_cdr3 <- NA_character_; seu$TRB_cdr3 <- NA_character_; seu$clone_id <- NA_character_
  seu$TRA_cdr3[hits] <- as.character(pair$TRA_cdr3[idx])
  seu$TRB_cdr3[hits] <- as.character(pair$TRB_cdr3[idx])
  seu$clone_id[hits] <- as.character(pair$clone_id[idx])
  seu$clone_size <- ifelse(is.na(seu$clone_id), 0, ave(seu$clone_id, seu$clone_id, FUN=length))
  seu
}
seu4 <- attach_pair(seu4, p4)
seu8 <- attach_pair(seu8, p8)

## --- Merge CD4 + CD8 ---
# 1) Integration features + prep
features <- SelectIntegrationFeatures(object.list = list(seu4, seu8), nfeatures = 3000)
seu4 <- PrepSCTIntegration(object.list = list(seu4), anchor.features = features)[[1]]
seu8 <- PrepSCTIntegration(object.list = list(seu8), anchor.features = features)[[1]]

# 2) Find anchors & integrate (SCT method)
anchors <- FindIntegrationAnchors(object.list = list(seu4, seu8),
                                  normalization.method = "SCT",
                                  anchor.features = features)
seu <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
# 3) DR on integrated assay
DefaultAssay(seu) <- "integrated"
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30, k.param = 30, nn.method = "annoy", annoy.metric = "cosine")
seu <- RunUMAP(seu, dims = 1:30, umap.method = "uwot", n.neighbors = 30, min.dist = 0.3, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.6)

sp(DimPlot(seu, group.by="celltype")+ggtitle("Merged CD4+CD8"), "umap_merged_celltype.png")

## Save objects & tables for next steps
saveRDS(seu, file.path(out_dir, "merged_seurat.rds"))
saveRDS(list(CD4=c4, CD8=c8), file.path(out_dir, "contigs_raw.rds"))
readr::write_tsv(p4, file.path(tab_dir,"pair_CD4.tsv"))
readr::write_tsv(p8, file.path(tab_dir,"pair_CD8.tsv"))

