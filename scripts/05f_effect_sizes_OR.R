## ========== 05f_effect_sizes_OR.R ==========
## Purpose: compute odds ratios and CIs for COVID-like enrichment by lineage and cluster.

# lineage
or_lineage <- lapply(split(seu@meta.data, seu$celltype), function(df){
  tab <- table(df$covid_like != "none")
  m <- matrix(c(tab["TRUE"], tab["FALSE"],
                sum(seu$covid_like!="none")-tab["TRUE"],
                sum(seu$covid_like=="none")-tab["FALSE"]), 2, byrow=TRUE)
  ft <- fisher.test(m)
  data.frame(level=unique(df$celltype), OR=unname(ft$estimate), CI_lo=ft$conf.int[1], CI_hi=ft$conf.int[2], p=ft$p.value)
}) |> dplyr::bind_rows()
readr::write_tsv(or_lineage, file.path(tab_dir, "covid_like_lineage_OR.tsv"))

# cluster
or_cluster <- lapply(split(seu@meta.data, seu$seurat_clusters), function(df){
  tab <- table(df$covid_like != "none")
  m <- matrix(c(tab["TRUE"], tab["FALSE"],
                sum(seu$covid_like!="none")-tab["TRUE"],
                sum(seu$covid_like=="none")-tab["FALSE"]), 2, byrow=TRUE)
  ft <- fisher.test(m)
  data.frame(cluster=unique(df$seurat_clusters), OR=unname(ft$estimate), CI_lo=ft$conf.int[1], CI_hi=ft$conf.int[2], p=ft$p.value)
}) |> dplyr::bind_rows() |> dplyr::arrange(dplyr::desc(OR))
readr::write_tsv(or_cluster, file.path(tab_dir, "covid_like_cluster_OR.tsv"))
