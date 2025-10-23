## ========== 08_cluster_DA_fisher.R ==========
## Cluster-level DA (COVID-like vs Other) with Fisher tests

setwd("/home/isadeghi/projects/covid_kids")   # <- adjust if needed
source("scripts/00_utils.R")

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(readr); library(ggplot2); library(tibble); library(forcats)
})

seu <- readRDS(file.path(out_dir, "merged_seurat_final.rds"))
stopifnot(all(c("covid_like","celltype","seurat_clusters") %in% colnames(seu@meta.data)))

# Define groups from covid_like
covid_like <- seu@meta.data$covid_like
group <- ifelse(covid_like %in% c("low","medium"), "COVID-like", "Other")  # <- CHOICE
lineage <- seu@meta.data$celltype

md <- seu@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  mutate(group   = factor(group, levels = c("Other","COVID-like")),
         lineage = factor(lineage, levels = c("CD4","CD8")),
         cluster = factor(seurat_clusters))

fisher_one <- function(df){
  tab <- table(df$cluster, df$group)  # rows = clusters, cols = groups
  rows <- rownames(tab)
  res <- lapply(seq_len(nrow(tab)), function(i){
    r <- tab[i,]
    mat <- matrix(c(r["COVID-like"], r["Other"],
                    sum(tab[,"COVID-like"]) - r["COVID-like"],
                    sum(tab[,"Other"])     - r["Other"]),
                  nrow=2, byrow=TRUE,
                  dimnames = list(c("in_cluster","out_cluster"),
                                  c("COVID-like","Other")))
    ft <- suppressWarnings(fisher.test(mat))
    tibble(cluster = rows[i],
           n_covid = as.integer(r["COVID-like"]),
           n_other = as.integer(r["Other"]),
           or = unname(ft$estimate),
           ci_lo = ft$conf.int[1], ci_hi = ft$conf.int[2],
           p = ft$p.value)
  }) %>% bind_rows() %>%
    mutate(q = p.adjust(p, "BH"),
           log2OR = log2(or))
  res
}

## A) All T cells combined
res_all <- fisher_one(md)
write_tsv(res_all, file.path(tab_dir, "DA_fisher_clusters_ALL.tsv"))

## B) By lineage (recommended)
res_lin <- md %>% group_by(lineage) %>% group_modify(~ fisher_one(.x)) %>% ungroup()
write_tsv(res_lin, file.path(tab_dir, "DA_fisher_clusters_byLineage.tsv"))

## Quick effects plots
for (lin in levels(md$lineage)) {
  d <- filter(res_lin, lineage == lin) %>% arrange(desc(log2OR))
  p <- ggplot(d, aes(x = fct_reorder(cluster, log2OR), y = log2OR)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_point(aes(size = n_covid + n_other, color = q < 0.05)) +
    scale_color_manual(values = c(`TRUE`="tomato", `FALSE`="grey40"), name="FDR < 0.05") +
    scale_size_continuous(name = "Cells in cluster") +
    coord_flip() +
    labs(title = paste0("Differential abundance (", lin, "): COVID-like vs Other"),
         x = "Seurat cluster", y = "log2(OR)") +
    theme_bw()
  ggsave(plot = p, file.path(fig_dir, sprintf("DA_fisher_by_cluster_%s.png", lin)), width = 7,height =  6)
}
