## ========== 05j_trbvj_delta_heatmap.R ==========
## Purpose: TRBV–TRBJ Δ fraction heatmap (covid - other), normalized per TRBV row.

ct <- vj_df %>% dplyr::count(group, TRB_v, TRB_j, name="n")
keep_pairs <- ct %>% dplyr::group_by(TRB_v, TRB_j) %>% dplyr::summarise(n_tot=sum(n), .groups="drop") %>% dplyr::filter(n_tot>=10) %>% dplyr::select(TRB_v, TRB_j)

cov <- ct %>% dplyr::filter(group=="covid") %>% dplyr::right_join(keep_pairs, by=c("TRB_v","TRB_j")) %>%
  tidyr::replace_na(list(n=0)) %>% dplyr::group_by(TRB_v) %>% dplyr::mutate(frac=n/sum(n)) %>% dplyr::ungroup() %>% dplyr::rename(frac_covid=frac)
oth <- ct %>% dplyr::filter(group=="other") %>% dplyr::right_join(keep_pairs, by=c("TRB_v","TRB_j")) %>%
  tidyr::replace_na(list(n=0)) %>% dplyr::group_by(TRB_v) %>% dplyr::mutate(frac=n/sum(n)) %>% dplyr::ungroup() %>% dplyr::rename(frac_other=frac)

Delta <- cov %>% dplyr::inner_join(oth, by=c("TRB_v","TRB_j")) %>% dplyr::transmute(TRB_v, TRB_j, delta=frac_covid-frac_other) %>%
  tidyr::pivot_wider(names_from=TRB_j, values_from=delta, values_fill=0)
mat <- as.matrix(Delta[,-1, drop=FALSE]); rownames(mat) <- Delta$TRB_v
if (any(duplicated(rownames(mat)))) { grp <- factor(rownames(mat)); mat <- rowsum(mat, grp) / as.vector(table(grp)) }

row_scale <- function(M){ mu <- rowMeans(M, na.rm=TRUE); sdv <- apply(M,1,sd,na.rm=TRUE); Z <- sweep(M,1,mu,"-"); Z <- sweep(Z,1,ifelse(sdv>0,sdv,1),"/"); Z[!is.finite(Z)] <- 0; Z }
mat_z <- row_scale(mat)
col_fun <- circlize::colorRamp2(c(min(mat_z), 0, max(mat_z)), c("#043915", "#FFF8E8", "#FF3F7F"))

ht <- ComplexHeatmap::Heatmap(mat_z, name="Δ frac", col=col_fun, cluster_rows=TRUE, cluster_columns=TRUE,
                              row_title="TRBV", column_title="TRBJ", column_title_side='bottom')
png(file.path(fig_dir, "TRBVJ_delta_fraction_covid_vs_other.png"), width=1400, height=1100, res=160); draw(ht); dev.off()
