# scTCR: COVID‑linked scTCR repertoire

This repository contains analysis code and figures for our COVID‑focused TCR/Transcriptome study in pediatric liver samples. It integrates **Seurat** RNA data with **scRepertoire** TCR annotations and adds a set of *Scirpy‑inspired* repertoire summaries (spectratypes, clonotype overlap, V/J skew, and clonotype‑level marker profiles).

> **Key outputs**: clonotype–epitope mapping against VDJdb, confidence‑tiered “COVID‑hit” labels per cell, clonotype expansion statistics, V/J enrichment and V–J chord diagrams, β‑CDR3 similarity networks, and differential expression for CD8 COVID‑hit clones.

---

## 🔧 Quick start

```r
# clone the repo and set working directory to the project root
renv::activate()           # optional but recommended
renv::restore()            # installs pinned packages

# project folders expected
# data/, figs/, tables/, rds/ (or adjust paths in 00_utils.R)

# run the main analysis
source("scripts/00_utils.R")
source("scripts/03_covid_mapping_reporting.R")  # map to VDJdb, assign covid_like tiers, baseline figs/tables
source("scripts/05_covid_focused.R")            # downstream COVID‑focused analyses (CD8‑centric)
```

If you prefer chunked modules, see `scripts/05_covid_*_*.R` variants (analysis is identical, split for easier debugging).

---

## 📁 Repository layout (what matters)

```
scripts/
  00_utils.R                         # paths (out_dir, fig_dir, tab_dir), gg saver `sp()`, common helpers
  03_covid_mapping_reporting.R       # VDJdb harmonization → hits; tiering: low/medium/high; summary figs/tables
  05_covid_focused.R                 # full COVID analyses (DE, V/J skew, chord diagrams, networks, spectratypes)
  05_covid_focus_01_labels.R         # (chunk) build covid_hit/covid_like + best epitope per barcode
  05_covid_focus_02_epitopes.R       # (chunk) epitope mix plots + UMAP by top epitopes
  05_covid_focus_03_expansion.R      # (chunk) clone size violin/box (log1p), lineage/cluster enrichments
  05_covid_focus_04_vj_skew.R        # (chunk) TRBV/TRBJ enrichment, Δfraction heatmap, V–J chord diagrams
  05_covid_focus_05_de_markers.R     # (chunk) DE for medium CD8 in dominant cluster + volcano
  05_covid_focus_06_topCT_markers.R  # (chunk) DotPlot panel across top COVID clonotypes; heatmap of z‑scores
  05_covid_focus_07_similarity.R     # (chunk) β‑CDR3 lev≤1 graph + components
  05_covid_focus_08_spectratype.R    # (chunk) TRA/TRB spectratypes + ECDF w/ stats
  helpers_vj_chord.R                 # `make_vj_chord()` and `build_trb_table()` utilities

figs/                                # all figures saved by scripts (png/pdf)
tables/                              # all tables saved by scripts (tsv)
rds/ or out/                         # serialized Seurat objects
```

> The repo contains additional utilities and notebooks; the core analysis hinges on the scripts above. Adjust paths inside `00_utils.R` as needed.

---

## 🧪 Data inputs & object expectations

- **Merged Seurat** object with TCR metadata: `rds/merged_seurat_final.rds`
  - `meta.data` columns required by scripts:  
    `orig.ident`, `celltype` (levels include `"CD8"` and optionally `"CD4"`),  
    `seurat_clusters`, `TRA_cdr3`, `TRB_cdr3`, `CTaa` (paired AA clonotype), `clone_id`.
- **scRepertoire pair tables** (barcode + V/J calls): `tables/pair_CD4.tsv`, `tables/pair_CD8.tsv`
- **VDJdb export** (tab‑delimited). The scripts harmonize to lowercase headers and select:  
  `gene`, `cdr3`, `epitope`, `mhc class`, `mhc a`, `mhc b`, `method`, `score`.

The mapper builds per‑cell matches with priority **TRB exact > TRA exact > TRB fuzzy (Levenshtein ≤1)**, split by lineage (MHCI→CD8, MHCII→CD4). We then derive confidence tiers:

- `low`: best single‑chain hit present (usually TRB exact).  
- `medium`: paired match (TRA **and** TRB) to *any* epitope.  
- `high`: paired TRA/TRB to the **same** epitope.

The chosen tier is controlled with `CONF_LEVEL` in the scripts.

---

## 📊 Main analyses & figures

Below is a compact index of what each plot/table represents and where it’s produced (file → script). You can reproduce everything by running the scripts in order.

### 1) Mapping & labeling
- **`tables/covid_hits_CD4_CD8.tsv`** — all per‑cell matches with chain/match‑type/epitope (**03**).  
- **`tables/covid_best_per_cell_CD4_CD8.tsv`** — single best epitope per cell (**03**).  
- **`tables/covid_paired_any_CD4_CD8.tsv`**, **`tables/covid_paired_same_epitope_CD4_CD8.tsv`** — pairing flags (**03**).  
- **`figs/umap_covid_like_confidence_CD4_CD8.png`** — UMAP colored by `covid_like` (**03**).

### 2) Expansion & distribution
- **Clone size violin/box**: `figs/clone_size_covid_like_colored.png` (**05/expansion**).  
  *Shows larger clones within COVID‑linked cells (log1p scale), Wilcoxon p on console.*  
- **Lineage/cluster enrichment**: `tables/covid_like_lineage_OR.tsv`, `covid_like_cluster_OR.tsv` (**05/expansion**).

### 3) Epitopes
- **Mix by tier & lineage**: `figs/covid_epitope_mix.png` (**05/epitopes**).  
- **UMAP by top epitopes**: `figs/covid_<tier>_umap_top_epitopes.png` (**05/epitopes**).

### 4) V/J skew & pair structure
- **TRBV/TRBJ enrichment barplots**: `figs/TRBV_enrichment_covid.png`, `TRBJ_enrichment_covid.png` (**05/vj_skew**).  
- **Δ fraction heatmap (COVID–Other)**: `figs/TRBVJ_delta_fraction_covid_vs_other.png` (**05/vj_skew**).  
- **Manual usage heatmaps**: `figs/vj_heatmap_TRB_manual.png`, `vj_usage_TRB_V_manual.png`, `vj_usage_TRB_J_manual.png`.  
- **V–J chord diagrams**:  
  `figs/vj_chord_TRB_all.png`, `vj_chord_TRB_CD8.png`, `vj_chord_TRB_CD8_covid.png`, `vj_chord_TRB_CD8_other.png` (**helpers_vj_chord.R**).

### 5) Clonotype marker biology
- **CD8 medium vs none (dominant cluster) volcano**: `figs/volcano_medium_CD8_cluster*.png` (**05/de_markers**).  
- **Canonical marker DotPlot** across top CTs: `figs/covid_<tier>_topCT_dotplot.png` (**05/topCT_markers**).  
- **Top‑CT z‑score heatmap**: `figs/topCT_marker_heatmap.png` (**05/topCT_markers**).

### 6) Repertoire structure
- **β‑CDR3 lev≤1 similarity network**: `tables/tcr_beta_lev1_communities.tsv` and PNG (**05/similarity**).  
- **Spectratype** (TRA/TRB, COVID vs other; lineage split):  
  `figs/scirpy_like_spectratype.png`, `figs/spectratype_trb_lineage_covid.png`.  
  ECDF panels with Wilcoxon p: `figs/ecdf_trb_len_by_lineage_covid_status.png` (**05/spectratype**).

---

## 🔬 Methods (concise)

- **Preprocessing**: Seurat `SCTransform(vst.flavor="v2")`, PCA, UMAP (dims 1–30), neighbors/clusters (`FindNeighbors/FindClusters`).  
- **TCR integration**: `scRepertoire::combineExpression()` matched barcodes (suffix harmonized), producing `TRA_cdr3`, `TRB_cdr3`, `CTaa`, V/J calls.  
- **VDJdb harmonization**: lower‑case headers; keep human records; split by MHC class (MHCI→CD8, MHCII→CD4); exact AA matching by chain, optional TRB fuzzy (Levenshtein ≤1, identical length). Priority: TRB exact > TRA exact > TRB fuzzy.  
- **Tiering**: “low/medium/high” per above; stored as `covid_like` and `covid_hit`.  
- **Statistics**: Fisher tests for enrichment (lineage/cluster), Wilcoxon for expansion, KS/Wilcoxon + Cliff’s delta for spectratypes. Multiple testing by BH FDR.  
- **V/J skew**: per‑feature 2×2 Fisher (COVID vs other) with ORs; Δ‑fraction heatmaps by row‑normalized TRBV; chord diagrams via `circlize::chordDiagram`.  
- **Differential expression**: `FindMarkers` (Wilcoxon; logFC≥0.25; min.pct≥0.1) within‑cluster, within‑lineage controls to reduce state confounding; volcano with labeled top genes.  
- **Visualization**: ggplot2, ComplexHeatmap, ggraph (FR layout for networks).

See inline comments in each script for the exact function calls and guards.

---

## 🧩 Reproducibility tips

- Ensure barcodes between Seurat and scRepertoire use the same suffix convention (we strip `"_1"`, `"_2"` before joining).  
- For *CD8‑centric* analysis, set `CONF_LEVEL="medium"` (paired) in `05_*` to prioritize higher specificity.  
- If V/J tables don’t draw chords (too sparse), lower `min_pair_count` in `make_vj_chord()`.

---

## 📦 Dependencies

R ≥ 4.2 with: `Seurat`, `scRepertoire`, `dplyr`, `tidyr`, `ggplot2`, `ggrepel`, `stringdist`, `FNN`,  
`ComplexHeatmap`, `circlize`, `ggraph`, `igraph`, `vegan`, `ggridges`, `scales`, `tibble` (and `renv` for lockfile).
Install via:

```r
install.packages(c("dplyr","tidyr","ggplot2","ggrepel","stringdist","FNN",
                   "Seurat","tibble","scales"))
BiocManager::install(c("ComplexHeatmap","circlize","igraph","ggraph","vegan"))
# scRepertoire from CRAN or GitHub as needed
```

---

## 🗂️ Outputs to expect (by default paths)

- `figs/` — all PNG/PDF figures listed above.  
- `tables/` — TSVs with hits, enrichment, DE markers, communities, etc.  
- `out/merged_seurat_final.rds` — updated object with `covid_hit`, `covid_like`, `covid_epitope` slots.

---

## ✍️ Citation / context

This repository accompanies our spatial single‑cell study of pediatric acute hepatitis potentially linked to prior SARS‑CoV‑2 exposure, with an emphasis on **CD8 T‑cell** repertoires. Please also see the related clinical/pathology paper for study background and specimen context.

---

## 🙋 FAQ

- **My `combineExpression()` says “<1% barcodes match.”**  
  Harmonize suffixes: in our scripts we add `base_bc = sub("_\\d+$","", barcode)` before joining.

- **DotPlot shows only “OtherCOVID.” Where are the CT groups?**  
  Make sure `topN` was built (COVID‑hit cells present) and identities were assigned. Check `table(Idents(seu))` for `CT:*` labels.

- **Chord diagram has single color.**  
  Pass a named `grid.col` vector (see `helpers_vj_chord.R`) and sufficient pair counts.

---

## 📜 License
MIT (see `LICENSE`).
