# Tpex RNA-seq Analysis

**Author:** Dr. Iman Sadeghi Dehcheshmeh  
**Repository:** [Tpex_RNAseq](https://github.com/isadeghi87/Tpex_RNAseq)  
**Last Updated:** October 2025

---

## 🧬 Project Overview

This repository contains a complete analysis workflow for **RNA-seq profiling of Tpex cells** under **AZM (Azithromycin)** vs **DMSO** treatment.  
The project integrates **differential expression**, **weighted gene co-expression network analysis (WGCNA)**, and **multi-database functional enrichment** (GO, Reactome, KEGG, WikiPathways, and MSigDB Hallmark).

---

## 🧠 Biological Goal

To identify transcriptional modules and gene networks modulated by **Azithromycin** treatment, with emphasis on:
- Immune regulation and T-cell exhaustion (Tpex-related genes)
- Translational and metabolic reprogramming
- Network hubs and functional pathways driving phenotypic response

---

## 📁 Repository Structure

```
Tpex_RNAseq/
├── data/
├── results/
│   ├── wgcna_allgenes/
│   ├── enrich_gprofiler/
│   ├── deg_analysis/
└── codes/
```

---

## 🧩 Key Analyses

### Preprocessing
- Top 8 000 genes by MAD
- Filtered for non-zero variance and quality

### Differential Expression
- Model: `~ Condition + Donor`
- Cutoffs: `|log2FC| > 0.5`, `p < 0.05`

### Weighted Gene Co-Expression Network
- Power selection via scale-free topology
- Module–trait correlation (AZM vs DMSO, Donor)
- HubScore = z(|kMEₒwn| + |GS| + kWithin)

### Functional Enrichment
- g:Profiler: GO, Reactome, KEGG, WikiPathways
- Hallmark: MSigDB H collection
- Integrated Up/Down dotplots (PDF)

### Network Visualization
- Cytoscape exports for per-module networks
- Global hub graph (`|bicor| ≥ 0.8`)

---

## ⚙️ Dependencies

| Package | Version |
|----------|----------|
| WGCNA | ≥ 1.73 |
| limma, DESeq2 | ≥ 3.56 |
| clusterProfiler, org.Hs.eg.db | ≥ 4.8 |
| gprofiler2 | ≥ 0.2 |
| msigdbr | ≥ 7.5 |
| igraph | ≥ 2.0 |
| ggplot2, ggrepel, pheatmap | ≥ 3.4 |

---

## ▶️ How to Run

```bash
git clone https://github.com/isadeghi87/Tpex_RNAseq.git
cd Tpex_RNAseq/codes

Rscript st1_preprocessing.R
Rscript st2_differential_expression.R
Rscript st3_wgcna_network.R
Rscript st4_enrichment_analysis.R
Rscript st5_visualization.R
```

---

## 🧠 Summary

- **Blue module** : immune/lysosomal upregulation by AZM  
- **Green module** : translational downregulation  
- **Hub genes** : IL7R, CXCR4, TSC22D3, FOXP3  
- **Hallmarks** : IL2-STAT5 signaling, apoptosis, oxidative stress

---

## 🧾 Citation

> Sadeghi Dehcheshmeh I., *et al.* (2025).  
> **Integrated transcriptomic network analysis reveals Azithromycin-induced transcriptional reprogramming in Tpex cells.**

---

## 🛠️ License

MIT License — see [LICENSE](LICENSE)

---

## 📬 Contact

**Dr. Iman Sadeghi Dehcheshmeh**  
Senior Bioinformatician | DKFZ Heidelberg  
📧 iman.sadeghi87@gmail.com  
🌐 [LinkedIn](https://www.linkedin.com/in/iman-sadeghi-dehcheshmeh/)
