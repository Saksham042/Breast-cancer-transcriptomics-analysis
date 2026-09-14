# Breast Cancer Gene Expression Analysis (GSE42568)

Differential gene expression, pathway enrichment, and protein-protein
interaction network analysis of a public breast cancer tumor-vs-normal
transcriptomics dataset.

## Dataset
- **Source:** NCBI GEO, accession [GSE42568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568)
- **Platform:** Affymetrix Human Genome U133 Plus 2.0 Array
- **Samples:** 104 breast tumor, 17 normal breast tissue

## Workflow
1. **Differential expression** — GEO2R (limma). A stricter cutoff (adj. p-value < 0.05, **|logFC| > 3**, vs. >1 in the companion projects below) was used because the standard threshold returned 4,156 genes — too many for interpretable analysis, reflecting the larger effect sizes typical of tumor-vs-normal comparisons. Result: **412 unique DEGs** (117 up, 295 down).
2. **Functional enrichment** — `clusterProfiler`: KEGG (9 pathways) and GO Biological Process (480 terms, simplified to 179)
3. **Pathway visualization** — `pathview`, logFC mapped onto the KEGG Regulation of Lipolysis in Adipocytes pathway (hsa04923)
4. **PPI network** — STRING-db → Cytoscape, node degree analysis to identify hub genes (330 nodes, 1,425 edges — the densest of the three projects)

## Key finding — and an important methodological lesson
The dominant enriched signal (both KEGG and GO) is an **adipose-tissue / lipid-metabolism signature** (PPARG, LEP, CD36, FABP4, ADIPOQ, LIPE, PLIN1/4, CIDEC/CIDEA). This most likely reflects **reduced adipocyte content in tumor tissue** relative to normal breast tissue (which is predominantly fat) — a **tissue-composition confound**, not necessarily tumor-cell-intrinsic regulation. A secondary, genuinely tumor-relevant signal was also identified: **CDH1** (E-cadherin, the network's dominant hub, degree = 61 — a well-known breast cancer tumor suppressor) alongside **ERBB2 (HER2)**, **IGF1**, and **EZH2**, all established breast cancer genes. Distinguishing these two signals is the central analytical takeaway of this project — full discussion in the report.

## Repo contents

| File | Description |
|---|---|
| `Breast_Cancer_DEG_Report.pdf` | Full write-up: methods, results, figures, discussion, limitations |
| `Breast_Cancer_DEG_Report.docx` | Same report, editable Word format |
| `analysis_pipeline.R` | Full R pipeline: gene ID conversion → KEGG/GO enrichment → pathview → hub gene calc, incl. cutoff sensitivity check |
| `breast_cancer_deg_list_full.csv` | Full DEG table — 412 genes with logFC, adj.P.Val, direction, etc. |
| `figures/kegg_barplot.png` | All 9 significant KEGG pathways |
| `figures/go_dotplot.png` / `go_barplot.png` | Top 20 non-redundant enriched GO terms |
| `figures/pathway_diagram_hsa04923.png` | KEGG lipolysis-regulation pathway with logFC overlay |
| `figures/ppi_network_cytoscape.png` | STRING/Cytoscape interaction network, styled by hub degree |
| `figures/string_interactions.tsv` | Raw STRING interaction edge list |

## Note on cutoff choice
Unlike the companion projects below (which use the field-standard |logFC| > 1), this project uses **|logFC| > 3**. This was a deliberate choice, tested against multiple thresholds (see `analysis_pipeline.R`), to keep the gene list analytically tractable given the much larger transcriptional differences typical of tumor-vs-normal-tissue comparisons. One consequence, discussed in the report: this stricter cutoff selectively favored the most dramatic tissue-composition-driven genes (adipose loss) over more moderate tumor-intrinsic signaling changes (e.g., cell-cycle genes), which is noted explicitly as a limitation rather than overlooked.

## Tools used
R (clusterProfiler, org.Hs.eg.db, enrichplot, pathview) · GEO2R · STRING-db · Cytoscape

---
*Companion projects, same pipeline applied to different diseases:*
- *[asthma-transcriptomics-analysis](https://github.com/Saksham042/asthma-transcriptomics-analysis)*
- *[Alzheimers-transcriptomics-analysis](https://github.com/Saksham042/Alzheimers-transcriptomics-analysis)*
