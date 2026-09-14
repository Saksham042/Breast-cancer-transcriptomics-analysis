## ===================================================================
## Breast Cancer GEO Dataset — DEG, Pathway Enrichment
## & PPI Network Analysis
## Dataset: GSE42568 (GEO, Affymetrix HG-U133 Plus 2.0)
##          104 tumor / 17 normal breast tissue
## ===================================================================
##
## NOTE ON REPRODUCIBILITY:
## The differential expression step (GEO2R group assignment + DEG
## export) was performed manually through NCBI's web-based GEO2R tool.
## The steps from DEG list onward are exactly what was run in this
## project.
##
## NOTE ON CUTOFF: unlike Projects 1-2 (|logFC| > 1), this project
## uses |logFC| > 3, because the standard cutoff returned 4,156 genes
## -- too many for interpretable enrichment/network analysis, and
## reflective of the substantially larger effect sizes typical of
## tumor-vs-normal-tissue comparisons. See report Section 2.2 for the
## full justification and the cutoff-comparison table.
## ===================================================================


## ---- Step 0: Install required packages (run once) -----------------
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db",
                        "enrichplot", "DOSE", "pathview"))


## ---- Step 1: Load & filter DEG list (from GEO2R) -------------------
## GEO2R export columns for this platform: ID, adj.P.Val, P.Value, t,
## B, logFC, Gene.symbol, Gene.title
## Cutoff applied: adj.P.Val < 0.05 AND |logFC| > 3

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

deg <- read.csv("breast_cancer_deg_list_full.csv")
head(deg)
str(deg)
nrow(deg)               # 412
table(deg$Direction)    # 117 Upregulated, 295 Downregulated

## Cutoff sensitivity check (optional, for reference):
## |logFC| > 1   -> 4156 genes (too many)
## |logFC| > 1.5 -> 2157 genes
## |logFC| > 2   -> 1176 genes
## |logFC| > 2.5 -> 694 genes
## |logFC| > 3   -> 412 genes  <- used in this project
## |logFC| > 3.5 -> 246 genes
## |logFC| > 4   -> 146 genes


## ---- Step 2: Convert gene symbols to Entrez IDs --------------------
gene_list <- deg$Gene.symbol
gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
length(gene_list)   # 412

gene_conversion <- bitr(gene_list,
                         fromType = "SYMBOL",
                         toType   = "ENTREZID",
                         OrgDb    = org.Hs.eg.db)

nrow(gene_conversion)   # 371/412 mapped successfully (90%)
entrez_genes <- gene_conversion$ENTREZID


## ---- Step 3: KEGG pathway enrichment --------------------------------
kegg_result <- enrichKEGG(gene = entrez_genes,
                           organism = "hsa",
                           pvalueCutoff = 0.05)

nrow(as.data.frame(kegg_result))   # 9 significant pathways
head(as.data.frame(kegg_result))

## Check which genes are driving these KEGG hits
kegg_genes_entrez <- unique(unlist(strsplit(as.data.frame(kegg_result)$geneID, "/")))
kegg_gene_symbols <- bitr(kegg_genes_entrez,
                           fromType = "ENTREZID",
                           toType = "SYMBOL",
                           OrgDb = org.Hs.eg.db)
kegg_gene_symbols
## -> dominated by adipocyte genes (PLIN1/4, FABP4/5, ADIPOQ, LEP, LIPE,
##    CIDEC, CIDEA) plus genuine cancer genes (ERBB2, ERBB3, IGF1, PRLR,
##    GHR, TWIST2) and ECM genes (COL11A1, ITGA7, SDC1, GPC3)
## See report Section 3.2 for full interpretation of the adipose
## tissue-composition confound.

dir.create("figures", showWarnings = FALSE)
png("figures/kegg_barplot.png", width = 1000, height = 800)
barplot(kegg_result, showCategory = 9)
dev.off()


## ---- Step 4: GO Biological Process enrichment -----------------------
go_result <- enrichGO(gene = entrez_genes,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       readable = TRUE)

nrow(as.data.frame(go_result))   # 480 significant terms (largest of all 3 projects)

go_simplified <- simplify(go_result, cutoff = 0.7,
                           by = "p.adjust", select_fun = min)
nrow(as.data.frame(go_simplified))   # 179 non-redundant terms

png("figures/go_dotplot.png", width = 1000, height = 900)
dotplot(go_simplified, showCategory = 20)
dev.off()

png("figures/go_barplot.png", width = 1000, height = 900)
barplot(go_simplified, showCategory = 20)
dev.off()


## ---- Step 5: Pathway-level visualization (fold-change overlay) -----
library(pathview)

fc_values <- deg$logFC
names(fc_values) <- gene_conversion$ENTREZID[match(deg$Gene.symbol,
                                                     gene_conversion$SYMBOL)]
fc_values <- fc_values[!is.na(names(fc_values))]

## Regulation of lipolysis in adipocytes - top KEGG hit by p-value
pathview(gene.data = fc_values,
          pathway.id = "hsa04923",
          species = "hsa")
## Produces: hsa04923.pathview.png in the working directory
## Result: every gene in the pathway is downregulated (green) -
## interpreted as adipocyte tissue-composition loss, not active
## suppression. See report Section 3.4.


## ---- Step 6: Export gene list for STRING ---------------------------
## Take deg$Gene.symbol (412 genes) to string-db.org:
##   STRING -> "Multiple proteins" -> paste gene list ->
##   organism: Homo sapiens -> Search
## Export interactions as TSV ("string_interactions.tsv"),
## included in figures/ for reference. 371 genes mapped to STRING IDs;
## network returned 330 connected/annotated nodes, 1425 edges.


## ---- Step 7: Hub gene identification (degree centrality) ------------
## Network was imported into Cytoscape from string_interactions.tsv,
## analyzed via Tools -> NetworkAnalyzer -> Analyze Network
## (undirected), styled by degree (node size + color), Prefuse Force
## Directed layout weighted by combined_score.
##
## Equivalent degree calculation done in R for reference:

interactions <- read.delim("figures/string_interactions.tsv",
                            stringsAsFactors = FALSE)

degree_table <- c(interactions$node1, interactions$node2) |>
  table() |>
  sort(decreasing = TRUE)

head(degree_table, 10)
## Top hubs in this project: CDH1 (61), PPARG (54), LEP (44),
## ERBB2 (42), IGF1 (40), CD36 (39) - CDH1 (E-cadherin) is a
## genuinely significant breast cancer tumor-suppressor gene.

## ===================================================================
## End of pipeline. See Breast_Cancer_DEG_Report.pdf for full
## write-up, discussion, and biological interpretation, including the
## adipose-tissue-confound discussion central to this project.
## ===================================================================
