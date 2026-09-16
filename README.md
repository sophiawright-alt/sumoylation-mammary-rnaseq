# SUMOylation and Mammary Tumor-Initiating Cells: An RNA-seq Analysis

An RNA-seq differential expression analysis investigating how loss of SUMOylation affects gene expression in mouse mammary cells, using public data from a study on cancer stem cell maintenance.

## Background

SUMOylation is a regulatory system where cells attach a small protein tag (SUMO) onto other proteins to control their activity, location, or stability. The gene *Ube2i* is the single essential enzyme required for this process. Prior published work found that SUMOylation maintains "tumor-initiating cells" — stem-cell-like cancer cells capable of starting new tumors — in mammary tissue, and that blocking it selectively damages this population without harming normal stem cells.

This analysis compares:
- **mCherry** - control cells, SUMOylation intact
- **Cre** - cells with *Ube2i* knocked out via Cre recombinase, SUMOylation disabled

## Dataset

- **Source**: [GEO accession GSE242485](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242485)
- **Organism**: *Mus musculus*
- **Design**: 3 mCherry replicates vs. 3 Cre replicates, raw RNA-seq counts
- **Associated paper**: Li et al., *Developmental Cell* 2025, PMID 40315856

## Methods

1. Loaded raw counts and cleaned the data — renamed an unlabeled gene column, removed 4 rows with non-standard numeric placeholder IDs instead of real gene symbols, and rounded fractional counts to whole numbers as required by DESeq2.
2. Ran differential expression analysis with **DESeq2** (R/Bioconductor), testing all genes with sufficient read counts.
3. Used a significance threshold of adjusted p-value (padj) < 0.05.

## Results

- 16,255 of 24,417 genes had sufficient counts to test.
- **835 genes significantly upregulated**, **563 significantly downregulated** in Cre vs. mCherry.
- PCA showed condition alone explained 69% of all variance in the dataset, with samples cleanly separating by group.
- The volcano plot showed confident, significant hits in both directions.
- A heatmap of the top 30 genes clustered samples cleanly into two groups by condition.

![PCA plot](https://github.com/sophiawright-alt/ube2i-knockout-rnaseq/blob/main/mCherry%20vs%20Cre/figures/pca_plot.png?raw=true)
![Volcano plot](https://github.com/sophiawright-alt/ube2i-knockout-rnaseq/blob/main/mCherry%20vs%20Cre/figures/volcano_plot.png?raw=true)
![Heatmap](https://github.com/sophiawright-alt/ube2i-knockout-rnaseq/blob/main/mCherry%20vs%20Cre/figures/heatmap.png?raw=true)

## Key finding

The three most statistically significant genes overall — **Pkp1**, **Sbsn**, and **Ivl** — are all established markers of terminal epithelial differentiation, and all were strongly upregulated when SUMOylation was removed. This is consistent with the published mechanism
