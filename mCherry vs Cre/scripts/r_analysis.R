library(openxlsx)
library(DESeq2)
library(ggplot2)
library(pheatmap)

# Load raw counts
counts <- read.xlsx("data/GSE242485_Ube2iF_F_mCherry_vs_Cre_RNAseq_Counts.xlsx")
colnames(counts)[1] <- "gene"

# Remove genes with non-standard numeric IDs (44986, 44987) instead of real gene symbols
counts <- counts[!counts$gene %in% c("44986", "44987"), ]


# Reshape into DESeq2's expected format: gene names as row names, whole-number counts
count_data <- as.data.frame(counts)
rownames(count_data) <- count_data$gene
count_data$gene <- NULL
count_data <- round(count_data)

# Build sample info table (which columns are Cre vs mCherry)
condition <- ifelse(grepl("Ad-Cre", colnames(count_data)), "Cre", "mCherry")
sample_info <- data.frame(condition = condition, row.names = colnames(count_data))
sample_info$condition <- factor(sample_info$condition, levels = c("mCherry", "Cre"))

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
summary(res)

# PCA plot
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")
ggsave("figures/pca_plot.png", width = 6, height = 5)

# Volcano plot
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey70", "firebrick"), na.value = "grey70") +
  theme_minimal() +
  labs(title = "SUMOylation knockout vs control",
       x = "log2 fold change",
       y = "-log10 adjusted p-value")
ggsave("figures/volcano_plot.png", width = 6, height = 5)

# Heatmap of top 30 genes
top_genes <- head(order(res$padj), 30)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)

short_names <- paste0(sample_info$condition, "_", gsub(".*Replicate ", "", colnames(count_data)))
colnames(mat) <- short_names

annotation_col <- data.frame(condition = sample_info$condition)
rownames(annotation_col) <- short_names

pheatmap(mat, annotation_col = annotation_col, show_rownames = TRUE,
         filename = "figures/heatmap.png")

# Top hits table
top_hits <- as.data.frame(res[order(res$padj), ])
top_hits$gene <- rownames(top_hits)
top_hits <- top_hits[, c("gene", "log2FoldChange", "padj")]
head(top_hits, 15)
