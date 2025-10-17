# DESeq2 pipeline for your featureCounts + metadata
# Assumes files "count2.txt" and "metadata.csv" are in your working directory.

# 0. Packages (install if needed)
if(!requireNamespace("DESeq2", quietly=TRUE)) install.packages("BiocManager") & BiocManager::install("DESeq2")
if(!requireNamespace("apeglm", quietly=TRUE)) BiocManager::install("apeglm")   # for LFC shrinkage (optional)
if(!requireNamespace("pheatmap", quietly=TRUE)) install.packages("pheatmap")
if(!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")

library(DESeq2)
library(pheatmap)
library(ggplot2)

# 1. Read files
fc <- read.table("count2.txt", header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.csv("metadata.csv", stringsAsFactors = FALSE)

# 2. Extract count columns and set rownames
# featureCounts standard: first 6 columns are annotation, counts start at column 7
countData <- fc[, 7:ncol(fc)]
rownames(countData) <- fc$Geneid

# 3. Simplify sample names to a shorter form that will be used in metadata too
# This regex removes everything up to last slash and drops '.fastq...' suffixes.
clean_names <- function(x) gsub("\\.fastq.*|.*\\/", "", x)
colnames(countData) <- clean_names(colnames(countData))
meta$Sample <- clean_names(meta$Sample)

# 4. Ensure metadata rows match count columns and are in the same order
# Show diagnostic info:
cat("count columns:\n"); print(colnames(countData))
cat("metadata samples:\n"); print(meta$Sample)

# Reorder metadata to match count columns. If metadata missing samples, show message.
missing_in_meta <- setdiff(colnames(countData), meta$Sample)
missing_in_counts <- setdiff(meta$Sample, colnames(countData))
if(length(missing_in_meta) > 0) {
  stop("ERROR: these samples are present in counts but missing in metadata:\n", paste(missing_in_meta, collapse=", "))
}
if(length(missing_in_counts) > 0) {
  warning("Warning: these samples are present in metadata but not in counts:\n", paste(missing_in_counts, collapse=", "))
}

meta <- meta[match(colnames(countData), meta$Sample), , drop = FALSE]
rownames(meta) <- meta$Sample

# 5. Convert count matrix to numeric and ensure no non-numeric entries
countData[] <- lapply(countData, function(x) as.numeric(as.character(x)))
# If NAs were introduced by coercion, replace with 0 (rare) but report count
na_count <- sum(is.na(as.matrix(countData)))
if(na_count > 0) {
  warning(sprintf("Found %d NA(s) in count matrix after coercion; replacing with 0.", na_count))
  countData[is.na(countData)] <- 0
}

# Ensure non-negative integers (round if decimals present)
if(any(as.matrix(countData) < 0)) stop("Negative counts found - remove/inspect them before proceeding.")
if(!all(as.matrix(countData) == round(as.matrix(countData)))) {
  warning("Some counts are non-integers; rounding them to integers.")
  countData <- round(countData)
}

# 6. Basic filtering (optional but recommended): remove genes with all zero counts
keep_genes <- rowSums(countData) > 0
cat(sprintf("Filtering: keeping %d genes (out of %d) with sum>0\n", sum(keep_genes), nrow(countData)))
countData <- countData[keep_genes, ]

# 7. Prepare design: ensure your condition column is a factor
# This example uses the column named "Infection" in your metadata (as shown earlier).
# If your column is named differently, change "Infection" below.
if(!"Infection" %in% colnames(meta)) stop("Metadata does not have a column named 'Infection'. Edit script to use the correct column.")
meta$Infection <- factor(meta$Infection)
# Optionally set the reference level:
# meta$Infection <- relevel(meta$Infection, ref = "Control")

# 8. Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = as.matrix(countData),
                              colData = meta,
                              design = ~ Infection)

# 9. Pre-filter rows with low counts (prevents zero-inflation; optional)
dds <- dds[rowSums(counts(dds)) >= 10, ]   # keep genes with at least 10 reads across all samples

# 10. Run DESeq
dds <- DESeq(dds)
final_res <- results(dds) #final DESeq results
head(final_res) #view final results

# 11. let's look at our differentially expressed genes
plot(x = final_res$log2FoldChange, 
     y = -log10(final_res$padj),
     cex = 0.25,
     pch = 19, 
     col = 'grey',
     ylim = c(0,20),
     ylab = 'Adjusted P-Value',
     xlab = 'Log2 FC')

abline(v = c(-2, 2), h = -log10(0.05), lwd = 0.5, lty = 2)

#where are the upregulated genes
upregulated <- subset(final_res, padj < 0.05 & log2FoldChange > 2)
points(upregulated$log2FoldChange,
       y = -log10(upregulated$padj), 
       cex = 0.35,
       pch = 19,
       col = 'salmon')

# downregulated genes
downregulated <- subset(final_res, padj < 0.05 & log2FoldChange < -2)
points(downregulated$log2FoldChange,
       y = -log10(downregulated$padj), 
       cex = 0.35,
       pch = 19,
       col = 'lightblue')

mtext('A simple volcano')

#we can merge the two to do a clean and less memory efficient heatmap
degs <- rbind(countData[rownames(upregulated),], 
              countData[rownames(downregulated),])
pheatmap(degs, 
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         scale = 'row',
         show_colnames = T)

#what are the genes that are upregulated
rownames(upregulated)
rownames(downregulated)

#exporting the files
write.csv(upregulated, 'upregulated.csv')
write.csv(downregulated, 'downregulated.csv')
write.csv(countData, 'raw_counts.csv')


# 12. Variance stabilizing transform (VST) for QC / PCA
vst <- vst(dds, blind = TRUE)
library(ggplot2)

# 13. PCA plot (first two PCs)
pcaData <- plotPCA(vst, intgroup = "Infection", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = Infection, label = name)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of samples (VST)") +
  theme_minimal()

# 14. Sample distance heatmap
sampleDists <- dist(t(assay(vst)))
sampleDistMat <- as.matrix(sampleDists)
rownames(sampleDistMat) <- colnames(vst)
colnames(sampleDistMat) <- colnames(vst)
pheatmap(sampleDistMat, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists, main = "Sample distances (VST)")

# 15. Get results for the primary contrast (DESeq2 default uses the last level of factor as numerator)
# If you want a specific contrast, use results(dds, contrast = c("Infection","level1","level2"))
res <- results(dds)       # default comparison
summary(res)

# 16. Log2 fold change shrinkage (use apeglm if available)
if(requireNamespace("apeglm", quietly = TRUE)) {
  resLFC <- lfcShrink(dds, coef=2, type="apeglm")   # coef number may differ; check resultsNames(dds)
} else {
  message("apeglm not available — using 'normal' shrinkage")
  resLFC <- lfcShrink(dds, coef=2, type="normal")
}

# NOTE: Check resultsNames(dds) and change coef index accordingly if needed:
cat("Available coefficients:\n"); print(resultsNames(dds))

# 17. Basic MA plot and Volcano-like plot
plotMA(resLFC, ylim = c(-5, 5), main = "MA plot (LFC-shrunk)")
# Volcano (simple)
resDF <- as.data.frame(resLFC)
resDF$gene <- rownames(resDF)
resDF$padj[is.na(resDF$padj)] <- 1
ggplot(resDF, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  xlab("log2 fold change (shrunken)") + ylab("-log10 adjusted p-value") +
  ggtitle("Volcano (shrunken LFC)")

# 18. Export results to CSV
write.csv(as.data.frame(resLFC), file = "deseq2_results_lfc_shrunk.csv", row.names = TRUE)
write.csv(as.data.frame(res), file = "deseq2_results_raw.csv", row.names = TRUE)
cat("Wrote: deseq2_results_lfc_shrunk.csv and deseq2_results_raw.csv\n")

# 19. Save normalized counts (VST or rlog)
normCounts <- assay(vst)
write.csv(normCounts, file = "vst_normalized_counts.csv", row.names = TRUE)
cat("Wrote: vst_normalized_counts.csv\n")
