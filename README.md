# Transcriptomic Profiling of *Staphylococcus aureus* During Acute vs Chronic Phases of Periprosthetic Joint Infection (PJI)

**Author:** Ifeoma Nwafor  
**Project Type:** RNA-seq Differential Expression Analysis  
**Dataset:** PRJNA867318  

---

## 🩺 Background and Rationale

Periprosthetic joint infections (PJIs) are among the most devastating complications of orthopedic implants. They increase morbidity, prolong hospital stays, and often require costly revision surgeries. *Staphylococcus aureus*—particularly methicillin-resistant strains (MRSA)—is a leading cause of PJIs.

One hallmark of *S. aureus* pathogenesis is its ability to switch between distinct phenotypes during infection:

- **Acute phase:** Bacteria exhibit an aggressive planktonic growth mode, expressing virulence factors such as toxins, adhesins, and immune evasion genes.  
- **Chronic phase:** Bacteria transition into a biofilm-like state, downregulating overt virulence genes and upregulating stress response, metabolic rewiring, and persistence pathways.

This switch underlies the resilience of chronic PJIs to antibiotics and immune clearance.

---

## 🎯 Project Goals

1. **Preprocessing and Quality Control**
   - Perform read trimming, alignment to the *S. aureus* reference genome, and QC.
   - Generate count matrices for acute and chronic isolates.

2. **Differential Gene Expression (DEG) Analysis**
   - Compare gene expression between acute and chronic PJI isolates.
   - Identify key virulence, metabolic, and biofilm-associated genes.

3. **Functional Enrichment and Pathway Mapping**
   - Conduct GO/KEGG enrichment analyses to highlight biological processes linked to persistence and immune evasion.

4. **Visualization and Reporting**
   - Generate PCA plots, volcano plots, and heatmaps of differentially expressed genes.
   - Summarize molecular adaptations distinguishing acute and chronic infections.

---

## 🧫 Dataset Information

| Accession Number | Infection State |
|------------------|----------------|
| SRR20959676 | Chronic |
| SRR20959677 | Chronic |
| SRR20959678 | Chronic |
| SRR20959679 | Chronic |
| SRR20959680 | Acute |
| SRR20959681 | Acute |
| SRR20959682 | Acute |

Data were obtained from **BioProject [PRJNA867318](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA867318)** and processed using **featureCounts** to generate gene-level count matrices.

---

## ⚙️ Workflow Overview

### 1. Quality Control and Alignment
Performed using:
- **FastQC** → read quality
- **FastP** → adapter/low-quality trimming
- **STAR** → alignment to *S. aureus* reference genome
- **featureCounts** → generation of raw count matrix

### 2. Differential Gene Expression with DESeq2

```r
library(DESeq2)

# Load featureCounts data and metadata
fc <- read.table("count2.txt", header = TRUE, sep = "\t", check.names = FALSE)
meta <- read.csv("metadata.csv", stringsAsFactors = FALSE)

# Extract count matrix
countData <- fc[, 7:ncol(fc)]
rownames(countData) <- fc$Geneid

# Simplify names and align metadata
clean_names <- function(x) gsub("\\.fastq.*|.*\\/", "", x)
colnames(countData) <- clean_names(colnames(countData))
meta$Sample <- clean_names(meta$Sample)
meta <- meta[match(colnames(countData), meta$Sample), ]
rownames(meta) <- meta$Sample

# Prepare DESeq2 dataset
meta$Infection <- factor(meta$Infection, levels = c("acute", "chronic"))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = meta, design = ~ Infection)

# Filter low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)
summary(res)

# Visualization and QC
vst <- vst(dds)
plotPCA(vst, intgroup = "Infection")

# Sample Distance Heatmap
library(pheatmap)
sampleDists <- dist(t(assay(vst)))
pheatmap(as.matrix(sampleDists), main = "Sample-to-Sample Distances")

# Volcano plot
library(ggplot2)
resDF <- as.data.frame(results(dds))
resDF$padj[is.na(resDF$padj)] <- 1

ggplot(resDF, aes(x = log2FoldChange, y = -log10(padj), color = padj < 0.05)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  xlab("Log2 Fold Change") + ylab("-Log10 Adjusted P-Value") +
  ggtitle("Volcano Plot: Chronic vs Acute")

let's look at our differentially expressed genes
plot(x = final_res$log2FoldChange, 
     y = -log10(final_res$padj),
     cex = 0.25,
     pch = 19, 
     col = 'grey',
     ylim = c(0,20),
     ylab = 'Adjusted P-Value',
     xlab = 'Log2 FC')

abline(v = c(-2, 2), h = -log10(0.5), lwd = 0.5, lty = 2)

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

mtext('A simple volcano plot')

# merge the two to do a clean and less memory efficient heatmap
degs <- rbind(countData[rownames(upregulated),], 
              countData[rownames(downregulated),])
pheatmap(degs, 
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         scale = 'row',
         show_colnames = T)

# upregulated and downregulated genes
rownames(upregulated)
rownames(downregulated)

# exporting the files
write.csv(upregulated, 'upregulated.csv')
write.csv(downregulated, 'downregulated.csv')
write.csv(countData, 'raw_counts.csv')

# Variance stabilizing transform (VST) for QC / PCA
vst <- vst(dds, blind = TRUE)
library(ggplot2)

# PCA plot (first two PCs)
pcaData <- plotPCA(vst, intgroup = "Infection", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = Infection, label = name)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of samples (VST)") +
  theme_minimal()

# Sample distance heatmap
sampleDists <- dist(t(assay(vst)))
sampleDistMat <- as.matrix(sampleDists)
rownames(sampleDistMat) <- colnames(vst)
colnames(sampleDistMat) <- colnames(vst)
pheatmap(sampleDistMat, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists, main = "Sample distances (VST)")

# Get results for the primary contrast (DESeq2 default uses the last level of factor as numerator)
# If you want a specific contrast, use results(dds, contrast = c("Infection","level1","level2"))
res <- results(dds)       # default comparison
summary(res)

# Basic MA plot and Volcano-like plot
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

# Export results to CSV
write.csv(as.data.frame(resLFC), file = "deseq2_results_lfc_shrunk.csv", row.names = TRUE)
write.csv(as.data.frame(res), file = "deseq2_results_raw.csv", row.names = TRUE)
cat("Wrote: deseq2_results_lfc_shrunk.csv and deseq2_results_raw.csv\n")

```
-----

## 🧠 Functional Enrichment Analysis (using ShinyGO)

Functional enrichment analysis of differentially expressed genes (DEGs) was conducted using ShinyGO v0.80 with _Staphylococcus aureus_ subsp. aureus NCTC 8325 as the reference genome. Out of 21 DEGs, 13 (62%) were successfully mapped to annotated S. aureus STRING database entries. Enrichment analysis revealed several significantly overrepresented pathways (FDR < 0.05), primarily associated with lipid metabolism and fatty acid degradation.

Terms such as lipid oxidation, fatty acid β-oxidation, lipid catabolic process, and AMP-binding and ABC transporter-like proteins were enriched (FDR = 7.0 × 10⁻³ to 1.7 × 10⁻²). These processes indicate metabolic reprogramming of _S. aureus_ during infection, likely reflecting enhanced utilization of host fatty acids and energy conservation strategies under stress conditions. Such metabolic adaptations are consistent with the bacterial transition to a chronic, biofilm-associated phenotype, which supports persistence and antibiotic tolerance in periprosthetic joint infections.

## 🧾 Conclusion
This study provides transcriptomic insight into the adaptive strategies employed by _Staphylococcus aureus_ during acute and chronic phases of periprosthetic joint infection. Differential expression and enrichment analyses revealed a strong association between chronic infection and genes involved in lipid oxidation and fatty acid catabolism. These findings highlight a metabolic transition that favors energy conservation, stress tolerance, and persistence within the host. Although only a subset of genes were functionally annotated, the enrichment of lipid metabolic pathways highlights the role of metabolic reprogramming in promoting chronic infection phenotypes. These results contribute to our understanding of _S. aureus_ pathogenesis and may inform future therapeutic approaches targeting metabolic pathways essential for bacterial persistence.
