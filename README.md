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
- **Trimmomatic** → adapter/low-quality trimming
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

```
-----
# Upregulated genes 

