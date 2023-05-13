# setwd("~/kdata")

install.packages("glue")
library(glue)

R.utils::setOption("clusterProfiler.download.method","auto")

#--------import the packages------
library(dplyr)
library(tidyverse)
library(GEOquery)
library(xml2)
library(DESeq2)
library(sva)
library(org.Hs.eg.db)

install.packages("enrichR")
library(enrichR)

install.packages("R.utils")
library(R.utils)

install.packages("bigmemory")
install.packages("bigalgebra")
library(bigmemory)
library(bigalgebra)

install.packages("magrittr")
library(magrittr)

install.packages("ff")
library(ff)

install.packages("VennDiagram")
library(VennDiagram)

install.packages("ggplot2")
library(ggplot2)

library(KEGGREST)
gene_list <- keggList("hsa")

install.packages('clusterProfiler')
library(clusterProfiler)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install("GEOquery")
BiocManager::install("DESeq2")
BiocManager::install("sva")
BiocManager::install("org.Hs.eg.db")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("KEGGREST")

if (!requireNamespace("ReactomePA", quietly = TRUE)) {
  install.packages("ReactomePA")
}

library(ReactomePA)

install.packages("stringr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rsbml")
library(rsbml)


#--------read in the data------

aa_counts_ff <- read.table.ffdf(file = "GSE215841_CCCR_counts_AA.txt", header = TRUE, sep = "\t", VERBOSE = TRUE)
ea_counts_ff <- read.table.ffdf(file = "GSE215841_CCCR_counts_EA.txt", header = TRUE, sep = "\t", VERBOSE = TRUE)
nl_counts_ff <- read.table.ffdf(file = "GSE215841_CCCR_counts_NL.txt", header = TRUE, sep = "\t", VERBOSE = TRUE)
sscl_counts_ff <- read.table.ffdf(file = "GSE215841_CCCR_counts_SSCL.txt", header = TRUE, sep = "\t", VERBOSE = TRUE)

# Aggregate counts for each dataset

aggregate_counts <- function(counts_df) {
  counts_agg <- aggregate(. ~ Symbol, data = counts_df, sum)
  return(counts_agg)
}

aa_counts_agg <- aggregate_counts(aa_counts_ff)
ea_counts_agg <- aggregate_counts(ea_counts_ff)
nl_counts_agg <- aggregate_counts(nl_counts_ff)
sscl_counts_agg <- aggregate_counts(sscl_counts_ff)

#-------merge the aggregated count tables using gene names/IDs as keys-------

all_counts_agg <- merge(aa_counts_agg, ea_counts_agg, by = "Symbol", all = TRUE)
all_counts_agg <- merge(all_counts_agg, nl_counts_agg, by = "Symbol", all = TRUE)
all_counts_agg <- merge(all_counts_agg, sscl_counts_agg, by = "Symbol", all = TRUE)

all_counts_agg[is.na(all_counts_agg)] <- 0

#----Convert all_counts_agg to a regular data frame----

all_counts <- as.data.frame(all_counts_agg)

#-----Prepare the sample information and the count matrix---
sample_info <- data.frame(
  sample = colnames(all_counts)[-1],
  race = factor(gsub(".*(AA|EA).*", "\\1", colnames(all_counts)[-1])),
  condition = factor(gsub(".*(NL|SSCL).*", "\\1", colnames(all_counts)[-1])),
  row.names = colnames(all_counts)[-1]
)

sample_info$group <- factor(paste0(sample_info$race, "_", sample_info$condition))

count_matrix <- as.matrix(all_counts[,-1])
rownames(count_matrix) <- all_counts$Symbol

table(sample_info$race, sample_info$condition)

# Create the DESeqDataSet object with the updated sample_info
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ group
)

#--------Preprocessing---------
#Pre-filter low count genes:

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#--------Normalization-----------

dds <- DESeq(dds)
norm_counts <- counts(dds, normalized = TRUE)

#-------Differential expression analysis--------

# AA-NL vs EA-NL
dds_aanl_eanl <- dds[ , dds$group %in% c("AA_NL", "EA_NL")]
dds_aanl_eanl <- DESeq(dds_aanl_eanl)
res_aanl_eanl <- results(dds_aanl_eanl)

# AA-SScL vs AA-NL
dds_aasscl_aanl <- dds[ , dds$group %in% c("AA_SSCL", "AA_NL")]
dds_aasscl_aanl <- DESeq(dds_aasscl_aanl)
res_aasscl_aanl <- results(dds_aasscl_aanl)

# EA-SScL vs EA-NL
dds_easscl_eanl <- dds[ , dds$group %in% c("EA_SSCL", "EA_NL")]
dds_easscl_eanl <- DESeq(dds_easscl_eanl)
res_easscl_eanl <- results(dds_easscl_eanl)

# AA-SScL vs EA-SScL
dds_aasscl_easscl <- dds[, dds$group %in% c("AA_SSCL", "EA_SSCL")]
dds_aasscl_easscl <- DESeq(dds_aasscl_easscl)
res_aasscl_easscl <- results(dds_aasscl_easscl)

DESeq2::summary(res_aanl_eanl)
DESeq2::summary(res_aasscl_aanl)
DESeq2::summary(res_easscl_eanl)
DESeq2::summary(res_aasscl_easscl)

#------Perform variance stabilizing transformation------
vst_counts <- vst(dds)

# Calculate PCA
pca_data <- plotPCA(vst_counts, intgroup = c("race", "condition"), returnData = TRUE)

# Create PCA plot
ggplot(pca_data, aes(PC1, PC2, color = race, shape = condition)) +
  geom_point(size = 3) +
  theme_bw() +
  xlab(paste0("PC1: ", round(100 * attr(pca_data, "percentVar")[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(100 * attr(pca_data, "percentVar")[2], 2), "% variance"))


#-------visualizations for the differential expression analysis-------

# Function to create volcano plot
create_volcano_plot <- function(res, title) {
  res_df <- as.data.frame(res)
  res_df$group <- factor(ifelse(res_df$padj < 0.05 & (res_df$log2FoldChange > 1 | res_df$log2FoldChange < -1), "Significant", "Not Significant"))
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = group)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    labs(title = title, x = "log2 Fold Change", y = "-log10 Adjusted p-value") +
    theme(legend.position = "bottom") +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "blue"))
}

# Create volcano plots
volcano_plot_aanl_eanl <- create_volcano_plot(res_aanl_eanl, "AA-NL vs EA-NL")
volcano_plot_aasscl_aanl <- create_volcano_plot(res_aasscl_aanl, "AA-SScL vs AA-NL")
volcano_plot_easscl_eanl <- create_volcano_plot(res_easscl_eanl, "EA-SScL vs EA-NL")
volcano_plot_aasscl_easscl <- create_volcano_plot(res_aasscl_easscl, "AA-SScL vs EA-SScL")

# Display volcano plots
volcano_plot_aanl_eanl
volcano_plot_aasscl_aanl
volcano_plot_easscl_eanl
volcano_plot_aasscl_easscl


#------comparing the numerical results o identify the differences between the datasets----
# Extract log2 fold change values and adjusted p-values from DESeqResults objects
extract_numerical_results <- function(res) {
  res_df <- as.data.frame(res)
  return(res_df[, c("log2FoldChange", "padj")])
}

res_aanl_eanl_num <- extract_numerical_results(res_aanl_eanl)
res_aasscl_aanl_num <- extract_numerical_results(res_aasscl_aanl)
res_easscl_eanl_num <- extract_numerical_results(res_easscl_eanl)

# Add a column to indicate the comparison group
res_aanl_eanl_num$comparison <- "AA-NL vs EA-NL"
res_aasscl_aanl_num$comparison <- "AA-SScL vs AA-NL"
res_easscl_eanl_num$comparison <- "EA-SScL vs EA-NL"

# Merge the data frames
all_comparisons <- rbind(res_aanl_eanl_num, res_aasscl_aanl_num, res_easscl_eanl_num)

# Plot the distribution of log2 fold change values using a boxplot
ggplot(all_comparisons, aes(x = comparison, y = log2FoldChange, fill = comparison)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribution of log2 Fold Change Values", y = "log2 Fold Change") +
  theme(legend.position = "none")

# Plot the distribution of -log10 adjusted p-values using a boxplot
ggplot(all_comparisons, aes(x = comparison, y = -log10(padj), fill = comparison)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribution of -log10 Adjusted p-values", y = "-log10 Adjusted p-value") +
  theme(legend.position = "none")



#-----Identify significant differentially expressed genes (DEGs)------

# Adjusted p-value threshold
alpha <- 0.05

# Log2 fold change threshold
lfc_threshold <- 1

# Find significant DEGs for each comparison
sig_DEGs_aanl_eanl <- rownames(res_aanl_eanl)[res_aanl_eanl$padj < alpha & abs(res_aanl_eanl$log2FoldChange) > lfc_threshold]
sig_DEGs_aasscl_aanl <- rownames(res_aasscl_aanl)[res_aasscl_aanl$padj < alpha & abs(res_aasscl_aanl$log2FoldChange) > lfc_threshold]
sig_DEGs_easscl_eanl <- rownames(res_easscl_eanl)[res_easscl_eanl$padj < alpha & abs(res_easscl_eanl$log2FoldChange) > lfc_threshold]

#------visualize the number of significant differentially expressed genes (DEGs) for each comparison-----

# Create a data frame with the counts of significant DEGs for each comparison
sig_DEGs_counts <- data.frame(
  Comparison = c("AA-NL vs EA-NL", "AA-SScL vs AA-NL", "EA-SScL vs EA-NL"),
  DEGs = c(length(sig_DEGs_aanl_eanl), length(sig_DEGs_aasscl_aanl), length(sig_DEGs_easscl_eanl))
)

# Create a bar plot to visualize the counts of significant DEGs
ggplot(sig_DEGs_counts, aes(x = Comparison, y = DEGs, fill = Comparison)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Significant Differentially Expressed Genes", x = "Comparison", y = "Number of DEGs") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("AA-NL vs EA-NL" = "dodgerblue", "AA-SScL vs AA-NL" = "darkorange", "EA-SScL vs EA-NL" = "purple"))


#-----Pathway analysis-------

sig_DEGs_entrez_aanl_eanl <- mapIds(org.Hs.eg.db, keys = sig_DEGs_aanl_eanl, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
sig_DEGs_entrez_aasscl_aanl <- mapIds(org.Hs.eg.db, keys = sig_DEGs_aasscl_aanl, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
sig_DEGs_entrez_easscl_eanl <- mapIds(org.Hs.eg.db, keys = sig_DEGs_easscl_eanl, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

#Run Gene Ontology (GO) enrichment analysis

go_aanl_eanl <- enrichGO(gene = sig_DEGs_entrez_aanl_eanl, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
go_aasscl_aanl <- enrichGO(gene = sig_DEGs_entrez_aasscl_aanl, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
go_easscl_eanl <- enrichGO(gene = sig_DEGs_entrez_easscl_eanl, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)

#Visualize the enriched GO terms using dot plots

# Run KEGG pathway enrichment analysis
kegg_aanl_eanl <- enrichKEGG(gene = sig_DEGs_entrez_aanl_eanl, organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.05)
kegg_aasscl_aanl <- enrichKEGG(gene = sig_DEGs_entrez_aasscl_aanl, organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.05)
kegg_easscl_eanl <- enrichKEGG(gene = sig_DEGs_entrez_easscl_eanl, organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.05)

dotplot(go_aanl_eanl, showCategory = 10) + labs(title = "GO Enrichment: AA-NL vs EA-NL") + 
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

dotplot(go_aasscl_aanl, showCategory = 10) + labs(title = "GO Enrichment: AA-SScL vs AA-NL") + 
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

dotplot(go_easscl_eanl, showCategory = 10) + labs(title = "GO Enrichment: EA-SScL vs EA-NL") + 
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

# Run KEGG pathway enrichment analysis with "hsapiens" as the organism parameter
kegg_aanl_eanl <- enrichKEGG(gene = sig_DEGs_entrez_aanl_eanl, organism = "hsapiens", pAdjustMethod = "BH", qvalueCutoff = 0.05)
kegg_aasscl_aanl <- enrichKEGG(gene = sig_DEGs_entrez_aasscl_aanl, organism = "hsapiens", pAdjustMethod = "BH", qvalueCutoff = 0.05)
kegg_easscl_eanl <- enrichKEGG(gene = sig_DEGs_entrez_easscl_eanl, organism = "hsapiens", pAdjustMethod = "BH", qvalueCutoff = 0.05)


#-----Run Reactome pathway enrichment analysis-------
reactome_aanl_eanl <- enrichPathway(gene = sig_DEGs_entrez_aanl_eanl, organism = "human", readable = TRUE, pAdjustMethod = "BH", qvalueCutoff = 0.05)
reactome_aasscl_aanl <- enrichPathway(gene = sig_DEGs_entrez_aasscl_aanl, organism = "human", readable = TRUE, pAdjustMethod = "BH", qvalueCutoff = 0.05)
reactome_easscl_eanl <- enrichPathway(gene = sig_DEGs_entrez_easscl_eanl, organism = "human", readable = TRUE, pAdjustMethod = "BH", qvalueCutoff = 0.05)

dotplot(reactome_aanl_eanl, showCategory = 10) + labs(title = "Reactome Enrichment: AA-NL vs EA-NL") +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

dotplot(reactome_aasscl_aanl, showCategory = 10) + labs(title = "Reactome Enrichment: AA-SScL vs AA-NL") +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

dotplot(reactome_easscl_eanl, showCategory = 10) + labs(title = "Reactome Enrichment: EA-SScL vs EA-NL") +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))


#----- Check if any genes are mappable in KEGG----
all_genes <- bitr_kegg(c(sig_DEGs_entrez_aanl_eanl, sig_DEGs_entrez_aasscl_aanl, sig_DEGs_entrez_easscl_eanl), fromType = "entrezgene", toType = "kegg", organism = "hsa")

mappable_genes_aanl_eanl <- sig_DEGs_entrez_aanl_eanl %in% names(hsa_genes)
mappable_genes_aasscl_aanl <- sig_DEGs_entrez_aasscl_aanl %in% names(hsa_genes)
mappable_genes_easscl_eanl <- sig_DEGs_entrez_easscl_eanl %in% names(hsa_genes)

print("Mappable genes in sig_DEGs_entrez_aanl_eanl:")
print(sum(mappable_genes_aanl_eanl))

print("Mappable genes in sig_DEGs_entrez_aasscl_aanl:")
print(sum(mappable_genes_aasscl_aanl))

print("Mappable genes in sig_DEGs_entrez_easscl_eanl:")
print(sum(mappable_genes_easscl_eanl))




#------Identify potential biomarkers--------

# Function to get the top DEGs by adjusted p-value and log2 fold change
get_top_DEGs <- function(res, n = 10) {
  res_sorted <- res[order(res$padj, decreasing = FALSE), ]
  top_DEGs <- rownames(res_sorted)[1:n]
  return(top_DEGs)
}

# Get the top 10 DEGs for each comparison
top_DEGs_aanl_eanl <- get_top_DEGs(res_aanl_eanl)
top_DEGs_aasscl_aanl <- get_top_DEGs(res_aasscl_aanl)
top_DEGs_easscl_eanl <- get_top_DEGs(res_easscl_eanl)

# Function to extract genes from enriched terms
get_genes_from_enriched_terms <- function(enriched_terms) {
  enriched_genes <- unlist(lapply(enriched_terms@result$geneID, function(x) strsplit(x, "/")[[1]]))
  return(enriched_genes)
}

# Get the genes involved in enriched GO terms and KEGG pathways
enriched_genes_kegg_aanl_eanl <- get_genes_from_enriched_terms(kegg_aanl_eanl)
enriched_genes_kegg_aasscl_aanl <- get_genes_from_enriched_terms(kegg_aasscl_aanl)
enriched_genes_kegg_easscl_eanl <- get_genes_from_enriched_terms(kegg_easscl_eanl)

# Function to find potential biomarkers
find_potential_biomarkers <- function(top_DEGs, enriched_genes) {
  potential_biomarkers <- intersect(top_DEGs, enriched_genes)
  return(potential_biomarkers)
}

# Find potential biomarkers for each comparison
potential_biomarkers_aanl_eanl <- find_potential_biomarkers(top_DEGs_aanl_eanl, c(enriched_genes_go_aanl_eanl, enriched_genes_kegg_aanl_eanl))
potential_biomarkers_aasscl_aanl <- find_potential_biomarkers(top_DEGs_aasscl_aanl, c(enriched_genes_go_aasscl_aanl, enriched_genes_kegg_aasscl_aanl))
potential_biomarkers_easscl_eanl <- find_potential_biomarkers(top_DEGs_easscl_eanl, c(enriched_genes_go_easscl_eanl,enriched_genes_kegg_easscl_eanl))
        
                                                                              
#-------- Download the Recon3D model in SBML format-----

Recon3D <- read_xml("Recon3D.xml")









