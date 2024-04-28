# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("sva")

#library(PROPER)
# library(limma)
library(ggplot2)
library(Matrix)
library(rhdf5)
library(sva)
#library(qvalue)
#library(fdrtool)
#library(DSS)
#library(DESeq2) 
# source("limma-voom.R")

library(DESeq2)

# file_path = "/storage10/shuchen/JerberData/alltime_fiaj_12100_5000.mtx"
# meta_path = "/storage10/shuchen/JerberData/alltime_fiaj_12100_5000.h5"
# h5ls(meta_path)

file_path = "C:/Nabil Awan/Nabil_UW-Madison/SEMESTERS_Courses_TAships/Spring 2024/BMI CS 776/Spring 2024/Project/Single cell data/alltime_fiaj_12100_5000.mtx"
meta_path = "C:/Nabil Awan/Nabil_UW-Madison/SEMESTERS_Courses_TAships/Spring 2024/BMI CS 776/Spring 2024/Project/Single cell data/alltime_fiaj_12100_5000.h5"
h5ls(meta_path)

data <- t(readMM(file_path))
cell_type <- h5read(meta_path, name="obs/celltype")
sample_id <- h5read(meta_path, name="obs/sample_id")
pool_id <- h5read(meta_path, name="obs/pool_id")
time_point <- h5read(meta_path, name="obs/time_point")
treatment <- h5read(meta_path, name="obs/treatment")
gene <- h5read(meta_path, name="var")

length(gene$index)

get_categories <- function(meta) {
  names = meta$codes
  for (i in 1:length(meta$categories)) {
    name = meta$categories[i]
    names[meta$codes==(i-1)] = name
  }
  return(names)
}

cell_type <- get_categories(cell_type)
sample_id <- get_categories(sample_id)
pool_id <- get_categories(pool_id)
time_point <- get_categories(time_point)
treatment <- get_categories(treatment)

data_pool <- ComBat(data, pool_id)
data_sample <- ComBat(data, sample_id)

unique(treatment)
# 'NONE''ROT'

min(data_pool[,which(time_point=="D52")])
# -1.86106558705981

temp_data = data_pool[,which(time_point=="D52")] - min(data_pool[,which(time_point=="D52")])
dim(temp_data)
# [1] 5000 7299
class(temp_data)
# [1] "matrix" "array"

# > dim(data_pool)
# [1]  5000 12100
# > class(data_pool)
# [1] "matrix" "array"

# > class(data_pool[,which(time_point=="D52")])
# [1] "matrix" "array" 
# > dim(data_pool[,which(time_point=="D52")])
# [1] 5000 7299


unique(as.factor(treatment[which(time_point=="D52")]))
# [1] NONE ROT 
# Levels: NONE ROT




# Load DESeq2 library
library(DESeq2)

# > dds_D52 <- DESeqDataSetFromMatrix(countData = temp_data,
#                                     +                                   colData = data.frame(treatment = trt),
#                                     +                                   design = ~ treatment)
# Error in DESeqDataSet(se, design = design, ignoreRank) : 
#   some values in assay are not integers

# Convert count data to integers
temp_data <- round(temp_data)

# Create DESeqDataSet object
dds_D52 <- DESeqDataSetFromMatrix(countData = temp_data,
                                  colData = data.frame(treatment = trt),
                                  design = ~ treatment)
# converting counts to integer mode

# Preprocess data
dds <- DESeq(dds_D52)

# > dds <- DESeq(dds_D52)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# Error in estimateDispersionsFit(object, fitType = fitType, quiet = quiet) : 
#   all gene-wise dispersion estimates are within 2 orders of magnitude
# from the minimum value, and so the standard curve fitting techniques will not work.
# One can instead use the gene-wise estimates as final estimates:
#   dds <- estimateDispersionsGeneEst(dds)
# dispersions(dds) <- mcols(dds)$dispGeneEst
# ...then continue with testing using nbinomWaldTest or nbinomLRT

# Estimate size factors
dds_D52 <- estimateSizeFactors(dds_D52)

# Estimate gene-wise dispersions using gene-wise estimates
dds_D52 <- estimateDispersionsGeneEst(dds_D52)

# Set the gene-wise estimates as final estimates
dispersions(dds_D52) <- mcols(dds_D52)$dispGeneEst

# Continue with testing using nbinomWaldTest or nbinomLRT
dds_D52 <- nbinomWaldTest(dds_D52)


# Perform differential expression analysis
res <- results(dds_D52)
dim(res)
# > dim(res)
# [1] 5000    6

write.csv(res, file="results_dds_D52.csv", row.names = F)

head(res)
# > head(res)
# log2 fold change (MLE): treatment ROT vs NONE 
# Wald test p-value: treatment ROT vs NONE 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat    pvalue      padj
# <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
#   1   2.23191     0.01024270 0.0228456  0.448345  0.653904  0.784137
# 2   2.26770     0.01686872 0.0226642  0.744288  0.456702  0.619107
# 3   2.05572    -0.00193479 0.0238057 -0.081274  0.935224        NA
# 4   2.30555    -0.04168231 0.0224829 -1.853959  0.063745  0.150707
# 5   2.17165    -0.02788269 0.0231641 -1.203704  0.228704  0.389537
# 6   2.10574    -0.00591979 0.0235212 -0.251679  0.801289        NA

# res$log2FoldChange: Log2 fold changes.
# res$pvalue: Raw p-values.
# res$padj: Adjusted p-values (if you perform multiple testing correction).

# Using raw statistic values from DESeq2 to correct for the adjusted p-values
library(fdrtool)
out.results <- fdrtool(res$stat, statistic = "normal", plot = FALSE)
head(out.results)

res$pval_fdrtool <- out.results$pval # Note: this is still unadjusted p-value
res$qval_fdrtool <- out.results$qval
res$lfdr_fdrtool <- out.results$lfdr
head(res)
dim(res)

hist(res$pvalue, breaks=100, xlim=c(0.0, 1.0), ylim=c(0, 400), xlab="p-values", main="Histogram of raw p-values (DESeq2)")
hist(res$padj, breaks=100, xlim=c(0.0, 1.0), ylim=c(0, 400), xlab="Adjusted p-values", main="Histogram of adjusted p-values (DESeq2)")
hist(res$pval_fdrtool, breaks=100, xlim=c(0.0, 1.0), ylim=c(0, 400), xlab="Adjusted p-values (fdrtool)", main="Histogram of unadjusted p-values (DESeq2)")
hist(res$qval_fdrtool, breaks=100, xlim=c(0.0, 1.0), ylim=c(0, 400), xlab="q-values", main="Histogram of q-values (DESeq2)")
hist(res$lfdr_fdrtool, breaks=100, xlim=c(0.0, 1.0), ylim=c(0, 400), xlab="Local FDR", main="Histogram of local FDR (DESeq2)")

sum(res$pvalue < 0.05)
sum(res$padj < 0.05, na.rm = T)
sum(res$pval_fdrtool < 0.05)
sum(res$qval_fdrtool < 0.05)
sum(res$lfdr_fdrtool < 0.05)
# > sum(res$pvalue < 0.05)
# [1] 1050
# > sum(res$padj < 0.05, na.rm = T)
# [1] 808
# > sum(res$pval_fdrtool < 0.05)
# [1] 1172
# > sum(res$qval_fdrtool < 0.05)
# [1] 869
# > sum(res$lfdr_fdrtool < 0.05)
# [1] 639


# install.packages("ggplot2")
# install.packages("ggrepel")
# install.packages("gtable", dependencies = TRUE)

library(ggplot2)
library(ggrepel)

library(ggplot2)
library(gtable)

# Sample data (replace with your actual data)
# data <- res[complete.cases(res), ]  # Filter out NA values
# data$log10padj <- -log10(data$padj)
data <- res
data$log10qval_fdrtool <- -log10(data$qval_fdrtool)


# Define colors based on adjusted p-value significance
# data$Significance <- ifelse(data$padj < 0.05, ifelse(data$log2FoldChange > 0, "Upregulated", "Downregulated"), "Not significant")
data$Significance <- ifelse(data$qval_fdrtool < 0.05, ifelse(data$log2FoldChange > 0, "Upregulated", "Downregulated"), "Not significant")


dim(data)
# > dim(data)
# [1] 5000   12
# When padj was used:
# > dim(data)
# [1] 2577    8

sum(data$qval_fdrtool < 0.05 & data$log2FoldChange > 0)
# [1] 176
# When padj was used:
# sum(data$padj < 0.05 & data$log2FoldChange > 0)
# [1] 163

sum(data$qval_fdrtool < 0.05 & data$log2FoldChange < 0)
# [1] 693
# When padj was used:
# sum(data$padj < 0.05 & data$log2FoldChange < 0)
# [1] 645

sum(data$qval_fdrtool > 0.05)
# [1] 4131
# When padj was used:
# sum(data$padj > 0.05)
# [1] 1769


table(data$Significance)
# > table(data$Significance)
# 
# Downregulated Not significant     Upregulated 
# 693            4131             176
# When padj was used:
# > table(data$Significance)
# 
# Downregulated Not significant     Upregulated 
# 645            1769             163 

### Volcano plot when padj (from DESeq2 output) was used

# library(ggplot2)
# 
# # Create a volcano plot
# myvolcanoplot <- ggplot(data = data, aes(x = log2FoldChange, y = log10padj, col = Significance)) +
#   geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
#   geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
#   geom_point(size = 2) +
#   scale_color_manual(values = c("red", "blue", "grey"),
#                      labels = c("Upregulated", "Downregulated", "Not significant")) +
#   coord_cartesian(ylim = c(0, 25), xlim = c(-0.2, 0.2)) +
#   labs(color = 'Significance',
#        x = expression("log"[2]*" Fold Change"),
#        y = expression("-log"[10]*" Adjusted p-value")) +
#   ggtitle('Volcano plot (DESeq2)') +
#   theme_minimal()
# 
# print(myvolcanoplot)
# 
# 
# ggplot(data = data, aes(x = log2FoldChange, y = log10padj, col = Significance)) +
#   geom_point(size = 2)
# 
# library(ggplot2)
# 
# # Create a volcano plot
# myvolcanoplot <- ggplot(data = data, aes(x = log2FoldChange, y = log10padj, col = Significance)) +
#   geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
#   geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
#   geom_point(size = 2) +
#   scale_color_manual(values = c("blue", "grey", "red"),  # Change color order
#                      breaks = c("Downregulated", "Not significant", "Upregulated"),  # Specify order of colors
#                      labels = c("Downregulated", "Not significant", "Upregulated")) +  # Adjust labels
#   coord_cartesian(ylim = c(0, 25), xlim = c(-0.2, 0.2)) +
#   labs(color = 'Significance',
#        x = expression("log"[2]*" Fold Change"),
#        y = expression("-log"[10]*" Adjusted p-value")) +
#   ggtitle('Volcano plot (DESeq2)') +
#   theme_minimal()
# 
# print(myvolcanoplot)


### Volcano plot when qval_fdrtool (from fdrtool output using raw stat from DESeq2) was used

library(ggplot2)

ggplot(data = data, aes(x = log2FoldChange, y = log10qval_fdrtool, col = Significance)) +
  geom_point(size = 2)

library(ggplot2)

# Create a volcano plot
myvolcanoplot <- ggplot(data = data, aes(x = log2FoldChange, y = log10qval_fdrtool, col = Significance)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"),  # Change color order
                     breaks = c("Downregulated", "Not significant", "Upregulated"),  # Specify order of colors
                     labels = c("Downregulated (n=693)", "Not significant (n=4131)", "Upregulated (n=176)")) +  # Adjust labels
  coord_cartesian(ylim = c(0, 15), xlim = c(-0.4, 0.4)) +
  labs(color = 'Significance',
       x = expression("log"[2]*" Fold Change"),
       y = expression("-log"[10]*" q-value")) +
  ggtitle('Volcano plot (DESeq2)') +
  theme_minimal()

print(myvolcanoplot)


head(data)
data$gene <- gene$index

write.csv(data, file="C:/Nabil Awan/Nabil_UW-Madison/SEMESTERS_Courses_TAships/Spring 2024/BMI CS 776/Spring 2024/Project/Single cell data/results_D52_DESeq2_fdrtool.csv", row.names = F)

data_qval0.05 <- subset(data, qval_fdrtool<0.05)
write.csv(data_qval0.05, file="C:/Nabil Awan/Nabil_UW-Madison/SEMESTERS_Courses_TAships/Spring 2024/BMI CS 776/Spring 2024/Project/Single cell data/results_D52_sig_genes_DESeq2_fdrtool.csv", row.names = F)
dim(data_qval0.05)
# [1] 869  13

######################################################
### Analysis for DA vs. non-DA cells
######################################################

trt <- as.factor(treatment[which(time_point=="D52")])
class(trt)
# [1] "factor"

table(trt)
# > table(trt)
# trt
# NONE  ROT 
# 3677 3622

# length(cell_type)
# [1] 12100

is.DA <- cell_type[which(time_point=="D52")]=="DA"
length(is.DA)

table(is.DA)
# is.DA
# FALSE  TRUE 
# 6865   434

table(trt, is.DA)
# > table(trt, is.DA)
# is.DA
# trt    FALSE TRUE
# NONE  3372  305
# ROT   3493  129

length(which(is.DA==T))
# > length(which(is.DA==T))
# [1] 434

dim(temp_data)
# [1] 5000 7299

DA_temp_data <- temp_data[,which(is.DA==T)]
dim(DA_temp_data)
# > dim(DA_temp_data)
# [1] 5000  434
nonDA_temp_data <- temp_data[,which(is.DA==F)]
dim(nonDA_temp_data)
# > dim(nonDA_temp_data)
# [1] 5000 6865

######################################################

### Analysis with the DA data

library(DESeq2)

# # Convert count data to integers
# temp_data <- round(temp_data)

trt <- as.factor(treatment[which(time_point=="D52")])
length(trt)
# [1] 7299
trt_DA <- as.factor(treatment[which(is.DA==T)])
length(trt_DA)
# [1] 434

# Create DESeqDataSet object
dds_D52_DA <- DESeqDataSetFromMatrix(countData = DA_temp_data,
                                  colData = data.frame(treatment = trt_DA),
                                  design = ~ treatment)
# converting counts to integer mode

# Preprocess data
dds_DA <- DESeq(dds_D52_DA)
# > dds_DA
# Error: object 'dds_DA' not found

# Estimate size factors
dds_D52_DA <- estimateSizeFactors(dds_D52_DA)

# Estimate gene-wise dispersions using gene-wise estimates
dds_D52_DA <- estimateDispersionsGeneEst(dds_D52_DA)

# Set the gene-wise estimates as final estimates
dispersions(dds_D52_DA) <- mcols(dds_D52_DA)$dispGeneEst

# Continue with testing using nbinomWaldTest or nbinomLRT
dds_D52_DA <- nbinomWaldTest(dds_D52_DA)

# Perform differential expression analysis
res_DA <- results(dds_D52_DA)
dim(res_DA)
# > dim(res_DA)
# [1] 5000    6

library(fdrtool)
out.results_DA <- fdrtool(res_DA$stat, statistic = "normal", plot = FALSE)
head(out.results_DA)

res_DA$pval_fdrtool <- out.results_DA$pval # Note: this is still unadjusted p-value
res_DA$qval_fdrtool <- out.results_DA$qval
res_DA$lfdr_fdrtool <- out.results_DA$lfdr
head(res_DA)
dim(res_DA)

sum(res_DA$pvalue < 0.05)
sum(res_DA$padj < 0.05, na.rm = T)
sum(res_DA$pval_fdrtool < 0.05)
sum(res_DA$qval_fdrtool < 0.05)
sum(res_DA$lfdr_fdrtool < 0.05)
# > sum(res_DA$pvalue < 0.05)
# [1] 0
# > sum(res_DA$padj < 0.05, na.rm = T)
# [1] 0
# > sum(res_DA$pval_fdrtool < 0.05)
# [1] 788
# > sum(res_DA$qval_fdrtool < 0.05)
# [1] 331
# > sum(res_DA$lfdr_fdrtool < 0.05)
# [1] 179

res_DA$log10qval_fdrtool <- -log10(res_DA$qval_fdrtool)

res_DA$Significance <- ifelse(res_DA$qval_fdrtool < 0.05, ifelse(res_DA$log2FoldChange > 0, "Upregulated", "Downregulated"), "Not significant")

table(res_DA$Significance)
# > table(res_DA$Significance)
# 
# Downregulated Not significant     Upregulated 
# 40            4669             291

res_DA$gene <- gene$index

write.csv(res_DA, file="C:/Nabil Awan/Nabil_UW-Madison/SEMESTERS_Courses_TAships/Spring 2024/BMI CS 776/Spring 2024/Project/Single cell data/results_DA_DESeq2_fdrtool.csv", row.names = F)

res_DA_qval0.05 <- subset(res_DA, qval_fdrtool<0.05)
write.csv(res_DA_qval0.05, file="C:/Nabil Awan/Nabil_UW-Madison/SEMESTERS_Courses_TAships/Spring 2024/BMI CS 776/Spring 2024/Project/Single cell data/results_DA_sig_genes_DESeq2_fdrtool.csv", row.names = F)
dim(res_DA_qval0.05)
# [1] 331  12


# Volcano plot for DA cells

library(ggplot2)

# Create a volcano plot
myvolcanoplot <- ggplot(data = res_DA, aes(x = log2FoldChange, y = log10qval_fdrtool, col = Significance)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"),  # Change color order
                     breaks = c("Downregulated", "Not significant", "Upregulated"),  # Specify order of colors
                     labels = c("Downregulated (n=40)", "Not significant (n=4669)", "Upregulated (n=291)")) +  # Adjust labels
  coord_cartesian(ylim = c(0, 15), xlim = c(-0.4, 0.4)) +
  labs(color = 'Significance',
       x = expression("log"[2]*" Fold Change"),
       y = expression("-log"[10]*" q-value")) +
  ggtitle('Volcano plot for DA cells (DESeq2)') +
  theme_minimal()

print(myvolcanoplot)


######################################################

### Analysis with the non-DA data

library(DESeq2)

# # Convert count data to integers
# temp_data <- round(temp_data)

trt <- as.factor(treatment[which(time_point=="D52")])
length(trt)
# [1] 7299
trt_nonDA <- as.factor(treatment[which(is.DA==F)])
length(trt_nonDA)
# [1] 6865

# Create DESeqDataSet object
dds_D52_nonDA <- DESeqDataSetFromMatrix(countData = nonDA_temp_data,
                                     colData = data.frame(treatment = trt_nonDA),
                                     design = ~ treatment)
# converting counts to integer mode

# Preprocess data
dds_nonDA <- DESeq(dds_D52_nonDA)
# > dds_nonDA
# Error: object 'dds_nonDA' not found

# Estimate size factors
dds_D52_nonDA <- estimateSizeFactors(dds_D52_nonDA)

# Estimate gene-wise dispersions using gene-wise estimates
dds_D52_nonDA <- estimateDispersionsGeneEst(dds_D52_nonDA)

# Set the gene-wise estimates as final estimates
dispersions(dds_D52_nonDA) <- mcols(dds_D52_nonDA)$dispGeneEst

# Continue with testing using nbinomWaldTest or nbinomLRT
dds_D52_nonDA <- nbinomWaldTest(dds_D52_nonDA)

# Perform differential expression analysis
res_nonDA <- results(dds_D52_nonDA)
dim(res_nonDA)
# > dim(res_nonDA)
# [1] 5000    6

library(fdrtool)
out.results_nonDA <- fdrtool(res_nonDA$stat, statistic = "normal", plot = FALSE)
head(out.results_nonDA)

res_nonDA$pval_fdrtool <- out.results_nonDA$pval # Note: this is still unadjusted p-value
res_nonDA$qval_fdrtool <- out.results_nonDA$qval
res_nonDA$lfdr_fdrtool <- out.results_nonDA$lfdr
head(res_nonDA)
dim(res_nonDA)

sum(res_nonDA$pvalue < 0.05)
sum(res_nonDA$padj < 0.05, na.rm = T)
sum(res_nonDA$pval_fdrtool < 0.05)
sum(res_nonDA$qval_fdrtool < 0.05)
sum(res_nonDA$lfdr_fdrtool < 0.05)
# > sum(res_nonDA$pvalue < 0.05)
# [1] 183
# > sum(res_nonDA$padj < 0.05, na.rm = T)
# [1] 66
# > sum(res_nonDA$pval_fdrtool < 0.05)
# [1] 1044
# > sum(res_nonDA$qval_fdrtool < 0.05)
# [1] 741
# > sum(res_nonDA$lfdr_fdrtool < 0.05)
# [1] 539

res_nonDA$log10qval_fdrtool <- -log10(res_nonDA$qval_fdrtool)

res_nonDA$Significance <- ifelse(res_nonDA$qval_fdrtool < 0.05, ifelse(res_nonDA$log2FoldChange > 0, "Upregulated", "Downregulated"), "Not significant")

table(res_nonDA$Significance)
# > table(res_nonDA$Significance)
# 
# Downregulated Not significant     Upregulated 
# 125            4259             616

res_nonDA$gene <- gene$index

write.csv(res_nonDA, file="C:/Nabil Awan/Nabil_UW-Madison/SEMESTERS_Courses_TAships/Spring 2024/BMI CS 776/Spring 2024/Project/Single cell data/results_nonDA_DESeq2_fdrtool.csv", row.names = F)

res_nonDA_qval0.05 <- subset(res_nonDA, qval_fdrtool<0.05)
write.csv(res_nonDA_qval0.05, file="C:/Nabil Awan/Nabil_UW-Madison/SEMESTERS_Courses_TAships/Spring 2024/BMI CS 776/Spring 2024/Project/Single cell data/results_nonDA_sig_genes_DESeq2_fdrtool.csv", row.names = F)
dim(res_nonDA_qval0.05)
# [1] 741  12

# Volcano plot for non-DA cells

library(ggplot2)

# Create a volcano plot
myvolcanoplot <- ggplot(data = res_nonDA, aes(x = log2FoldChange, y = log10qval_fdrtool, col = Significance)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"),  # Change color order
                     breaks = c("Downregulated", "Not significant", "Upregulated"),  # Specify order of colors
                     labels = c("Downregulated (n=125)", "Not significant (n=4259)", "Upregulated (n=616)")) +  # Adjust labels
  coord_cartesian(ylim = c(0, 15), xlim = c(-0.4, 0.4)) +
  labs(color = 'Significance',
       x = expression("log"[2]*" Fold Change"),
       y = expression("-log"[10]*" q-value")) +
  ggtitle('Volcano plot for non-DA cells (DESeq2)') +
  theme_minimal()

print(myvolcanoplot)

######################################################

### Comparing the genes that were significant in DA vs. non-DA cells

# Number of common significant genes 
length(intersect(res_nonDA_qval0.05$gene, res_DA_qval0.05$gene))
# [1] 150

# Number of common significant genes unique to DA cells
length(setdiff(res_DA_qval0.05$gene, res_nonDA_qval0.05$gene))
# [1] 181

# Number of common significant genes unique to non-DA cells
length(setdiff(res_nonDA_qval0.05$gene, res_DA_qval0.05$gene))
# [1] 591

# Venn diagram

# Install and load the VennDiagram package if you haven't already
# install.packages("VennDiagram")
library(VennDiagram)

# Define the two vectors
nonDA_genes <- res_nonDA_qval0.05$gene
DA_genes <- res_DA_qval0.05$gene

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = list(nonDA = nonDA_genes, DA = DA_genes),
  category.names = c("Non-DA cells", "DA cells"),
  filename = NULL
)

# Plot the Venn diagram
grid.draw(venn.plot)

## A fancier version

# Create the Venn diagram with smaller circles
venn.plot <- venn.diagram(
  x = list(nonDA = nonDA_genes, DA = DA_genes),
  category.names = c("Non-DA cells", "DA cells"),
  filename = NULL,
  fill = c("blue", "red"),
  alpha = 0.5,
  cat.col = c("blue", "red"),
  cat.fontface = "bold",
  cat.fontfamily = "serif",
  cat.cex = 1.5,
  cex = 1.5,
  fontfamily = "serif",
  scaled = TRUE, # Set scaled to TRUE to make circles smaller
  margin = 0.05 # Set margin between circles
)

# Plot the Venn diagram
grid.draw(venn.plot)



######################################################

