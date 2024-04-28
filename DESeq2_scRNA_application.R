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

unique(cell_type)
# > unique(cell_type)
# [1] "P_FPP"   "NB"      "FPP"     "Epen1"   "U_Neur1" "DA"      "U_Neur2" "Sert"    "Astro"   "Epen2"  
# [11] "U_Neur3"

unique(pool_id)
# > unique(pool_id)
# [1] "pool4" "pool5"

unique(sample_id)
# > unique(sample_id)
# [1] "5245STDY7487301" "5245STDY7487302" "5245STDY7520697" "5245STDY7520698" "5245STDY7520699"
# [6] "5245STDY7520700" "5245STDY7544154" "5245STDY7544155" "5245STDY7544156" "5245STDY7544157"
# [11] "5245STDY7544158" "5245STDY7544159" "5245STDY7544160" "5245STDY7544161" "5245STDY7577237"
# [16] "5245STDY7577238" "5245STDY7577239" "5245STDY7577240" "5245STDY7603454" "5245STDY7603455"
# [21] "5245STDY7603456" "5245STDY7603457" "5245STDY7631339" "5245STDY7631340" "5245STDY7631341"
# [26] "5245STDY7631342" "5245STDY7631343" "5245STDY7631344" "5245STDY7631345" "5245STDY7631346"

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

cell_type_temp_data <- cell_type[which(time_point=="D52")]
length(cell_type_temp_data)
# [1] 7299

table(cell_type_temp_data)
# > table(cell_type_temp_data)
# cell_type_temp_data
# Astro      DA   Epen1   Epen2     FPP   P_FPP    Sert U_Neur1 U_Neur3 
# 2937     434    1590      56    1085     720       9     425      43 

# pseudo-bulk option:
# 7299-434
# DA: 434, mat dim: 5000*434
# nonDA: 434, mat dim: 5000*434
# 6865/434
# Each bin has 15 cell types: average of expression for these random 15 cell types for each gene

# Options: 1) equal bin operation, 2) pseudo bulk on sample ID, 3) downsampling

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

### Heatmap of temp_data by DA cells and treatment group combination

dim(temp_data)

length(trt)
table(trt)

# > length(trt)
# [1] 7299
# > table(trt)
# trt
# NONE  ROT 
# 3677 3622 

length(is.DA)
table(is.DA)
# > length(is.DA)
# [1] 7299
# > table(is.DA)
# is.DA
# FALSE  TRUE 
# 6865   434 


# Install and load the gplots package if you haven't already
install.packages("gplots")
library(gplots)

# Generate sample data (replace this with your actual data)
# set.seed(123)
# temp_data <- matrix(rpois(5000*7299, lambda = 10), nrow = 5000, ncol = 7299)
# trt <- sample(c("NONE", "ROT"), size = 7299, replace = TRUE)
# is.DA <- sample(c(FALSE, TRUE), size = 7299, replace = TRUE)
# 
# # Define group labels
# group_labels <- ifelse(is.DA, "DA", "non-DA")  # Cell type labels
# group_labels <- paste(group_labels, trt, sep = "-")  # Combination of cell type and treatment
# 
# # Create a data frame with the group labels as column names
# colnames(temp_data) <- group_labels
# 
# # Convert row names to gene names (replace "genes" with your actual gene names)
# rownames(temp_data) <- paste("gene", 1:5000, sep = "")
# 
# # Generate heatmap
# heatmap.2(
#   temp_data,
#   Rowv = FALSE,
#   Colv = FALSE,
#   dendrogram = "none",
#   scale = "row",  # Scale rows (genes)
#   trace = "none",
#   col = colorRampPalette(c("blue", "white", "red"))(100),  # Define heatmap colors
#   key = TRUE,  # Include color key
#   keysize = 1,  # Set size of color key
#   margins = c(10, 10),  # Set margins around the heatmap
#   main = "Heatmap of Gene Expression",
#   xlab = "Samples",
#   ylab = "Genes"
# )

# 
# # Define group labels
# group_labels <- ifelse(is.DA, "DA", "non-DA")  # Cell type labels
# group_labels <- paste(group_labels, trt, sep = "-")  # Combination of cell type and treatment
# 
# # Create a data frame with the group labels as column names
# colnames(temp_data) <- group_labels
# 
# # Convert row names to gene names (replace "genes" with your actual gene names)
# rownames(temp_data) <- gene$index
# 
# # Generate heatmap
# heatmap.2(
#   temp_data,
#   Rowv = FALSE,
#   Colv = FALSE,
#   dendrogram = "none",
#   scale = "row",  # Scale rows (genes)
#   trace = "none",
#   col = colorRampPalette(c("blue", "white", "red"))(100),  # Define heatmap colors
#   key = TRUE,  # Include color key
#   keysize = 1,  # Set size of color key
#   margins = c(10, 10),  # Set margins around the heatmap
#   main = "Heatmap of Gene Expression",
#   xlab = "Cells",
#   ylab = "Genes"
# )

######################################################

### Comparing results after downsampling non-DA cells

# down_sample for non-DA 
set.seed(1234)
nonDA_down = sample(1:ncol(nonDA_temp_data), ncol(DA_temp_data))
# data_nonDA <- data_nonDA[,nonDA_down]
# meta_nonDA <- meta_nonDA[nonDA_down]
# dim(data_nonDA)

library(DESeq2)

# # Convert count data to integers
# temp_data <- round(temp_data)

trt <- as.factor(treatment[which(time_point=="D52")])
length(trt)
# [1] 7299
trt_nonDA <- as.factor(treatment[which(is.DA==F)])
length(trt_nonDA)
# [1] 6865
trt_nonDA_down <- as.factor(treatment[nonDA_down])
length(trt_nonDA_down)
# [1] 434

nonDA_temp_data_down <- nonDA_temp_data[,nonDA_down]
dim(nonDA_temp_data_down)
# [1] 5000  434

# Create DESeqDataSet object
dds_D52_nonDA_down <- DESeqDataSetFromMatrix(countData = nonDA_temp_data_down,
                                        colData = data.frame(treatment = trt_nonDA_down),
                                        design = ~ treatment)
# converting counts to integer mode

# # Preprocess data
# dds_nonDA_down <- DESeq(dds_D52_nonDA_down)
# # > dds_nonDA
# # Error: object 'dds_nonDA' not found

# Estimate size factors
dds_D52_nonDA_down <- estimateSizeFactors(dds_D52_nonDA_down)

# Estimate gene-wise dispersions using gene-wise estimates
dds_D52_nonDA_down <- estimateDispersionsGeneEst(dds_D52_nonDA_down)

# Set the gene-wise estimates as final estimates
dispersions(dds_D52_nonDA_down) <- mcols(dds_D52_nonDA_down)$dispGeneEst

# Continue with testing using nbinomWaldTest or nbinomLRT
dds_D52_nonDA_down <- nbinomWaldTest(dds_D52_nonDA_down)

# Perform differential expression analysis
res_nonDA_down <- results(dds_D52_nonDA_down)
dim(res_nonDA_down)
# > dim(res_nonDA_down)
# [1] 5000    6

library(fdrtool)
out.results_nonDA_down <- fdrtool(res_nonDA_down$stat, statistic = "normal", plot = FALSE)
head(out.results_nonDA_down)

res_nonDA_down$pval_fdrtool <- out.results_nonDA_down$pval # Note: this is still unadjusted p-value
res_nonDA_down$qval_fdrtool <- out.results_nonDA_down$qval
res_nonDA_down$lfdr_fdrtool <- out.results_nonDA_down$lfdr
head(res_nonDA_down)
dim(res_nonDA_down)

sum(res_nonDA_down$pvalue < 0.05)
sum(res_nonDA_down$padj < 0.05, na.rm = T)
sum(res_nonDA_down$pval_fdrtool < 0.05)
sum(res_nonDA_down$qval_fdrtool < 0.05)
sum(res_nonDA_down$lfdr_fdrtool < 0.05)
# > dim(res_nonDA_down)
# [1] 5000    9
# > sum(res_nonDA_down$pvalue < 0.05)
# [1] 11
# > sum(res_nonDA_down$padj < 0.05, na.rm = T)
# [1] 0
# > sum(res_nonDA_down$pval_fdrtool < 0.05)
# [1] 865
# > sum(res_nonDA_down$qval_fdrtool < 0.05)
# [1] 459
# > sum(res_nonDA_down$lfdr_fdrtool < 0.05)
# [1] 281

res_nonDA_down$log10qval_fdrtool <- -log10(res_nonDA_down$qval_fdrtool)

res_nonDA_down$Significance <- ifelse(res_nonDA_down$qval_fdrtool < 0.05, ifelse(res_nonDA_down$log2FoldChange > 0, "Upregulated", "Downregulated"), "Not significant")

table(res_nonDA_down$Significance)
# > table(res_nonDA_down$Significance)
# 
# Downregulated Not significant     Upregulated 
# 100            4541             359 

res_nonDA_down$gene <- gene$index

write.csv(res_nonDA_down, file="C:/Nabil Awan/Nabil_UW-Madison/SEMESTERS_Courses_TAships/Spring 2024/BMI CS 776/Spring 2024/Project/Single cell data/results_nonDA_down_DESeq2_fdrtool.csv", row.names = F)

res_nonDA_qval0.05_down <- subset(res_nonDA_down, qval_fdrtool<0.05)
write.csv(res_nonDA_qval0.05_down, file="C:/Nabil Awan/Nabil_UW-Madison/SEMESTERS_Courses_TAships/Spring 2024/BMI CS 776/Spring 2024/Project/Single cell data/results_nonDA_down_sig_genes_DESeq2_fdrtool.csv", row.names = F)
dim(res_nonDA_qval0.05_down)
# [1] 459  12

# Volcano plot for non-DA cells

library(ggplot2)

# Create a volcano plot
myvolcanoplot <- ggplot(data = res_nonDA_down, aes(x = log2FoldChange, y = log10qval_fdrtool, col = Significance)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"),  # Change color order
                     breaks = c("Downregulated", "Not significant", "Upregulated"),  # Specify order of colors
                     labels = c("Downregulated (n=100)", "Not significant (n=4541)", "Upregulated (n=359)")) +  # Adjust labels
  coord_cartesian(ylim = c(0, 15), xlim = c(-0.4, 0.4)) +
  labs(color = 'Significance',
       x = expression("log"[2]*" Fold Change"),
       y = expression("-log"[10]*" q-value")) +
  ggtitle('Volcano plot for non-DA cells (DESeq2, downsampled)') +
  theme_minimal()

print(myvolcanoplot)


### Venn diagram

library(VennDiagram)

res_DA_qval0.05 <- read.csv("C:/Nabil Awan/Nabil_UW-Madison/SEMESTERS_Courses_TAships/Spring 2024/BMI CS 776/Spring 2024/Project/Single cell data/results_DA_sig_genes_DESeq2_fdrtool.csv", header = T)

# Define the two vectors
nonDA_genes <- res_nonDA_qval0.05_down$gene
DA_genes <- res_DA_qval0.05$gene



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

### Comparison of number of DEGs

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


# > sum(res_nonDA_down$pvalue < 0.05)
# [1] 11
# > sum(res_nonDA_down$padj < 0.05, na.rm = T)
# [1] 0
# > sum(res_nonDA_down$pval_fdrtool < 0.05)
# [1] 865
# > sum(res_nonDA_down$qval_fdrtool < 0.05)
# [1] 459
# > sum(res_nonDA_down$lfdr_fdrtool < 0.05)
# [1] 281

######################################################

### Using pseudo bulk data

# a metadata dataframe for the variables that we are interested in
meta = data.frame(sample_id=as.factor(sample_id), cell_type=as.factor(cell_type), time_point=as.factor(time_point), treatment=as.factor(treatment))
meta

# pseudobulk metadata with 200 bulks in total
bulk_meta = unique(meta)
bulk_meta

# Aggregate counts by groups defined by bulk samples
aggr_data <- aggregate(t(data_sample), by = list(meta$sample_id, meta$cell_type, meta$time_point, meta$treatment), FUN = mean)
aggr_data = t(aggr_data[,-c(1:4)])
aggr_data

aggr_data_pos = aggr_data - min(aggr_data)
data_pos = data_sample - min(data_sample)
aggr_data_pos

# change bulk metadata to characters for future use
bulk_meta <- as.data.frame(lapply(bulk_meta, function(x) as.character(x)))
bulk_meta

aggr_data_pos = aggr_data - min(aggr_data)
data_pos = data_sample - min(data_sample)
aggr_data_pos

# aggregated data
data_DA <- aggr_data_pos[, which(bulk_meta$time_point=="D52" & bulk_meta$cell_type=="DA")]
meta_DA <- bulk_meta[which(bulk_meta$time_point=="D52" & bulk_meta$cell_type=="DA"),]$treatment
dim(data_DA)
# > dim(data_DA)
# [1] 5000   16
data_nonDA <- aggr_data_pos[, which(bulk_meta$time_point=="D52" & bulk_meta$cell_type!="DA")]
meta_nonDA <- bulk_meta[which(bulk_meta$time_point=="D52" & bulk_meta$cell_type!="DA"),]$treatment
dim(data_nonDA)
# > dim(data_nonDA)
# [1] 5000  113

# single cell data
data_DA <- data_pos[,which(time_point=="D52" & cell_type=="DA")]
meta_DA <- treatment[which(time_point=="D52" & cell_type=="DA")]
dim(data_DA)
# > dim(data_DA)
# [1] 5000  434
data_nonDA <- data_pos[,which(time_point=="D52" & cell_type!="DA")]
meta_nonDA <- treatment[which(time_point=="D52" & cell_type!="DA")]
dim(data_nonDA)
# > dim(data_nonDA)
# [1] 5000 6865

# down_sample for non-DA 
set.seed(1234)
nonDA_down = sample(1:ncol(data_nonDA), ncol(data_DA))
data_nonDA <- data_nonDA[,nonDA_down]
meta_nonDA <- meta_nonDA[nonDA_down]
dim(data_nonDA)

library(DESeq2)

# # Convert count data to integers
# temp_data <- round(temp_data)

# trt <- as.factor(treatment[which(time_point=="D52")])
# length(trt)
# # [1] 7299
# trt_nonDA <- as.factor(treatment[which(is.DA==F)])
# length(trt_nonDA)
# # [1] 6865
trt_nonDA_down <- as.factor(treatment[nonDA_down])
length(trt_nonDA_down)
# [1] 434

trt_nonDA_bulk_down <- as.factor(meta_nonDA)
length(trt_nonDA_bulk_down)
# [1] 434

dim(data_nonDA)

nonDA_temp_data_down <- nonDA_temp_data[,nonDA_down]
dim(nonDA_temp_data_down)
# [1] 5000  434

data_nonDA <- round(data_nonDA)

# Create DESeqDataSet object
dds_D52_nonDA_bulk_down <- DESeqDataSetFromMatrix(countData = data_nonDA,
                                             colData = data.frame(treatment = trt_nonDA_bulk_down),
                                             design = ~ treatment)
# converting counts to integer mode

# # Preprocess data
# dds_nonDA_down <- DESeq(dds_D52_nonDA_down)
# # > dds_nonDA
# # Error: object 'dds_nonDA' not found

# Estimate size factors
dds_D52_nonDA_bulk_down <- estimateSizeFactors(dds_D52_nonDA_bulk_down)

# Estimate gene-wise dispersions using gene-wise estimates
dds_D52_nonDA_bulk_down <- estimateDispersionsGeneEst(dds_D52_nonDA_bulk_down)

# Set the gene-wise estimates as final estimates
dispersions(dds_D52_nonDA_bulk_down) <- mcols(dds_D52_nonDA_bulk_down)$dispGeneEst

# Continue with testing using nbinomWaldTest or nbinomLRT
dds_D52_nonDA_bulk_down <- nbinomWaldTest(dds_D52_nonDA_bulk_down)

# Perform differential expression analysis
res_nonDA_bulk_down <- results(dds_D52_nonDA_bulk_down)
dim(res_nonDA_bulk_down)
# > dim(res_nonDA_down)
# [1] 5000    6

library(fdrtool)
out.results_nonDA_bulk_down <- fdrtool(res_nonDA_bulk_down$stat, statistic = "normal", plot = FALSE)
head(out.results_nonDA_bulk_down)

res_nonDA_bulk_down$pval_fdrtool <- out.results_nonDA_bulk_down$pval # Note: this is still unadjusted p-value
res_nonDA_bulk_down$qval_fdrtool <- out.results_nonDA_bulk_down$qval
res_nonDA_bulk_down$lfdr_fdrtool <- out.results_nonDA_bulk_down$lfdr
head(res_nonDA_bulk_down)
dim(res_nonDA_bulk_down)

sum(res_nonDA_bulk_down$pvalue < 0.05)
sum(res_nonDA_bulk_down$padj < 0.05, na.rm = T)
sum(res_nonDA_bulk_down$pval_fdrtool < 0.05)
sum(res_nonDA_bulk_down$qval_fdrtool < 0.05)
sum(res_nonDA_bulk_down$lfdr_fdrtool < 0.05)
# > sum(res_nonDA_bulk_down$pvalue < 0.05)
# [1] 2
# > sum(res_nonDA_bulk_down$padj < 0.05, na.rm = T)
# [1] 0
# > sum(res_nonDA_bulk_down$pval_fdrtool < 0.05)
# [1] 960
# > sum(res_nonDA_bulk_down$qval_fdrtool < 0.05)
# [1] 543
# > sum(res_nonDA_bulk_down$lfdr_fdrtool < 0.05)
# [1] 338

res_nonDA_bulk_down$log10qval_fdrtool <- -log10(res_nonDA_bulk_down$qval_fdrtool)

res_nonDA_bulk_down$Significance <- ifelse(res_nonDA_bulk_down$qval_fdrtool < 0.05, ifelse(res_nonDA_bulk_down$log2FoldChange > 0, "Upregulated", "Downregulated"), "Not significant")

table(res_nonDA_bulk_down$Significance)
# > table(res_nonDA_bulk_down$Significance)
# 
# Downregulated Not significant     Upregulated 
# 447            4457              96 

res_nonDA_bulk_down$gene <- gene$index

write.csv(res_nonDA_bulk_down, file="C:/Nabil Awan/Nabil_UW-Madison/SEMESTERS_Courses_TAships/Spring 2024/BMI CS 776/Spring 2024/Project/Single cell data/results_nonDA_bulk_down_DESeq2_fdrtool.csv", row.names = F)

res_nonDA_qval0.05_bulk_down <- subset(res_nonDA_bulk_down, qval_fdrtool<0.05)
write.csv(res_nonDA_qval0.05_bulk_down, file="C:/Nabil Awan/Nabil_UW-Madison/SEMESTERS_Courses_TAships/Spring 2024/BMI CS 776/Spring 2024/Project/Single cell data/results_nonDA_bulk_down_sig_genes_DESeq2_fdrtool.csv", row.names = F)
dim(res_nonDA_qval0.05_bulk_down)
# [1] 543  12

# Volcano plot for non-DA cells

library(ggplot2)

# Create a volcano plot
myvolcanoplot <- ggplot(data = res_nonDA_bulk_down, aes(x = log2FoldChange, y = log10qval_fdrtool, col = Significance)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"),  # Change color order
                     breaks = c("Downregulated", "Not significant", "Upregulated"),  # Specify order of colors
                     labels = c("Downregulated (n=447)", "Not significant (n=4457)", "Upregulated (n=96)")) +  # Adjust labels
  coord_cartesian(ylim = c(0, 15), xlim = c(-0.4, 0.4)) +
  labs(color = 'Significance',
       x = expression("log"[2]*" Fold Change"),
       y = expression("-log"[10]*" q-value")) +
  ggtitle('Volcano plot for non-DA cells (DESeq2, downsampled)') +
  theme_minimal()

print(myvolcanoplot)


### Venn diagram

library(VennDiagram)

res_DA_qval0.05 <- read.csv("C:/Nabil Awan/Nabil_UW-Madison/SEMESTERS_Courses_TAships/Spring 2024/BMI CS 776/Spring 2024/Project/Single cell data/results_DA_sig_genes_DESeq2_fdrtool.csv", header = T)

# Define the two vectors
nonDA_genes <- res_nonDA_qval0.05_bulk_down$gene
DA_genes <- res_DA_qval0.05$gene



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
