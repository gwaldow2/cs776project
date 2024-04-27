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

trt <- as.factor(treatment[which(time_point=="D52")])
class(trt)
# [1] "factor"

table(trt)
# > table(trt)
# trt
# NONE  ROT 
# 3677 3622

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

hist(res$pvalue, breaks=100, xlim=c(0.0, 1.0), ylim=c(0, 400), xlab="p-values", main="Histogram of p-values (DESeq2)")
hist(res$padj, breaks=100, xlim=c(0.0, 1.0), ylim=c(0, 400), xlab="Adjusted p-values", main="Histogram of adjusted p-values (DESeq2)")

sum(res$pvalue < 0.05)
sum(res$padj < 0.05, na.rm = T)

# > sum(res$pvalue < 0.05)
# [1] 1050
# > sum(res$padj < 0.05)
# [1] NA
# > sum(res$pvalue < 0.05)
# [1] 1050
# > sum(res$padj < 0.05, na.rm = T)
# [1] 808

# install.packages("ggplot2")
# install.packages("ggrepel")
# install.packages("gtable", dependencies = TRUE)

library(ggplot2)
library(ggrepel)

library(ggplot2)
library(gtable)

# Sample data (replace with your actual data)
data <- res[complete.cases(res), ]  # Filter out NA values
data$log10padj <- -log10(data$padj)

# Define colors based on adjusted p-value significance
data$color <- ifelse(data$padj < 0.05, ifelse(data$log2FoldChange > 0, "red", "blue"), "grey")

table(data$color)
# > table(data$color)
# 
# blue grey  red 
# 645 1769  163 

### Another version

library(ggplot2)

# Create a volcano plot
myvolcanoplot <- ggplot(data = data, aes(x = log2FoldChange, y = log10padj, col = color)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("red", "blue", "grey"),
                     labels = c("Upregulated", "Downregulated", "Not significant")) +
  coord_cartesian(ylim = c(0, 25), xlim = c(-0.2, 0.2)) +
  labs(color = 'Significance',
       x = expression("log"[2]*" Fold Change"),
       y = expression("-log"[10]*" Adjusted p-value")) +
  ggtitle('Volcano plot (DESeq2)') +
  theme_minimal()

print(myvolcanoplot)


