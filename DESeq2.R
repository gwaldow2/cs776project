################################################################################################

### Simulating data 

library(PROPER)
?RNAseq.SimOptions.2grp
?simRNAseq

# Simulation options

sim.opts.Cheung = RNAseq.SimOptions.2grp(ngenes = 20000, p.DE=0.05,
                                         lOD="cheung", lBaselineExpr="cheung", sim.seed=123)

sim.opts.Bottomly = RNAseq.SimOptions.2grp(ngenes = 20000, p.DE=0.05,
                                           lOD="bottomly", lBaselineExpr="bottomly", sim.seed=123)

sim.opts.MAQC = RNAseq.SimOptions.2grp(ngenes = 20000, p.DE=0.05,
                                       lOD="maqc", lBaselineExpr="maqc", sim.seed=123)

# Simulated data

sim.data.Cheung <- simRNAseq(sim.opts.Cheung, n1=10, n2=10)
sim.data.Bottomly <- simRNAseq(sim.opts.Bottomly, n1=3, n2=3)
sim.data.MAQC <- simRNAseq(sim.opts.MAQC, n1=3, n2=3)


################################################################################################

### Installing and reading DESeq2 and related packages

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
# 
# install.packages("magrittr")  # For pipes %>%

# install.packages("caret", dependencies = TRUE)
# detach("package:cli", unload = TRUE)
# library(caret)

library(DESeq2)
library(magrittr)
# library(caret)

################################################################################################

### Results for sim.data.Cheung

################################################################################################

### Prepare DESeq2 analysis and calculate the confusion matrix:

# Create a data frame for the experimental design
design <- data.frame(condition = as.factor(c(rep("Condition1", 10), rep("Condition2", 10))))

# Create DESeqDataSet objects
# Assuming you have already loaded the necessary packages, generated the simulated count data 'sim.data.Cheung',
# and created the experimental design 'design'

# Create the DESeqDataSet object with the experimental design
dds_Cheung_mat <- DESeqDataSetFromMatrix(countData = sim.data.Cheung$counts, colData = design, design = ~ condition)

# Run DESeq2 analysis
dds_Cheung <- DESeq(dds_Cheung_mat)

# Get DESeq2 results
results_Cheung <- results(dds_Cheung)

str(results_Cheung)
dim(results_Cheung)
head(results_Cheung)

# Checking the distribution of p-values

hist(results_Cheung$padj)
hist(results_Cheung$pvalue)

sum(is.na(results_Cheung$pvalue))
sum(is.na(results_Cheung$padj))
# > sum(is.na(results_Cheung$pvalue))
# [1] 15690
# > sum(is.na(results_Cheung$padj))
# [1] 16919

# sum(results_Cheung$padj < 0.05, na.rm = T)
# sum(results_Cheung$pvalue < 0.05, na.rm = T)
# > sum(results_Cheung$padj < 0.05, na.rm = T)
# [1] 91
# > sum(results_Cheung$pvalue < 0.05, na.rm = T)
# [1] 261

?DESeq2::results

names(results_Cheung)
head(results_Cheung)

summary(results_Cheung$stat)

hist(results_Cheung$stat)

sum(is.na(results_Cheung$stat))
# > sum(is.na(results_Cheung$stat))
# [1] 15653

# The stat column in this case represents the Wald statistic, which is calculated as the log2 fold change divided by the standard error of the log2 fold change (lfcSE).
# Mathematically, it is calculated as stat = log2FoldChange / lfcSE.

### Now, let's calibrate the p-values using fdrtool

# Using raw statistic values from DSS
library(fdrtool)
complete.results_Cheung <- na.omit(results_Cheung)
out.results_Cheung <- fdrtool(complete.results_Cheung$stat, statistic = "normal", plot = FALSE)

hist(out.results_Cheung$pval)
hist(out.results_Cheung$lfdr)
hist(out.results_Cheung$qval)

sum(out.results_Cheung$lfdr < 0.05)
sum(out.results_Cheung$qval < 0.05)

sum(results_Cheung$padj < 0.05, na.rm = T)

# > sum(out.results_Cheung$lfdr < 0.05)
# [1] 83
# > 
#   > sum(results_Cheung$padj < 0.05, na.rm = T)
# [1] 91




results_Cheung$Gene_ID <- 1:20000

# Assuming 'results_Cheung' is already generated from DESeq2 analysis
# and 'sim.data.Cheung$DEid' contains the IDs of DE genes

# Extract DE genes based on log2 fold change threshold (e.g., abs(log2FoldChange) > 1)
DE_genes_Cheung <- subset(results_Cheung, abs(log2FoldChange) > 1 & padj < 0.05)
dim(DE_genes_Cheung)
DE_genes_Cheung <- subset(results_Cheung, padj < 0.05)
dim(DE_genes_Cheung)

# Check the structure of DE_genes_Cheung to see if 'gene_id' is a column
str(DE_genes_Cheung)

results_Cheung$predicted <- ifelse(abs(results_Cheung$log2FoldChange) > 1 & results_Cheung$padj < 0.05, "DE", "Not DE")
table(results_Cheung$predicted)
sum(is.na(results_Cheung$predicted))

# # Create a vector indicating whether each gene is DE or not based on its ID
# DE_genes_Cheung$predicted <- ifelse(DE_genes_Cheung$Gene_ID %in% sim.data.Cheung$DEid, "DE", "Not DE")
# 
# DE_genes_Cheung$actual <- ifelse(DE_genes_Cheung$baseMean > 10, "DE", "Not DE")  # Example threshold for baseMean

results_Cheung$actual <- ifelse(results_Cheung$Gene_ID %in% sim.data.Cheung$DEid, "DE", "Not DE")
table(results_Cheung$actual)

# Confusion matrix
confusion_matrix_Cheung <- table(results_Cheung$actual, results_Cheung$predicted)


# Function to calculate performance metrics
calculate_metrics <- function(confusion_matrix) {
  TP <- confusion_matrix["DE", "DE"]
  FP <- confusion_matrix["Not DE", "DE"]
  TN <- confusion_matrix["Not DE", "Not DE"]
  FN <- confusion_matrix["DE", "Not DE"]
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  accuracy <- (TP + TN) / sum(confusion_matrix)
  balanced_accuracy <- (sensitivity + specificity) / 2
  NPV <- TN / (TN + FN)
  PPV <- TP / (TP + FP)
  precision <- TP / (TP + FP)
  recall <- sensitivity
  
  metrics <- c(sensitivity = sensitivity, specificity = specificity, accuracy = accuracy,
               balanced_accuracy = balanced_accuracy, NPV = NPV, PPV = PPV,
               precision = precision, recall = recall)
  
  return(metrics)
}

# Calculate metrics for the confusion matrix
metrics_Cheung <- calculate_metrics(confusion_matrix_Cheung)
metrics_Cheung

################################################################################################

################################################################################################

### Results for sim.data.Bottomly

################################################################################################

### Prepare DESeq2 analysis and calculate the confusion matrix:

# Create a data frame for the experimental design
design <- data.frame(condition = as.factor(c(rep("Condition1", 3), rep("Condition2", 3))))

# Create DESeqDataSet objects
# Assuming you have already loaded the necessary packages, generated the simulated count data 'sim.data.Bottomly',
# and created the experimental design 'design'

# Create the DESeqDataSet object with the experimental design
dds_Bottomly <- DESeqDataSetFromMatrix(countData = sim.data.Bottomly$counts, colData = design, design = ~ condition)

# Run DESeq2 analysis
dds_Bottomly <- DESeq(dds_Bottomly)

# Get DESeq2 results
results_Bottomly <- results(dds_Bottomly)

str(results_Bottomly)
dim(results_Bottomly)
head(results_Bottomly)

results_Bottomly$Gene_ID <- 1:20000

# Assuming 'results_Bottomly' is already generated from DESeq2 analysis
# and 'sim.data.Bottomly$DEid' contains the IDs of DE genes

# Extract DE genes based on log2 fold change threshold (e.g., abs(log2FoldChange) > 1)
DE_genes_Bottomly <- subset(results_Bottomly, abs(log2FoldChange) > 1 & padj < 0.05)
dim(DE_genes_Bottomly)
DE_genes_Bottomly <- subset(results_Bottomly, padj < 0.05)
dim(DE_genes_Bottomly)

# Check the structure of DE_genes_Bottomly to see if 'gene_id' is a column
str(DE_genes_Bottomly)

results_Bottomly$predicted <- ifelse(abs(results_Bottomly$log2FoldChange) > 1 & results_Bottomly$padj < 0.05, "DE", "Not DE")
table(results_Bottomly$predicted)
sum(is.na(results_Bottomly$predicted))

# # Create a vector indicating whether each gene is DE or not based on its ID
# DE_genes_Bottomly$predicted <- ifelse(DE_genes_Bottomly$Gene_ID %in% sim.data.Bottomly$DEid, "DE", "Not DE")
# 
# DE_genes_Bottomly$actual <- ifelse(DE_genes_Bottomly$baseMean > 10, "DE", "Not DE")  # Example threshold for baseMean

results_Bottomly$actual <- ifelse(results_Bottomly$Gene_ID %in% sim.data.Bottomly$DEid, "DE", "Not DE")
table(results_Bottomly$actual)

# Confusion matrix
confusion_matrix_Bottomly <- table(results_Bottomly$actual, results_Bottomly$predicted)


# Function to calculate performance metrics
calculate_metrics <- function(confusion_matrix) {
  TP <- confusion_matrix["DE", "DE"]
  FP <- confusion_matrix["Not DE", "DE"]
  TN <- confusion_matrix["Not DE", "Not DE"]
  FN <- confusion_matrix["DE", "Not DE"]
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  accuracy <- (TP + TN) / sum(confusion_matrix)
  balanced_accuracy <- (sensitivity + specificity) / 2
  NPV <- TN / (TN + FN)
  PPV <- TP / (TP + FP)
  precision <- TP / (TP + FP)
  recall <- sensitivity
  
  metrics <- c(sensitivity = sensitivity, specificity = specificity, accuracy = accuracy,
               balanced_accuracy = balanced_accuracy, NPV = NPV, PPV = PPV,
               precision = precision, recall = recall)
  
  return(metrics)
}

# Calculate metrics for the confusion matrix
metrics_Bottomly <- calculate_metrics(confusion_matrix_Bottomly)
metrics_Bottomly

################################################################################################

################################################################################################

### Results for sim.data.MAQC

################################################################################################

### Prepare DESeq2 analysis and calculate the confusion matrix:

# Create a data frame for the experimental design
design <- data.frame(condition = as.factor(c(rep("Condition1", 3), rep("Condition2", 3))))

# Create DESeqDataSet objects
# Assuming you have already loaded the necessary packages, generated the simulated count data 'sim.data.MAQC',
# and created the experimental design 'design'

# Create the DESeqDataSet object with the experimental design
dds_MAQC <- DESeqDataSetFromMatrix(countData = sim.data.MAQC$counts, colData = design, design = ~ condition)

# Run DESeq2 analysis
dds_MAQC <- DESeq(dds_MAQC)

# Get DESeq2 results
results_MAQC <- results(dds_MAQC)

str(results_MAQC)
dim(results_MAQC)
head(results_MAQC)

results_MAQC$Gene_ID <- 1:20000

# Assuming 'results_MAQC' is already generated from DESeq2 analysis
# and 'sim.data.MAQC$DEid' contains the IDs of DE genes

# Extract DE genes based on log2 fold change threshold (e.g., abs(log2FoldChange) > 1)
DE_genes_MAQC <- subset(results_MAQC, abs(log2FoldChange) > 1 & padj < 0.05)
dim(DE_genes_MAQC)
DE_genes_MAQC <- subset(results_MAQC, padj < 0.05)
dim(DE_genes_MAQC)

# Check the structure of DE_genes_MAQC to see if 'gene_id' is a column
str(DE_genes_MAQC)

results_MAQC$predicted <- ifelse(abs(results_MAQC$log2FoldChange) > 1 & results_MAQC$padj < 0.05, "DE", "Not DE")
table(results_MAQC$predicted)
sum(is.na(results_MAQC$predicted))

# # Create a vector indicating whether each gene is DE or not based on its ID
# DE_genes_MAQC$predicted <- ifelse(DE_genes_MAQC$Gene_ID %in% sim.data.MAQC$DEid, "DE", "Not DE")
# 
# DE_genes_MAQC$actual <- ifelse(DE_genes_MAQC$baseMean > 10, "DE", "Not DE")  # Example threshold for baseMean

results_MAQC$actual <- ifelse(results_MAQC$Gene_ID %in% sim.data.MAQC$DEid, "DE", "Not DE")
table(results_MAQC$actual)

# Confusion matrix
confusion_matrix_MAQC <- table(results_MAQC$actual, results_MAQC$predicted)


# Function to calculate performance metrics
calculate_metrics <- function(confusion_matrix) {
  TP <- confusion_matrix["DE", "DE"]
  FP <- confusion_matrix["Not DE", "DE"]
  TN <- confusion_matrix["Not DE", "Not DE"]
  FN <- confusion_matrix["DE", "Not DE"]
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  accuracy <- (TP + TN) / sum(confusion_matrix)
  balanced_accuracy <- (sensitivity + specificity) / 2
  NPV <- TN / (TN + FN)
  PPV <- TP / (TP + FP)
  precision <- TP / (TP + FP)
  recall <- sensitivity
  
  metrics <- c(sensitivity = sensitivity, specificity = specificity, accuracy = accuracy,
               balanced_accuracy = balanced_accuracy, NPV = NPV, PPV = PPV,
               precision = precision, recall = recall)
  
  return(metrics)
}

# Calculate metrics for the confusion matrix
metrics_MAQC <- calculate_metrics(confusion_matrix_MAQC)
metrics_MAQC

################################################################################################

### Creating the data frame

metrics_Cheung_rounded <- round(metrics_Cheung, digits = 3)
metrics_Bottomly_rounded <- round(metrics_Bottomly, digits = 3)
metrics_MAQC_rounded <- round(metrics_MAQC, digits = 3)
sim.results.DESeq2 <- data.frame(metrics_Cheung_rounded, metrics_Bottomly_rounded, metrics_MAQC_rounded)
sim.results.DESeq2

################################################################################################


################################################################################################


################################################################################################


################################################################################################
