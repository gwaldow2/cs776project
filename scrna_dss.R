  # Function to calculate performance metrics
library(ggplot2)
library(edgeR)
library(DSS)
library(BiocManager)

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

library(DSS)
library(BiocManager)

# Define the differential expression function
DE_dss <- function(counts, group, contrfml, p_thres=0.05, lfc_thres=1, 
                   cutoff=1, batch=NULL, thres_type="adjp", adj_method = "BH"){
  
  
  design <- matrix((as.integer(trt) - 1), ncol=1)
  # Initialize DSS SeqCountSet
  seqData <- newSeqCountSet(counts, design)  # DSS expects zero-based groups
  print(head(seqData))
  # Normalize and estimate dispersion
  seqData <- estNormFactors(seqData)
  seqData <- estDispersion(seqData)
  
  # Handle batch effects if necessary
  #if (!is.null(batch) && !is.null(batch$batch)) {
  #  seqData <- batchCorrect(seqData, batch$batch)
  #}
  
  # Specify model matrix and make contrasts
  #if (!is.null(batch) && !is.null(batch$batch)) {
  #  model_matrix <- model.matrix(~ group + batch$batch)
  #} else {
  #  model_matrix <- model.matrix(~ group)
  #}
  model_matrix <- model.matrix(~ group - 1)
  #contrast_matrix <- matrix(c(-1, 1), ncol = 1)
  #colnames(contrast_matrix) <- "groupROT-groupNONE"
  #rownames(contrast_matrix) <- colnames(model_matrix)
  print(model_matrix)
  print(dim(model_matrix))
  #print(contrast_matrix)
  #print(dim(contrast_matrix))
  
  sampleA <- which(group == "NONE")
  sampleB <- which(group == "ROT")
  
  # Perform Wald test with the specified samples
  result <- waldTest(seqData, sampleA = sampleA, sampleB = sampleB, equal.var = FALSE)
  print(str(result))
  print(head(result))
  print(unique(result$muA))
  
  if ("muA" %in% names(result) && "muB" %in% names(result)) {
    result$logFC <- log2(result$muB / result$muA)
  }
  # Adjust p-values and apply cutoffs
  result$padj <- p.adjust(result$pval, method = adj_method)
  de_results <- result[result$padj < p_thres & abs(result$logFC) >= lfc_thres, ]
  result$logFC
  adj_method
  # Generate volcano plot
  volcano_plot <- ggplot(de_results, aes(x=logFC, y=-log10(padj))) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(p_thres), linetype = "dashed") +
    labs(title = "Volcano Plot", x = "Log Fold Change", y = "-log10 adjusted P-value")
  
  return(list(top_table = de_results, volcano_plot = volcano_plot))
}


# Load necessary libraries
library(Matrix)
library(rhdf5)
library(sva)  # For batch effect correction using ComBat
library(DSS)


# Load data and metadata
file_path = "alltime_fiaj_12100_5000.mtx"
meta_path = "alltime_fiaj_12100_5000.h5"
h5ls(meta_path)

data <- t(readMM(file_path))  # Transpose to get genes in rows and samples in columns
cell_type <- h5read(meta_path, name="obs/celltype")
sample_id <- h5read(meta_path, name="obs/sample_id")
pool_id <- h5read(meta_path, name="obs/pool_id")
time_point <- h5read(meta_path, name="obs/time_point")
treatment <- h5read(meta_path, name="obs/treatment")
gene <- h5read(meta_path, name="var")

# Decode metadata categories
get_categories <- function(meta) {
  names = meta$codes
  for (i in 1:length(meta$categories)) {
    name = meta$categories[i]
    names[meta$codes == (i-1)] = name
  }
  return(names)
}

cell_type <- get_categories(cell_type)
sample_id <- get_categories(sample_id)
pool_id <- get_categories(pool_id)
time_point <- get_categories(time_point)
treatment <- get_categories(treatment)

# Apply ComBat for batch correction
data_pool <- ComBat(data, pool_id)
data_sample <- ComBat(data, sample_id)

# Filter data to a specific time point and normalize
temp_data = data_pool[, which(time_point == "D52")] - min(data_pool[, which(time_point == "D52")])
trt <- as.factor(treatment[which(time_point=="D52")])
temp_data <- round(temp_data)

# Assign column names to temp_data based on the sequence of trt
colnames(temp_data) <- as.character(seq_along(trt))


# Perform differential expression analysis with DSS
res <- DE_dss(counts=temp_data, group=trt, 
              contrfml = "groupROT-groupNONE", cutoff=1.0, thres_type="adjp",
              p_thres=0.05, lfc_thres=0, adj_method = "BH")

# Another DE analysis comparing different time points
#res2 <- DE_dss(counts=temp_data, group=trt, 
#               contrfml = "groupROT-groupNONE", cutoff=1.0, thres_type="adjp",
#               p_thres=0.05, lfc_thres=0, adj_method = "BH")

# Visualize results: Volcano plots, histograms of p-values, and adjusted p-values
print(res$volcano_plot)
#print(res2$volcano_plot)

hist(res$toptable$P.Value, breaks=100, main="Histogram of P-values")
#hist(res$toptable$adj.P.Val, breaks=100, main="Histogram of Adjusted P-values")

# Additional analysis based on results
sum(res$toptable$P.Value < 0.05)  # Count of significant p-values
