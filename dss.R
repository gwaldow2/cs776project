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

################################################################################################
### Simulating data 

library(PROPER)
?RNAseq.SimOptions.2grp
?simRNAseq

# Simulation options

#library(DSS)
#library(DESeq2) 

#dds_dss <- DESeqDataSet_to_DSS(dds)
#dss_results <- DSS::testDifferentialExpression(dds_dss, method = "MLE")

sim.opts.Cheung = RNAseq.SimOptions.2grp(ngenes = 20000, p.DE=0.05,
                                         lOD="cheung", lBaselineExpr="cheung", sim.seed=123)

sim.opts.Bottomly = RNAseq.SimOptions.2grp(ngenes = 20000, p.DE=0.05,
                                           lOD="bottomly", lBaselineExpr="bottomly", sim.seed=123)

sim.opts.MAQC = RNAseq.SimOptions.2grp(ngenes = 20000, p.DE=0.05,
                                       lOD="maqc", lBaselineExpr="maqc", sim.seed=123)

# Simulated data

sim.data.Cheung <- simRNAseq(sim.opts.Cheung, n1=3, n2=3)
sim.data.Bottomly <- simRNAseq(sim.opts.Bottomly, n1=3, n2=3)
sim.data.MAQC <- simRNAseq(sim.opts.MAQC, n1=3, n2=3)
summary(sim.data.Cheung)
sim.data$counts
#######################################################################
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DSS")
library(DSS)

# Combine datasets into a list for easy iteration
sim.data.list <- list(
  Cheung = sim.data.Cheung,
  Bottomly = sim.data.Bottomly,
  MAQC = sim.data.MAQC
)

# Initialize a list to store results from all datasets
all_results_list <- list()

# Loop through each dataset
for (data_name in names(sim.data.list)) {
  sim.data <- sim.data.list[[data_name]]
  
  counts = sim.data$counts
  true_DE_genes = sim.data$DEid
  designs = c(rep(0, 3), rep(1, 3))
  designs
  #design <- data.frame(condition = as.factor(c(rep("Condition1", 3), rep("Condition2", 3))))
  seqData = newSeqCountSet(counts, designs)
  
  seqData = estNormFactors(seqData)
  seqData = estDispersion(seqData)
  
  result = waldTest(seqData, 0, 1) # diff expression analysis
  
  # parameter grid search:
  pval_cutoffs <- c(0.05, 0.01, 0.001)
  lfc_thresholds <- c(1, 1.5, 2)
  results_list <- list()
  
  for (pval_cutoff in pval_cutoffs) {
    for (lfc_threshold in lfc_thresholds) {
      
      # Perform DE analysis
      result$padj <- p.adjust(result$pval, method = "BH")  # Adjust p-values using Benjamin-Hochberg 
      predicted_DE_genes <- result[result$padj < pval_cutoff & abs(result$lfc) > lfc_threshold, ]
      rownames(predicted_DE_genes) <- predicted_DE_genes$geneIndex
      
      # Generate confusion matrix
      TP_ids <- intersect(rownames(predicted_DE_genes), true_DE_genes)
      FP_ids <- setdiff(rownames(predicted_DE_genes), true_DE_genes)
      FN_ids <- setdiff(true_DE_genes, rownames(predicted_DE_genes))
      TN_ids <- setdiff(rownames(result), union(true_DE_genes, rownames(predicted_DE_genes)))
      
      confusion_matrix <- matrix(c(length(TP_ids), length(FP_ids),
                                   length(FN_ids), length(TN_ids)),
                                 nrow = 2, byrow = TRUE,
                                 dimnames = list(Predicted = c("DE", "Not DE"),
                                                 Actual = c("DE", "Not DE")))
      
      # Calculate performance metrics
      metrics <- calculate_metrics(confusion_matrix)
      
      # Store results
      results_list[[paste("pval", pval_cutoff, "lfc", lfc_threshold, sep="_")]] <- metrics
      
      hist(predicted_DE_genes$pval)
      hist(predicted_DE_genes$padj)
    }
  }
  
  # Store results
  all_results_list[[data_name]] <- results_list
}

# all_results_list now contains the results for all datasets
all_results_list


### METHOD FOR USE FOR scRNA-Seq data ###

library(DSS)
library(BiocManager)

# Define the differential expression function
DE_dss <- function(counts, group, contrfml, p_thres=0.05, lfc_thres=1, 
                   cutoff=1, batch=NULL, thres_type="adjp", adj_method = "BH"){
  # Check input assumptions
  if (!is.factor(group)) {
    stop("Group must be a factor.")
  }
  
  # Initialize DSS SeqCountSet
  seqData <- newSeqCountSet(counts, as.integer(group) - 1)  # DSS expects zero-based groups
  
  # Normalize and estimate dispersion
  seqData <- estNormFactors(seqData)
  seqData <- estDispersion(seqData, method="pooled")
  
  # Handle batch effects if necessary
  if (!is.null(batch) && !is.null(batch$batch)) {
    seqData <- batchCorrect(seqData, batch$batch)
  }
  
  # Specify model matrix and make contrasts
  model_matrix <- model.matrix(~ group + batch)
  contrast_matrix <- makeContrasts(contrasts=contrfml, levels=model_matrix)
  
  # Perform Wald test
  result <- waldTest(seqData, contrast_matrix)
  
  # Adjust p-values and apply cutoffs
  result$padj <- p.adjust(result$pval, method = adj_method)
  de_results <- result[result$padj < p_thres & abs(result$logFC) >= lfc_thres, ]
  
  # Generate volcano plot
  volcano_plot <- ggplot(de_results, aes(x=logFC, y=-log10(padj))) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(p_thres), linetype = "dashed") +
    labs(title = "Volcano Plot", x = "Log Fold Change", y = "-log10 adjusted P-value")
  
  return(list(top_table = de_results, volcano_plot = volcano_plot))
}

