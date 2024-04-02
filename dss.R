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

sim.data.Cheung <- simRNAseq(sim.opts.Cheung, n1=3, n2=3)
sim.data.Bottomly <- simRNAseq(sim.opts.Bottomly, n1=3, n2=3)
sim.data.MAQC <- simRNAseq(sim.opts.MAQC, n1=3, n2=3)

#######################################################################
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DSS")

library(DSS)

counts = sim.data.Cheung$counts
designs = c(rep(0, 3), rep(1, 3))
#designs <- data.frame(condition = as.factor(c(rep("Condition1", 3), rep("Condition2", 3))))
designs
seqData = newSeqCountSet(counts, designs)

seqData = estNormFactors(seqData)
seqData = estDispersion(seqData)

result = waldTest(seqData, 0, 1) # diff expression analysis
head(result, 5)

# Adjust p-values for multiple testing, using the Benjamini-Hochberg method as an example
result$padj <- p.adjust(result$pval, method = "BH")

# Filter for DE genes using the adjusted p-values and log2 fold changes
predicted_DE_genes <- result[result$padj < 0.05 & abs(result$lfc) > 1, ]

rownames(predicted_DE_genes) <- predicted_DE_genes$geneIndex
#result$predicted <- ifelse(rownames(result) %in% sim.data.Cheung$DEid, "DE", "Not DE")
#table(result$predicted)
predicted_DE_genes
#table(result$predicted)
#sum(is.na(result$predicted))
true_DE_genes <- sim.data.Cheung$DEid

TP_ids <- intersect(rownames(predicted_DE_genes), true_DE_genes)
FP_ids <- setdiff(rownames(predicted_DE_genes), true_DE_genes)
FN_ids <- setdiff(true_DE_genes, rownames(predicted_DE_genes))
TN_ids <- setdiff(rownames(result), union(true_DE_genes, rownames(predicted_DE_genes)))

confusion_matrix <- matrix(c(length(TP_ids), length(FP_ids),
                             length(FN_ids), length(TN_ids)),
                           nrow = 2, byrow = TRUE,
                           dimnames = list(Predicted = c("DE", "Not DE"),
                                           Actual = c("DE", "Not DE")))
confusion_matrix
# Performance Metrics
calculate_metrics <- function(confusion_matrix) {
  TP <- confusion_matrix["DE", "DE"]
  FP <- confusion_matrix["Not DE", "DE"]
  TN <- confusion_matrix["Not DE", "Not DE"]
  FN <- confusion_matrix["DE", "Not DE"]
  
  # Calculate metrics based on TP, FP, TN, FN
  sensitivity = TP / (TP + FN)
  specificity = TN / (TN + FP)
  precision = TP / (TP + FP)
  recall = TP / (TP + FN)  # same as sensitivity
  F1_score = 2 * (precision * recall) / (precision + recall)
  
  metrics = list(sensitivity = sensitivity, specificity = specificity, precision = precision,
                 recall = recall, F1_score = F1_score)
  return(metrics)
}

metrics <- calculate_metrics(confusion_matrix)
metrics
