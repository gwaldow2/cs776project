library(edgeR)
#library(ggplot2)

#' Apply DE analysis using limma-voom
#'
#' @param count count matrix as a data.frame with genes by samples
#' @param group a factor object as the sample information, with the length = number of samples
#' @param contrfml a string like "A-B" to make contrast between A and B
#' @param p_thres a p value threshold for significance
#' @param lfc_thres a log FC threshold for significance
#' @param cutoff cutoff for filtering low expressed genes
#' @param batch batch
#' @param thres_type string "adjp" or "p", using adjusted p value or p value for significance
#' @param adj_method the FDR methods
#' @return a list: 1. toptab: topTable for the result, 2. de_genes: DE gene names, 3 volcano_plot
#' @export

DE_limmavoom <- function(counts, group, contrfml, p_thres=0.05, lfc_thres=1, 
                         cutoff=1, batch=NULL, thres_type="adjp", adj_method = "BH"){
  # 1. data input
  # Create DGEList object
  # gene_names <- rownames(counts)
  d0 <- DGEList(counts)
  
  # 2. Preprocessing
  # Calculate normalization factors
  d0 <- calcNormFactors(d0)
  
  # Filter low-expressed genes
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  if (length(drop) != 0) {
    d <- d0[-drop,] 
  } else {
    d <- d0
  }
  
  print(paste0(nrow(d), " out of ", nrow(d0), " genes left after filtering low-expressed genes"))
  
  # 3. Voom transformation and calculation of variance weights
  # Specify the model to be fitted.  We do this before using voom since voom uses variances of the model residuals (observed - fitted)
  if (is.null(batch)){
    mm <- model.matrix(~0 + group)
  } else {
    mm <- model.matrix(~0 + group + batch)
  }
  
  # The above specifies a model where each coefficient corresponds to a group mean
  
  # Voom
  # - Counts are transformed to log2 counts per million reads (CPM), where "per million reads" is defined based on the normalization factors we calculated earlier
  # - A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated
  # - A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression (see red line in plot above)
  # - The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs.
  y <- voom(d, mm, plot = FALSE)
  
  # 4. Fitting linear models in limma
  # lmFit fits a linear model using weighted least squares for each gene:
  fit <- lmFit(y, mm)
  
  # Comparisons between groups (log fold-changes) are obtained as _contrasts_ of these fitted linear models:
  contr <- limma::makeContrasts(contrasts=contrfml, levels = colnames(coef(fit)))
  print(contr)
  
  # Estimate contrast for each gene
  tmp <- contrasts.fit(fit, contr)
  
  # Empirical Bayes smoothing of standard errors (shrinks standard errors that are much larger or smaller than those from other genes towards the average standard error)
  tmp <- eBayes(tmp)
  
  # What genes are most differentially expressed?
  top.table <- topTable(tmp, sort.by = "P", adjust.method = adj_method, n = Inf)
  top.table$genes <- rownames(top.table)
  # head(top.table, 20)
  # - logFC: log2 fold change of I5.9/I5.6 
  # - AveExpr: Average expression across all samples, in log2 CPM
  # - t: logFC divided by its standard error 
  # - P.Value: Raw p-value (based on t) from test that logFC differs from 0
  # - adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
  # - B: log-odds that gene is DE (arguably less useful than the other columns)
  
  # DE gene index by adjusted p value
  # Filter significant genes
  if (thres_type == "adjp") {
    de_genes <- rownames(top.table)[top.table$adj.P.Val < p_thres]
    significant_genes <- subset(top.table, adj.P.Val < p_thres & abs(logFC) > lfc_thres)
  }
  else if (thres_type == "p") {
    de_genes <- rownames(top.table)[top.table$P.Value < p_thres]
    significant_genes <- subset(top.table, P.Value < p_thres & abs(logFC) > lfc_thres)
  }
  else {
    print(paste0("Wrong thres_type", thres_type, ". It should be adjp or p."))
    return(NULL)
  }
  
  # make volcano plot
  volcano_plot <- volcano_plot <- ggplot(top.table, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(color = "black", alpha = 0.5) +
    geom_text(data = significant_genes, aes(label = genes), hjust = -0.1, vjust = -0.1, size = 3) + # Add text labels for significant genes
    geom_vline(xintercept = 0, linetype = "dashed", color = "blue") + # Add a vertical line at logFC = 0
    geom_hline(yintercept = p_thres, linetype = "dashed", color = "red") + # Add a horizontal line at significance threshold
    labs(x = "logFC", y = "-log10(adj.P.Val)", title = "Volcano Plot") +
    theme_minimal()

  
  return(list(toptable=top.table, de_genes = de_genes, volcano_plot=volcano_plot))
}

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