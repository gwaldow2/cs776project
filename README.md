This is the repository of scripts used for the cs776 project:
"Comparison of popular methods for differential gene expression analysis using simulated and real RNA-seq data"
The purpose of each script is described here:
--------------------------------------------------
DESeq2.R :                    Provides method for performing preprocessing and DESeq2 on given data.
DESeq2_scRNA_application.R :  Applies the DESeq2 method to the Jerber2021 dataset.
PROPER_data_generation.R :    File to generate the PROPER datasets.
dss.R :                       Applies dss method to simulated PROPER datasets. 
jerber_py.ipynb :             Pythonic Jupyter Notebook script to convert h5ad data into a format that limma-voom and DESeq2 can work with.
limma-voom.R :                Provides method for performing preprocessing and limma-voom on given data.
limma-voom.Rmd :              Applies the limma-voom method from limma-voom.R to simulated PROPER datasets.
limma-voom_jerber.ipynb :     R Jupyter Notebook applies the limma-voom method from limma-voom.R to Jerber2021 dataset.
scrna_dss.R :                 [Depricated] Script to apply dss method to Jerber2021 data.
