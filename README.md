# CRISPR Barcoding Project

Welcome to the CRISPR Barcoding Project repository. This project focuses on the development and implementation of a novel CRISPR barcoding method to investigate dynamic transitions in single cells readout.

In this manuscript, we have developed a novel CRISPR barcoding method that  allows dynamic editing for over 12 days, and enables the profiling of barcodes and transcriptomes at the single-cell level. We believe this method marks a substantial improvement over current methodologies, expanding the potential applications of this technology. 

This repository contains the supplementary tables in the Sup_mat folder and the R code used for developing the project in the Code folder:
- 1-script_QC.R --> QC of seurat object
- 2-script_merge.R --> combine the data with the Barcode Library
- 3-script_Plots.R --> some plots for the manuscript.
- DEA_intersect_fisher.R --> Perform the Differential Expression Analysis 
- Filter_UMIs.R --> to apply the filters of quality control of UMIs
- function-readtable.R --> to extract the information about the CRISPR-Barcodes

