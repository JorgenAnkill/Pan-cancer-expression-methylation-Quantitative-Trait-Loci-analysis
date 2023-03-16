# A comprehensive analysis of expression-methylation Quantitative Trait Loci across 33 cancer types provides novel insight into the role of DNA methylation on tumor progression

***Jørgen Ankill, Zhi Zhao, Xavier Tekpli, Vessela N. Kristensen, Anthony Mathelier, Thomas Fleischer***

Contact: Jørgen Ankill, jorgen.ankill2@rr-research.no

#
**Description:**

This repository contains scripts related to the paper "An integrative pan-cancer analysis of DNA methylation and gene expression identifies a common link between aberrant enhancer methylation and proliferation in cancer". The code show how the pan-cancer expression-methylation Quantitative Trait Loci analysis was performed and the sources for the TCGA data.

In brief, we impute the missing methylation values in the TCGA-PANCAN dataset using knn.impute. After this we identify tumor samples with matching DNA methylation and gene expression data available. We then split the samples in each cancer type into two groups for the discovery and validation TCGA dataset. The samples for the discovery dataset is then combined for all cancer types, as for the validation dataset. We then perform a correlation analysis between DNA methylation (IQR>0.1) and gene expression (IQR>0) in the discovery dataset. Only correlations with an absolute correlation coefficient >0.5 and a bonferroni corrected p-value<0.05 are kept. These associations are then reanalyzed in the validation dataset. 

DNA methylation, gene expression, ATAC-seq and clinical data was downloaded from the Xenabrowser: https://xenabrowser.net/datapages/?cohort=GDC%20Pan-Cancer%20(PANCAN)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443 8th of May 2022 [1]. 

#
**References:**
  1. Goldman, M., et al., The UCSC Xena platform for public and private cancer genomics data visualization and interpretation. bioRxiv, 2019: p. 326470.
