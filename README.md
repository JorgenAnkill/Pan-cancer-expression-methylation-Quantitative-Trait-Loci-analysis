# Integrative pan-cancer analysis reveals a common architecture of dysregulated transcriptional networks characterized by loss of enhancer methylation

* *Jørgen Ankill, Zhi Zhao, Xavier Tekpli, Elin H. Kure, Vessela N. Kristensen, Anthony Mathelier, Thomas Fleischer

Contact: Jørgen Ankill, jorgen.ankill2@rr-research.no

#
**Description:**

This repository contains scripts related to the paper "An integrative pan-cancer analysis of DNA methylation and gene expression identifies a common link between aberrant enhancer methylation and proliferation in cancer". 

The code and associated data is a minimal working example for pan-cancer expression-methylation Quantitative Trait Loci (emQTL) analysis.

This code performs correlation analysis between two matrices, metMat and exprMat, which contain normalized expression data and normalized CpG methylation data from the same samples, respectively. It correlates all rows in these matrices, i.e., all genes to all CpGs, to identify significant CpG-gene associations known as emQTLs. The resulting associations are then grouped by bipartite network analysis using COmplex Network Desrciption Of Regulators (CONDOR) into emQTL communities. 

This example includes 100 random samples, 1000 random CpGs, and 100 random genes sourced from The Cancer Genome Atlas pan-cancer (TCGA-PANCAN) data, downloaded from the Xena browser. For more details on the data source, refer to: Goldman M et al. [1]. Data Source: https://xenabrowser.net/datapages/?cohort=GDC%20Pan-Cancer%20(PANCAN)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443.

**References:**
  1. Goldman, M., et al., The UCSC Xena platform for public and private cancer genomics data visualization and interpretation. bioRxiv, 2019: p. 326470.
