###################################################################################

# The code presented shows the work flow of the pan-cancer expression-methylation 
# Quantitative Trait Loci (emQTL) analysis:

# Title:
# Integrative pan-cancer analysis reveals a common architecture of dysregulated 
# transcriptional networks characterized by loss of enhancer methylation

# Authors:
# Jørgen Ankill, Zhi Zhao, Xavier Tekpli, Elin H. Kure, Vessela N. Kristensen, 
# Anthony Mathelier, Thomas Fleischer.

# Contact:
# Jørgen Ankill (2023-September)
# jorgen.ankill2@rr-research.no

# This code performs correlation analysis between two matrices, metMat and exprMat,
# which contain normalized expression data and normalized CpG methylation data
# from the same samples, respectively. It correlates all rows in these matrices,
# i.e., all genes to all CpGs, to identify significant CpG-gene associations known
# as emQTLs. The resulting associations are then grouped by bipartite network
# analysis using COmplex Network Desrciption Of Regulators into emQTL communities. 

# This example includes 100 random samples, 1000 random CpGs, and 100 random genes
# sourced from The Cancer Genome Atlas pan-cancer (TCGA-PANCAN) data, downloaded
# from the Xena browser. For more details on the data source, refer to:
# Goldman M et al. "The UCSC Xena platform for public and private cancer genomics
# data visualization and interpretation." bioRxiv 2019:326470.
# Data Source: https://xenabrowser.net/datapages/?cohort=GDC%20Pan-Cancer%20(PANCAN)
# &removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443.

###################################################################################
# -- Dependencies
###################################################################################

library(fst)
library(condor)
library(igraph)
library(pkgload)
library(devtools)

###################################################################################
# -- Perform pan-cancer emQTL analysis
###################################################################################

# Load data
mm <- read.table("HumanMethylation450_TCGA-PANCAN_100x1000.txt", sep = "\t")
me <- read.table("HTSeq-FPKM-UQ_TCGA-PANCAN_100x100.txt", sep = "\t")

# Match columns
mm <- mm[, match(colnames(me), colnames(mm))]

# Remove genes with no variance in expression level across the samples (IQR > 0.1)
iqr <- apply(me, 1, IQR)
me <- me[iqr > 0, ]

# Remove CpGs with little or no variance in methylation level across the samples (IQR >= 0.1)
iqr <- apply(mm, 1, IQR)
mm <- mm[iqr >= 0.1, ]

# Correlate DNA methylation and gene expression by Pearson correlation
r <- cor(t(me), t(mm), method = "pearson")

# Calculate p-values from correlation coefficients
n <- ncol(me)
t <- r * sqrt(n - 2)
t <- t / sqrt(1 - r^2)
df <- n - 2
p <- 2 * pt(abs(t), df, lower.tail = FALSE)

# Filter CpGs on the sex chromosomes
probeinfo <- read.table("Probeinfo.txt", sep = "\t")
probeinfo <- probeinfo[!probeinfo$Chr %in% c("X", "Y"), ]

# Filter correlation data based on probeinfo
r <- r[, colnames(r) %in% probeinfo$Probe]
p <- p[, colnames(p) %in% probeinfo$Probe]

# Deine significant cutoff
m <- p
pcut <- 0.05/(as.numeric(nrow(m))*as.numeric(ncol(m)))

# Define correlation cutoff
rcut <- 0.5

# Apply p-value and correlation cutoffs
p[p>pcut] <- NA
p[p<pcut] <- 0
tmp <- p + r

# Calculate absolute values and apply correlation cutoff
tmp <- abs(tmp)
tmp[tmp<rcut] <- NA
tmp[!is.na(tmp)] <- 1

# Calculate row and column sums
rowsums <- rowSums(tmp, na.rm = TRUE)
colsums <- colSums(tmp, na.rm = TRUE)

# Filter data
tmp <- tmp[!rowsums == 0, !colsums == 0]

# Filter the original p-values and associations data
m <- m[rownames(m) %in% rownames(tmp), colnames(m) %in% colnames(tmp)]

# Filter rows and columns in the 'r' matrix based on common row and column names with 'm'
r <- r[rownames(r) %in% rownames(m), colnames(r) %in% colnames(m)]

# Set elements in 'tmp' matrix with a value of 1 to 0
tmp[tmp == 1] <- 0

# Add the 'tmp' matrix values to the 'r' matrix
r <- r + tmp

# Reshape the 'r' matrix into a long format using the 'melt' function
r <- reshape2::melt(r)

# Remove rows with missing ('NA') values
r <- r[!is.na(r$value),]

# Create a new column 'emQTL_ID' by concatenating 'Var2' and 'Var1' columns
r$emQTL_ID <- paste0(r$Var2, "_", r$Var1)

# Reshape the 'tmp' matrix into a long format using the 'melt' function
tmp <- reshape2::melt(tmp)

# Remove rows with missing ('NA') values
tmp <- tmp[!is.na(tmp$value),]

# Create a new column 'emQTL_ID' by concatenating 'Var2' and 'Var1' columns
tmp$emQTL_ID <- paste0(tmp$Var2, "_", tmp$Var1)

# Merge 'tmp' and 'r' data frames based on the 'emQTL_ID' column
tmp <- merge(tmp, r, by = "emQTL_ID")

# Select and retain specific columns in the merged data frame
tmp <- tmp[, colnames(tmp) %in% c("emQTL_ID", "Var1.x", "Var2.x", "value.y")]

# Rename the columns for clarity
colnames(tmp) <- c("emQTL_ID", "Gene", "Probe", "rvalue")

###################################################################################
# -- emQTL community detection by bipartite network analysis using CONDOR 
###################################################################################

# Load the list of validated emQTLs after r-value filtering
elist <- tmp

# Filter emQTLs with negative r-values
elist <- elist[elist$rvalue < 0, ]
elist$rvalue <- abs(elist$rvalue)

# Select relevant columns
elist <- elist[, c("Probe", "Gene", "rvalue")]
reds <- elist$Probe
blues <- elist$Gene

# Create a CONDOR object
condor.object <- create.condor.object(elist)

# Set a random seed
set.seed(42)

# Cluster emQTLs using CONDOR
condor.object <- condor.cluster(condor.object, cs.method = "LCS")

# Make dataframes with CpGs and genes with community annotations
geneinfo <- condor.object$blue.memb
probeinfo <- condor.object$red.memb

###############################################################################
sessionInfo()































