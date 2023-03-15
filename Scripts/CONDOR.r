######################################################################################################################################################################################################################################################################
# -- CONDOR
######################################################################################################################################################################################################################################################################5

#install.packages("pkgload",dependencies=TRUE)
#install.packages("devtools",dependencies=TRUE)
#devtools::install_github("jplatig/condor",dependencies=TRUE)

# Parameters
output_dir <- paste0("/my_directory/")

lapply(c("devtools","pkgload","igraph","condor","fst"),require,character.only=TRUE)
elist <- read_fst("/my_directory/List of validated emQTLs with r and p-values.fst")
elist <- elist[elist$rvalue<0,]
elist$rvalue <- abs(elist$rvalue)

elist <- elist[,c("Probe","Gene","rvalue")]

# Create a CONDOR object
condor.object <- create.condor.object(elist)
names(condor.object)

set.seed(42)
condor.object <- condor.cluster(condor.object,cs.method="LCS")
print(condor.object$red.memb)
print(condor.object$blue.memb)

table(condor.object$blue.memb$com)
table(condor.object$red.memb$com)

geneinfo <- condor.object$blue.memb # File with genes with annotations of which community they belongs to
probeinfo <- condor.object$red.memb # File with CpGs with annotations of which community they belongs to