######################################################################################################################################################################################################################################################################
# -- DNA methylation data for PANCAN-TCGA
######################################################################################################################################################################################################################################################################

#Transpose 450K data for each cancer type in bash
cd ..
cd .. 
cd "my_directory"

ls *\450K.txt | xargs -I% ./fixme.sh %

# Impute missing methylation data in R

# Load required packages
lapply(c("readr","svMisc","impute"),require,character.only=TRUE)

setwd("/my_directory")
temp <- list.files(pattern="*450K.txt_transposed.txt")

for(i in 1:length(temp)){
temporary <- read_delim(temp[i],delim="\t",col_names=TRUE); 
dat <- temporary

dat <- data.frame(dat)
rownames(dat) <- dat$REF
dat <- dat[,-1]
colnames(dat) <- gsub(".","-",colnames(dat),fixed=TRUE)

# -- Remove probes with more than 50% missing values
dat <- dat[rowSums(is.na(dat))<=sum(ncol(dat)/2),]

# Impute missing methylation values
if(ncol(dat)<=10){m_impute <- impute.knn(as.matrix(dat),k=3,maxp=1500,rng.seed=362436069)}
if(ncol(dat)>10){m_impute <- impute.knn(as.matrix(dat),k=10,maxp=1500,rng.seed=362436069)}

saveRDS(m_impute,paste("/my_directory/",strsplit(temp[i],"_")[[1]][1],".rds",sep=""))

progress(i,max.value=length(temp))
Sys.sleep(0.01)
if(i == length(temp)) cat("Done!\n")}