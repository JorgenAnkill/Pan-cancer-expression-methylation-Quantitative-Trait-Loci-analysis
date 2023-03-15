######################################################################################################################################################################################################################################################################
# -- Divide the tumor samples from TCGA into a discovery and validation dataset
######################################################################################################################################################################################################################################################################

lapply(c("svMisc"),require,character.only=TRUE)

setwd("/my_directory/") # Directory for the imputed methylation data from the script "Script for imputing missing methylation values part 2"
temp <- list.files(pattern="*.rds")

sampleinfo <- readRDS("/my_directory/GDC-PANCAN.basic_phenotype.rds")   # Info about the TCGA samples obtained from: https://xenabrowser.net/datapages/?dataset=GDC-PANCAN.basic_phenotype.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

discovery.samples <- vector()
validation.samples <- vector()
common.probes <- vector()

set.seed(42)
for(i in 1:length(temp)){
temporary <- readRDS(temp[i]); dat <- temporary

# Find samples with expression and methylation data
met.samples <- colnames(dat)
expr.samples <- readLines("/my_directory/GDC-PANCAN.htseq_fpkm-uq.tsv",n=1) # Gene expression data from TCGA obtained from: https://xenabrowser.net/datapages/?dataset=GDC-PANCAN.htseq_fpkm-uq.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
expr.samples <- unlist(strsplit(expr.samples,"\t"))[-1]
met.expr.samples <- intersect(met.samples,expr.samples)

# Obtain tumor samples from TCGA
tmp <- sampleinfo
tmp <- tmp[tmp$program=="TCGA" & tmp$sample_type%in%c("Primary Tumor","Primary Blood Derived Cancer - Peripheral Blood","Primary Blood Derived Cancer - Bone Marrow") & !is.na(tmp$Cancer_type) & tmp$Cancer_type !="" & tmp$Cancer_type==strsplit(temp[i],".rds")[[1]][1],]
tmp <- tmp[rownames(tmp)%in%intersect(met.expr.samples,rownames(tmp)),]
tmp$short_id <- substring(rownames(tmp),1,12)

# Remove duplicates by only keeping first vial
if(nrow(tmp[duplicated(tmp$short_id),])>=1){

dup1 <- tmp[duplicated(tmp$short_id),]
dup <- tmp[tmp$short_id%in%dup1$short_id,]
dup$vial <- sub(".*-01","",rownames(dup))
tmp <- tmp[!rownames(tmp)%in%rownames(dup),]

for(g in 1:nrow(dup1)){
te <- dup[dup$short_id%in%dup1$short_id[g],]
te$vial[te$vial=="A"] <- 1 
te$vial[te$vial=="B"] <- 2
te$vial[te$vial=="C"] <- 3
te$vial[te$vial=="D"] <- 4
te$vial[te$vial=="E"] <- 5
te$vial[te$vial=="F"] <- 6
te$vial[te$vial=="G"] <- 7
te$vial[te$vial=="H"] <- 8
te$vial[te$vial=="I"] <- 9
te$vial[te$vial=="J"] <- 10
te$vial[te$vial=="K"] <- 11
te$vial[te$vial=="L"] <- 12
te$vial[te$vial=="M"] <- 13
te$vial[te$vial=="N"] <- 14
te$vial[te$vial=="O"] <- 15
te$vial[te$vial=="P"] <- 16
te$vial[te$vial=="Q"] <- 17
te$vial[te$vial=="R"] <- 18
te$vial[te$vial=="S"] <- 19
te$vial[te$vial=="T"] <- 20
te$vial[te$vial=="U"] <- 21
te$vial[te$vial=="V"] <- 22
te$vial[te$vial=="W"] <- 23
te$vial[te$vial=="X"] <- 24
te$vial[te$vial=="Y"] <- 25
te$vial[te$vial=="Z"] <- 26
te <- te[order(te$vial,decreasing=FALSE),]
te <- te[1,!colnames(te)%in%"vial",drop=FALSE]
tmp <- rbind(tmp,te)
}}

if(!(nrow(tmp)/2)%%1==0){
random.samples.discovery <- sample(rownames(tmp),((nrow(tmp)-1)/2)+1,replace=FALSE)
random.samples.validation <- sample(rownames(tmp)[!rownames(tmp)%in%random.samples.discovery],(nrow(tmp)-1)/2,replace=FALSE)}

if((nrow(tmp)/2)%%1==0){
random.samples.discovery <- sample(rownames(tmp),nrow(tmp)/2,replace=FALSE)
random.samples.validation <- sample(rownames(tmp)[!rownames(tmp)%in%random.samples.discovery],nrow(tmp)/2,replace=FALSE)}

discovery.samples <- append(discovery.samples,random.samples.discovery)
validation.samples <- append(validation.samples,random.samples.validation)

dat <- dat[,colnames(dat)%in%c(random.samples.discovery,random.samples.validation)]

if(i==1){common.probes <- append(common.probes,rownames(dat))}
if(!i==1){common.probes <- intersect(common.probes,rownames(dat))}

progress(i,max.value=length(temp))
Sys.sleep(0.01)
if(i == length(temp)) cat("Done!\n")}

write.table(discovery.samples,paste("/my_directory/Discovery_cohort_samples.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(validation.samples,paste("/my_directory/Validation_cohort_samples.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(common.probes,paste("/my_directory/Common probes for all cancer types.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)

######################################################################################################################################################################################################################################################################
# -- Make a master file with methylation data for the validation and discovery dataset
######################################################################################################################################################################################################################################################################

setwd("/my_directory/") # Directory for the imputed methylation data from the script "Script for imputing missing methylation values part 2"
temp <- list.files(pattern="*.rds")

#Discovery
template <- readRDS(temp[1])
disc <- template[rownames(template)%in%common.probes,colnames(template)%in%discovery.samples]
val <- template[rownames(template)%in%common.probes,colnames(template)%in%validation.samples]

for(i in 2:length(temp)){
temporary <- readRDS(temp[i])

#Discovery
disc.tmp <- temporary[rownames(temporary)%in%common.probes,colnames(temporary)%in%discovery.samples]
disc.tmp <- disc.tmp[match(rownames(disc),rownames(disc.tmp)),]
disc <- cbind(disc,disc.tmp)

#Validation
val.tmp <- temporary[rownames(temporary)%in%common.probes,colnames(temporary)%in%validation.samples]
val.tmp <- val.tmp[match(rownames(val),rownames(val.tmp)),]
val <- cbind(val,val.tmp)

print(paste(i," of 33 completed!",sep=""))}

saveRDS(disc,paste("/my_directory/450K data discovery cohort.rds",sep=""))
saveRDS(val,paste("/my_directory/450K data validation cohort.rds",sep=""))

######################################################################################################################################################################################################################################################################
# -- Make a master file for DNA methylation data from normal samples
######################################################################################################################################################################################################################################################################

lapply(c("svMisc"),require,character.only=TRUE)

setwd("/my_directory/") # Directory for the imputed methylation data from the script "Script for imputing missing methylation values part 2"
temp <- list.files(pattern="*.rds")

sampleinfo <- readRDS("/my_directory/GDC-PANCAN.basic_phenotype.rds")   # Info about the TCGA samples obtained from: https://xenabrowser.net/datapages/?dataset=GDC-PANCAN.basic_phenotype.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
sampleinfo <- sampleinfo[sampleinfo$program%in%"TCGA" & sampleinfo$sample_type%in%c("Solid Tissue Normal","Blood Derived Normal","Buccal Cell Normal","Bone Marrow Normal"),]

temporary <- readRDS(temp[2]); dat <- temporary
template <- dat[,colnames(dat)%in%rownames(sampleinfo)]

for(i in 1:length(temp)){

temporary <- readRDS(temp[i]); dat <- temporary
if(!ncol(dat[,colnames(dat)%in%rownames(sampleinfo)])==0){
if(!i==2){dat <- dat[,colnames(dat)%in%rownames(sampleinfo)]
template <- transform(merge(template,dat,by="row.names",all=TRUE),row.names=Row.names,Row.names=NULL)}}

progress(i,max.value=length(temp))
Sys.sleep(0.01)
if(i == length(temp)) cat("Done!\n")
}

saveRDS(template,"/my_directory/PANCAN.Normal.sample.data.rds")

######################################################################################################################################################################################################################################################################
# -- Preparation of the sampleinfo file
######################################################################################################################################################################################################################################################################

sampleinfo <- readRDS("/my_directory/GDC-PANCAN.basic_phenotype.rds")   # Info about the TCGA samples obtained from: https://xenabrowser.net/datapages/?dataset=GDC-PANCAN.basic_phenotype.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
sampleinfo <- sampleinfo[sampleinfo$program=="TCGA",colnames(sampleinfo)%in%c("sample_type","Cancer_type","ER")]

# Add InfiniumPurify tumor purities to sampleinfo
setwd("/my_directory/")
temp <- list.files(pattern="*.txt")

template <- read.table(temp[1],header=TRUE)
if(ncol(template)==2){template$Purity_ABSOLUTE <- "-"}

for(i in 2:length(temp)){
temporary <- read.table(temp[i],header=TRUE); dat <- temporary
if(ncol(temporary)==2){temporary$Purity_ABSOLUTE <- "-"}
template <- rbind(template,temporary)}

template$short_ID <- substring(template$SampleName,1,16)

# Average duplicates
non.dup <- template[!duplicated(template$short_ID),]

dup <- template[duplicated(template$short_ID),]$short_ID

tmp <- template[template$short_ID%in%dup[1],]
tmp$Purity_InfiniumPurify <- mean(tmp$Purity_InfiniumPurify)
tmp$Purity_ABSOLUTE <- mean(as.numeric(tmp$Purity_ABSOLUTE))
tmp <- tmp[1,,drop=FALSE]

for(i in 2:length(dup)){
te <- template[template$short_ID%in%dup[i],];te
te$Purity_InfiniumPurify <- mean(te$Purity_InfiniumPurify)
if(!sum(te$Purity_ABSOLUTE=="-")>1){te$Purity_ABSOLUTE <- mean(as.numeric(te$Purity_ABSOLUTE))}
te <- te[1,,drop=FALSE]
tmp <- rbind(te,tmp)}

template <- rbind(tmp,non.dup)

sampleinfo <- merge(sampleinfo,template,by.x="row.names",by.y="short_ID",all=TRUE)
sampleinfo$Purity_ABSOLUTE[sampleinfo$Purity_ABSOLUTE=="-"] <- NA
sampleinfo$Purity_InfiniumPurify[sampleinfo$Purity_InfiniumPurify=="-"] <- NA

saveRDS(sampleinfo,"/my_directory/Sampleinfo_TCGA.rds")
sampleinfo <- readRDS("/my_directory/Sampleinfo_TCGA.rds")

# Make an overview of all cancer types
tmp <- data.frame(c("LAML","ACC","LUSC","CESC","SKCM","PCPG","COAD","HNSC","ESCA","UCS","READ","ER+","ER-","STAD","THCA","THYM","PAAD","KIRC","LUAD","MESO","TGCT","PRAD","KICH","LIHC","DLBC","UCEC","KIRP","BLCA","LGG","GBM","SARC","UVM","CHOL","OV","BRCA"))
colnames(tmp) <- "Cancer_type";rownames(tmp) <- tmp$Cancer_type
tmp$Description <- tmp$Cancer_type

tmp$Description[tmp$Description=="LAML"] <- "Acute Myeloid Leukemia"
tmp$Description[tmp$Description=="ACC"] <- "Adrenocortical Cancer"
tmp$Description[tmp$Description=="LUSC"] <- "Lung Squamous Cell Carcinoma"
tmp$Description[tmp$Description=="CESC"] <- "Cervical Cancer"
tmp$Description[tmp$Description=="SKCM"] <- "Melanoma"
tmp$Description[tmp$Description=="PCPG"] <- "Pheochromocytoma & Paraganglioma"
tmp$Description[tmp$Description=="COAD"] <- "Colon Cancer"
tmp$Description[tmp$Description=="HNSC"] <- "Head and Neck Cancer"
tmp$Description[tmp$Description=="ESCA"] <- "Esophageal Cancer"
tmp$Description[tmp$Description=="UCS"] <- "Uterine Carcinosarcoma"
tmp$Description[tmp$Description=="READ"] <- "Rectal Cancer"
tmp$Description[tmp$Description=="STAD"] <- "Stomach Cancer"
tmp$Description[tmp$Description=="THCA"] <- "Thyroid Cancer"
tmp$Description[tmp$Description=="THYM"] <- "Thymoma"
tmp$Description[tmp$Description=="PAAD"] <- "Pancreatic Cancer"
tmp$Description[tmp$Description=="KIRC"] <- "Kidney Clear Cell Carcinoma"
tmp$Description[tmp$Description=="LUAD"] <- "Lung Adenocarcinoma"
tmp$Description[tmp$Description=="MESO"] <- "Mesothelioma"
tmp$Description[tmp$Description=="TGCT"] <- "Testicular Cancer"
tmp$Description[tmp$Description=="PRAD"] <- "Prostate Cancer"
tmp$Description[tmp$Description=="KICH"] <- "Kidney Chromophobe"
tmp$Description[tmp$Description=="LIHC"] <- "Liver Cancer"
tmp$Description[tmp$Description=="DLBC"] <- "Large B-cell Lymphoma"
tmp$Description[tmp$Description=="UCEC"] <- "Endometroid Cancer"
tmp$Description[tmp$Description=="KIRP"] <- "Kidney Papillary Cell Carcinoma"
tmp$Description[tmp$Description=="BLCA"] <- "Bladder Cancer"
tmp$Description[tmp$Description=="LGG"] <- "Lower Grade Glioma"
tmp$Description[tmp$Description=="GBM"] <- "Glioblastoma"
tmp$Description[tmp$Description=="SARC"] <- "Sarcoma"
tmp$Description[tmp$Description=="UVM"] <- "Ocular melanomas"
tmp$Description[tmp$Description=="CHOL"] <- "Bile Duct Cancer"
tmp$Description[tmp$Description=="OV"] <- "Ovarian Cancer"
tmp$Description[tmp$Description=="BRCA"] <- "Breast cancer"

te <- tmp
te$Cancer_type <- paste(te$Cancer_type,"_Normal",sep="")
tmp <- rbind(tmp,te)

saveRDS(tmp,"/my_directory/TCGA cancer types overview.rds")

######################################################################################################################################################################################################################################################################
# -- Pan-cancer emQTL -- Discovery cohort
######################################################################################################################################################################################################################################################################

mm <- readRDS("/my_directory/450K data discovery cohort.rds")
me <- readRDS("/my_directory/PANCAN.Discovery.gene_expression.protein.coding.rds")
mm <- mm[,match(colnames(me),colnames(mm))]

#Remove genes with no variance in expression level across the samples(IQR>0.1)
iqr <- apply(me,1,IQR)
me  <- me[iqr>0,]; nrow(me)

#Remove CpGs with little or no variance in methylation level across the samples(IQR>0.1)
iqr <- apply(mm,1,IQR)
mm  <- mm[iqr>=0.1,]; nrow(mm)

#Calculate the correlation coefficients
r <- cor(t(me),t(mm))
saveRDS(r,"/my_directory/emQTL_PANCAN_discovery_r_uncorrected.rds")

#Calculate p-values from correlation coefficients
n <- ncol(me)
t <- r*sqrt(n-2)
t = t/sqrt(1-r^2)
df <- n-2
p <- 2*pt(abs(t), df, lower.tail=F)
saveRDS(p,"/my_directory/emQTL_PANCAN_discovery_p_uncorrected.rds")

# Filter CpGs on the sex chromosomes
probeinfo <- read.table("/my_directory/Probeinfo2017.txt",header=T,sep="\t",na.strings="NA",row.names=NULL,stringsAsFactors=F)
probeinfo <- probeinfo[!probeinfo$Chr%in%c("X","Y"),]

r <- r[,colnames(r)%in%probeinfo$Probe]
p <- p[,colnames(p)%in%probeinfo$Probe]

#Keep all CpGs and genes with at least one significant CpG-gene association
m <- p
pcut <- 0.05/(as.numeric(nrow(m))*as.numeric(ncol(m)))
saveRDS(pcut,"/my_directory/emQTL_PANCAN_discovery_pcut.rds")

rcut <- 0.5
tmp.r <- r
tmp.p <- p

tmp.p[tmp.p>pcut] <- NA
tmp.p[tmp.p<pcut] <- 0; nrow(tmp.p)==nrow(tmp.r) ; sum(rownames(tmp.p)==rownames(tmp.r))==nrow(tmp.p)

tmp <- tmp.p + tmp.r

tmp <- abs(tmp)
tmp[tmp<rcut] <- NA
tmp[!is.na(tmp)] <- 1

rowsums <- rowSums(tmp,na.rm=TRUE)
colsums <- colSums(tmp,na.rm=TRUE)

tmp <- tmp[!rowsums==0,!colsums==0]; dim(tmp); sum(!is.na(tmp))

write.table(rownames(tmp),"/my_directory/Unvalidated_significant_genes_PANCAN_discovery.txt",sep="\t")
write.table(colnames(tmp),"/my_directory/Unvalidated_significant_CpGs_PANCAN_discovery.txt",sep="\t")

m <- m[rownames(m)%in%rownames(tmp),colnames(m)%in%colnames(tmp)]
saveRDS(m,"/my_directory/emQTL_PANCAN_p_values.rds")
saveRDS(tmp,"/my_directory/emQTL_PANCAN_associations_binary.rds")

######################################################################################################################################################################################################################################################################
# -- Pan-cancer emQTL -- Validation cohort (Reanalyzing the significant associations in the validation cohort)
######################################################################################################################################################################################################################################################################

mm <- readRDS("/my_directory/450K data validation cohort.rds")  # 450k methylation data
me <- readRDS("/my_directory/PANCAN.Validation.gene_expression.protein.coding.rds") # Expression data of protein coding genes
mm <- mm[,match(colnames(me),colnames(mm))] # Find samples with both methylation and expression data available

# Remove genes with no variance in expression level across the samples(IQR>0)
iqr <- apply(me,1,IQR)
me  <- me[iqr>0,]

# Remove CpGs with little or no variance in methylation level across the samples(IQR>0)
iqr <- apply(mm,1,IQR)
mm  <- mm[iqr>0,]

genes <- read.table("/my_directory/Unvalidated_significant_genes_PANCAN_discovery.txt",sep="\t")
probes <- read.table("/my_directory/Unvalidated_significant_CpGs_PANCAN_discovery.txt",sep="\t")

me <- me[rownames(me)%in%genes$x,]
mm <- mm[rownames(mm)%in%probes$x,]

#Calculate the correlation coefficients
r <- cor(t(me),t(mm))
saveRDS(r,"/my_directory/emQTL_PANCAN_validation_r.rds")

# Adjust p-values for the number of tests - reanalyzis of significant correlations from emQTL analysis in discovery dataset
n_tests <- readRDS("/my_directory/emQTL_PANCAN_associations_binary.rds")
n_tests <- n_tests[rownames(n_tests)%in%rownames(me),colnames(n_tests)%in%rownames(mm)]
n_tests <- sum(!is.na(n_tests)); n_tests

#Calculate p-values from correlation coefficients
n <- ncol(me)
t <- r*sqrt(n-2)
t = t/sqrt(1-r^2)
df <- n-2
p <- 2*pt(abs(t), df, lower.tail=F)

#Keep all CpGs and genes with at least one significant CpG-gene association
m <- p
saveRDS(p,"/my_directory/emQTL_PANCAN_validation_p_uncorrected.rds")

pcut <- 0.05/(n_tests)
saveRDS(pcut,"/my_directory/emQTL_PANCAN_validation_pcut.rds")
pcut <- readRDS("/my_directory/emQTL_PANCAN_validation_pcut.rds")

rcut <- 0.5

tmp.r <- r
tmp.p <- p

tmp.p[tmp.p>pcut] <- NA
tmp.p[tmp.p<pcut] <- 0

tmp <- tmp.p + tmp.r

tmp <- abs(tmp)
tmp[tmp<rcut] <- NA
tmp[!is.na(tmp)] <- 1

# emQTLs
disc.assos <- readRDS("/my_directory/emQTL_PANCAN_associations_binary.rds"); disc <- disc.assos

disc <- disc[rownames(disc)%in%rownames(tmp),colnames(disc)%in%colnames(tmp)];rownames(disc)==rownames(tmp); colnames(tmp)%in%colnames(disc)

temp <- tmp + disc
temp[temp==1] <- NA
temp[temp==2] <- 1

rowsums <- rowSums(temp,na.rm=TRUE)
colsums <- colSums(temp,na.rm=TRUE)

# Validated associations
tmp <- temp[!rowsums==0,!colsums==0]; dim(tmp);sum(!is.na(tmp))

#Make new dataframes consisting only of CpGs and genes with validated associations
tmp.r <- r[rownames(tmp.r)%in%rownames(tmp),colnames(tmp.r)%in%colnames(tmp)]
tmp.p <- p[rownames(tmp.p)%in%rownames(tmp),colnames(tmp.p)%in%colnames(tmp)]

saveRDS(tmp.p,"/my_directory/emQTL_PANCAN_p.rds")
saveRDS(tmp.r,"/my_directory/emQTL_PANCAN_validation_r.rds")
saveRDS(tmp,"/my_directory/emQTL_PANCAN_validated_associations_binary.rds")

write.table(colnames(tmp),"/my_directory/validated_CpGs.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)  # CpGs with validated gene associations
write.table(rownames(tmp),"/my_directory/validated_genes.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE) # Genes with validated CpG associations
