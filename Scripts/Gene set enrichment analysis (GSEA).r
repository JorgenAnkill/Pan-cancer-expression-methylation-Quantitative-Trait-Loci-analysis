######################################################################################################################################################################################################################################################################
# -- GSEA using gene sets from MSigDB
######################################################################################################################################################################################################################################################################

GSEA <- function(genes,geneset.collection,pcut,geneset.threshold,version,background,minGenes,include,plot.output,plot.output.name,top.genesets.to.show,top.genesets.to.show.in.plot,plot.title){
if(missing(minGenes)){minGenes <- 0} 
if(missing(pcut)){pcut <- 0.05}
if(missing(geneset.threshold)){geneset.threshold <- 50}
if(missing(geneset.collection)){geneset.collection <- c("H","C1","C2","C4","C5","C6","C7","C8")}
if(missing(version)){version <- "v7.4"}
if(missing(plot.output)){plot.output <- FALSE}
if(missing(top.genesets.to.show.in.plot)){top.genesets.to.show.in.plot <- 10}
if(missing(plot.title)){plot.title <- ""}

if(version=="v7.0"){genesets <- readRDS("/my_directory/GSEA/Molecular signatures database/v7.0/msigdb_v7.0_GMTs/MSigDB_GSVA_data.rds")}
if(version=="v7.1"){genesets <- readRDS("/my_directory/GSEA/Molecular signatures database/v7.1/msigdb_v7.1_GMTs/MSigDB_GSVA_data.rds")}
if(version=="v7.4"){genesets <- readRDS("/my_directory/GSEA/MSigDB genesets/v7.4/msigdb_v7.4_GMTs/MSigDB_GSVA_data.rds")}
if(version=="v7.5.1"){genesets <- readRDS("/my_directory/GSEA/Molecular signatures database/v7.5.1/msigdb_v7.5.1_GMTs/MSigDB_GSVA_data.rds")}

# -- With unspecified background
if(missing(background)){
total.genes <- length(unique(genesets$Gene))
genesets <- genesets[genesets$GeneSetCollection%in%geneset.collection,]
genes.length <- sum(unique(genesets$Gene)%in%genes)

if(!missing(include)){genesets <- genesets[genesets$GeneSetName%in%genesets[genesets$GeneSetName%in%unique(grep(paste(include,collapse="|"),genesets$GeneSetName,value=TRUE)),]$GeneSetName,]} #HP, GOBP, GOMF, GOCC

res <- data.frame(matrix(nrow=length(unique(genesets$GeneSetName)),ncol=5))
colnames(res) <- c("GeneSetName","GenesInGeneset","GenesInOverlap","p.value","FDR.q.value")
res$GeneSetName <- unique(genesets$GeneSetName)

dup <- genesets[!duplicated(genesets$GeneSetName),]
res$GenesInGeneset <- dup$GenesInGeneset

out <- split(genesets,f=genesets$GeneSetName); out <- out[match(res$GeneSetName,names(out))]
vec <- vector()
for(i in 1:nrow(res)){vec[i] <- length(intersect(out[[i]]$Gene,genes))}

res$GenesInOverlap <- vec

# N: background total (#total balls in urn)	
# q: test success (#white balls drawn)			   
# m: bacground success (#white balls in urn)	
# n=N-m: background failure (#black balls in urn)	
# k: test total (#total balls drawn)	 			  

pvals <- vector()
fe <- vector()
N <- total.genes
k <- genes.length
for(i in 1:nrow(res)){
q <- res$GenesInOverlap[i]
m <- res$GenesInGeneset[i]
n <- N-m
pvals[i] <- phyper(q-1,m,n,k,lower.tail=FALSE)
fe[i] <- q*N/(m*k)}

res$p.value <- pvals
res$FE <- fe
res$FDR.q.value <- p.adjust(pvals,method="fdr")
res <- res[order(res$p.value,decreasing=FALSE),]
res <- res[res$GenesInGeneset>=geneset.threshold,]
result <- res[res$FDR.q.value<=pcut & res$GenesInOverlap>=minGenes,]
print(noquote("GSEA completed"))
}

# -- With specified background
if(!missing(background)){
total.genes <- length(unique(background))
genes.length <- sum(unique(background)%in%genes)

genesets <- genesets[genesets$GeneSetCollection%in%geneset.collection,]

res <- data.frame(matrix(nrow=length(unique(genesets$GeneSetName)),ncol=5))
colnames(res) <- c("GeneSetName","GenesInGeneset","GenesInOverlap","p.value","FDR.q.value")
res$GeneSetName <- unique(genesets$GeneSetName)

vecgenesetdim <- NULL
for(currgeneset in res$GeneSetName){
  vecgenesetdim <- c(vecgenesetdim, 
                     sum(unique(background)%in%genesets$Gene[which(genesets$GeneSetName==currgeneset)]) )
}
res$GenesInGeneset <- vecgenesetdim

out <- split(genesets,f=genesets$GeneSetName); out <- out[match(res$GeneSetName,names(out))]
vec <- vector()
for(i in 1:nrow(res)){vec[i] <- length(intersect(intersect(out[[i]]$Gene,background),genes))}

res$GenesInOverlap <- vec

pvals <- vector()
fe <- vector()
N <- total.genes
k <- genes.length
for(i in 1:nrow(res)){
q <- res$GenesInOverlap[i]
m <- res$GenesInGeneset[i]
n <- N-m
pvals[i] <- phyper(q-1,m,n,k,lower.tail=FALSE)
fe[i] <- q*N/(m*k)}

res$p.value <- pvals
res$FE <- fe
res$FDR.q.value <- p.adjust(pvals,method="fdr")
res <- res[order(res$p.value,decreasing=FALSE),]
res <- res[res$GenesInGeneset>=geneset.threshold,]
result <- res[res$FDR.q.value<=pcut & res$GenesInOverlap>=minGenes,]
print(noquote("GSEA completed"))
}

if(sum(plot.output)==1){

    # -- Plot biclusters	
    png(file = paste0("GSEA barplot ",plot.title,".png"),width=2000,height=600)
    par(mar=c(8,75,8,6))

    res <- result
    res$Gene_Ontology <- NA
    res$Gene_Ontology <- lapply(strsplit(res$GeneSetName,split="_"),function(l) l[[1]])

    annot <- strsplit(res$GeneSetName,split="_")
    for(i in 1:length(strsplit(res$GeneSetName,split="_"))){
    annot[[i]] <- annot[[i]][-1]
    res$GeneSetName[i] <- paste(annot[[i]],collapse="_")}

    # Annotation of GO terms
    res$Gene_Ontology[res$Gene_Ontology=="GOBP"] <- "GO (Biological process)" ; res$Gene_Ontology[res$Gene_Ontology=="GOMF"] <- "GO (Molecular function)" ; res$Gene_Ontology[res$Gene_Ontology=="GOCC"] <- "GO (Cellular component)" ;res$Gene_Ontology[res$Gene_Ontology=="HALLMARK"] <- "H (Hallmark geneset)" 

    if(nrow(res)>=top.genesets.to.show.in.plot){res <- res[1:top.genesets.to.show.in.plot,]}
    res <- res[order(res$p.value,decreasing=TRUE),]

    res$GeneSetName <- factor(res$GeneSetName,levels=unique(res$GeneSetName),ordered=TRUE)
    res$FDR.q.value <- -log10(res$FDR.q.value)

    # -- Add colors by GO term
    mycols <- data.frame(matrix(nrow=nrow(res),ncol=1));colnames(mycols) <- "Color";rownames(mycols) <- rownames(res)
    for(i in 1:top.genesets.to.show.in.plot){if(res$Gene_Ontology[i]=="GO (Biological process)"){mycols$Color[i] <- "orange"}; if(res$Gene_Ontology[i]=="H (Hallmark geneset)"){mycols$Color[i] <- "red2"};if(res$Gene_Ontology[i]=="GO (Molecular function)"){mycols$Color[i] <- "limegreen"};if(res$Gene_Ontology[i]=="GO (Cellular component)"){mycols$Color[i] <- "royalblue"}}

    x1 <- barplot(res$FDR.q.value,
        horiz	= TRUE,					
        names.arg=res$GeneSetName,
        las=1,
        col=mycols$Color,
        cex.lab=4,
        cex.names=2.6,
        cex.axis=2.8,
        xlim=c(0,max(res$FDR.q.value)*5/3))
    abline(v=0)
    title(paste0(plot.title),adj=0,cex.main=4)
    text(x=0.4,y=x1,paste0("FE=",signif(res$FE,digits=2)),pos=4,cex=2.5)
    text(cex=3,x=max(res$FDR.q.value)*4/3/2,y=-2,expression(paste("-log"["10"]*"(",italic("p"),"-value)")),srt=0,xpd=TRUE,pos=3)
    legend("topright",c("(H)  Hallmark genesets","GO (Biological process)","GO (Molecular function)","GO (Cellular component)"),fill=c("red2","orange","limegreen","royalblue"),cex=1.8)
   print(noquote("Plot completed, saved at work directory if not defined"))
   dev.off()
   }
return(result)
}
