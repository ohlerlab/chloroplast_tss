library(here)
setwd(here("/CSAR-Style/"))
library(DESeq2)
library("pheatmap")
require("gplots")
require(ggplot2)
library(memoise)

DESeq <- memoise(DESeq2::DESeq)

#Using parameters rather than magic numbers here
#this at least is easy to change
SAMPLESCORETHRESHOLD <- 10
FDR_THRESHOLD <- 0.01
CONTROLSCORETHRESHOLD <- 5
FCTHRESHvalues <- c(1,1.5,2)
# 
for(FDR_THRESHOLD in c(0.05,0.01)){
for(CONTROLSCORETHRESHOLD in c(5,10)){
#   
#This creates a lot of files in a lot of places with a lot of functions
#moving into a new working directory is probably the best option
setwd(here("/CSAR-Style/"))
getwd()
outdir <- paste0('run','fdr_',FDR_THRESHOLD,'_pseudocount_',CONTROLSCORETHRESHOLD)
dir.create(showWarnings = F,outdir) 
setwd(outdir)
library(tidyverse)
if('beds-Step2' %in% list.files(outdir)){
  list.files(outdir,full=T,rec=T)%>%file.remove
}

  
ids<-list.files(pattern=".counts$","../counts")
res<-read.table(paste("../counts/",ids[1],sep=""),header=T)
for(id in ids[-1]){
  # print(id)
  temp<-read.table(paste("../counts/",id,sep=""),header=T)[,7]
  res<-cbind(res,temp)
}


mid<-as.integer((res$Start+res$End)/2)
rrna<-sapply(1:length(res$Strand),function(x){if(res$Strand[x]=="+"){rrna<-c(101012:102502,104691:107500,107599:107701,107949:108069)}
  else{rrna<-c(130580:130700,130948:131050,131149:133958,136147:137637)};
  is.element(mid[x],rrna)})
res<-res[!rrna,]
coor<-res[,1:6]
counts<-res[,7:length(res[1,])]
colnames(counts)<-ids;rownames(counts)<-res$Geneid
ww<-which(rowSums(counts[,grep("sample",colnames(counts))]>SAMPLESCORETHRESHOLD)>0)
#coor<-coor[ww,];counts<-counts[ww,]#This I didn't use for Julia
for(i in grep("control",colnames(counts))){
  counts[counts[,i]<CONTROLSCORETHRESHOLD,i]<-CONTROLSCORETHRESHOLD
} ##THis I used for Julia
write.csv(counts,file="counts.csv")

####Annotate
mdis=5000
gff<-read.table(here("ext_data/Araport11_GFF3_genes_transposons.201606.gtf"),sep="\t")
gff<-gff[gff$V3=="gene" & gff$V1=="ChrC",]
dis<-lapply(coor$Start[coor$Strand=="+"],function(x){temp<- x-gff$V4[gff$V7=="+"];temp[abs(temp)<1000]})
# hist(unlist(dis),breaks=1000)
disp<-lapply(coor$Start[coor$Strand=="+"],function(x){temp<- x-gff$V4[gff$V7=="+"];names(temp)<-gff$V9[gff$V7=="+"];temp[temp==max(temp[(temp)< 0 & temp> -mdis ])]})
names(disp)<-coor$Geneid[coor$Strand=="+"]
disn<-lapply(coor$Start[coor$Strand=="-"],function(x){temp<- gff$V5[gff$V7=="-"]-x;names(temp)<-gff$V9[gff$V7=="-"];temp[temp==max(temp[(temp)< 0 & temp> -mdis ])]})
names(disn)<-coor$Geneid[coor$Strand=="-"]
temp<-sapply(as.character(coor$Geneid),function(x){temp<-unlist(c(disp[x],disn[x]));c(names(temp)," transcript_id None;")[1]})
target<-unlist(lapply(strsplit(temp,"transcript_id "),function(x)x[2]))
coor$target<-sub(";","",target)
coor$dis<-sapply(as.character(coor$Geneid),function(x){temp<-unlist(c(disp[x],disn[x]));c(temp,mdis)[1]})
##

geno<-as.factor(unlist(lapply(strsplit(ids,"_"),function(x)x[1])))
dev<-as.factor(unlist(lapply(strsplit(ids,"_"),function(x)x[2])))
rep<-as.factor(unlist(lapply(strsplit(ids,"control|sample"),function(x)x[2])))
cond<-as.factor(substr(unlist(lapply(strsplit(ids,"_|\\."),function(x)x[3])),1,3))

colData<-data.frame(cond=factor(cond),rep=as.factor(rep),geno=geno,dev=dev);rownames(colData)<-ids
save(counts,coor,file=file.path("Counts.RD"))


samp<-paste(geno,dev,sep="_")

for(i in unique(samp)){
  print(i)
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts[,samp==i]),colData = colData[samp==i,],design = ~ cond+rep)
  message(digest::digest(dds))
  
  dds1 <- DESeq(dds,fitType="local")
  
  res<-as.data.frame(results(dds1,contrast=c("cond","sam","con"),altHypothesis="greater"))
  res$pos<-coor$Start;res$pos2<-coor$End;res$strand<-coor$Strand;res$target<-coor$target;res$dis<-coor$dis
  #Julia wanted to filter for the number of control samples bigger than the sample
  #If I an interpretting this correctly then this occurs at the level of DESeq - there's no explicit setting
  #for how many samples are less than or greater
  #I'll therefore just use a different FDR threshold
  deseq_sig <- (res$padj<FDR_THRESHOLD & !is.na(res$padj) )
  sig<-res[deseq_sig,7:10];
  sig$pos<-sig$pos## It is already in0-based coordinates
  sig$chr<-"ChrC";sig$point<-".";sig<-sig[,c(5,1,2,4,6,3)]
  dir.create(file.path(("beds-Step2/")))
  
  write.table(sig,file=file.path(paste("beds-Step2/",i,"-SigFDR.bed",sep="")),sep="\t",col.names=F,row.names=F,quote=F)
  write.csv(res,file=file.path(paste(i,"-Step2.csv",sep="")))
  
}



allbeds =   Sys.glob(file.path(paste("beds-Step2/","*","-SigFDR.bed",sep="")))%>%setNames(.,.)
names(allbeds)%<>%basename%>%str_replace('\\-SigFD.*','')
library(rtracklayer)
allbeddf<-allbeds%>%map(import)%>%GRangesList%>%unlist%>%unique
names(allbeddf)<-NULL
allbeddf$score = 0
allbeddf%>%export('beds-Step2/merge4julia.bed')

################
library(DESeq2)
library("pheatmap")
ids<-list.files(pattern="-Step2.csv")
res<-read.csv(ids[1],header=T)[,8:11]
resf<-read.csv(ids[1],header=T)[,8:11]
for(id in ids){
  print(id)
  temp<-read.csv(id,header=T)[,7]
  res<-cbind(res,temp)
  tempf<-read.csv(id,header=T)[,3]
  resf<-cbind(resf,tempf)
};
res[is.na(res)]<-1;colnames(res)[5:13]<-ids;
#res1<-res[rowSums(res[,5:13]<0.05)>0,];rownames(res1)<-paste(res1$pos,res1$target,sep="-")
#adding explicit fold change threshold here 
FCTHRESH=1
for(FCTHRESH in FCTHRESHvalues){
  
  deseq_sig <- (res$padj<FDR_THRESHOLD & !is.na(res$padj) )
    
  res1<-res[rowSums(res[,5:13]<FDR_THRESHOLD & resf[,5:13]>FCTHRESH )>0,];
  rownames(res1)<-paste(res1$pos,res1$target,sep="-")#added nov 2019
  
  #and writing to appropriate file
  write.csv(res1,file.path(paste0("Targets-FDR-lFC",FCTHRESH,".csv")))
 
  stopifnot(nrow(res1) <= length(allbeddf))
  
  
  res2<-res1[,5:13]
  res2[res2<0.05]<-0
  res2[res2>=0.05]<-1
 # pheatmap((res2),scale="none", cluster_rows=T, show_rownames=T,cluster_cols=T,fontsize_row=5)
}

ids<-list.files(pattern="-Step2.csv")
res<-read.csv(ids[1],header=T)[,8:11]
for(id in ids){
  print(id)
  temp<-read.csv(id,header=T)[,3]
  res<-cbind(res,temp)
};res[is.na(res)]<- -5;colnames(res)[5:13]<-ids
res2<-res[,5:13]
#pheatmap((res2),scale="none", cluster_rows=T, show_rownames=T,cluster_cols=T,fontsize_row=5)

##All
colData<-data.frame(cond=factor(cond),rep=as.factor(rep),geno=geno,dev=dev);rownames(colData)<-colnames(counts)
dds <- DESeqDataSetFromMatrix(countData = counts,colData = colData,design = ~ cond+rep+geno+dev)
dds1 <- DESeq(dds,fitType="local")
vsd <- varianceStabilizingTransformation(dds1,fitType="local")
plotfile<-file.path("PCA-vsdNomrPV.pdf")
pdf(plotfile)
pcaData <- plotPCA(vsd, intgroup=c("cond", "geno"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
print(ggplot(pcaData, aes(PC1, PC2, color=cond, shape=geno)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed())
message(normalizePath(plotfile))
dev.off()
write.csv(cbind(coor,as.data.frame(assay(vsd))),file=file.path("ALLGenesVSDNormPV.csv"))

fin<-read.csv("ALLGenesVSDNormPV.csv");rownames(fin)<-fin$X;fin$X<-NULL

#pheatmap(as.matrix(fin), cluster_rows=T, show_rownames=F,cluster_cols=T,scale="none",main="scale by ")#, annotation_row =colData)	


ids<-colnames(fin)
samples<-ids[grep("sample",ids)]
controls<-sub("sample","control",samples)

res<-as.data.frame(lapply(1:length(samples),function(x){print(x);fin[,samples[x]]-fin[,controls[x]]}))
colnames(res)<-samples;rownames(res)<-rownames(fin)
temp<-strsplit(colnames(res),"_|.bam")
genotype<-unlist(lapply(temp,function(x)x[1]))
dev<-unlist(lapply(temp,function(x)x[2]))
cond0<-unlist(lapply(temp,function(x)x[3]));cond<-cond0
cond[grep("sample",cond)]<-"sample";cond[grep("control",cond)]<-"control"
colData<-data.frame(cond,dev,genotype);rownames(colData)<-colnames(res)

# for(i in 1){for(i in 1){for(i in 1){
plotfile <- file.path("Log2RatioDeseq2-heatmapPV.pdf")
pdf(plotfile,h=14,w=14)
pheatmap(res,scale="none", cluster_rows=T, show_rownames=T,cluster_cols=T, annotation_col=colData)
pheatmap(res[,dev=="early"],scale="none", cluster_rows=T, show_rownames=T,cluster_cols=T, annotation_col=colData[dev=="early",])
pheatmap(res[,dev=="late"],scale="none", cluster_rows=T, show_rownames=T,cluster_cols=T, annotation_col=colData[dev=="late",])
sampleDists <- cor(res)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(res)
colnames(sampleDistMatrix) <- colnames(res)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=dist(sampleDists),
         clustering_distance_cols=dist(sampleDists),
         main="Pearson Correlation",annotation_col=colData)
pca <- prcomp(res,center = TRUE, scale = TRUE) 
plot(pca[[2]][,1:2],lwd=.001)
points(pca[[2]][genotype=="clb19"  & dev=="late",1:2],col="red",lwd=4)
points(pca[[2]][genotype=="wt"  & dev=="late",1:2],col="blue",lwd=4)
points(pca[[2]][genotype=="ptac2"  & dev=="late",1:2],col="green",lwd=4)
points(pca[[2]][genotype=="rpotmp"  & dev=="late",1:2],col="brown",lwd=4)
points(pca[[2]][genotype=="rpotp"  & dev=="late",1:2],col="yellow",lwd=4)

points(pca[[2]][genotype=="clb19" & dev=="early",1:2],col="red",lwd=4,pch=24)
points(pca[[2]][genotype=="wt" & dev=="early",1:2],col="blue",lwd=4,pch=24)
points(pca[[2]][genotype=="ptac2" & dev=="early",1:2],col="green",lwd=4,pch=24)
points(pca[[2]][genotype=="rpotmp" & dev=="early",1:2],col="brown",lwd=4,pch=24)
points(pca[[2]][genotype=="rpotp" & dev=="early",1:2],col="yellow",lwd=4,pch=24)
legend("topright",legend=c("clb19","wt","ptac2","rpotmp","rpotp"),col=c("red","blue","green","brown","yellow"),lwd=3)
dev.off()
message(normalizePath(plotfile))
# 
# pdf(plotfile);for(i in 1)qplot(1) ; dev.off()
plotfile%>%file.info%>%.$size%>%`>`(4e3)

write.csv(res,file.path("Log2RatioDeseq2-heatmapPV.csv"))

}}