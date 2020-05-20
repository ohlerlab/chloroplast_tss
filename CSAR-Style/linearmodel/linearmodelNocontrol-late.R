require(DESeq2)
require("pheatmap")
require(RColorBrewer)
FC=2
library(here)
setwd(here("/CSAR-Style/"))
library(DESeq2)
library("pheatmap")
require("gplots")
require(ggplot2)
library(memoise)
pdf <- grDevices::pdf
DESeq <- memoise(DESeq2::DESeq)

for(FDR_THRESHOLD in c(0.05)){
for(CONTROLSCORETHRESHOLD in c(5)){
for(FC in c(1.5,2)){
#   
#This creates a lot of files in a lot of places with a lot of functions
#moving into a new working directory is probably the best option
setwd(here("/CSAR-Style/"))
getwd()
outdir <- paste0('run','fdr_',FDR_THRESHOLD,'_pseudocount_',CONTROLSCORETHRESHOLD)
dir.create(showWarnings = F,outdir) 
setwd(outdir)
library(tidyverse)

res0<-read.csv(paste("Targets-FDR-lFC",FC,".csv",sep=""))
fdr<-read.csv(paste("Targets-FDR-lFC",FC,".csv",sep=""));rownames(fdr)<-fdr$X;fdr$X<-NULL

message("foobarbar")
load("Counts.RD");
counts<-counts[,-grep("early|control|clb19",colnames(counts))];ids<-colnames(counts)
geno<-as.factor(unlist(lapply(strsplit(ids,"_"),function(x)x[1])))
dev<-as.factor(unlist(lapply(strsplit(ids,"_"),function(x)x[2])))
rep<-as.factor(unlist(lapply(strsplit(ids,"control|sample"),function(x)x[2])))
cond<-as.factor(substr(unlist(lapply(strsplit(ids,"_|\\."),function(x)x[3])),1,3))
samp<-paste(geno,dev,sep="_")
colData<-data.frame(cond=factor(cond),rep=as.factor(rep),geno=geno,dev=dev);rownames(colData)<-ids
dds <- DESeqDataSetFromMatrix(countData = counts,colData = colData,design = ~ geno)
dds1 <- DESeq(dds,fitType="local")
vsd <- varianceStabilizingTransformation(dds1,fitType="local")
rtp<-results(dds1,contrast=c("geno","rpotp","wt"))
rmp<-results(dds1,contrast=c("geno","rpotmp","wt"))
ptac2<-results(dds1,contrast=c("geno","ptac2","wt"))
pv0<- -log10(data.frame(rpotp=rtp[,c(5)],rpotmp=rmp[,c(5)],ptac2=ptac2[,c(5)]));rownames(pv0)<-rownames(vsd)
fc<-data.frame(rtp[,c(2)],rmp[,c(2)],ptac2[,c(2)]);rownames(fc)<-rownames(vsd)
#pv<-pv0*fc/abs(fc)
sigtss<-which(is.element(paste(coor$Start,coor$target,sep="-"),res0$X))
exp<-as.data.frame(pv0[sigtss,]);
exp1<-as.data.frame(fc[sigtss,]);
res1=exp1

order<-sapply(rownames(res1),function(x)grep(x,fdr$pos))
fdr1<-fdr[order,]


message("foobar")
if(FALSE){
####
data<-assay(vsd)
data<-data[names(order),]##Only significants
colnames(data)<-sub(".bam.counts","Normalizedcounts.BedGraph",colnames(data))
for(i in 1:13){
dat<-data.frame("ChrC",as.integer(rownames(data))-1,as.integer(rownames(data)),data[,i])
con<-file(paste("bedGraphs/",colnames(data)[i],sep=""),open="w+")
writeLines('track type=bedGraph name="BedGraph Format" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20',con)
write.table(dat,file=con,sep="\t",col.names=F,quote=F,row.names=F,append=T)
close(con)
}
}
####




temp<-strsplit(colnames(res1)," ")
genotype<-unlist(lapply(temp,function(x)x[2]))
dev<-unlist(lapply(temp,function(x)x[1]))
colData<-data.frame(dev,genotype);rownames(colData)<-colnames(res1)
anno<-read.csv(here("ext_data/PlastidAnnotation.csv"),stringsAsFactors=F)
symbol<-sapply(as.character(fdr1$target),function(x)c(anno$symbol[anno$TAIR_id==x],NA)[1])
type<-sapply(as.character(fdr1$target),function(x)c(anno$type[anno$TAIR_id==x],NA)[1])
rownames(res1)<-paste(fdr1$pos,fdr1$strand,symbol)
type[is.na(type)]<-"No_assigned"
rowData<-data.frame(type);rownames(rowData)<-rownames(res1)
ann_colors=list(type=c(other="light pink",Photosynthesis="green",Ribosome="burlywood4",RNA_polymerase="dark red",rRNA="red",tRNA="dark blue",No_assigned="yellow"))

colnames(res1)<-substr(colnames(res1),1,3)


library(fpc)
message("foo")
pdf(paste("LateFCvsWT-Kmeans-lFC",FC,".pdf",sep=""))
for(nc in 2:10){
par(mfrow=c(4,3))
cluster=kmeans(res1[,1:3],nc,iter.max=10^6,nstart=100)$cluster
 table(cluster)
plotcluster(res1[,1:3],cluster,main=paste("psbA is cluster", as.integer(cluster["1515 - PSBA"])))
boxplot(res1[,1:3],main="All");abline(v=0,col="green",lwd=3);abline(h=0,lwd=3,col="green");abline(h=0,lwd=3,col="green")
for (i in 1:nc){
        boxplot(res1[cluster==i,1:3],main=paste("cluster",i,"\n n=",sum(cluster==i)));
        abline(h=0,lwd=3,col="green")}
        write.csv(data.frame(res1,kmers=cluster),file=paste("LateFCvsWT-Kmeans",nc,"-FC",FC,".csv",sep=""))
}
dev.off()


message("ffffff")
pdf(paste("LateFCvsWT-pamk-lFC",FC,".pdf",sep=""))
cluster1=pamk(res1[,1:3],3:15)[[1]][[3]]
 table(cluster1)
 par(mfrow=c(3,3))
plotcluster(res1[,1:3],cluster1,main=paste("psbA is cluster", as.integer(cluster1["1515 - PSBA"])))
boxplot(res1[,1:3],main="All");abline(h=0,lwd=3,col="green")
for (i in 1:max(cluster1)){boxplot(res1[cluster1==i,1:3],main=paste("cluster",i,"\n n=",sum(cluster1==i)));abline(h=0,lwd=3,col="green")}
write.csv(data.frame(res1,kmers=cluster1),file=paste("LateFCvsWT-pamK-",FC,".csv",sep=""))
dev.off()






message("peh")
pdf(paste0("LateFCvsWT-Mclust-lFC",FC,".pdf"))
library(mclust)
cluster=Mclust(res1[,1:3],G=3:9)$classification
 table(cluster)
 par(mfrow=c(3,3))
plotcluster(res1[,1:3],cluster,main=paste("psbA is cluster", as.integer(cluster["1515 - PSBA"])))
boxplot(res1[,1:3],main="All")
for (i in unique(cluster))boxplot(res1[cluster==i,1:3],main=paste("cluster",i,"\n n=",sum(cluster==i)))
dev.off()

write.csv(data.frame(res1,kmers=cluster,pamk=cluster1),file="LateFCvsWT-2Classifiers.csv")

write.csv(res1,file="LateFCvsWT-all.csv")
res2=res1;res2[exp>  -log10(0.05)]<-0
write.csv(res2,file="LateFCvsWT-onlysig.csv")
}
}
}
