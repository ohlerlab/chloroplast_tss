setwd("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment")
library(DESeq2)
library(Rsamtools)
library(foreach)
library(doParallel)
require(rtracklayer)
library("pheatmap")
library("RColorBrewer")
library("vsn")

run<-function(id,flag=0){
what <- c("flag", "strand","rname", "pos","mapq","qwidth")
param <- ScanBamParam( what = what)
bam <- scanBam(paste("bamFiles/original/",id,sep=""),param=param)[[1]];
bam$pos[bam$flag==16]<-bam$pos[bam$flag==16]+bam$qwidth[bam$flag==16]-1 ##1-based coordinates
bam$pos[bam$flag==0]<-bam$pos[bam$flag==0] ##1-based coordinates
pos<-bam$pos[bam$rname=="ChrC" & bam$flag==flag]
print(table(bam$flag))

tt=as.data.frame(table(pos));colnames(tt)<-c("pos",id);tt$pos<-as.integer(as.character(tt$pos))
return(tt)
}


ids<-list.files(pattern=".bam$","bamFiles/original")
#ids<-ids[grep("wt",ids)]
##flag 4 unmapped; 16 reverse; 0 foward
fin0<-run(ids[1],flag=0)
for(id in ids[-1]){
print(id)
temp<-run(id,flag=0)
fin0<-merge(fin0,temp,by="pos",all=T)
}
fin0[is.na(fin0)]<-0; rownames(fin0)<-paste(fin0$pos,"+");fin0$pos<-NULL
fin0[fin0<10]<-10
fin10<-fin0[rowSums(fin0[,grep("sample",colnames(fin0))]>10)>2,]

fin0<-run(ids[1],flag=16)
for(id in ids[-1]){
print(id)
temp<-run(id,flag=16)
fin0<-merge(fin0,temp,by="pos",all=T)
}
fin0[is.na(fin0)]<-0; rownames(fin0)<-paste(fin0$pos,"-");fin0$pos<-NULL
fin20<-fin0[rowSums(fin0>10)>2,]
fin1<-rbind(fin10,fin20)
fin1[fin1<10]<-10
setwd("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment/QC")
###QC
sort(round(colSums(fin1)/1000))[1:3]
# rpotp_late_control3.bam   clb19_late_sample4.bam clb19_early_control2.bam
#                       2                        3                       16


temp<-strsplit(colnames(fin1),"_|.bam")
genotype<-unlist(lapply(temp,function(x)x[1]))
dev<-unlist(lapply(temp,function(x)x[2]))
cond0<-unlist(lapply(temp,function(x)x[3]));cond<-cond0
rep<-sub("sample|control","",cond0)
cond[grep("sample",cond)]<-"sample";cond[grep("control",cond)]<-"control"
fin1<-fin1[rowSums(fin1[,cond=="sample"]>10)>2,]

colData<-data.frame(cond,dev,genotype,rep);rownames(colData)<-colnames(fin1)
dds <- DESeqDataSetFromMatrix(countData = as.matrix(fin1),colData = colData,design = ~ cond+genotype+cond+rep)
dds1 <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds1)
save(fin1,colData,dds1,vsd,genotype,dev,cond0,cond,file="DESeq2-Analysis.RD")

setwd("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment/QC")
library(DESeq2)
library(Rsamtools)
library(foreach)
library(doParallel)
require(rtracklayer)
library("pheatmap")
library("RColorBrewer")
library("vsn")
load(file="DESeq2-Analysis.RD")
fin<-as.data.frame(assays(vsd))[,-c(1:2)];rownames(fin)<-rownames(fin1)

nr<-colSums(fin)/10000
colData<-data.frame(colData,coverage=nr)
pdf("heatmapPerGenotype.pdf")
for(i in 1:5 ){
rs<-rowSums(fin[,cond=="sample" &genotype==unique(genotype)[i]])
pheatmap(fin[ rs>10,genotype==unique(genotype)[i]  ], cluster_rows=T, show_rownames=FALSE,cluster_cols=T, annotation_col=colData,main=unique(genotype)[i])

sampleDists <- (cor((fin[ rs>10,genotype==unique(genotype)[i]  ])))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(fin[ rs>10,genotype==unique(genotype)[i]  ])
colnames(sampleDistMatrix) <- colnames(fin[ rs>10,genotype==unique(genotype)[i]  ])
pheatmap(sampleDistMatrix,
         clustering_distance_rows=dist(sampleDists),
         clustering_distance_cols=dist(sampleDists),
         main="Pearson Correlation",annotation_col=colData)

}
dev.off()


ids<-colnames(fin1)
samples<-ids[grep("sample",ids)]
controls<-sub("sample","control",samples)

res<-as.data.frame(lapply(1:length(samples),function(x){print(x);fin[,samples[x]]-fin[,controls[x]]}))
colnames(res)<-samples;rownames(res)<-rownames(fin1)
temp<-strsplit(colnames(res),"_|.bam")
genotype<-unlist(lapply(temp,function(x)x[1]))
dev<-unlist(lapply(temp,function(x)x[2]))
cond0<-unlist(lapply(temp,function(x)x[3]));cond<-cond0
cond[grep("sample",cond)]<-"sample";cond[grep("control",cond)]<-"control"
colData<-data.frame(cond,dev,genotype);rownames(colData)<-colnames(res)

res1<-res[rowSums(res>2)>4,]
#res1[res1<0]<-0
pdf("Log2RatioDeseq2-heatmap.pdf")
pheatmap(res1,scale="none", cluster_rows=T, show_rownames=T,cluster_cols=T, annotation_col=colData)
dev.off()
write.csv(res1,"Log2RatioDeseq2-heatmap.csv")

#####RPKM BAD
fin2<-t(t(fin1)/colSums(fin1))
ids<-colnames(fin2)
samples<-ids[grep("sample",ids)]
controls<-sub("sample","control",samples)
res<-as.data.frame(lapply(1:length(samples),function(x){print(x);fin2[,samples[x]]/fin2[,controls[x]]}))
colnames(res)<-samples
res1<-res[rowSums(res>2)>2,]
pheatmap(log(res1+1),scale="none", cluster_rows=T, show_rownames=FALSE,cluster_cols=T)
