cd /scratch/AG_Ohler/jmuino/Julia/Tex-experiment/nutResolution/beds
rm merged*.bed
cat *Sigpv05.bed >> merged.bed
sort  -k1,1 -k2,2n -i merged.bed > mergeds.bed
mergeBed -i mergeds.bed -iobuf 5G -d 1 -s  >SigPV05-merged.bed

setwd("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment/nutResolution/beds")

q1<-read.table("SigPV05-merged.bed");q1$V1<-as.character(q1$V1)
gtf<-data.frame(q1$V1,"DESEQ2","gene",q1$V2,q1$V3,"1000",q1$V4,".",paste("gene_id \"",q1$V2,"\"; transcript_id \"",q1$V2,"\";",sep=""))
gtf[,9]<-as.character(gtf[,9]);gtf[,7]<-as.character(gtf[,7])
write.table(gtf,file="SigPV-merged.gtf",quote=F,row.names=F,col.names=F,sep="\t")

setwd("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment/nutResolution/")
library(DESeq2)
library("pheatmap")
require("gplots")
 require(ggplot2)

ids<-list.files(pattern=".counts$","countsPV")
res<-read.table(paste("countsPV/",ids[1],sep=""),header=T)
for(id in ids[-1]){
print(id)
temp<-read.table(paste("countsPV/",id,sep=""),header=T)[,7]
res<-cbind(res,temp)
}
coor<-res[,1:6]
counts<-res[,7:length(res[1,])]
colnames(counts)<-ids;rownames(counts)<-res$Geneid
for(i in grep("control",colnames(counts))){counts[counts[,i]<2,i]<-2}

####Annotate
mdis=5000
gff<-read.table("/scratch/AG_Ohler/jmuino/genomes/TAIR9/Araport11_GFF3_genes_transposons.201606.gtf",sep="\t")
gff<-gff[gff$V3=="gene" & gff$V1=="ChrC",]
dis<-lapply(coor$Start[coor$Strand=="+"],function(x){temp<- x-gff$V4[gff$V7=="+"];temp[abs(temp)<1000]})
hist(unlist(dis),breaks=1000)
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

samp<-paste(geno,dev,sep="_")
colData<-data.frame(cond=factor(cond),rep=as.factor(rep),geno=geno,dev=dev);rownames(colData)<-ids
for(i in unique(samp)){
print(i)
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts[,samp==i]),colData = colData[samp==i,],design = ~ cond+rep)
dds1 <- DESeq(dds,fitType="local")
res<-as.data.frame(results(dds1,contrast=c("cond","sam","con"),altHypothesis="greater"))
res$pos<-coor$Start;res$pos2<-coor$End;res$strand<-coor$Strand;res$target<-coor$target;res$dis<-coor$dis
sig<-res[res$padj<0.05 & !is.na(res$padj),7:10];sig$pos<-sig$pos## It is already in0-based coordinates
sig$chr<-"ChrC";sig$point<-".";sig<-sig[,c(5,1,2,4,6,3)]
write.table(sig,file=paste("beds-Step2/",i,"-SigFDR.bed",sep=""),sep="\t",col.names=F,row.names=F,quote=F)
write.csv(res,file=paste(i,"-Step2.csv",sep=""))
print(res[res$pos==54617,])
}

ids<-list.files(pattern="-Step2.csv")
res<-read.csv(ids[1],header=T)[,8:11]
for(id in ids){
print(id)
temp<-read.csv(id,header=T)[,7]
res<-cbind(res,temp)
};res[is.na(res)]<-1;colnames(res)[5:13]<-ids
res1<-res[rowSums(res[,5:13]<0.05)>0,];rownames(res1)<-paste(res1$pos,res1$target,sep="-")
write.csv(res1,"Targets-FDR.csv")
res2<-res1[,5:13]
res2[res2<0.05]<-0
res2[res2>=0.05]<-1
pheatmap((res2),scale="none", cluster_rows=T, show_rownames=T,cluster_cols=T,fontsize_row=5)


ids<-list.files(pattern="-Step2.csv")
res<-read.csv(ids[1],header=T)[,8:11]
for(id in ids){
print(id)
temp<-read.csv(id,header=T)[,3]
res<-cbind(res,temp)
};res[is.na(res)]<- -5;colnames(res)[5:13]<-ids
res2<-res[,5:13]
pheatmap((res2),scale="none", cluster_rows=T, show_rownames=T,cluster_cols=T,fontsize_row=5)

##All
colData<-data.frame(cond=factor(cond),rep=as.factor(rep),geno=geno,dev=dev);rownames(colData)<-colnames(counts)
dds <- DESeqDataSetFromMatrix(countData = counts,colData = colData,design = ~ cond+rep+geno+dev)
dds1 <- DESeq(dds,fitType="local")
vsd <- varianceStabilizingTransformation(dds1,fitType="local")
pdf("PCA-vsdNomrPV.pdf")
pcaData <- plotPCA(vsd, intgroup=c("cond", "geno"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=cond, shape=geno)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
  dev.off()
write.csv(as.data.frame(assay(vsd)),file="ALLGenesVSDNormPV.csv")

 fin<-read.csv("ALLGenesVSDNormPV.csv");rownames(fin)<-fin$X;fin$X<-NULL
pheatmap(as.matrix(fin), cluster_rows=T, show_rownames=F,cluster_cols=T,scale="none",main="scale by ")#, annotation_row =colData)	


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

#res1<-res[rowSums(res>2)>4,]
#res1[res1<0]<-0
pdf("Log2RatioDeseq2-heatmapPV.pdf")
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
write.csv(res,"Log2RatioDeseq2-heatmapPV.csv")

