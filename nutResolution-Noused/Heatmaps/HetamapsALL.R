
##All
setwd("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment/CSAR-Style/Heatmaps")
require("pheatmap")
fin<-read.csv("../ALLGenesVSDNormPV.csv");rownames(fin)<-fin$X;fin$X<-NULL
anno<-read.csv("/scratch/AG_Ohler/jmuino/CRPC-Marie/PlastidAnnotation.csv",stringsAsFactors=F)
symbol<-sapply(as.character(fin$target),function(x)c(anno$symbol[anno$TAIR.id==x],NA)[1])
type<-sapply(as.character(fin$target),function(x)c(anno$type[anno$TAIR.id==x],NA)[1])
rownames(fin)<-paste(fin$Start,fin$Strand,symbol)
names(type)<-paste(fin$Start,fin$strand,symbol)
fin<-fin[,9:62]

ids<-colnames(fin)
samples<-ids[grep("sample",ids)]
controls<-sub("sample","control",samples)
res<-as.data.frame(lapply(1:length(samples),function(x){print(x);(fin[,samples[x]])-fin[,controls[x]]}))
colnames(res)<-samples;rownames(res)<-rownames(fin)
temp<-strsplit(colnames(res),"_|.bam")
genotype<-unlist(lapply(temp,function(x)x[1]))
dev<-unlist(lapply(temp,function(x)x[2]))
cond0<-unlist(lapply(temp,function(x)x[3]));cond<-cond0
cond[grep("sample",cond)]<-"sample";cond[grep("control",cond)]<-"control"
colData<-data.frame(cond,dev,genotype);rownames(colData)<-colnames(res)


ids<-paste(dev,genotype);
res1<-as.data.frame(lapply(unique(ids),function(x){rowMeans(res[,ids==x])}))
pv<-as.data.frame(lapply(unique(ids),function(x){temp<-(res[,ids==x]); sapply(1:length(res[,1]),function(q)t.test(temp[q,],alternative="greater",mu=0)$p.value)}))
colnames(pv)<-unique(ids);rownames(pv)<-rownames(res)
colnames(res1)<-unique(ids)
temp<-strsplit(colnames(res1)," ")
genotype<-unlist(lapply(temp,function(x)x[2]))
dev<-unlist(lapply(temp,function(x)x[1]))
ww= which(rowSums(pv<0.1,na.rm=T)>0)
res2<-res1[ww,];pv2<-pv[ww,]
#res2[pv2>0.05]<-0
sig<-pv2;sig[pv2>0.05 | is.na(pv2)]<-1;sig[pv2<0.05]<-0;sig<-sig==0;sig<-data.frame(sig);rownames(sig)<-rownames(res2)
colData<-data.frame(dev,genotype);rownames(colData)<-colnames(res2)
pheatmap(as.matrix(res2), cluster_rows=T, show_rownames=F,cluster_cols=T,scale="none",annotation_col=colData)#,annotation_row=sig)





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

