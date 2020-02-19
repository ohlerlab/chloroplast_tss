
##All
setwd("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment/CSAR-Style/Heatmaps")
require(DESeq2)
require("pheatmap")
require(RColorBrewer)
res1<-read.csv("../Targets-FDR.csv")
#res1<-res1[res1$wt_early.Step2.csv<0.05 | res1$wt_late.Step2.csv<0.05,]#adding 5 sept

load("../Counts.RD");counts<-counts[,grep("sample",colnames(counts))];ids<-colnames(counts)
geno<-as.factor(unlist(lapply(strsplit(ids,"_"),function(x)x[1])))
dev<-as.factor(unlist(lapply(strsplit(ids,"_"),function(x)x[2])))
rep<-as.factor(unlist(lapply(strsplit(ids,"control|sample"),function(x)x[2])))
cond<-as.factor(substr(unlist(lapply(strsplit(ids,"_|\\."),function(x)x[3])),1,3))
samp<-paste(geno,dev,sep="_")
colData<-data.frame(cond=factor(cond),rep=as.factor(rep),geno=geno,dev=dev);rownames(colData)<-ids
dds <- DESeqDataSetFromMatrix(countData = counts,colData = colData,design = ~ geno+dev)
dds1 <- DESeq(dds,fitType="local")
vsd <- varianceStabilizingTransformation(dds1,fitType="local")
sigtss<-which(is.element(paste(coor$Start,coor$target,sep="-"),res1$X))
exp<-as.data.frame(assay(vsd[sigtss]))

ids<-paste(dev,geno);
res1<-as.data.frame(lapply(unique(ids),function(x){rowMeans(exp[,ids==x])}))
colnames(res1)<-unique(ids)
temp<-strsplit(colnames(res1)," ")
genotype<-unlist(lapply(temp,function(x)x[2]))
dev<-unlist(lapply(temp,function(x)x[1]))
colData<-data.frame(dev,genotype);rownames(colData)<-colnames(res1)
fdr<-read.csv("../Targets-FDR.csv");rownames(fdr)<-fdr$X;fdr$X<-NULL
order<-sapply(rownames(res1),function(x)which(fdr$pos==x))
fdr1<-fdr[order,]
anno<-read.csv("/scratch/AG_Ohler/jmuino/CRPC-Marie/PlastidAnnotation.csv",stringsAsFactors=F)
symbol<-sapply(as.character(fdr1$target),function(x)c(anno$symbol[anno$TAIR_id==x],NA)[1])
type<-sapply(as.character(fdr1$target),function(x)c(anno$type[anno$TAIR_id==x],NA)[1])
rownames(res1)<-paste(fdr1$pos,fdr1$strand,symbol)
type[is.na(type)]<-"No_assigned"
rowData<-data.frame(type);rownames(rowData)<-rownames(res1)

require(tsne)
#qq=(tsne(res1,perplexity=5,max_iter=10000))
plot(qq)

pca=prcomp(t(res1))
plot(pca[[2]][,1:2])
points(pca[[2]][dev=="early",1:2],col="red")
text(pca[[2]][,1:2], label=samp)
ann_colors=list(type=c(other="light pink",Photosynthesis="green",Ribosome="burlywood4",RNA_polymerase="dark red",rRNA="red",tRNA="dark blue",No_assigned="yellow"))

avg=res1$'early wt'
res2=res1-avg
pp<-kmeans(res2,5);
pp1<-sort(pp$cluster)
cluster=pp1[rownames(res2)]
rowData1<-data.frame(type,cluster);rownames(rowData1)<-rownames(res1)
pheatmap(as.matrix(res2[names(pp1),]), cluster_rows=T, show_rownames=T,cluster_cols=T,scale="none",annotation_col=colData,fontsize_row=3,
annotation_row=rowData1,annotation_colors=ann_colors)

pdf(paste("JuliaHeatmaps-NormalizedCountsDistance.pdf"))
for(dd in c("euclidean", "maximum", "manhattan", "canberra","binary","minkowski","correlation")){
pheatmap(as.matrix((res1-res1$'early wt')[-which(names(res1)=="early wt")]), cluster_rows=T, show_rownames=T,cluster_cols=T,scale="none",annotation_col=colData,fontsize_row=3,
annotation_row=rowData,annotation_colors=ann_colors,main=dd,clustering_distance_rows = dd,clustering_distance_cols = dd)
}
dev.off()

pdf(paste("Effec2EachMtype-OnlyWTrefEarlyWt.pdf"))
res2<-res1-res1$'early wt'
for(x in unique(type)){
plot(density(res2[type==x,1],bw=.5),ylim=c(0,1),lwd=4,main=x)
for(i in 2:9)lines(density(res2[type==x,i],bw=.5),col=i,lwd=4)
abline(v=0)
legend("topright",legend=colnames(res2),col=1:9,lwd=3)
}
dev.off()






pdf("JuliaHeatmaps-NormalizedCountsOnlyWT.pdf")
tf<-list(c(Sig="blue",NoSig="grey"))
ann_colors=list(type=c(other="light pink",Photosynthesis="green",Ribosome="burlywood4",RNA_polymerase="dark red",rRNA="red",tRNA="dark blue",No_assigned="yellow"))
#pheatmap(as.matrix(res1), cluster_rows=T, show_rownames=T,cluster_cols=T,scale="none",annotation_col=colData,fontsize_row=3,annotation_row=rowData,annotation_colors=ann_colors)
#pheatmap(as.matrix(res1), cluster_rows=T, show_rownames=T,cluster_cols=T,scale="row",annotation_col=colData,fontsize_row=3,annotation_row=rowData,annotation_colors=ann_colors,main="Standarized by row")
avg<-res1$'late wt'
med<-sapply(1:length(res1[,1]),function(x)median(unlist(res1[x,])))
pheatmap(as.matrix(log2(res1/avg)), cluster_rows=T, show_rownames=T,cluster_cols=T,scale="none",annotation_col=colData,fontsize_row=3,annotation_row=rowData,annotation_colors=ann_colors
,main="log2 ratio value/avg value in late WT", color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),breaks=unique(c(seq(-3,0,length.out=51),seq(0,2,length.out=50))))
pheatmap(as.matrix(res1-avg), cluster_rows=T, show_rownames=T,cluster_cols=T,scale="none",annotation_col=colData,fontsize_row=3,annotation_row=rowData,annotation_colors=ann_colors
,main="Difference: value - avg value in late WT",breaks=unique(c(seq(-4,0,length.out=51),seq(0,6,length.out=50))))
dev.off()


#for(dd in c("euclidean", "maximum", "manhattan", "canberra","binary","minkowski","correlation")){
dd="euclidean"
pdf(paste("JuliaHeatmaps-NormalizedCountsOnlyWTDistance-Early-",dd,".pdf"))
tf<-list(c(Sig="blue",NoSig="grey"))
ann_colors=list(type=c(other="light pink",Photosynthesis="green",Ribosome="burlywood4",RNA_polymerase="dark red",rRNA="red",tRNA="dark blue",No_assigned="yellow"))
w=grep("early",colnames(res1))
pheatmap(as.matrix(res1)[,w], cluster_rows=T, show_rownames=T,cluster_cols=T,scale="none",annotation_col=colData,fontsize_row=3,annotation_row=rowData,annotation_colors=ann_colors,clustering_distance_rows = dd,clustering_distance_cols = dd)
pheatmap(as.matrix(res1)[,w], cluster_rows=T, show_rownames=T,cluster_cols=T,scale="row",annotation_col=colData,fontsize_row=3,annotation_row=rowData,annotation_colors=ann_colors,main="Standarized by row",clustering_distance_rows = dd,clustering_distance_cols = dd)
avg<-rowMeans(res1[,w])
med<-sapply(1:length(res1[,1]),function(x)median(unlist(res1[x,w])))
pheatmap(as.matrix(log2(res1[,w]/avg)), cluster_rows=T, show_rownames=T,cluster_cols=T,scale="none",annotation_col=colData,fontsize_row=3,annotation_row=rowData,annotation_colors=ann_colors
,main="log2 ratio value/avg value across all conditions", color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),breaks=unique(c(seq(-2,0,length.out=51),seq(0,1,length.out=50))),clustering_distance_rows = dd,clustering_distance_cols = dd)
pheatmap(as.matrix(res1[,w]-avg), cluster_rows=T, show_rownames=T,cluster_cols=T,scale="none",annotation_col=colData,fontsize_row=3,annotation_row=rowData,annotation_colors=ann_colors
,main="Difference: value - avg value across all conditions",breaks=seq(-4,4,length.out=100),clustering_distance_rows = dd,clustering_distance_cols = dd)
dev.off()

pdf(paste("JuliaHeatmaps-NormalizedCountsOnlyWTDistance-Late-",dd,".pdf"))
tf<-list(c(Sig="blue",NoSig="grey"))
ann_colors=list(type=c(other="light pink",Photosynthesis="green",Ribosome="burlywood4",RNA_polymerase="dark red",rRNA="red",tRNA="dark blue",No_assigned="yellow"))
w=grep("late",colnames(res1))
pheatmap(as.matrix(res1)[,w], cluster_rows=T, show_rownames=T,cluster_cols=T,scale="none",annotation_col=colData,fontsize_row=3,annotation_row=rowData,annotation_colors=ann_colors,clustering_distance_rows = dd,clustering_distance_cols = dd)
pheatmap(as.matrix(res1)[,w], cluster_rows=T, show_rownames=T,cluster_cols=T,scale="row",annotation_col=colData,fontsize_row=3,annotation_row=rowData,annotation_colors=ann_colors,main="Standarized by row",clustering_distance_rows = dd,clustering_distance_cols = dd)
avg<-rowMeans(res1[,w])
med<-sapply(1:length(res1[,1]),function(x)median(unlist(res1[x,w])))
pheatmap(as.matrix(log2(res1[,w]/avg)), cluster_rows=T, show_rownames=T,cluster_cols=T,scale="none",annotation_col=colData,fontsize_row=3,annotation_row=rowData,annotation_colors=ann_colors
,main="log2 ratio value/avg value across all conditions", color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),breaks=unique(c(seq(-2,0,length.out=51),seq(0,1,length.out=50))),clustering_distance_rows = dd,clustering_distance_cols = dd)
pheatmap(as.matrix(res1[,w]-avg), cluster_rows=T, show_rownames=T,cluster_cols=T,scale="none",annotation_col=colData,fontsize_row=3,annotation_row=rowData,annotation_colors=ann_colors
,main="Difference: value - avg value across all conditions",breaks=seq(-4,4,length.out=100),clustering_distance_rows = dd,clustering_distance_cols = dd)
dev.off()
#}
