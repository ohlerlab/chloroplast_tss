setwd("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment/nutResolution")
library(DESeq2)
library("pheatmap")
##Significants
ids<-list.files(pattern="-Step2.csv")
res<-read.csv(ids[1],header=T)[,8:11]
for(id in ids){
print(id)
temp<-read.csv(id,header=T)[,7]
res<-cbind(res,temp)
};res[is.na(res)]<-1;colnames(res)[5:13]<-ids
res1<-res[rowSums(res[,5:13]<0.05)>0,];rownames(res1)<-paste(res1$pos,res1$target,sep="-")
write.csv(res1,"Targets-FDR.csv")


ids<-list.files(pattern="-Step2.csv")
resf<-read.csv(ids[1],header=T)[,8:11]
for(id in ids){
print(id)
temp<-read.csv(id,header=T)[,3]
resf<-cbind(resf,temp)
};colnames(resf)[5:13]<-ids
res1<-resf[rowSums(res[,5:13]<0.05)>0,];rownames(res1)<-paste(res1$pos,res1$target,sep="-")
write.csv(res1,"Targets-log2FC.csv")

setwd("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment/nutResolution/Heatmaps")
library("pheatmap")


res1<-read.csv("../Targets-log2FC.csv");rownames(res1)<-res1$X;res1$X<-NULL
fdr1<-read.csv("../Targets-FDR.csv");rownames(fdr1)<-fdr1$X;fdr1$X<-NULL
table(rownames(res1)==rownames(res1))
anno<-read.csv("/scratch/AG_Ohler/jmuino/CRPC-Marie/PlastidAnnotation.csv",stringsAsFactors=F)
symbol<-sapply(as.character(res1$target),function(x)c(anno$symbol[anno$TAIR_id==x],NA)[1])
type<-sapply(as.character(res1$target),function(x)c(anno$type[anno$TAIR_id==x],NA)[1])
res2<-res1[,5:13];fdr2<-fdr1[,5:13]
colnames(res2)<-sub(".Step2.csv","",colnames(res2))
colnames(fdr2)<-sub(".Step2.csv","",colnames(fdr2))
rownames(res2)<-paste(res1$pos,res1$strand,symbol)
rownames(fdr2)<-paste(res1$pos,res1$strand,symbol)
type[is.na(type)]<-"No_assigned"
fdr2[fdr2<0.05]<-0;fdr2[fdr2>0.05]<-1
fdr2[fdr2=="0"]<-"Sig";fdr2[fdr2=="1"]<-"NoSig"
rowData<-data.frame(type,fdr2);rownames(rowData)<-rownames(res2)
rowData<-rowData[,10:1]


pdf("JuliaHeatmaps-log2FC.pdf")
tf<-list(c(Sig="blue",NoSig="grey"))
ann_colors=list(type=c(other="light pink",Photosynthesis="green",Ribosome="burlywood4",RNA_polymerase="dark red",rRNA="red",tRNA="dark blue",No_assigned="yellow"))
ann_colors<-c(ann_colors,rep(tf,9));names(ann_colors)<-c("type",colnames(fdr2))
pheatmap(res2,scale="none", cluster_rows=T, show_rownames=T,cluster_cols=T,fontsize_row=3,annotation_row=rowData,annotation_colors=ann_colors)
pheatmap((res2[,grep("late",colnames(res2))]),scale="none", cluster_rows=T,fontsize_row=3,annotation_row=rowData[,c(grep("late",colnames(rowData)),10)],annotation_colors=ann_colors)
pheatmap((res2[,grep("early",colnames(res2))]),scale="none", cluster_rows=T,fontsize_row=3,annotation_row=rowData[,c(grep("early",colnames(rowData)),10)],annotation_colors=ann_colors)
pheatmap((res2[,grep("wt",colnames(res2))]),scale="none", cluster_rows=T,fontsize_row=3,annotation_row=rowData[,c(grep("wt",colnames(rowData)),10)],annotation_colors=ann_colors)
dev.off()

pdf("JuliaHeatmaps-log2FC-Standarized.pdf")
tf<-list(c(Sig="blue",NoSig="grey"))
ann_colors=list(type=c(other="light pink",Photosynthesis="green",Ribosome="burlywood4",RNA_polymerase="dark red",rRNA="red",tRNA="dark blue",No_assigned="yellow"))
ann_colors<-c(ann_colors,rep(tf,9));names(ann_colors)<-c("type",colnames(fdr2))
pheatmap(res2,scale="row", cluster_rows=T, show_rownames=T,cluster_cols=T,fontsize_row=3,annotation_row=rowData,annotation_colors=ann_colors)
pheatmap((res2[,grep("late",colnames(res2))]),scale="row", cluster_rows=T,fontsize_row=3,annotation_row=rowData[,c(grep("late",colnames(rowData)),10)],annotation_colors=ann_colors)
pheatmap((res2[,grep("early",colnames(res2))]),scale="row", cluster_rows=T,fontsize_row=3,annotation_row=rowData[,c(grep("early",colnames(rowData)),10)],annotation_colors=ann_colors)
pheatmap((res2[,grep("wt",colnames(res2))]),scale="row", cluster_rows=T,fontsize_row=3,annotation_row=rowData[,c(grep("wt",colnames(rowData)),10)],annotation_colors=ann_colors)
dev.off()