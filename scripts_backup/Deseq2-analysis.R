setwd("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment")
library(DESeq2)
require(Biostrings)
require(ImpulseDE)
library("BiocParallel")
library("pheatmap")
require("gplots")
files<-list.files("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment/counts",pattern="counts$")
ids<-unlist(lapply(strsplit(files,".counts$"),function(x)x[1]))

res<-read.table(paste("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment/counts/",files[1],sep=""),sep="\t",skip=2)$V7
gene<-read.table(paste("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment/counts/",files[1],sep=""),sep="\t",skip=2)$V1
for (i in files[-1]){
temp<-read.table(paste("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment/counts/",i,sep=""),sep="\t",skip=2)
print(table(gene==temp$V1))
res<-cbind(res,temp$V7)
}
rownames(res)<-gene
colnames(res)<-ids

#resn0<-res[substr(gene,3,3)!="C" & substr(gene,3,3)!="M",]
range(round(colSums(res)/10^6,2))


#aa<-unique(substr(names(readAAStringSet("/scratch/AG_Ohler/jmuino/genomes/TAIR9/Araport11_genes.201606.pep.fasta")),1,9))
#resn0<-resn0[is.element(rownames(resn0),aa),]
#plot(hclust(as.dist(1-cor(resn0[rowSums(resn0)>100,],use="complete.obs"))),main="raw data",cex=1,xlab="Using 1-Pearson cor as distance")

resn<-res[rowSums(res>5)> 1,]
write.csv(as.data.frame(resn),file="CountData.csv")

geno<-as.factor(unlist(lapply(strsplit(ids,"_"),function(x)x[1])))
rep<-as.factor(unlist(lapply(strsplit(ids,"control|sample"),function(x)x[2])))
cond<-as.factor(unlist(lapply(strsplit(ids,"_|1|2|3"),function(x)x[2])))

colData<-data.frame(cond=factor(cond),rep=as.factor(rep),geno=geno)
dds <- DESeqDataSetFromMatrix(countData = resn,colData = colData,design = ~ cond+rep+geno)
dds1 <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds1)
write.csv(as.data.frame(assay(vsd)),file="ALLGenesVSDNorm.csv")

 res<-read.csv("ALLGenesVSDNorm.csv")



######
for(i in unique(geno)){
dds <- DESeqDataSetFromMatrix(countData = resn[,geno==i],colData = colData[geno==i,],design = ~ cond+rep)
dds1 <- DESeq(dds)
res1<-results(dds1,contrast=c("cond","sample","control"),altHypothesis="greater")
write.csv(res1,file=paste(i,"_Deseq2.csv",sep="\t"))
sig<-res1[res1$padj<0.05 & !is.na(res1$padj),]
q1<-read.table("merged-1b2w1minl0.bed")
q1$V5<-".";q1$V6<-q1$V4;q1$V4<-"."
q1<-q1[is.element(q1$V2,rownames(sig)),]
write.table(q1,file=paste("merged-1b2w1minl0Sig",i,".bed",sep=""),quote=F,sep="\t",col.names=F,row.names=F)
}


mm<-read.csv(paste(unique(geno)[1],"_Deseq2.csv",sep="\t"))
ids<-as.character(mm$X[mm$padj<0.05])
mm<-mm[,c(1,3)]
for(i in unique(geno)[-1]){
temp<-read.csv(paste(i,"_Deseq2.csv",sep="\t"))

ids<-c(ids,as.character(temp$X[temp$padj<0.05]))
mm<-cbind(mm,temp[,3])
}
rownames(mm)<-mm$X;mm$X<-NULL
colnames(mm)<-unique(geno)
mm1<-mm[is.element(rownames(mm),unique(as.character(ids))),]
mm1[mm1<0 | is.na(mm1)]<-0

m1<-unique(as.character(read.table("motif1.gff")$V1))
m2<-unique(as.character(read.table("motif2.gff")$V1))
m3<-unique(as.character(read.table("motif3.gff")$V1))
m<-data.frame(c(m1,m2,m3),c(rep("m1",length(m1)),rep("m2",length(m2)),rep("m3",length(m3))))
id1<-sapply(as.character(rownames(mm1)),function(x)c(as.character(m[as.character(m[,1])==x,2]),NA)[1])
df=data.frame(motif=id1);rownames(df)<-rownames(mm1)
pdf("Heatmap.pdf")
pheatmap(as.matrix(mm1), cluster_rows=T, show_rownames=F,cluster_cols=T,scale="row",main="scale by row", annotation_row =df)	
pheatmap(as.matrix(mm1), cluster_rows=T, show_rownames=F,cluster_cols=T,scale="none",main="scale by none")	
dev.off()


ids<-c();for( i in unique(geno)){
mm<-read.csv(paste(i,"_Deseq2.csv",sep="\t"))
ids<-c(ids,as.character(mm$X[mm$padj<0.05 & !is.na(mm$padj) & mm$log2FoldChange>2]))
}

q1<-read.table("merged-1b2w1minl0.bed")
q1$V5<-".";q1$V6<-q1$V4;q1$V4<-"."
q1<-q1[is.element(as.character(q1$V2),ids),]
q1$V10<-round(q1$V2+q1$V3)/2
q1<-q1[q1$V1=="ChrC",]
gtf<-data.frame(q1$V1,"CSAR","gene",q1$V10-100,q1$V10+100,"1000",q1$V6,".",paste("gene_id \"",q1$V2,"\"; transcript_id \"",q1$V2,"\";",sep=""))
gtf[,9]<-as.character(gtf[,9]);gtf[,7]<-as.character(gtf[,7])
write.table(gtf,file="merged-1b2w1minl0Siglog2FC2.gtf",quote=F,sep="\t",col.names=F,row.names=F)


#gffread -w merged-1b2w1minl0Siglog2FC2.fa  -g /scratch/AG_Ohler/jmuino/genomes/TAIR9/chr-TAIR9.fas merged-1b2w1minl0Siglog2FC2.gtf
####################

pdf("Normalization.pdf")
plot(log2(resn[,4]/resn["spike-inJose",4]),log2(resn[,5]/resn["spike-inJose",5]),xlab="sample-r1",ylab="sample-r2",main="by spike-inJulia");lines(-1000:10^6,-1000:10^6)
plot(log2(resn[,4]/sum(resn[,4])),log2(resn[,5]/sum(resn[-104,5])),xlab="sample-r1",ylab="sample-r2",main="by coverage");lines(-1000:10^6,-1000:10^6)
plot(assay(vsd)[,4],assay(vsd)[,5],xlab="sample-r1",ylab="sample-r2",main="by DESEQ2");lines(-1000:10^6,-1000:10^6)
resc3<-assay(vsd)
resn1<-resn+1
resc1<-log2(resn1/resn1["spike-inJose",])
resc2<-resn1
resc2<-log2(resc2/colSums(resc2))
plot(resc1[,4]-resc1[,1],resc1[,5]-resc1[,2],main="by spike-inJose");lines(-1000:10^6,-1000:10^6)
plot(resc2[,4]-resc2[,1],resc2[,5]-resc2[,2],main="by coverage");lines(-1000:10^6,-1000:10^6)
plot(resc3[,4]-resc3[,1],resc3[,5]-resc3[,2],main="by DESEQ2");lines(-1000:10^6,-1000:10^6)
points(resc3["spike-inJose",4]-resc3["spike-inJose",1],resc3["spike-inJose",5]-resc3["spike-inJose",2],main="by DESEQ2",col="red",lwd=3);lines(-1000:10^6,-1000:10^6)
dev.off()



###Hetamaps for Julia/Tex-experiment

res<-read.csv("ALLGenesVSDNorm.csv")

w1<-res$wt_sample1-res$wt_control1
w2<-res$wt_sample2-res$wt_control2
w3<-res$wt_sample3-res$wt_control3
ww<-sapply(1:473,function(x)mean(c(w1[x],w2[x],w3[x])))

p1<-res$pep_sample1-res$pep_control1
p2<-res$pep_sample2-res$pep_control2
ww<-sapply(1:473,function(x)mean(c(w1[x],w2[x],w3[x])))

rm1<-res$rpotmp_sample1-res$rpotmp_control1
rm2<-res$rpotmp_sample2-res$rpotmp_control2


r1<-res$rpotp_sample1-res$rpotp_control1
r2<-res$rpotp_sample2-res$rpotp_control2
ww<-sapply(1:473,function(x)mean(c(w1[x],w2[x],w3[x])))

fin<-data.frame(id=res$X,pep=rowMeans(data.frame(p1,p2)),rpotmp=rowMeans(data.frame(rm1,rm2)),rpotp=rowMeans(data.frame(r1,r2)))

ids<-c();for( i in unique(geno)){
mm<-read.csv(paste(i,"_Deseq2.csv",sep="\t"))
ids<-c(ids,as.character(mm$X[mm$padj<0.05 & !is.na(mm$padj) & (mm$log2FoldChange)>2]))
}

an<-read.csv("annotationTSS.csv")
fin1<-fin[is.element(fin$id,ids),]
id1<-sapply(as.character(fin1$id),function(x)c(as.character(an$gene[as.character(an$id)==x]),NA)[1])
id2<-paste(id1,names(id1),sep="_")
rownames(fin1)<-id2
an1<-read.csv("/scratch/AG_Ohler/jmuino/genomes/TAIR9/ArabidopsisPlastidAnnotation.csv")
type<-as.character(sapply(as.character(id1),function(x)c(as.character(an1$type[as.character(an1$TAIR.id)==x]),NA)[1]))
type[is.na(type) | type=="other" |type==" " | type=="unknown"]<-"Other/NA"
df=data.frame(type);rownames(df)<-id2

m1<-unique(as.character(read.table("motif1.gff")$V1))
m2<-unique(as.character(read.table("motif2.gff")$V1))
m3<-unique(as.character(read.table("motif3.gff")$V1))
m<-data.frame(c(m1,m2,m3),c(rep("m1",length(m1)),rep("m2",length(m2)),rep("m3",length(m3))))
fin1<-fin[is.element(fin$id,ids),]
id1<-sapply(as.character(fin1$id),function(x)c(as.character(m[as.character(m[,1])==x,2]),NA)[1])
#id1[is.na(id1)]<-"None"
df=data.frame(motif=id1);rownames(df)<-fin1$id
rownames(fin1)<-fin1$id

colors<-list("white","green","red","blue","brown","black");names(colors)<-unique(type)
pdf("HeatmapMotif.pdf")
pheatmap(as.matrix(fin1[,-1]), cluster_rows=T, show_rownames=F,cluster_cols=T,scale="row", annotation_row =df)
dev.off()

