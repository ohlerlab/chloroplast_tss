setwd("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment")
library(DESeq2)
library(Rsamtools)
library(foreach)
library(doParallel)
require(rtracklayer)

run<-function(id,flag=0){
what <- c("flag", "strand","rname", "pos","mapq","qwidth")
param <- ScanBamParam( what = what)
bam <- scanBam(paste("bamFiles/",id,sep=""),param=param)[[1]];
bam$pos[bam$flag==16]<-bam$pos[bam$flag==16]+bam$qwidth[bam$flag==16]-1 ##1-based coordinates
bam$pos[bam$flag==0]<-bam$pos[bam$flag==0] ##1-based coordinates
if(flag==0){rrna<-c(101012:102502,104691:107500,107599:107701,107949:108069)}else{rrna<-c(130580:130700,130948:131050,131149:133958,136147:137637)}
pos<-bam$pos[bam$rname=="ChrC" & bam$flag==flag & !is.element(bam$pos,rrna) ]
print(table(bam$flag))

tt=as.data.frame(table(pos));colnames(tt)<-c("pos",id);tt$pos<-as.integer(as.character(tt$pos))
return(tt)
}
tag="wt_late"
tag="wt_early"
tag="clb19_early"
ids<-list.files(pattern=".bam$","bamFiles")
tags<-unique(unlist(lapply(strsplit(ids,"_"),function(x)paste(x[[1]],x[[2]],sep="_"))))
filter=5
backg=0
for(tag in tags){
print(tag)
ids<-list.files(pattern=".bam$","bamFiles")
ids<-ids[grep(tag,ids)]

fin0<-run(ids[1],flag=0)##It is strand specific
for(id in ids[-1]){
print(id)
temp<-run(id,flag=0)
fin0<-merge(fin0,temp,by="pos",all=T)
}
fin0[is.na(fin0)]<-0; rownames(fin0)<-fin0$pos;fin0$strand<-"+"
fin10<-fin0

fin0<-run(ids[1],flag=16)
for(id in ids[-1]){
print(id)
temp<-run(id,flag=16)
fin0<-merge(fin0,temp,by="pos",all=T)
}
fin0[is.na(fin0)]<-0; rownames(fin0)<-fin0$pos;fin0$strand<-"-"
fin1<-rbind(fin10,fin0)
for(i in grep("control",colnames(fin1))){fin1[fin1[,i]<backg,i]<-backg}
fin1<-fin1[rowSums(fin1[,grep("sample",colnames(fin1))]>filter)>0,]

geno<-as.factor(unlist(lapply(strsplit(ids,"_"),function(x)x[1])))
dev<-as.factor(unlist(lapply(strsplit(ids,"_"),function(x)x[2])))
rep<-as.factor(unlist(lapply(strsplit(ids,"control|sample"),function(x)x[2])))
cond<-as.factor(substr(unlist(lapply(strsplit(ids,"_|\\."),function(x)x[3])),1,3))

colData<-data.frame(cond=cond,rep=rep)
dds <- DESeqDataSetFromMatrix(countData = as.matrix(fin1[,grep("sample|control",colnames(fin1))]),colData = colData,design = ~ cond+rep)
dds1 <- DESeq(dds,fitType="local")
res<-as.data.frame(results(dds1,contrast=c("cond","sam","con"),altHypothesis="greater"))
vsd <- varianceStabilizingTransformation(dds1,fitType="local")
fin<-as.data.frame(assays(vsd));
res$pos<-fin1$pos;res$strand<-fin1$strand
#write.csv(res,file=paste("nutResolution/",tag,"_res.csv",sep=""))
sig<-res[res$pv<0.05 & !is.na(res$padj),7:8];sig$pos<-sig$pos-1## to get in 0-based coordinates
sig$pos2<-sig$pos+1;sig$chr<-"ChrC";sig$point<-".";sig<-sig[,c(4,1,3,5,5,2)]
write.table(sig,file=paste("nutResolution/beds/",tag,"-Sigpv05.bed",sep=""),sep="\t",col.names=F,row.names=F,quote=F)
sig<-res[res$pv<0.1 & !is.na(res$padj),7:8];sig$pos<-sig$pos-1## to get in 0-based coordinates
sig$pos2<-sig$pos+1;sig$chr<-"ChrC";sig$point<-".";sig<-sig[,c(4,1,3,5,5,2)]
write.table(sig,file=paste("nutResolution/beds/",tag,"-Sigpv1.bed",sep=""),sep="\t",col.names=F,row.names=F,quote=F)
sig<-res[res$pv<0.5 & !is.na(res$padj),7:8];sig$pos<-sig$pos-1## to get in 0-based coordinates
sig$pos2<-sig$pos+1;sig$chr<-"ChrC";sig$point<-".";sig<-sig[,c(4,1,3,5,5,2)]
write.table(sig,file=paste("nutResolution/beds/",tag,"-Sigpv5.bed",sep=""),sep="\t",col.names=F,row.names=F,quote=F)

#sig<-res[res$padj<0.05 & !is.na(res$padj),7:8];sig$pos<-sig$pos-1## to get in 0-based coordinates
#sig$pos2<-sig$pos+1
#sig$chr<-"ChrC";sig$point<-".";sig<-sig[,c(4,1,3,5,5,2)]
#write.table(sig,file=paste("nutResolution/beds/",tag,"-SigFDR.bed",sep=""),sep="\t",col.names=F,row.names=F,quote=F)

#write.table(data.frame(as.integer(rownames(vsd)),res[,2]),file=file,append=TRUE,sep="\t",col.names=FALSE,row.names=FALSE,quote=F);

file=paste("nutResolution/wigs/log2FC-",tag,".wig",sep="")
write(paste(sep='',"variableStep  chrom=ChrC  span=",1),file=file,1,append=F,sep="\t");
write.table(data.frame(res$pos,res$log2FoldChange),file=file,append=TRUE,sep="\t",col.names=FALSE,row.names=FALSE,quote=F);

file=paste("nutResolution/wigs/log10FDR-",tag,".wig",sep="")
write(paste(sep='',"variableStep  chrom=ChrC  span=",1),file=file,1,append=F,sep="\t");
pv<- -log10(res$padj)*10
write.table(data.frame(res$pos[!is.na(pv)],pv[!is.na(pv)]),file=file,append=TRUE,sep="\t",col.names=FALSE,row.names=FALSE,quote=F);
for(i in 3:dim(fin)[2] ){
print(i)
file=paste("nutResolution/wigs/",sub(".bam","",colnames(fin)[i]),".wig",sep="")
write(paste(sep='',"variableStep  chrom=ChrC  span=",1),file=file,1,append=F,sep="\t");
write.table(data.frame(fin1$pos,fin[,i]),file=file,append=TRUE,sep="\t",col.names=FALSE,row.names=FALSE,quote=F);
}
}

