setwd("/fast/AG_Ohler/jmuino/Julia/Tex-experiment/CSAR-Style/linearmodel")
require(DESeq2)
require("pheatmap")
require(RColorBrewer)
FC=1.5
res0<-read.csv(paste("../Targets-FDR-lFC",FC,".csv",sep=""))
#res1<-res1[res1$wt_early.Step2.csv<0.05 | res1$wt_late.Step2.csv<0.05,]#adding 5 sept

load("../Counts.RD");
counts<-counts[,-grep("late|control",colnames(counts))];ids<-colnames(counts)
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
clb<-results(dds1,contrast=c("geno","clb19","wt"))
pv0<- -log10(data.frame(rpotp=rtp[,c(5)],rpotmp=rmp[,c(5)],clb19=clb[,c(5)]));rownames(pv0)<-rownames(vsd)
fc<-data.frame(rtp[,c(2)],rmp[,c(2)],clb[,c(2)]);rownames(fc)<-rownames(vsd)
#pv<-pv0*fc/abs(fc)
sigtss<-which(is.element(paste(coor$Start,coor$target,sep="-"),res0$X))
exp<-as.data.frame(pv0[sigtss,]);
exp1<-as.data.frame(fc[sigtss,]);
res1=exp1

fdr<-read.csv(paste("../Targets-FDR-lFC",FC,".csv",sep=""));rownames(fdr)<-fdr$X;fdr$X<-NULL
order<-sapply(rownames(res1),function(x)grep(x,fdr$pos))
fdr1<-fdr[order,]


####
if(F){
data<-assay(vsd)
data<-data[names(order),]##Only significants
colnames(data)<-sub(".bam.counts","Normalizedcounts.BedGraph",colnames(data))
for(i in 1:12){
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
anno<-read.csv("/scratch/AG_Ohler/jmuino/CRPC-Marie/PlastidAnnotation.csv",stringsAsFactors=F)
symbol<-sapply(as.character(fdr1$target),function(x)c(anno$symbol[anno$TAIR_id==x],NA)[1])
type<-sapply(as.character(fdr1$target),function(x)c(anno$type[anno$TAIR_id==x],NA)[1])
rownames(res1)<-paste(fdr1$pos,fdr1$strand,symbol)
type[is.na(type)]<-"No_assigned"
rowData<-data.frame(type);rownames(rowData)<-rownames(res1)
ann_colors=list(type=c(other="light pink",Photosynthesis="green",Ribosome="burlywood4",RNA_polymerase="dark red",rRNA="red",tRNA="dark blue",No_assigned="yellow"))

colnames(res1)<-substr(colnames(res1),1,3)


library(fpc)




pdf(paste("EarlyFCvsWT-Kmeans-lFC",FC,".pdf",sep=""))
for(nc in 2:10){
par(mfrow=c(4,3))
cluster=kmeans(res1[,1:3],nc,iter.max=10^6,nstart=100)$cluster
 table(cluster)
plotcluster(res1[,1:3],cluster,main=paste("psbA is cluster", as.integer(cluster["1515 - PSBA"])))
boxplot(res1[,1:3],main="All");abline(v=0,col="green",lwd=3);abline(h=0,lwd=3,col="green");abline(h=0,lwd=3,col="green")
for (i in 1:nc){boxplot(res1[cluster==i,1:3],main=paste("cluster",i,"\n n=",sum(cluster==i)));abline(h=0,lwd=3,col="green")}
write.csv(data.frame(res1,kmers=cluster),file=paste("EarlyFCvsWT-Kmeans",nc,"-FC",FC,".csv",sep=""))
}
dev.off()


pdf(paste("EarlyFCvsWT-pamk-lFC",FC,".pdf",sep=""))
sig<- rowSums(exp> -log10(0.05),na.rm=T)>0
res2=res1#[sig,]
cluster1=pamk(res2[,1:3],3:15,scaling=F)[[1]][[3]]
 table(cluster1)
 par(mfrow=c(3,3))
plotcluster(res2[,1:3],cluster1,main=paste("psbA is cluster", as.integer(cluster1["1515 - PSBA"])))
boxplot(res2[,1:3],main="All");abline(h=0,lwd=3,col="green")
for (i in 1:max(cluster1)){boxplot(res2[cluster1==i,1:3],main=paste("cluster",i,"\n n=",sum(cluster1==i)));abline(h=0,lwd=3,col="green")}
write.csv(data.frame(res2,kmers=cluster1),file=paste("EarlyFCvsWT-pamk-FC",FC,".csv",sep=""))
dev.off()






pdf("EarlyFCvsWT-heatmap-lFC1.pdf")
out<-pheatmap(res1[,1:3],scale="none", cluster_rows=T, show_rownames=T,cluster_cols=T, annotation_row=rowData,annotation_colors=ann_colors,fontsize_row=3)
for(nc in 2:10){
par(mfrow=c(4,3))
cluster<-cutree(out$tree_row, k=nc)
 table(cluster)
plotcluster(res1[,1:3],cluster,main=paste("psbA is cluster", as.integer(cluster["1515 - PSBA"])))
boxplot(res1[,1:3],main="All");abline(v=0,col="green",lwd=3);abline(h=0,lwd=3,col="green");abline(h=0,lwd=3,col="green")
for (i in 1:nc){boxplot(res1[cluster==i,1:3],main=paste("cluster",i,"\n n=",sum(cluster==i)));abline(h=0,lwd=3,col="green")}
write.csv(data.frame(res1,kmers=cluster),file=paste("EarlyFCvsWT-Heatmap",nc,".csv",sep=""))
}
dev.off()

pdf("EarlyFCvsWT-Mclust-lFC1.pdf")
library(mclust)
cluster=Mclust(res1[,1:3],G=3:9)$classification
 table(cluster)
 par(mfrow=c(3,3))
plotcluster(res1[,1:3],cluster,main=paste("psbA is cluster", as.integer(cluster["1515 - PSBA"])))
boxplot(res1[,1:3],main="All")
for (i in 1:4)boxplot(res1[cluster==i,1:3],main=paste("cluster",i,"\n n=",sum(cluster==i)))
dev.off()

write.csv(data.frame(res1,kmers=cluster,pamk=cluster1),file="EarlyFCvsWT-2Classifiers.csv")

write.csv(res1,file="EarlyFCvsWT-all.csv")
res2=res1;res2[exp> -log10(exp)]<-0
write.csv(res2,file="EarlyFCvsWT-onlysig.csv")





##linear model
res0<-cbind(res1[,5:8]-res1[,1:4],-res1[,9:11]+res1[,12:14],-res1[,15:17]+res1[,18:20],-res1[,21:23]+res1[,24:26])
res0$tss<-rownames(res0)

res<-res0[,c(1,14)];res$type=colnames(res0)[1];colnames(res)[1]<-"counts"
for(i in 2:13){
temp<-res0[,c(i,14)];temp$type=colnames(res0)[i];colnames(temp)[1]<-"counts"
res<-rbind(res,temp)
}
res$type<-as.factor(unlist(lapply(strsplit(res$type,"_"),function(x)x[1])))
res$tss<-as.factor(res$tss)
res$type<-relevel(res$type,ref="wt")
ll<-lm(counts~type,res[grep("100867",res$tss),])## in rpotmp should be down

 summary(ll)
Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  -0.2091     0.4143  -0.505   0.6259
typeptac2     1.5500     0.5481   2.828   0.0198 *
typerpotmp   -0.3111     0.5859  -0.531   0.6083
typerpotp     0.4865     0.5859   0.830   0.4278
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.7176 on 9 degrees of freedom
Multiple R-squared:  0.606,     Adjusted R-squared:  0.4747
F-statistic: 4.615 on 3 and 9 DF,  p-value: 0.03217



##
PEP dependent (genotypes clb19 and ptac2): psbA (1515), psaA (42039), atpB (54665)
NEP dependent (genotypes rpotmp and rpotp): ycf1 (129272), clpP (71934)
rpotmp dependent: rrn16 (100891)
w<-res1[grep("129272|71934|1515|42039|54665|100891",rownames(res1)),]
colnames(w)<-sub("..bam.counts","",colnames(w))
pheatmap(as.matrix(w),scale="row")

ll<-lapply(unique(res$tss),function(x){lm(counts~type,res[res$tss==x,])})
names(ll)<-unique(res$tss)
summary(ll[[grep("1515",names(ll))]])
ll1<-t(as.data.frame(lapply(ll,function(x) {temp<-summary(x)[[4]]; temp1<-temp[,1]; temp1[temp[,4]>0.05]<-0;temp1 })))
rownames(ll1)<-names(ll);res1<-ll1[,2:4]
write.csv(res1,file="LMresults-LateNoCLBvsControl-0isnosignificant.csv")
ll1<-t(as.data.frame(lapply(ll,function(x) {temp<-summary(x)[[4]]; temp1<-temp[,1]; temp1 })))
rownames(ll1)<-names(ll);res1<-ll1[,2:4]
write.csv(res1,file="LMresults-LateNoCLBvsControl.csv")
