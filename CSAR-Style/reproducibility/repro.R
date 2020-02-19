setwd("/fast/AG_Ohler/jmuino/Julia/Tex-experiment/CSAR-Style/reproducibility")
require(pheatmap)
ids<-list.files(pattern=".counts$","../counts")
res<-read.table(paste("../counts/",ids[1],sep=""),header=T)
for(id in ids[-1]){
print(id)
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
cs<-counts[,grep("sample",colnames(counts))]
colnames(cs)<-sub("_sample..bam.counts","",colnames(cs))

cc<-cor(log2(cs+1))
hist(log10(colSums(cs)),breaks=10)
pdf("HeatmapCorCandidateRegions.pdf")
pheatmap(as.matrix(cc))
dev.off()
plot(log10(colSums(cs)),rowMeans(cc))


ww<-which(rowSums(counts[,grep("sample",colnames(counts))]>10)>0)
#coor<-coor[ww,];counts<-counts[ww,]#This I didn't use for Julia
for(i in grep("control",colnames(counts))){counts[counts[,i]<10,i]<-10} 
