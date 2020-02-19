setwd("/fast/AG_Ohler/jmuino/Julia/Tex-experiment/motif")

require(Biostrings)
q<-read.csv("../CSAR-Style/linearmodel/EarlyFCvsWT-pamk.csv",stringsAsFactors=F)
seq<-readDNAStringSet("/scratch/AG_Ohler/jmuino/genomes/TAIR9/chr-TAIR9.fas")[[7]]
temp<-strsplit(q$X," ")
q$pos<-as.integer(unlist(lapply(temp,function(x)x[[1]])))
q$strand<-(unlist(lapply(temp,function(x)x[[2]])))

bss<-DNAStringSet(seq,q$pos-25,q$pos+25)
bss[q$strand=="-"]<-reverseComplement(bss[q$strand=="-"])
names(bss)<-q$X
#writeXStringSet(bss,file="AllTSSs.fasta")
for(i in 1:6){
#writeXStringSet(bss[q$kmers==i],file=paste("TSSs-cluster",i,".fasta",sep=""))
}

##kmer
pdf("Kmer.pdf")
for(cc in 1:6){
par(mfrow=c(1,1))
bss<-DNAStringSet(seq,q$pos-15,q$pos+15)
bss[q$strand=="-"]<-reverseComplement(bss[q$strand=="-"])
names(bss)<-q$X
wd=4

c3=oligonucleotideFrequency(bss[q$kmers==cc],width=wd,simplify.as="collapsed",as.prob=T)
all=oligonucleotideFrequency(bss[q$kmers!=cc],width=wd,simplify.as="collapsed",as.prob=T)
w=which(c3-all>0.006)
plot(all,c3,ylab="proportion of kmers in Cluster", xlab="Proportion of kmers in other TSSs",xlim=c(0,0.030),ylim=c(0,0.030),main=paste("Cluster",cc));
points(all[w],c3[w],lwd=3,col="red")
text(all[w],c3[w],pos=1,label=names(all[w]),col="red")
lines(0:1,0:1)

kmer=names(w)
pos<- -100:100
res<-rep(NA,length(pos))
for(i in 1:length(pos)){
bss<-DNAStringSet(seq,(q$pos+(-15+pos[i])),(q$pos+(15+pos[i])))
c3=oligonucleotideFrequency(bss[q$kmers==cc],width=wd,as.prob=F)
res[i]<-sum(rowSums(c3[,kmer])>0)/sum(rowSums(c3[,kmer])> -1)
}
plot(pos,res,xlab="relative position to TSS",ylab="Proportion of TSSs containing Kmer combination",main=paste("Cluster",cc));abline(v=0)
}
dev.off()

##Hetamap
res=t(as.data.frame(lapply(1:6,function(i)oligonucleotideFrequency(bss[q$kmers==i],width=wd,simplify.as="collapsed",as.prob=T))))
rownames(res)<-1:6

require(pheatmap)
pheatmap(res[1:5,colSums(res)>.05])
