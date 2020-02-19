setwd("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment")
library(CSAR)
library(Rsamtools)
library(foreach)
library(doParallel)

run<-function(name,rep,n=-1,b=2,w=1L){
bam <- scanBam(paste("bamFiles/",name,"_sample",rep,".bam",sep=""),param=param)[[1]];names(bam)<-c("chr","strand","pos","lengthRead","mq");bam2<-as.data.frame(bam);bam2$Nhits=1;
ip=bam2
bam <- scanBam(paste("bamFiles/",name,"_control",rep,".bam",sep=""),param=param)[[1]];names(bam)<-c("chr","strand","pos","lengthRead","mq");bam1<-as.data.frame(bam);bam1$Nhits=1;
control<-bam1
for(strand in c("Reverse", "Foward")){
chr=c("ChrC","ChrM")
chrL=c(154478,366924)
minlength=0##with this I eliminate the spike-in in the reverse strand

print(w)
nhitsS<-     mappedReads2Nhits(bam2[!is.na(bam2$pos)  & bam2$lengthRead>minlength,], "hitS", w = w, considerStrand = strand, uniquelyMapped = F, uniquePosition = F,chr=chr,chrL=chrL)
nhitsS$digits=0
score2wig(nhitsS,file=paste(name,"_sample",rep,"Strand",strand,n,"b",b,"w",w,"minl",minlength,".wig",sep=""),times=10000)
nhitsC<-     mappedReads2Nhits(bam1[!is.na(bam1$pos) & bam1$lengthRead>minlength ,], "hitC", w = w, considerStrand = strand, uniquelyMapped = F, uniquePosition = F,chr=chr,chrL=chrL)
nhitsC$digits=0
score2wig(nhitsC,file=paste(name,"_control",rep,"Strand",strand,n,"b",b,"w",w,"minl",minlength,".wig",sep=""),times=10000)

test3<-ChIPseqScore(control=nhitsC,sample=nhitsS,file="test3",norm= n,backg=b,times=1000000,test="Ratio")
score2wig(test3,file=paste("TestRatio-",name,"-r",rep,"Strand",strand,n,"b",b,"w",w,"minl",minlength,".wig",sep=""),times=1000000,t=1)
win<-sigWin(test3,g=3)
if(!is.null(win)){
foreach(i=1:150, .packages='CSAR')%dopar%{permutatedWinScores(nn=i,sample=bam2[!is.na(bam2$pos) & bam2$lengthRead>minlength ,],control=bam1[!is.na(bam1$pos) & bam1$lengthRead>minlength ,],fileOutput="per",
backg=b,chr = chr, chrL = chrL, w = w, considerStrand = strand, uniquelyMapped = F, uniquePosition = F,norm= -1,g=5)}
nulldist<-getPermutatedWinScores(file="per",nn=1:150)
save(nulldist,win,file=paste("Perm-",name,"-r",rep,strand,n,"b",b,"w",w,"minl",minlength,".RD",sep=""))
load(file=paste("Perm-",name,"-r",rep,strand,n,"b",b,"w",w,"minl",minlength,".RD",sep=""))
tr=getThreshold(winscores=values(win)$score,permutatedScores=nulldist,FDR=.05)
win1<-as.data.frame(win)[values(win)$score>tr$threshold,1:3]
win1$name<-".";win1$score<-"."
if(strand=="Foward"){win1$strand<-"+"}
if(strand=="Reverse"){win1$strand<-"-"}
write.table(win1,file=paste("TestRatio-",name,"-r",rep,"Strand",strand,n,"b",b,"w",w,"minl",minlength,".bed",sep=""),quote=F,sep="\t",col.names=F,row.names=F)
}
}
}


ids<-list.files(pattern=".bam","bamFiles")
ids<-unique(sub("..bam","",ids))
ids<-unique(sub("_sample","",ids))
ids<-unique(sub("_control","",ids))

cores=detectCores()
cl <- makeCluster(25) #not to overload your computer
registerDoParallel(cl)

##check ATCG00450.1 as control
##No QC
what <- c("qwidth", "strand","rname", "pos","mapq")
param <- ScanBamParam( what = what)

for(nam in ids){for(y in 1:2){print(nam);run(name=nam,rep=y)}}
run(name="wt",rep=3)

stopCluster(cl)

##merge with script at run,sh

q1<-read.table("merged-1b2w1minl0.bed");q1$V1<-as.character(q1$V1)
q1<-rbind(q1,c("spike-in",1,485,"+"))
gtf<-data.frame(q1$V1,"CSAR","gene",q1$V2,q1$V3,"1000",q1$V4,".",paste("gene_id \"",q1$V2,"\"; transcript_id \"",q1$V2,"\";",sep=""))
gtf[,9]<-as.character(gtf[,9]);gtf[,7]<-as.character(gtf[,7])
gtf[length(gtf[,1]),9]<-paste("gene_id \"spike-inJose\"; transcript_id \"spike-inJose\";",sep="")
write.table(gtf,file="merged-1b2w1minl0.gtf",quote=F,row.names=F,col.names=F,sep="\t")

q1<-read.table("merged-1b2w1minl0.bed");q1<-q1[q1$V1=="ChrC",]
gff<-read.table("/scratch/AG_Ohler/jmuino/genomes/TAIR9/Araport11_GFF3_genes_transposons.201606.gtf",sep="\t")
gff1<-gff[gff$V1=="ChrC" & gff$V3=="gene",]
disF<-sapply(q1$V3[q1$V4=="+"],function(x){temp<-gff1$V4[gff1$V7=="+"]-x;c(as.character(gff1$V9[gff1$V7=="+"])[temp==min(temp[temp>0 & temp <1000])],NA)[1]})
disR<-sapply(q1$V2[q1$V4=="-"],function(x){temp<-gff1$V5[gff1$V7=="-"]-x;c(as.character(gff1$V9[gff1$V7=="-"])[temp==max(temp[temp<0 & temp > -1000])],NA)[1]})
fin<-data.frame(id=c(q1$V2[q1$V4=="+"],q1$V2[q1$V4=="-"]),gene=c(disF,disR))
fin$gene<- unlist(lapply(strsplit(as.character(fin$gene),split=" |;"),function(x)x[2]))
write.csv(fin,file="annotationTSS1kb.csv")