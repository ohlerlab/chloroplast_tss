setwd("/scratch/AG_Ohler/jmuino/Julia/Tex-experiment/CSAR-Style")
library(CSAR)
library(Rsamtools)
library(foreach)
library(doParallel)
 "rpotmp_early_sample2"

run<-function(name,rep,n=-1,b=2,w=1L,nper=100){
bam<-NULL
bam <- scanBam(paste("../bamFiles/",name,".bam",sep=""),param=param)[[1]];names(bam)<-c("chr","strand","pos","lengthRead","mq");bam2<-as.data.frame(bam);bam2$Nhits=1;
bam2$pos<-bam2$pos-1##to put it in 0-based coordinates (as genome browser)
bam2$pos[bam2$strand=="-"]<-bam2$pos[bam2$strand=="-"]+bam2$lengthRead[bam2$strand=="-"]-1
if(is.null(bam))return(NA)
ip=bam2
name<-sub("sample","control",name)
bam <- scanBam(paste("../bamFiles/",name,".bam",sep=""),param=param)[[1]];names(bam)<-c("chr","strand","pos","lengthRead","mq");bam1<-as.data.frame(bam);bam1$Nhits=1;
bam1$pos<-bam1$pos-1##to put it in 0-based coordinates (as genome browser)
bam1$pos[bam1$strand=="-"]<-bam1$pos[bam1$strand=="-"]+bam1$lengthRead[bam1$strand=="-"]-1
control<-bam1
for(strand in c("Reverse", "Foward")){
chr=c("ChrC")#,"ChrM")
chrL=c(154478)#,366924)
minlength=0
name<-sub("control","",name)
print(w)
rrna= -1;#if(strand=="Foward"){rrna<-c(101012:102502,104691:107500,107599:107701,107949:108069)}else{rrna<-c(130580:130700,130948:131050,131149:133958,136147:137637)}
nhitsS<-     mappedReads2Nhits(bam2[!is.na(bam2$pos)  & bam2$lengthRead>minlength & !is.element(bam2$pos,rrna) ,], "hitS", w = w, considerStrand = strand, uniquelyMapped = F, uniquePosition = F,chr=chr,chrL=chrL)
nhitsS$digits=0
score2wig(nhitsS,file=paste(name,"_sample","Strand",strand,n,"b",b,"w",w,"minl",minlength,".wig",sep=""),times=10000)
nhitsC<-     mappedReads2Nhits(bam1[!is.na(bam1$pos) & bam1$lengthRead>minlength & !is.element(bam1$pos,rrna) ,], "hitC", w = w, considerStrand = strand, uniquelyMapped = F, uniquePosition = F,chr=chr,chrL=chrL)
nhitsC$digits=0
score2wig(nhitsC,file=paste(name,"_control","Strand",strand,n,"b",b,"w",w,"minl",minlength,".wig",sep=""),times=10000)

test3<-ChIPseqScore(control=nhitsC,sample=nhitsS,file="test3",norm= n,backg=b,times=1000000,test="Ratio")
score2wig(test3,file=paste("wig/TestRatio-",name,"Strand",strand,n,"b",b,"w",w,"minl",minlength,".wig",sep=""),times=1000000,t=1)
win<-sigWin(test3,g=3)
if(!is.null(win)){
foreach(i=1:nper, .packages='CSAR')%dopar%{permutatedWinScores(nn=i,sample=bam2[!is.na(bam2$pos) & bam2$lengthRead>minlength & !is.element(bam2$pos,rrna) ,],control=bam1[!is.na(bam1$pos) & bam1$lengthRead>minlength & !is.element(bam1$pos,rrna),],fileOutput="per",
backg=b,chr = chr, chrL = chrL, w = w, considerStrand = strand, uniquelyMapped = F, uniquePosition = F,norm= n,g=3)}
nulldist<-getPermutatedWinScores(file="per",nn=1:nper)
unlink("*permutatedWin")
save(nulldist,win,file=paste("Perm-",name,strand,n,"b",b,"w",w,"minl",minlength,".RD",sep=""))
load(file=paste("Perm-",name,strand,n,"b",b,"w",w,"minl",minlength,".RD",sep=""))
tr=getThreshold(winscores=values(win)$score,permutatedScores=nulldist,FDR=.05)
win1<-NULL;
win1<-as.data.frame(win)[values(win)$score>tr$threshold,c(1,6,6)]
print(head(win1))
if(dim(win1)[1]>0){
win1$posPeak<-win1$posPeak-5;win1$posPeak.1<-win1$posPeak.1+5;win1$posPeak[win1$posPeak<1]<-1
win1$name<-".";win1$score<-"."
if(strand=="Foward"){win1$strand<-"+"}
if(strand=="Reverse"){win1$strand<-"-"}
write.table(win1,file=paste("bed-Step1/TestRatio-",name,"Strand",strand,n,"b",b,"w",w,"minl",minlength,".bed",sep=""),quote=F,sep="\t",col.names=F,row.names=F)
}

}
}
gc()
return("done")
}

ids<-list.files(pattern=".bam$","../bamFiles")
ids<-unique(sub(".bam","",ids))
ids<-ids[grep("sample",ids)]
#run(name="wt_late_sample1",rep=1,nper=2)
cores=detectCores()
cl <- makeCluster(25) #not to overload your computer
registerDoParallel(cl)

what <- c("qwidth", "strand","rname", "pos","mapq")
param <- ScanBamParam( what = what)
for(nam in ids){print(nam);run(name=nam,rep=1)}
stopCluster(cl)
 
##merge with script at run,sh

q1<-read.table("bed-Step1/merged-1b2w1minl0.bed");q1$V1<-as.character(q1$V1)
#q1<-rbind(q1,c("spike-in",1,485,"+"))
gtf<-data.frame(q1$V1,"CSAR","gene",q1$V2,q1$V3,"1000",q1$V4,".",paste("gene_id \"",q1$V2,"\"; transcript_id \"",q1$V2,"\";",sep=""))
gtf[,9]<-as.character(gtf[,9]);gtf[,7]<-as.character(gtf[,7])
#gtf[length(gtf[,1]),9]<-paste("gene_id \"spike-inJose\"; transcript_id \"spike-inJose\";",sep="")
write.table(gtf,file="bed-Step1/merged-1b2w1minl0.gtf",quote=F,row.names=F,col.names=F,sep="\t")

q1<-read.table("merged-1b2w1minl0.bed");q1<-q1[q1$V1=="ChrC",]
gff<-read.table("/scratch/AG_Ohler/jmuino/genomes/TAIR9/Araport11_GFF3_genes_transposons.201606.gtf",sep="\t")
gff1<-gff[gff$V1=="ChrC" & gff$V3=="gene",]
disF<-sapply(q1$V3[q1$V4=="+"],function(x){temp<-gff1$V4[gff1$V7=="+"]-x;c(as.character(gff1$V9[gff1$V7=="+"])[temp==min(temp[temp>0 & temp <1000])],NA)[1]})
disR<-sapply(q1$V2[q1$V4=="-"],function(x){temp<-gff1$V5[gff1$V7=="-"]-x;c(as.character(gff1$V9[gff1$V7=="-"])[temp==max(temp[temp<0 & temp > -1000])],NA)[1]})
fin<-data.frame(id=c(q1$V2[q1$V4=="+"],q1$V2[q1$V4=="-"]),gene=c(disF,disR))
fin$gene<- unlist(lapply(strsplit(as.character(fin$gene),split=" |;"),function(x)x[2]))
write.csv(fin,file="annotationTSS1kb.csv")