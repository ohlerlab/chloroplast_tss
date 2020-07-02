library(here)
library(CSAR)
library(Rsamtools)
library(foreach)
library(doParallel)
 "rpotmp_early_sample2"
 setwd(here("CSAR-Style"))

 bamfolder = "/Users/dharnet/Dropbox/peak_data/"
 outfolder='bed-Step1/'
 stopifnot(file.exists(bamfolder))
 dir.create(outfolder,showWarn=FALSE)
 

 exportenv <- function(){
       fenv <- parent.frame()
       for(i in ls(envir=fenv)){
           assign(i,get(i,envir=fenv),envir=.GlobalEnv)
       }
   }
 
#BiocManager::install(c('doParallel','CSAR','GenomeInfoDb','foreach'))
#This is jose's function that runs CSAR on a single bam file 
run<-function(name,rep,n=-1,b=2,w=1L,nper=100){
  bam<-paste(bamfolder,name,".bam",sep="")
  stopifnot(file.exists(bam))
  #read in the bam data as a list
  bam <- scanBam(bam,param=param)[[1]];names(bam)<-c("chr","strand","pos","lengthRead","mq");
  tail(as.data.frame(bam))
  #copy this data for some reason

  bam2<-as.data.frame(bam);bam2$Nhits=1;
  #! filter out the bam lines that have no value for 'chr' (unmapped reads)
  bam2 <- bam2[!is.na(bam2[['chr']]),]
  bam2$pos<-bam2$pos-1##to put it in 0-based coordinates (as genome browser)
  bam2$pos[bam2$strand=="-"]<-bam2$pos[bam2$strand=="-"]+bam2$lengthRead[bam2$strand=="-"]-1

  if(is.null(bam))return(NA)
  ip=bam2
  #Now get the control bam for this sample
  name<-sub("sample","control",name)
  controlbam <- paste(bamfolder,name,".bam",sep="")
  bam <- scanBam(controlbam,param=param)[[1]];names(bam)<-c("chr","strand","pos","lengthRead","mq");
  bam1<-as.data.frame(bam);bam1$Nhits=1;
  #! filter out the bam lines that have no value for 'chr' (unmapped reads)
  bam1 <- bam1[!is.na(bam1[['chr']]),]
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
    #using 'bam2' which is the ip, create one of these objects CSAR wants describing read overlaps
    nhitsS<-     mappedReads2Nhits(
      bam2[!is.na(bam2$pos)  & bam2$lengthRead>minlength & !is.element(bam2$pos,rrna) ,], 
      "hitS", w = w, considerStrand = strand, uniquelyMapped = F, uniquePosition = F,chr=chr,chrL=chrL
    )
    nhitsS$digits=0
    score2wig(nhitsS,file=paste(name,"_sample","Strand",strand,n,"b",b,"w",w,"minl",minlength,".wig",sep=""),times=10000)
    #Do the same for controls
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
      write.table(win1,file=paste(outfolder,"TestRatio-",name,"Strand",strand,n,"b",b,"w",w,"minl",minlength,".bed",sep=""),quote=F,sep="\t",col.names=F,row.names=F)
      }
      
    }
  }
  gc()
  return("done")
  
}

ids<-list.files(pattern=".bam$","../bamFiles")
ids<-list.files(pattern=".bam$",bamfolder)
ids<-unique(sub(".bam","",ids))
ids<-ids[grep("sample",ids)]
#run(name="wt_late_sample1",rep=1,nper=2)
cores=detectCores()
cl <- makeCluster(4) #not to overload your computer
registerDoParallel(cl)

what <- c("qwidth", "strand","rname", "pos","mapq")
param <- ScanBamParam( what = what)
for(nam in ids){print(nam);run(name=nam,rep=1)}
stopCluster(cl)
 
# ##merge with script at run,sh
#   #run in shell
# cd /scratch/AG_Ohler/jmuino/Julia/Tex-experiment/CSAR-Style/bed-Step1
# rm merged*.bed
# cat *.bed >> merged.bed
# sort  -k1,1 -k2,2n -i merged.bed > mergeds.bed
#  mergeBed -i mergeds.bed -iobuf 5G -s  >merged-1b2w1minl0.bed

#  ##Look in DESEq2-analysis.R how to get fasta
# bedtools getfasta -s -fi /scratch/AG_Ohler/jmuino/genomes/TAIR9/chr-TAIR9.fas -bed merged-1b2w1minl0.bed > merged-1b2w1minl0.bed

system(paste0('cat ',paste(list.files(outfolder,full=T),collapse=' '),'>> merged.bed'))
system('sort  -k1,1 -k2,2n -i merged.bed > mergeds.bed')
if(system('mergeBed')==0)stop('need to install mergeBed with e.g. conda install bedtools ')
system('mergeBed -i mergeds.bed -iobuf 5G -s  >merged-1b2w1minl0.bed')

q1<-read.table(file.path(outfolder,"merged-1b2w1minl0.bed"));q1$V1<-as.character(q1$V1)
#q1<-rbind(q1,c("spike-in",1,485,"+"))
gtf<-data.frame(q1$V1,"CSAR","gene",q1$V2,q1$V3,"1000",q1$V4,".",paste("gene_id \"",q1$V2,"\"; transcript_id \"",q1$V2,"\";",sep=""))
gtf[,9]<-as.character(gtf[,9]);gtf[,7]<-as.character(gtf[,7])
#gtf[length(gtf[,1]),9]<-paste("gene_id \"spike-inJose\"; transcript_id \"spike-inJose\";",sep="")
write.table(gtf,file=file.path(outfolder,"merged-1b2w1minl0.gtf"),quote=F,row.names=F,col.names=F,sep="\t")

q1<-read.table("merged-1b2w1minl0.bed");q1<-q1[q1$V1=="ChrC",]
gff<-read.table("ext_data/Araport11_GFF3_genes_transposons.201606.gtf",sep="\t")
gff1<-gff[gff$V1=="ChrC" & gff$V3=="gene",]
disF<-sapply(q1$V3[q1$V4=="+"],function(x){temp<-gff1$V4[gff1$V7=="+"]-x;c(as.character(gff1$V9[gff1$V7=="+"])[temp==min(temp[temp>0 & temp <1000])],NA)[1]})
disR<-sapply(q1$V2[q1$V4=="-"],function(x){temp<-gff1$V5[gff1$V7=="-"]-x;c(as.character(gff1$V9[gff1$V7=="-"])[temp==max(temp[temp<0 & temp > -1000])],NA)[1]})
fin<-data.frame(id=c(q1$V2[q1$V4=="+"],q1$V2[q1$V4=="-"]),gene=c(disF,disR))
fin$gene<- unlist(lapply(strsplit(as.character(fin$gene),split=" |;"),function(x)x[2]))
write.csv(fin,file="annotationTSS1kb.csv")