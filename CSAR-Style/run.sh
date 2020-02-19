  #run in shell
cd /scratch/AG_Ohler/jmuino/Julia/Tex-experiment/CSAR-Style/bed-Step1
rm merged*.bed
cat *.bed >> merged.bed
sort  -k1,1 -k2,2n -i merged.bed > mergeds.bed
 mergeBed -i mergeds.bed -iobuf 5G -s  >merged-1b2w1minl0.bed

 ##Look in DESEq2-analysis.R how to get fasta
bedtools getfasta -s -fi /scratch/AG_Ohler/jmuino/genomes/TAIR9/chr-TAIR9.fas -bed merged-1b2w1minl0.bed > merged-1b2w1minl0.bed


