  #run in shell
cd /scratch/AG_Ohler/jmuino/Julia/Tex-experiment
cat TestRatio*.bed >> merged.bed

sort  -k1,1 -k2,2n -i merged.bed > mergeds.bed
 mergeBed -i mergeds.bed -iobuf 5G -s  >merged-1b2w1minl0.bed

 ##Look in DESEq2-analysis.R how to get fasta
bedtools getfasta -s -fi /scratch/AG_Ohler/jmuino/genomes/TAIR9/chr-TAIR9.fas -bed merged-1b2w1minl0.bed > merged-1b2w1minl0.bed

##Add spikein in file
[1] "pep"    "rpotmp" "rpotp"  "wt"
##
featureCounts bamFiles/wt_control1.bam -a merged-1b2w1minl0.gtf -o counts/wt_control1.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5
featureCounts bamFiles/wt_control2.bam -a merged-1b2w1minl0.gtf -o counts/wt_control2.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5
featureCounts bamFiles/wt_sample1.bam -a merged-1b2w1minl0.gtf -o counts/wt_sample1.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5
featureCounts bamFiles/wt_sample2.bam -a merged-1b2w1minl0.gtf -o counts/wt_sample2.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5
featureCounts bamFiles/wt_control3.bam -a merged-1b2w1minl0.gtf -o counts/wt_control3.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5
featureCounts bamFiles/wt_sample3.bam -a merged-1b2w1minl0.gtf -o counts/wt_sample3.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5

featureCounts bamFiles/pep_control1.bam -a merged-1b2w1minl0.gtf -o counts/pep_control1.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5
featureCounts bamFiles/pep_control2.bam -a merged-1b2w1minl0.gtf -o counts/pep_control2.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5
featureCounts bamFiles/pep_sample1.bam -a merged-1b2w1minl0.gtf -o counts/pep_sample1.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5
featureCounts bamFiles/pep_sample2.bam -a merged-1b2w1minl0.gtf -o counts/pep_sample2.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5

featureCounts bamFiles/rpotmp_control1.bam -a merged-1b2w1minl0.gtf -o counts/rpotmp_control1.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5
featureCounts bamFiles/rpotmp_control2.bam -a merged-1b2w1minl0.gtf -o counts/rpotmp_control2.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5
featureCounts bamFiles/rpotmp_sample1.bam -a merged-1b2w1minl0.gtf -o counts/rpotmp_sample1.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5
featureCounts bamFiles/rpotmp_sample2.bam -a merged-1b2w1minl0.gtf -o counts/rpotmp_sample2.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5

featureCounts bamFiles/rpotp_control1.bam -a merged-1b2w1minl0.gtf -o counts/rpotp_control1.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5
featureCounts bamFiles/rpotp_control2.bam -a merged-1b2w1minl0.gtf -o counts/rpotp_control2.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5
featureCounts bamFiles/rpotp_sample1.bam -a merged-1b2w1minl0.gtf -o counts/rpotp_sample1.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5
featureCounts bamFiles/rpotp_sample2.bam -a merged-1b2w1minl0.gtf -o counts/rpotp_sample2.counts -t transcript -g transcript_id -T 20 -s 1 -M --read2pos 5

