#!/bin/bash

#$ -t 1-62
#$ -N Julia
#$ -o /scratch/AG_Ohler/jmuino/Julia/Tex-experiment/bamFiles/logs
#$ -e /scratch/AG_Ohler/jmuino/Julia/Tex-experiment/bamFiles/logs

export PATH="/home/jmuinoa/.guix-profile/bin${PATH:+:}$PATH"
cd /scratch/AG_Ohler/jmuino/Julia/Tex-experiment/bamFiles


FilewithIDS=/scratch/AG_Ohler/jmuino/Julia/Tex-experiment/bamFiles/ids.txt
ID=$(awk "NR==$SGE_TASK_ID"  $FilewithIDS)

#samtools view -b -F 4 original/${ID} > ${ID} ## This parameter filter multimapped reads
#samtools index ${ID}
#samtools idxstats ${ID} > idxstats/${ID}.idxstats
bamCoverage -b ${ID} -o coverage/${ID}-F.bw -bs 1 --filterRNAstrand reverse --skipNonCoveredRegions --normalizeUsingRPKM --extendReads 1
bamCoverage -b ${ID} -o coverage/${ID}-R.bw -bs 1 --filterRNAstrand forward --skipNonCoveredRegions --normalizeUsingRPKM --extendReads 1
#mkdir ${ID}temp
#featureCounts ${ID} -a /scratch/AG_Ohler/jmuino/Julia/Tex-experiment/CSAR-Style/merged-1b2w1minl0.gtf -o counts/${ID}.counts -t gene -g transcript_id -T 20 -s 1 -M --read2pos 5
#rm -r ${ID}temp