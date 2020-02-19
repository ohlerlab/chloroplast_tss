#!/bin/bash

#$ -t 1-62
#$ -N Julia
#$ -o /scratch/AG_Ohler/jmuino/Julia/Tex-experiment/bamFiles/logs
#$ -e /scratch/AG_Ohler/jmuino/Julia/Tex-experiment/bamFiles/logs

export PATH="/home/jmuinoa/.guix-profile/bin${PATH:+:}$PATH"
cd /scratch/AG_Ohler/jmuino/Julia/Tex-experiment/CSAR-Style


FilewithIDS=/scratch/AG_Ohler/jmuino/Julia/Tex-experiment/bamFiles/ids.txt
ID=$(awk "NR==$SGE_TASK_ID"  $FilewithIDS)

mkdir ${ID}temp
featureCounts /scratch/AG_Ohler/jmuino/Julia/Tex-experiment/bamFiles/${ID} -a /scratch/AG_Ohler/jmuino/Julia/Tex-experiment/CSAR-Style/bed-Step1/merged-1b2w1minl0.gtf -o /scratch/AG_Ohler/jmuino/Julia/Tex-experiment/CSAR-Style/counts/${ID}.counts -t gene -g transcript_id -T 20 -s 1 -M --read2pos 5
rm -r ${ID}temp
