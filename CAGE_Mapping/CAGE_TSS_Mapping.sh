#!/bin/sh
#PBS -q normal
#PBS -N "CAGE_ctss"
#PBS -l select=1:ncpus=1:mem=50GB
#PBS -l walltime=8:00:00

module load samtools
module load bedtools2/2.24.0

Wdir=/path/to/inputs 
cd $Wdir


#install moirai programs as instructed in the readme.md, run bam_to_ctssbed.sh 

export PATH=/shares/powell/quan/CAGE/Software/moirai20140528/bin:$PATH
for file in *.gz.se.bam; do name1=`basename $file .fastq.gz.se.bam`; bam_to_ctssbed.sh 10 $file >${name1}.ctss_all_q10_cardio.bed; done

#format  bed files ready for CAGEr Analysis
for file in *.ctss_all_q10_cardio.bed; do name2=`basename $file .bed`; awk '{OFS="\t"}{print "chr"$1, $2, $6,$5}' $file >${name2}.cageRctss; echo $file; done

