#!/bin/sh
#PBS -N "map"
#PBS -l select=1:ncpus=8:mem=50GB
#PBS -l walltime=8:00:00

module load bwa/0.7.15
module load samtools

Wdir=/path/to/inputs 
cd $Wdir

#map 
for fastqfile in `ls fastq_files/*_merged.fastq.gz`; do
	bwa mem -t 8 Homo_sapiens.GRCh38.dna.primary_assembly.fa  $fastqfile > ${fastqfile}.se.sam
	samtools view -bS ${fastqfile}.se.sam > ${fastqfile}.se.bam
    samtools view -q10 -b ${fastqfile}.se.bam >${fastqfile}.se.q10.bam
	samtools sort ${fastqfile}.se.q10.bam>${fastqfile}.se.q10.sorted.bam
	samtools index ${fastqfile}.se.q10.sorted.bam
done 

#get mapping statistics all reads 
for file in *se.bam; do newname1=`basename $file .fastq.gz.se.bam`; echo $file >>${newname1}.mapped.qc; samtools flagstat ${file} >>${newname1}.mapped.qc;  done 
#get uniquely mapped CAGE tags with a minimal mapping quality 10
for file in *q10.sorted.bam; do newname2=`basename $file .fastq.gz.se.q10.bam`; echo $file >>${newname2}.mapped.qc; samtools flagstat ${file} >>${newname2}.mapped.qc;  done 
