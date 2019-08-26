#!/bin/sh
#PBS -N "CAGE_qc"
#PBS -l select=1:ncpus=8:mem=50GB
#PBS -l walltime=4:00:00


module load fastqc/0.11.4 
conda activate /your/env #with multiqc active


RUN_PATH=/path/to/inputs 
cd $RUN_PATH
for file in $(ls $RUN_PATH)
do
    SAMPLE=`basename $file`
    fastqc -t 8 ${SAMPLE} -o /path/to/outputs
done

multiqc .

#Note:  multiqc duplicated reads are applied for mapped files like bam files not single-end read fastq files                                                                                                                                     

