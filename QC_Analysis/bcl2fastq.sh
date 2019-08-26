#!/bin/sh
#PBS -N "Name"
#PBS -l select=1:ncpus=16:mem=32GB
#PBS -l walltime=6:00:00

wdir=/path/to/workingdir
cd $wdir

module load bcl2fastq/2.17
bcl2fastq --runfolder-dir ./path/to/BCL/Files  -p 16 --output-dir $wdir/fastq_files --no-lane-splitting
