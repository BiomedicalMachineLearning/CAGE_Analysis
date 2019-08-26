# This folder is for quality control of CAGE sequencing data
* Two main steps are shown below 


## Step 1: generate fastq files from Illumina BCL files 
* To generate fastq files and test the quality of the fastq files, use bcl2fastq to split samples using index sequences: barcode-mismatches=0


## Step 2: check the quality of the fastq files
* To generate fastq files and test the quality of the fastq files
	* Install multiqc (run this once)
	* Fastq quality check using [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), followed by [multiQC] (https://multiqc.info/docs/)
	* For setting up multiqc, use python virtual environment 
``` 
conda activate /your/env
pip install multiqc 
```
