# This folder is for finding CAGE enhancers and annotating them
1) use commands in "script" to transfer .BED or .BedGraph files into CAGEfightR-recognisable BigWig files (single nucleotide based);
2) run R script enhancer.r to find enhancers (this script will use all .bw files in the working directory as input, make sure no redundant files are there)
