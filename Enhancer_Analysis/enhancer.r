library(CAGEfightR)
genomeInfo <- SeqinfoForUCSCGenome("hg38")
files <- list.files()

# Quantify CAGE TSSs

# this script will read all .bw files in the working directory, make sure no redundant files are there
bw_plus <- BigWigFileList(files[grep('.plus.bw', files)])
bw_minus <- BigWigFileList(files[grep('.minus.bw', files)])

# users can customise how many characters to be taken as sample names
names(bw_plus) <- substr(files[grep('.plus.bw', files)],1,9)
names(bw_minus) <- substr(files[grep('.minus.bw', files)],1,9)
CTSSs <- quantifyCTSSs(plusStrand = bw_plus, minusStrand = bw_minus, genome = genomeInfo)

# Find enhancers
enhancers <- quickEnhancers(CTSSs)

# Annotate enhancers
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
enhancers <- assignTxType(enhancers, txModels=txdb)
enhancers <- subset(enhancers, txType %in% c("intergenic", "intron"))

# build bedGraph files
bedgraph <- matrix(,nrow(enhancers),4)
bedgraph[,1] <- sapply(rownames(enhancers), function(x) strsplit(x, ':')[[1]][1])
bedgraph[,2] <- sapply(rownames(enhancers), function(x) strsplit(strsplit(x, ':')[[1]][2], '-')[[1]][1])
bedgraph[,3] <- sapply(rownames(enhancers), function(x) strsplit(strsplit(x, ':')[[1]][2], '-')[[1]][2])
bedgraph[,4] <- rowData(enhancers)[,1]

write.table(bedgraph,'enhancers.bedGraph',sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
saveRDS(enhancers, 'enhancers.RDS')
