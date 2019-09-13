library(CAGEfightR)
genomeInfo <- SeqinfoForUCSCGenome("hg38")
files <- list.files()

# Quantify CAGE TSSs

# this script will read all .bw files in the working directory, make sure no redundant files are there
bw_plus <- files[grep('.pos.bw', files)]
bw_minus <- files[grep('.neg.bw', files)]

# users can customise how many characters to be taken as sample names
names(bw_plus) <- substr(bw_plus,5,8)
names(bw_minus) <- substr(bw_minus,5,8)
CTSSs <- quantifyCTSSs(plusStrand = bw_plus, minusStrand = bw_minus, genome = genomeInfo)

# Find enhancers
enhancers <- quickEnhancers(CTSSs)

# Annotate enhancers
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
enhancers <- assignTxType(enhancers, txModels=txdb)
enhancers <- subset(enhancers, txType %in% c("intergenic", "intron"))

saveRDS(enhancers, 'enhancers.RDS')
