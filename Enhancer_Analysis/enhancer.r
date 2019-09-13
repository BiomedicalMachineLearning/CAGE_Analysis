library(CAGEfightR)
genomeInfo <- SeqinfoForUCSCGenome("hg38")

# Quantify CAGE TSSs
bw_plus <- BigWigFileList(c('ESC_0980_S1_R1_001.ctss_all_q10_CSC.pos.bw', 
                            'ESC_1041_S2_R1_001.ctss_all_q10_CSC.pos.bw', 
                            'ESC_1067_S3_R1_001.ctss_all_q10_CSC.pos.bw', 
                            'ESC_1156_S4_R1_001.ctss_all_q10_CSC.pos.bw',
                            'ESC_1271_S5_R1_001.ctss_all_q10_CSC.pos.bw', 
                            'ESC_1369_S6_R1_001.ctss_all_q10_CSC.pos.bw',
                            'ESC_1412_S7_R1_001.ctss_all_q10_CSC.pos.bw', 
                            'ESC_1575_S8_R1_001.ctss_all_q10_CSC.pos.bw'))
bw_minus <- BigWigFileList(c('ESC_0980_S1_R1_001.ctss_all_q10_CSC.neg.bw', 
                             'ESC_1041_S2_R1_001.ctss_all_q10_CSC.neg.bw', 
                             'ESC_1067_S3_R1_001.ctss_all_q10_CSC.neg.bw', 
                             'ESC_1156_S4_R1_001.ctss_all_q10_CSC.neg.bw',
                             'ESC_1271_S5_R1_001.ctss_all_q10_CSC.neg.bw', 
                             'ESC_1369_S6_R1_001.ctss_all_q10_CSC.neg.bw',
                             'ESC_1412_S7_R1_001.ctss_all_q10_CSC.neg.bw', 
                             'ESC_1575_S8_R1_001.ctss_all_q10_CSC.neg.bw'))
names(bw_plus) <- c('0980','1041','1067','1156','1271','1369','1412','1575')
names(bw_minus) <- c('0980','1041','1067','1156','1271','1369','1412','1575')
CTSSs <- quantifyCTSSs(plusStrand = bw_plus, minusStrand = bw_minus, genome = genomeInfo)

# Find enhancers
enhancers <- quickEnhancers(CTSSs)

# Annotate enhancers
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
enhancers <- assignTxType(enhancers, txModels=txdb)
enhancers <- subset(enhancers, txType %in% c("intergenic", "intron"))

saveRDS(enhancers, 'enhancers.RDS')
