library("MultiAssayExperiment")
library("SummarizedExperiment")
library(CAGEr)
library("BSgenome.Hsapiens.UCSC.hg38")
name_prefix="K"
sub_string=".ctss_q10.cageRctss"

inputFiles=list.files("./", pattern=paste0(name_prefix,".*.cageRctss"))
ce <- CAGEexp( genomeName     = "BSgenome.Hsapiens.UCSC.hg38"
             , inputFiles     = inputFiles
             , inputFilesType = "ctss"
             , sampleLabels   = gsub(sub_string, "", inputFiles)
)

#show sample label 
sampleLabels(ce)

#read the files, define CTSS, and count reads for CTSS 
getCTSS(ce)

#extract raw tag couts 
tag_counts_SEobject <- CTSStagCountSE(ce)
tag_counts_df <- assay(tag_counts)

#CAGEr QC result 1: how many CTSS per samples 
NonZero_CTSS<- apply(as.matrix(tag_counts_df), 2, function(x){sum(x>0)})
hist(NonZero_CTSS[-1]) #plot total number of CAGE CTSS 

#get coordinates for CTSSs 
CTSScoordinatesGR(ce)

#get read counts for each CTSS in each sample 
CTSStagCountDF(ce)

#CAGEr QC result 2: how many reads per samples
librarySizes(ce)

#annotation 
load("../../Reference/Annotation_CAGEr_hg38_annot.RData")
annotateCTSS(ce, ranges=gff)

#Promoter QC result 3: ploting genomic distribution 
colData(ce)[,c("librarySizes", "promoter", "exon", "intron", "unknown")]
CTSScoordinatesGR(ce)
pdf("genomics_distribution.pdf")
plotAnnot(ce, "counts")
dev.off()

#Promoter QC result 4: ploting correlation between samples 
pdf("samples_correlation_raw_5tagsThres.pdf")
corr.m <- plotCorrelation2( ce, samples = "all"
                          , tagCountThreshold = 5, applyThresholdBoth = TRUE
                          , method = "pearson")
dev.off()

#optional QC: powerlay distribution of cage tags per promoters 
#Normalization
librarySizes(ce)

pdf("powerlaws_before_norm.pdf")
plotReverseCumulatives(ce, fitInRange = c(5, 1000), onePlot = TRUE)
dev.off()

normalizeTagCount(ce, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.2, T = 5*10^4)

pdf("powerlaws_after_norm.pdf")
plotReverseCumulatives(ce, values = "normalized", fitInRange = c(5, 1000), onePlot = TRUE)
dev.off()
