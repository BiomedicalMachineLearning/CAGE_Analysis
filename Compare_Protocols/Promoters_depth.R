
library("MultiAssayExperiment")
library("SummarizedExperiment")
library(CAGEr)
library("BSgenome.Hsapiens.UCSC.hg38")

wdir="/gpfs1/scratch/90days/uqhnguy8/CAGE/ESC/inhouse_CAGE"
name_prefix="ESC"
sub_string="_R1_001.ctss_all_q10_cardio.cageRctss"

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
tag_counts_df <- assay(tag_counts_SEobject)
#CAGEr QC result 1: how many CTSS per samples 
NonZero_CTSS<- apply(as.matrix(tag_counts_df), 2, function(x){sum(x>0)})

lib_size = librarySizes(ce)/1e6

CTSS_per_M = NonZero_CTSS[-1]/lib_size

CTSS_number <- data.frame("CTSS_count" = CTSS_per_M, "Samples" = sampleLabels(ce))


#-------repeat for commercial kit data----------#
path_com <- "/90days/uqhnguy8/CAGE/ESC/commercial_CAGE/"

inputFiles_com =list.files(path_com , pattern=paste0(name_prefix,".*.cageRctss"))

sub_string="_R1_001.ctss_all_q10_CSC.cageRctss"
ce_com <- CAGEexp( genomeName     = "BSgenome.Hsapiens.UCSC.hg38"
             , inputFiles     = paste0(path_com, inputFiles_com)
             , inputFilesType = "ctss"
             , sampleLabels   = gsub(sub_string, "", inputFiles_com)
)

#show sample label 
sampleLabels(ce_com)
#read the files, define CTSS, and count reads for CTSS 
getCTSS(ce_com)

#extract raw tag couts 
tag_counts_SEobject_com <- CTSStagCountSE(ce_com)
tag_counts_df_com <- assay(tag_counts_SEobject_com)
#CAGEr QC result 1: how many CTSS per samples 
NonZero_CTSS_com <- apply(as.matrix(tag_counts_df_com), 2, function(x){sum(x>0)})

lib_size_com = librarySizes(ce_com)/1e6

CTSS_per_M_com = NonZero_CTSS_com[-1]/lib_size_com

CTSS_number_com <- data.frame("CTSS_count_com" = CTSS_per_M_com, "Samples_com" = sampleLabels(ce_com))

#------------------------------------------------
#plotting 
library(ggplot2)
CTSS_number_both <- cbind(CTSS_number_com, CTSS_number)
CTSS_number_both_reshaped  <- melt(CTSS_number_both, ID=Samples)


pdf("TotalGenesPerSample.pdf")
ggplot(CTSS_number_both_reshaped, aes(Samples, value, fill=factor(variable))) + geom_bar(stat="identity", position = "dodge") + theme(axis.text.x = element_text(angle=90, size=16))
dev.off()

#-------------#----------------#------------------

