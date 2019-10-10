enhancers <- readRDS('enhancers.RDS')
# differential CTSSs analysis - day20 vs day20 replate
exp <- assay(enhancers)[, 10:15]
t=apply(exp, 1, FUN=function(x) sum(x==0))
countdata <- exp[t==0,]
coldata <- matrix(,ncol(countdata),2)
rownames(coldata) <- colnames(countdata)
coldata[,1] <- c('Ori','Ori','Ori','Replate','Replate','Replate')
coldata[,2] <- c('Day20','Day20','Day20','Day20','Day20','Day20')
colnames(coldata) <- c('condition','day')
coldata <- as.data.frame(coldata)
coldata$condition <- relevel(coldata$condition, ref="Ori")
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~condition)
dds <- DESeq(dds)
res <- results(dds, name="Day20_vs_Day20Replate", contrast=c("condition","Ori","Replate"))
resOrdered <- res[order(res$padj),]
res_nona <- resOrdered[!is.na(resOrdered$padj),]
write.table(res_nona[res_nona$padj < 0.05,], 'active_enhancers_Day20_vs_Day20Replate.csv', quote=FALSE)

# differential CTSSs analysis - day20 vs day20 KD replate
exp <- assay(enhancers)[, c(10:12,16:18)]
t=apply(exp, 1, FUN=function(x) sum(x==0))
countdata <- exp[t==0,]
coldata <- matrix(,ncol(countdata),2)
rownames(coldata) <- colnames(countdata)
coldata[,1] <- c('Ori','Ori','Ori','KDReplate','KDReplate','KDReplate')
coldata[,2] <- c('Day20','Day20','Day20','Day20','Day20','Day20')
colnames(coldata) <- c('condition','day')
coldata <- as.data.frame(coldata)
coldata$condition <- relevel(coldata$condition, ref="Ori")
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~condition)
dds <- DESeq(dds)
res <- results(dds, name="Day20_vs_Day20KDReplate", contrast=c("condition","Ori","KDReplate"))
resOrdered <- res[order(res$padj),]
res_nona <- resOrdered[!is.na(resOrdered$padj),]
write.table(res_nona[res_nona$padj < 0.05,], 'active_enhancers_Day20_vs_Day20KDReplate.csv', quote=FALSE)

# differential CTSSs analysis - day20 replate vs day20 KD replate
exp <- assay(enhancers)[, 13:18]
t=apply(exp, 1, FUN=function(x) sum(x==0))
countdata <- exp[t==0,]
coldata <- matrix(,ncol(countdata),2)
rownames(coldata) <- colnames(countdata)
coldata[,1] <- c('Replate','Replate','Replate','KDReplate','KDReplate','KDReplate')
coldata[,2] <- c('Day20','Day20','Day20','Day20','Day20','Day20')
colnames(coldata) <- c('condition','day')
coldata <- as.data.frame(coldata)
coldata$condition <- relevel(coldata$condition, ref="Replate")
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~condition)
dds <- DESeq(dds)
res <- results(dds, name="Day20Replate_vs_Day20KDReplate", contrast=c("condition","Replate","KD_Replate"))
resOrdered <- res[order(res$padj),]
res_nona <- resOrdered[!is.na(resOrdered$padj),]
write.table(res_nona[res_nona$padj < 0.05,], 'active_enhancers_Day20Replate_vs_Day20KDReplate.csv', quote=FALSE)

# differential CTSSs analysis - day 0,2,5,20
exp <- assay(enhancers)[, 1:12]
t=apply(exp, 1, FUN=function(x) sum(x==0))
countdata <- exp[t==0,]
coldata <- matrix(,ncol(countdata),2)
rownames(coldata) <- colnames(countdata)
coldata[,1] <- c('Ori','Ori','Ori','Ori','Ori','Ori','Ori','Ori','Ori','Ori','Ori','Ori')
coldata[,2] <- c('Day0','Day0','Day0','Day2','Day2','Day2','Day5','Day5','Day5','Day20','Day20','Day20')
colnames(coldata) <- c('condition','day')
coldata <- as.data.frame(coldata)
coldata$day <- relevel(coldata$day, ref="Day0")
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~day)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
res_nona <- resOrdered[!is.na(resOrdered$padj),]
write.table(res_nona[res_nona$padj < 0.05,], 'active_enhancers_Day0_vs_Day2_5_20.csv', quote=FALSE)

# differential CTSSs analysis - day 2 vs day 5
exp <- assay(enhancers)[, 4:9]
t=apply(exp, 1, FUN=function(x) sum(x==0))
countdata <- exp[t==0,]
coldata <- matrix(,ncol(countdata),2)
rownames(coldata) <- colnames(countdata)
coldata[,1] <- c('Ori','Ori','Ori','Ori','Ori','Ori')
coldata[,2] <- c('Day2','Day2','Day2','Day5','Day5','Day5')
colnames(coldata) <- c('condition','day')
coldata <- as.data.frame(coldata)
coldata$day <- relevel(coldata$day, ref="Day2")
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~day)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
res_nona <- resOrdered[!is.na(resOrdered$padj),]
write.table(res_nona[res_nona$padj < 0.05,], 'active_enhancers_Day2_vs_Day5.csv', quote=FALSE)

# differential CTSSs analysis - day 5 vs day 20
exp <- assay(enhancers)[, 7:12]
t=apply(exp, 1, FUN=function(x) sum(x==0))
countdata <- exp[t==0,]
coldata <- matrix(,ncol(countdata),2)
rownames(coldata) <- colnames(countdata)
coldata[,1] <- c('Ori','Ori','Ori','Ori','Ori','Ori')
coldata[,2] <- c('Day5','Day5','Day5','Day20','Day20','Day20')
colnames(coldata) <- c('condition','day')
coldata <- as.data.frame(coldata)
coldata$day <- relevel(coldata$day, ref="Day5")
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~day)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
res_nona <- resOrdered[!is.na(resOrdered$padj),]
write.table(res_nona[res_nona$padj < 0.05,], 'active_enhancers_Day5_vs_Day20.csv', quote=FALSE)

# differential CTSSs analysis - day 0 vs day 2_high_chiron
exp <- assay(enhancers)[, c(1:3,19:21)]
t=apply(exp, 1, FUN=function(x) sum(x==0))
countdata <- exp[t==0,]
coldata <- matrix(,ncol(countdata),2)
rownames(coldata) <- colnames(countdata)
coldata[,1] <- c('Ori','Ori','Ori','ori','Ori','Ori')
coldata[,2] <- c('Day0','Day0','Day0','Day2c','Day2c','Day2c')
colnames(coldata) <- c('condition','day')
coldata <- as.data.frame(coldata)
coldata$day <- relevel(coldata$day, ref="Day0")
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~day)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
res_nona <- resOrdered[!is.na(resOrdered$padj),]
write.table(res_nona[res_nona$padj < 0.05,], 'active_enhancers_Day0_vs_Day2chiron.csv', quote=FALSE)

# differential CTSSs analysis - day 2 vs day 2_high_chiron
exp <- assay(enhancers)[, c(4:6,19:21)]
t=apply(exp, 1, FUN=function(x) sum(x==0))
countdata <- exp[t==0,]
coldata <- matrix(,ncol(countdata),2)
rownames(coldata) <- colnames(countdata)
coldata[,1] <- c('Ori','Ori','Ori','Ori','Ori','Ori')
coldata[,2] <- c('Day2','Day2','Day2','Day2c','Day2c','Day2c')
colnames(coldata) <- c('condition','day')
coldata <- as.data.frame(coldata)
coldata$day <- relevel(coldata$day, ref="Day2")
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~day)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
res_nona <- resOrdered[!is.na(resOrdered$padj),]
write.table(res_nona[res_nona$padj < 0.05,], 'active_enhancers_Day2_vs_Day2chiron.csv', quote=FALSE)

active_enhancer_closest_genes <- function(sample1, sample2){
  # please make sure sample1 and sample2 are the same as those being used in the input files
  # get results of closest run by "bedtools closest -k2 -a 1.bed -b 2.bed"
  closest <- read.csv('enhancers_closest_genes.txt', FALSE, '\t')
  pos <- apply(closest[,1:3],1,paste,collapse=':')
  pos <- sub("^([^_]*:[^_]*):", "\\1-", pos)
  pos <- (gsub(" ","",pos))
  closest2 <- cbind(pos, closest)
  closest_significant = matrix(ncol = ncol(closest2))
  colnames(closest_significant) = colnames(closest2)
  # read in the active enhancers
  res_nona <- read.csv(paste0('active_enhancers_',sample1,'_vs_',sample2,'.csv'), sep=' ')
  for (i in 1:nrow(res_nona)) {
    item = rownames(res_nona)[i]
    closest_significant = rbind(closest_significant, closest2[closest2$pos == item,])
  }
  output=closest_significant[-1,c(1,5,6,7,8,10,22)]
  colnames(output) = c('enhancer','gene chr','gene start','gene end','gene id','strand','gene name')
  write.table(output, paste0('active_enhancers_closest_genes_',sample1,'_vs_',sample2,'.csv'),
              sep=',', quote=FALSE, row.names = FALSE)
  return(output)
}

closest_genes20_20RP <- active_enhancer_closest_genes('Day20','Day20Replate')
closest_genes20_20KDRP <- active_enhancer_closest_genes('Day20','Day20KDReplate')
closest_genes20_20KDRP <- active_enhancer_closest_genes('Day20Replate','Day20KDReplate')
closest_genes20_20KDRP <- active_enhancer_closest_genes('Day0','Day2_5_20')
closest_genes20_20KDRP <- active_enhancer_closest_genes('Day2','Day5')
closest_genes20_20KDRP <- active_enhancer_closest_genes('Day5','Day20')
closest_genes20_20KDRP <- active_enhancer_closest_genes('Day0','Day2Chiron')
closest_genes20_20KDRP <- active_enhancer_closest_genes('Day2','Day2Chiron')

