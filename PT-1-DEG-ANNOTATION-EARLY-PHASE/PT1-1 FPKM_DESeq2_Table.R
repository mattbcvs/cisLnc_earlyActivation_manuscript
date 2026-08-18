#Build sample table with FPKM/DESeq2 from RSEM output

library(tximport)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

samplenames <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/rsem_4timepoints(NovelTx)_NonMergedEnsID", pattern = "*genes.results", full.names = TRUE)[1:16]

#sample info copied/pasted from Amira (matches PCA characteristics later)
#she orderered by patient, then condition
actualnames <- c("1) Pt151 - 0h 2) Pt151 - 4h 3) Pt151 - 8h 4) Pt151 - 24h 
                 5) Pt157 - 0h 6) Pt157 - 4h 7) Pt157 - 8h 8) Pt157 - 24h 
                 9) Pt134 - 0h 10) Pt134 - 4h 11) Pt134 - 8h 12) Pt134 - 24h 
                 13) Pt2279 - 0h14) Pt2279 - 4h 15) Pt2279 - 8h 16) Pt2279 - 24h")
trial <- unlist(strsplit(actualnames, "[0-9]) "))[2:17]
actualnames <- strsplit(trial, "h")
actualnames <- sapply(actualnames, "[[", 1)

sampleInfo_Timecourse <- data.frame("SampleNo" = c(1,10:16,2:9),
                                    "FileName" = samplenames,
                                    "ActualNames" = actualnames[c(1,10:16,2:9)])

#sample names re-ordered into the same
sampleInfo_Timecourse <- sampleInfo_Timecourse[c(1,12,16,5,
                                   9,13,2,6,
                                   10,14,3,7,
                                   11,15,4,8),]

#generate counts table for downstream use
genes<-read.table(sampleInfo_Timecourse$FileName[1],
                  header=TRUE,sep="\t",stringsAsFactors = FALSE)[,1]

fpkm <-do.call(cbind, lapply(sampleInfo_Timecourse$FileName, function(fn){
  read.table(fn, header=TRUE, sep="\t", stringsAsFactors = FALSE)[,7]
  }
  ))

#fpkm <- TPM
fpkm <- data.frame(genes,fpkm,stringsAsFactors = FALSE)
colnames(fpkm)<-c("ENSEMBL",sampleInfo_Timecourse$ActualNames)

fpkm_list_ctrl <- as.list(as.data.frame(t(fpkm[,2:5])))
fpkm_list_pd <- as.list(as.data.frame(t(fpkm[,6:9])))
fpkm_list_il <- as.list(as.data.frame(t(fpkm[,10:13])))
fpkm_list_bo <- as.list(as.data.frame(t(fpkm[,14:17])))

fpkm_mean_ctrl <- sapply(fpkm_list_ctrl, mean)
fpkm_mean_pd <- sapply(fpkm_list_pd, mean)
fpkm_mean_il <- sapply(fpkm_list_il, mean)
fpkm_mean_bo <- sapply(fpkm_list_bo, mean)

std <- function(x) sd(x)/sqrt(length(x))

fpkm_se_basal <- sapply(fpkm_list_ctrl, std)
fpkm_se_chol <- sapply(fpkm_list_pd, std)
fpkm_se_mig <- sapply(fpkm_list_il, std)
fpkm_se_pdgf <- sapply(fpkm_list_bo, std)

fpkm_mean_treatment <- data.frame("Hour0_meanFPKM" = fpkm_mean_ctrl, "Hour0_seFPKM" = fpkm_se_basal,
                                  "Hour4_meanFPKM" = fpkm_mean_pd, "Hour4_seFPKM" = fpkm_se_chol,
                                  "Hour8_meanFPKM" = fpkm_mean_il, "Hour8_seFPKM" = fpkm_se_mig,
                                  "Hour24_meanFPKM" = fpkm_mean_bo, "Hour24_seFPKM" = fpkm_se_pdgf)

trial <- as.list(as.data.frame(t(fpkm_mean_treatment)))
fpkm_max_treatment <- as.numeric(sapply(trial, max))

#trial <- as.list(as.data.frame(t(fpkm_mean_treatment[,1:4])))
#fpkm_max_treatment <- as.numeric(sapply(trial, max))

fpkm <- cbind(fpkm, fpkm_mean_treatment, fpkm_max_treatment)
#fpkm4 <- cbind(fpkm, fpkm_max_treatment4)

#tximport, better use of isoform level info during DESeq2 norm see: https://support.bioconductor.org/p/94003/
samplenames <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/rsem_4timepoints(NovelTx)_NonMergedEnsID/", 
                          pattern = "*isoforms.results", full.names = TRUE)[1:16]

samplenames_order <- samplenames[c(1,12,16,5,
                                   9,13,2,6,
                                   10,14,3,7,
                                   11,15,4,8)]
iso_filenames <- samplenames_order

tx2gene <- read.csv(iso_filenames[1], sep = "\t", stringsAsFactors = F) #85369 transcripts

#build meta-data table
sampleInfo_Timecourse$condition <- c(rep("0hr",4),rep("4hr",4),rep("8hr",4),rep("24hr",4))
sampleInfo_Timecourse$patient <- rep(c("pt151","pt157","pt134","pt2279"),4)

sampleInfo_Timecourse$condition <- factor(sampleInfo_Timecourse$condition, levels(as.factor(sampleInfo_Timecourse$condition))[c(1,3,4,2)])
row.names(sampleInfo_Timecourse) <- sampleInfo_Timecourse$ActualNames

#import in format for DESeq2 input
txi <- tximport(iso_filenames, type = "rsem", tx2gene = tx2gene)

#LRT model
dds <- DESeqDataSetFromTximport(txi, sampleInfo_Timecourse, design=~patient+condition)
dds <-dds[rowSums(counts(dds))>10, ]
ddsLRT <- DESeq(dds, test = "LRT", reduced = ~patient)
resLRT <- results(ddsLRT)

#pca check
rld <- rlog(dds, blind=FALSE)

#ds2 normalised counts:
trial <- assay(rld)
#write.csv(trial, "ds2_counts_timecourse.csv")

pcaData <- plotPCA(rld, intgroup = c("patient", "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

colnames(pcaData)[c(8,9)] <- c("Timepoint", "Patient")

ggplot(pcaData, aes(x = PC1, y = PC2, color = Timepoint, shape = Patient
                    )) +
  geom_point(size =3) +
  xlab(paste0("\nPC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance\n")) +
  coord_fixed() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size =16.5),
        axis.title.x = element_text(size = 16.5),
        legend.text = element_text(size =15),
        legend.title = element_text(size = 16.5))
#clear batch effects (patient) but consistent change in terms of treatment (colour)

assay(rld) <- limma::removeBatchEffect(assay(rld), rld$patient)

pcaData <- plotPCA(rld, intgroup = c("patient", "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

colnames(pcaData)[c(8,9)] <- c("Timepoint", "Patient")

pcaData$Patient <- as.character(pcaData$Patient)
pcaData$Patient <- gsub("pt134", "Pt1", pcaData$Patient)
pcaData$Patient <- gsub("pt151", "Pt2", pcaData$Patient)
pcaData$Patient <- gsub("pt157", "Pt3", pcaData$Patient)
pcaData$Patient <- gsub("pt2279", "Pt4", pcaData$Patient)

ggplot(pcaData, aes(x = PC1, y = PC2, color = Timepoint, shape = Patient
)) +
  geom_point(size =4.5, alpha =0.7) +
  xlab(paste0("\nPC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance\n")) +
  #coord_fixed() +
  theme_minimal() +
  theme(text = element_text(size=24))


sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$patient, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
par(mar=c(10,10,10,10))
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#clearer with the PCA, doesn't cluster by treatment tho - hence need for batch term in design

LRTres <- as.data.frame(resLRT)#FC is 0-24, p is for any change between any timepoints

#pairwise comparisons
dds <- DESeqDataSetFromTximport(txi, sampleInfo_Timecourse, design=~patient+condition)
dds<-dds[rowSums(counts(dds))>10, ]
dds <- DESeq(dds)
#can access pairwise changes as follows, p values are adjusted within each pair however
#could argue this is insufficient if wanting to find point at which gene changes exp.
#p values readjustment trial:
res0_4 <- results(dds, contrast = c("condition", "4hr", "0hr"))
res4_8 <- results(dds, contrast = c("condition", "8hr", "4hr"))
res8_24 <- results(dds, contrast = c("condition", "24hr", "8hr"))
res0_8 <- results(dds, contrast = c("condition", "8hr", "0hr"))
res0_24 <- results(dds, contrast = c("condition", "24hr", "0hr"))
res4_24 <- results(dds, contrast = c("condition", "24hr", "4hr"))

resAllPairs <- as.data.frame(rbind(res0_4, res4_8, res8_24, res0_8, res0_24, res4_24))
resAllPairs$preadj <- p.adjust(resAllPairs$pvalue, method = "BH")
#recapitulate NA values from padj:
resAllPairs$preadj[is.na(resAllPairs$padj)] <- NA

#assemble table
resLRT <- as.data.frame(resLRT)#FC is 0-24, p is for any change between any timepoints

res0_4adjust <- as.data.frame(res0_4)
resLRT$LogFC_0_4 <- res0_4adjust$log2FoldChange
dim(res0_4adjust)#16886 expected
resLRT$preadj_0_4 <- resAllPairs$preadj[1:18001]

res4_8adjust <- as.data.frame(res4_8)
resLRT$LogFC_4_8 <- res4_8adjust$log2FoldChange
dim(res4_8adjust)#same again
dim(res4_8adjust)[1]*2
resLRT$preadj_4_8 <- resAllPairs$preadj[18002:36002]

res8_24adjust <- as.data.frame(res8_24)
resLRT$LogFC_8_24 <- res8_24adjust$log2FoldChange
resLRT$preadj_8_24 <- resAllPairs$preadj[36003:54003]

res0_8adjust <- as.data.frame(res0_8)
resLRT$LogFC_0_8 <- res0_8adjust$log2FoldChange
resLRT$preadj_0_8 <- resAllPairs$preadj[54004:72004]

res0_24adjust <- as.data.frame(res0_24)
resLRT$LogFC_0_24 <- res0_24adjust$log2FoldChange
resLRT$preadj_0_24 <- resAllPairs$preadj[72005:90005]

res4_24adjust <- as.data.frame(res4_24)
resLRT$LogFC_4_24 <- res4_24adjust$log2FoldChange
resLRT$preadj_4_24 <- (resAllPairs$preadj[90006:108006])

Timecourse_DE <- data.frame("ENSEMBL" = rownames(LRTres), resLRT[,6:18])

#take to PT1-2:
Timecourse_4Timepoints_DEnonDE <- merge(fpkm, Timecourse_DE, by = "ENSEMBL", all.x = T)
