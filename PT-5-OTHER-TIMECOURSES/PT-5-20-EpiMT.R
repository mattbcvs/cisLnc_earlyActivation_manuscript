#Checking for a general pattern in timecourse data
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)

#aims
#1 which timecourses show similar features in terms of
#a) no. DEGs over time
#b) gene sigs over time
#c) types of TFs controlling over time
#expecting the infection/immune related to do this but maybe others

#2 which timecourses have an early bias of lncRNAs, are candidate cis-acting amongst these?

#### set up ####
#option 1:
#obtain gene level counts and do a DEG analysis followed by genomic proximity
#FANTOM CAT repository (access via paper) has SOME timecourses including from the ER set:
#7x cell types, 9x treatments on gene level, all the human stuff basically
#Includes the SMC and some inflammatory, some mesenchymal stuff
#PhD work has gene level counts, and the aoSMC + FGF looks to be in there, check for all human ER timecourse samples:

#plan has some updates:
#a) now including infection timecourses - not ER but sim timeframe to SMC
#b) taking counts, aiming to re-run DE analysis
#c) taking only robust genes now

#Supp Table 18 has the library IDs across the timecourses used in that paper:
Supp18 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/nature21374-s2/SuppTab18.csv")

#Human Early response data, including some additional timecourses in similar timeframe to SVSMC but classed as "activation":
Supp18_ER <- filter(Supp18, grepl("MCF", series_name) | grepl("smooth", series_name) | grepl("LPS", series_name) |
                      grepl("ARPE", series_name) | grepl("VEGFC", series_name) | grepl("Saos", series_name) | 
                      grepl("mesenchymal", series_name) | grepl("inderpest", series_name) | grepl("influenza", series_name))

#number of comparisons per 13x potentially useful ERs
table(Supp18_ER$series_tag)

#all libraries across ERs:
Supp18_ER_libs <- unique(c(unlist(strsplit(Supp18_ER$ref_sample_ID, ",")), unlist(strsplit(Supp18_ER$qry_sample_ID, ","))))

Libraries_bigTable <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Libraries", header = F)
Libraries_bigTable <- as.character(Libraries_bigTable)

#are all libraries in this big table? yes (414):
sum(Supp18_ER_libs %in% Libraries_bigTable)

#subset table in Eddie (for firepower)
#write.table(Supp18_ER_libs, "Supp18_ER_libs.txt")

#subsetted fantom atlas, now diverging from previous by getting counts aiming to re-run DE analysis:
genes <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/genes_robust", 
                  header = T)
ER_counts <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/FANTOM_CAT.expression_ERtime.genes.lv3_robust.count.csv", 
                       header = T)
ER_counts <- cbind(genes, ER_counts)

#try and match comparison number and approach to SVSMC
#need a 4/8/24hr timepoints or as close as possible (24 and 2 between 0 and around 8-12 would probs do...)
#these three timecourses fail: MCF7, aoSMC, lymphangio

#some meta-data:
#Actual timepoints, rep info from Hon Table S1:
All_FANTOM_samples <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/All_FANTOM_sample.csv")
Timecourse_FANTOM_samples <- filter(All_FANTOM_samples, grepl(("day|[0-9]hr|[0-9]min"), sample_name))
Timecourse_FANTOM_samples <- filter(Timecourse_FANTOM_samples, grepl(("rep|donor"), sample_name))
Timecourse_FANTOM_samples$time <- as.character(sapply(strsplit(Timecourse_FANTOM_samples$sample_name, ","), function(x){
  x[grepl("day|hr|min", x)]
}))
Timecourse_FANTOM_samples$rep <- as.character(sapply(strsplit(Timecourse_FANTOM_samples$sample_name, ","), function(x){
  x[grepl("rep|donor", x)]
}))
Timecourse_FANTOM_samples <- Timecourse_FANTOM_samples[,c(1,3,18,19)]

Timecourse_FANTOM_samples$rep <- sapply(strsplit(Timecourse_FANTOM_samples$rep, " "), "[[", 2)

#add in biotype, use their definition for lncRNAs? match to 3PLAR later?
GeneBiotypes <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/nature21374-s2/table11.csv")
#simplify lncRNA class:
GeneBiotypes$CAT_geneClassII <- GeneBiotypes$CAT_geneClass
GeneBiotypes$CAT_geneClassII[
  grepl("lncRNA", GeneBiotypes$CAT_geneClassII)] <- "lncRNA"


#### QC on lib sizes ####

#obtain table of 0,4,8,24:
Supp18_epiMT <- filter(Supp18_ER, grepl("EMT", series_name))

#include some finer increments, 2hrs, can subset later if needed
epiMT_counts <- ER_counts[,c("geneID", unique(unlist(strsplit(Supp18_epiMT$ref_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "02hr00min")$qry_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "04hr00min")$qry_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "06hr00min")$qry_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "08hr00min")$qry_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "12hr00min")$qry_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "16hr00min")$qry_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "24hr00min")$qry_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "42hr00min")$qry_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "60hr00min")$qry_sample_ID, ",")))
)]

dim(epiMT_counts)
#note n=2 in 12h

epiMT_meta <- filter(Timecourse_FANTOM_samples, CAGE_lib_ID %in% colnames(epiMT_counts))
colnames(epiMT_counts)[-1] == epiMT_meta$CAGE_lib_ID
#create sample name
epiMT_meta$time <- gsub(" ","", epiMT_meta$time)
epiMT_meta$sample_nameII <- paste(epiMT_meta$time, epiMT_meta$rep, sep = "-")
epiMT_meta$sample_nameII <- gsub("00min-biol", "", epiMT_meta$sample_nameII)
rownames(epiMT_meta) <- epiMT_meta$sample_nameII
colnames(epiMT_counts)[2:(1+length(rownames(epiMT_meta)))] <- epiMT_meta$sample_nameII

epiMT_counts <- as.matrix(epiMT_counts[,2:(1+length(rownames(epiMT_meta)))])
rownames(epiMT_counts) <- genes$geneID

#essential basic QC step not to be overlooked, library size:
epiMT_libSizes <- data.frame("sample" = colnames(epiMT_counts), "libSize" = colSums(epiMT_counts))
ggplot(epiMT_libSizes) + aes(x = sample, y = libSize) +
  geom_bar(stat = "identity") + Seurat::RotatedAxis()
#only 1 24hr timepoint really...

summary(epiMT_libSizes$libSize/1000000)
#median 1.9mill counts
#big issue in 24hr and 6hr rep2


#
#### for epiMT, find DEGs ####

#obtain table of 0,4,8,24:
Supp18_epiMT <- filter(Supp18_ER, grepl("EMT", series_name))

#ideally should obtain table of 0,4,8,24 BUT

#depth issue in 24hrs, and 1 of 16hr
#there are only 2 12hr samples
epiMT_counts <- ER_counts[,c("geneID", unique(unlist(strsplit(Supp18_epiMT$ref_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "02hr00min")$qry_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "04hr00min")$qry_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "06hr00min")$qry_sample_ID, ",")))[-2],
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "08hr00min")$qry_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "12hr00min")$qry_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "16hr00min")$qry_sample_ID, ","))),
                             #unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "24hr00min")$qry_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "42hr00min")$qry_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "60hr00min")$qry_sample_ID, ",")))
)]

epiMT_meta <- filter(Timecourse_FANTOM_samples, CAGE_lib_ID %in% colnames(epiMT_counts))
colnames(epiMT_counts)[-1] == epiMT_meta$CAGE_lib_ID
#create sample name
epiMT_meta$time <- gsub(" ","", epiMT_meta$time)
epiMT_meta$sample_nameII <- paste(epiMT_meta$time, epiMT_meta$rep, sep = "-")
epiMT_meta$sample_nameII <- gsub("00min-biol", "", epiMT_meta$sample_nameII)
rownames(epiMT_meta) <- epiMT_meta$sample_nameII
colnames(epiMT_counts)[2:(1+length(rownames(epiMT_meta)))] <- epiMT_meta$sample_nameII

epiMT_counts <- as.matrix(epiMT_counts[,2:(1+length(rownames(epiMT_meta)))])
rownames(epiMT_counts) <- genes$geneID

library(DESeq2)
dds <- DESeqDataSetFromMatrix(epiMT_counts[,], 
                              epiMT_meta[,], design=~rep#+batch
                              +time)
dds <- dds[rowSums(counts(dds))>10, ]
dds <- DESeq(dds)

rld <- rlog(dds, blind=FALSE)

pcaData <- plotPCA(rld, intgroup = c("rep", "time"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

colnames(pcaData)
colnames(pcaData)[c(7,8)] <- c("Timepoint", "Replicate")

ggplot(pcaData, aes(x = PC1, y = PC2, color = Timepoint, shape = Replicate)) +
  geom_point(size =3) +
  xlab(paste0("\nPC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance\n")) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size =16.5),
        axis.title.x = element_text(size = 16.5),
        legend.text = element_text(size =15),
        legend.title = element_text(size = 16.5))
#perfect

#with this many timepoints, stick to a 0hr comparison
DEGs_2_0 <- results(dds, contrast=c("time","02hr00min","00hr00min"))
DEGs_4_0 <- results(dds, contrast=c("time","04hr00min","00hr00min"))
DEGs_6_0 <- results(dds, contrast=c("time","06hr00min","00hr00min"))
DEGs_8_0 <- results(dds, contrast=c("time","08hr00min","00hr00min"))
DEGs_12_0 <- results(dds, contrast=c("time","12hr00min","00hr00min"))
DEGs_16_0 <- results(dds, contrast=c("time","16hr00min","00hr00min"))
DEGs_42_0 <- results(dds, contrast=c("time","42hr00min","00hr00min"))
DEGs_60_0 <- results(dds, contrast=c("time","60hr00min","00hr00min"))

#identify DEGs WT v mut in either treatment:
epiMT_DEGs <- data.frame(geneID = rownames(DEGs_4_0), 
                         LFC_0_2 = DEGs_2_0$log2FoldChange, pa_0_2 = DEGs_2_0$padj,
                         LFC_0_4 = DEGs_4_0$log2FoldChange, pa_0_4 = DEGs_4_0$padj,
                         LFC_0_6 = DEGs_6_0$log2FoldChange, pa_0_6 = DEGs_6_0$padj,
                         LFC_0_8 = DEGs_8_0$log2FoldChange, pa_0_8 = DEGs_8_0$padj,
                         LFC_0_12 = DEGs_12_0$log2FoldChange, pa_0_12 = DEGs_12_0$padj,
                         LFC_0_16 = DEGs_16_0$log2FoldChange, pa_0_16 = DEGs_16_0$padj,
                         LFC_0_42 = DEGs_42_0$log2FoldChange, pa_0_42 = DEGs_42_0$padj,
                         LFC_0_60 = DEGs_60_0$log2FoldChange, pa_0_60 = DEGs_60_0$padj
                         )

#get a CPM table alongside:
trial <- lapply(as.data.frame(epiMT_counts), function(x){
  x/sum(x)*1000000
})
epiMT_CPM <- bind_cols(trial)
epiMT_CPM$Hour0_meanCPM <- rowMeans(epiMT_CPM[,1:3])
epiMT_CPM$Hour2_meanCPM <- rowMeans(epiMT_CPM[,4:6])
epiMT_CPM$Hour4_meanCPM <- rowMeans(epiMT_CPM[,7:9])
epiMT_CPM$Hour6_meanCPM <- rowMeans(epiMT_CPM[,10:11])
epiMT_CPM$Hour8_meanCPM <- rowMeans(epiMT_CPM[,12:14])
epiMT_CPM$Hour12_meanCPM <- rowMeans(epiMT_CPM[,15:16])
epiMT_CPM$Hour16_meanCPM <- rowMeans(epiMT_CPM[,17:19])
epiMT_CPM$Hour42_meanCPM <- rowMeans(epiMT_CPM[,20:22])
epiMT_CPM$Hour60_meanCPM <- rowMeans(epiMT_CPM[,23:25])

epiMT_CPM$geneID <- genes$geneID
colnames(epiMT_CPM)
epiMT_DEGs_CPM <- merge(epiMT_DEGs, epiMT_CPM[,], by = "geneID")

#count how many samples reach >1 (too far, no diff to DEGs)
#epiMT_CPM$Hour0_CPM1_no <-  apply(epiMT_CPM[,1:3], 1, function(x){sum(x >1)})
#epiMT_CPM$Hour4_CPM1_no <-  apply(epiMT_CPM[,4:6], 1, function(x){sum(x >1)})
#epiMT_CPM$Hour8_CPM1_no <-  apply(epiMT_CPM[,7:9], 1, function(x){sum(x >1)})
#epiMT_CPM$Hour24_CPM1_no <-  apply(epiMT_CPM[,1:3], 1, function(x){sum(x >1)})

epiMT_DEGs_DE <- filter(epiMT_DEGs_CPM, 
                        #to find DEGs in WT sham vs. cryo:
                        ((Hour0_meanCPM >1 | Hour2_meanCPM >1) & abs(LFC_0_2) > log2(1.5) & pa_0_2 <0.05) |
                          ((Hour0_meanCPM >1 | Hour4_meanCPM >1) & abs(LFC_0_4) > log2(1.5) & pa_0_4 <0.05) |
                          ((Hour0_meanCPM >1 | Hour6_meanCPM >1) & abs(LFC_0_6) > log2(1.5) & pa_0_6 <0.05) |
                          ((Hour0_meanCPM >1 | Hour8_meanCPM >1) & abs(LFC_0_8) > log2(1.5) & pa_0_8 <0.05) |
                          ((Hour0_meanCPM >1 | Hour12_meanCPM >1) & abs(LFC_0_12) > log2(1.5) & pa_0_12 <0.05) |
                          ((Hour0_meanCPM >1 | Hour16_meanCPM >1) & abs(LFC_0_16) > log2(1.5) & pa_0_16 <0.05) |
                          ((Hour0_meanCPM >1 | Hour42_meanCPM >1) & abs(LFC_0_42) > log2(1.5) & pa_0_42 <0.05) |
                          ((Hour0_meanCPM >1 | Hour60_meanCPM >1) & abs(LFC_0_60) > log2(1.5) & pa_0_60 <0.05))

dim(epiMT_DEGs_DE)[1] #no. DEGs total
sum(Biobase::rowMax(as.matrix(epiMT_CPM[,26:34])) >1)
8738/22268 #39%
#4 timepoint version for reference
4721/20272

#biotypes:
trial <- merge(epiMT_CPM, GeneBiotypes, by.x = "geneID", by.y = "CAT_geneID")
table(filter(trial, Biobase::rowMax(as.matrix(trial[,27:35])) >1)$CAT_geneClassII)

#temporal ordering:
#regulated within 2 hours:
epiMT_DEGs_DE_Upwithin_2 <- filter(epiMT_DEGs_DE, (Hour0_meanCPM >=1 | Hour2_meanCPM >=1) &
                                     (LFC_0_2 >= log2(1.5) & pa_0_2 <0.05))
epiMT_DEGs_DE_Downwithin_2 <- filter(epiMT_DEGs_DE, (Hour0_meanCPM >=1 | Hour2_meanCPM >=1) &
                                       (LFC_0_2 < -log2(1.5) & pa_0_2 <0.05))

#of remainder, regulated within 4 hours:
epiMT_DEGs_DE_Upwithin_4 <- filter(epiMT_DEGs_DE, !geneID %in% c(epiMT_DEGs_DE_Upwithin_2$geneID, epiMT_DEGs_DE_Downwithin_2$geneID),
                                   (Hour0_meanCPM >=1 | Hour4_meanCPM >=1) &
                                   (LFC_0_4 >= log2(1.5) & pa_0_4 <0.05))
epiMT_DEGs_DE_Downwithin_4 <- filter(epiMT_DEGs_DE, !geneID %in% c(epiMT_DEGs_DE_Upwithin_2$geneID, epiMT_DEGs_DE_Downwithin_2$geneID),
                                     (Hour0_meanCPM >=1 | Hour4_meanCPM >=1) &
                                     (LFC_0_4 < -log2(1.5) & pa_0_4 <0.05))

#of remainder, regulated within 6 hours:
epiMT_DEGs_DE_Upwithin_6 <- filter(epiMT_DEGs_DE, !geneID %in% c(epiMT_DEGs_DE_Upwithin_2$geneID, epiMT_DEGs_DE_Downwithin_2$geneID,
                                                                 epiMT_DEGs_DE_Upwithin_4$geneID, epiMT_DEGs_DE_Downwithin_4$geneID),
                                   (Hour0_meanCPM >=1 | Hour6_meanCPM >=1) &
                                     (LFC_0_6 >= log2(1.5) & pa_0_6 <0.05))
epiMT_DEGs_DE_Downwithin_6 <- filter(epiMT_DEGs_DE, !geneID %in% c(epiMT_DEGs_DE_Upwithin_2$geneID, epiMT_DEGs_DE_Downwithin_2$geneID,
                                                                   epiMT_DEGs_DE_Upwithin_4$geneID, epiMT_DEGs_DE_Downwithin_4$geneID),
                                     (Hour0_meanCPM >=1 | Hour6_meanCPM >=1) &
                                       (LFC_0_6 < -log2(1.5) & pa_0_6 <0.05))

#of remainder, regulated within 8 hours:
epiMT_DEGs_DE_Upwithin_8 <- filter(epiMT_DEGs_DE, !geneID %in% c(epiMT_DEGs_DE_Upwithin_2$geneID, epiMT_DEGs_DE_Downwithin_2$geneID,
                                                                 epiMT_DEGs_DE_Upwithin_4$geneID, epiMT_DEGs_DE_Downwithin_4$geneID,
                                                                 epiMT_DEGs_DE_Upwithin_6$geneID, epiMT_DEGs_DE_Downwithin_6$geneID),
                                   (Hour0_meanCPM >=1 | Hour8_meanCPM >=1) &
                                     (LFC_0_8 >= log2(1.5) & pa_0_8 <0.05))
epiMT_DEGs_DE_Downwithin_8 <- filter(epiMT_DEGs_DE, !geneID %in% c(epiMT_DEGs_DE_Upwithin_2$geneID, epiMT_DEGs_DE_Downwithin_2$geneID,
                                                                   epiMT_DEGs_DE_Upwithin_4$geneID, epiMT_DEGs_DE_Downwithin_4$geneID,
                                                                   epiMT_DEGs_DE_Upwithin_6$geneID, epiMT_DEGs_DE_Downwithin_6$geneID),
                                     (Hour0_meanCPM >=1 | Hour8_meanCPM >=1) &
                                       (LFC_0_8 < -log2(1.5) & pa_0_8 <0.05))

#of remainder, regulated within 12 hours:
epiMT_DEGs_DE_Upwithin_12 <- filter(epiMT_DEGs_DE, !geneID %in% c(epiMT_DEGs_DE_Upwithin_2$geneID, epiMT_DEGs_DE_Downwithin_2$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_4$geneID, epiMT_DEGs_DE_Downwithin_4$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_6$geneID, epiMT_DEGs_DE_Downwithin_6$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_8$geneID, epiMT_DEGs_DE_Downwithin_8$geneID),
                                   (Hour0_meanCPM >=1 | Hour12_meanCPM >=1) &
                                     (LFC_0_12 >= log2(1.5) & pa_0_12 <0.05))
epiMT_DEGs_DE_Downwithin_12 <- filter(epiMT_DEGs_DE, !geneID %in% c(epiMT_DEGs_DE_Upwithin_2$geneID, epiMT_DEGs_DE_Downwithin_2$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_4$geneID, epiMT_DEGs_DE_Downwithin_4$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_6$geneID, epiMT_DEGs_DE_Downwithin_6$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_8$geneID, epiMT_DEGs_DE_Downwithin_8$geneID),
                                     (Hour0_meanCPM >=1 | Hour12_meanCPM >=1) &
                                       (LFC_0_12 < -log2(1.5) & pa_0_12 <0.05))

#of remainder, regulated within 16 hours:
epiMT_DEGs_DE_Upwithin_16 <- filter(epiMT_DEGs_DE, !geneID %in% c(epiMT_DEGs_DE_Upwithin_2$geneID, epiMT_DEGs_DE_Downwithin_2$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_4$geneID, epiMT_DEGs_DE_Downwithin_4$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_6$geneID, epiMT_DEGs_DE_Downwithin_6$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_8$geneID, epiMT_DEGs_DE_Downwithin_8$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_12$geneID, epiMT_DEGs_DE_Downwithin_12$geneID),
                                   (Hour0_meanCPM >=1 | Hour16_meanCPM >=1) &
                                     (LFC_0_16 >= log2(1.5) & pa_0_16 <0.05))
epiMT_DEGs_DE_Downwithin_16 <- filter(epiMT_DEGs_DE, !geneID %in% c(epiMT_DEGs_DE_Upwithin_2$geneID, epiMT_DEGs_DE_Downwithin_2$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_4$geneID, epiMT_DEGs_DE_Downwithin_4$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_6$geneID, epiMT_DEGs_DE_Downwithin_6$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_8$geneID, epiMT_DEGs_DE_Downwithin_8$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_12$geneID, epiMT_DEGs_DE_Downwithin_12$geneID),
                                     (Hour0_meanCPM >=1 | Hour16_meanCPM >=1) &
                                       (LFC_0_16 < -log2(1.5) & pa_0_16 <0.05))

#of remainder, regulated within 42 hours:
epiMT_DEGs_DE_Upwithin_42 <- filter(epiMT_DEGs_DE, !geneID %in% c(epiMT_DEGs_DE_Upwithin_2$geneID, epiMT_DEGs_DE_Downwithin_2$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_4$geneID, epiMT_DEGs_DE_Downwithin_4$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_6$geneID, epiMT_DEGs_DE_Downwithin_6$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_8$geneID, epiMT_DEGs_DE_Downwithin_8$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_12$geneID, epiMT_DEGs_DE_Downwithin_12$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_16$geneID, epiMT_DEGs_DE_Downwithin_16$geneID),
                                   (Hour0_meanCPM >=1 | Hour42_meanCPM >=1) &
                                     (LFC_0_42 >= log2(1.5) & pa_0_42 <0.05))
epiMT_DEGs_DE_Downwithin_42 <- filter(epiMT_DEGs_DE, !geneID %in% c(epiMT_DEGs_DE_Upwithin_2$geneID, epiMT_DEGs_DE_Downwithin_2$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_4$geneID, epiMT_DEGs_DE_Downwithin_4$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_6$geneID, epiMT_DEGs_DE_Downwithin_6$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_8$geneID, epiMT_DEGs_DE_Downwithin_8$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_12$geneID, epiMT_DEGs_DE_Downwithin_12$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_16$geneID, epiMT_DEGs_DE_Downwithin_16$geneID),
                                     (Hour0_meanCPM >=1 | Hour42_meanCPM >=1) &
                                       (LFC_0_42 < -log2(1.5) & pa_0_42 <0.05))

#of remainder, regulated within 60 hours:
epiMT_DEGs_DE_Upwithin_60 <- filter(epiMT_DEGs_DE, !geneID %in% c(epiMT_DEGs_DE_Upwithin_2$geneID, epiMT_DEGs_DE_Downwithin_2$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_4$geneID, epiMT_DEGs_DE_Downwithin_4$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_6$geneID, epiMT_DEGs_DE_Downwithin_6$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_8$geneID, epiMT_DEGs_DE_Downwithin_8$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_12$geneID, epiMT_DEGs_DE_Downwithin_12$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_16$geneID, epiMT_DEGs_DE_Downwithin_16$geneID,
                                                                  epiMT_DEGs_DE_Upwithin_42$geneID, epiMT_DEGs_DE_Downwithin_42$geneID),
                                   (Hour0_meanCPM >=1 | Hour60_meanCPM >=1) &
                                     (LFC_0_60 >= log2(1.5) & pa_0_60 <0.05))
epiMT_DEGs_DE_Downwithin_60 <- filter(epiMT_DEGs_DE, !geneID %in% c(epiMT_DEGs_DE_Upwithin_2$geneID, epiMT_DEGs_DE_Downwithin_2$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_4$geneID, epiMT_DEGs_DE_Downwithin_4$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_6$geneID, epiMT_DEGs_DE_Downwithin_6$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_8$geneID, epiMT_DEGs_DE_Downwithin_8$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_12$geneID, epiMT_DEGs_DE_Downwithin_12$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_16$geneID, epiMT_DEGs_DE_Downwithin_16$geneID,
                                                                    epiMT_DEGs_DE_Upwithin_42$geneID, epiMT_DEGs_DE_Downwithin_42$geneID),
                                     (Hour0_meanCPM >=1 | Hour60_meanCPM >=1) &
                                       (LFC_0_60 < -log2(1.5) & pa_0_60 <0.05))


#for figures, assign each gene to wave, should have all the genes seperated by the above:
epiMT_DEGs_DE$RegulationStart[epiMT_DEGs_DE$geneID %in% epiMT_DEGs_DE_Upwithin_2$geneID] <- "Induced <2hrs"
epiMT_DEGs_DE$RegulationStart[epiMT_DEGs_DE$geneID %in% epiMT_DEGs_DE_Downwithin_2$geneID] <- "Repressed <2hrs"

epiMT_DEGs_DE$RegulationStart[epiMT_DEGs_DE$geneID %in% epiMT_DEGs_DE_Upwithin_4$geneID] <- "Induced 2-4hrs"
epiMT_DEGs_DE$RegulationStart[epiMT_DEGs_DE$geneID %in% epiMT_DEGs_DE_Downwithin_4$geneID] <- "Repressed 2-4hrs"

epiMT_DEGs_DE$RegulationStart[epiMT_DEGs_DE$geneID %in% epiMT_DEGs_DE_Upwithin_6$geneID] <- "Induced 4-6hrs"
epiMT_DEGs_DE$RegulationStart[epiMT_DEGs_DE$geneID %in% epiMT_DEGs_DE_Downwithin_6$geneID] <- "Repressed 4-6hrs"

epiMT_DEGs_DE$RegulationStart[epiMT_DEGs_DE$geneID %in% epiMT_DEGs_DE_Upwithin_8$geneID] <- "Induced 6-8hrs"
epiMT_DEGs_DE$RegulationStart[epiMT_DEGs_DE$geneID %in% epiMT_DEGs_DE_Downwithin_8$geneID] <- "Repressed 6-8hrs"

epiMT_DEGs_DE$RegulationStart[epiMT_DEGs_DE$geneID %in% epiMT_DEGs_DE_Upwithin_12$geneID] <- "Induced 8-12hrs"
epiMT_DEGs_DE$RegulationStart[epiMT_DEGs_DE$geneID %in% epiMT_DEGs_DE_Downwithin_12$geneID] <- "Repressed 8-12hrs"

epiMT_DEGs_DE$RegulationStart[epiMT_DEGs_DE$geneID %in% epiMT_DEGs_DE_Upwithin_16$geneID] <- "Induced 12-16hrs"
epiMT_DEGs_DE$RegulationStart[epiMT_DEGs_DE$geneID %in% epiMT_DEGs_DE_Downwithin_16$geneID] <- "Repressed 12-16hrs"

epiMT_DEGs_DE$RegulationStart[epiMT_DEGs_DE$geneID %in% epiMT_DEGs_DE_Upwithin_42$geneID] <- "Induced 16-42hrs"
epiMT_DEGs_DE$RegulationStart[epiMT_DEGs_DE$geneID %in% epiMT_DEGs_DE_Downwithin_42$geneID] <- "Repressed 16-42hrs"

epiMT_DEGs_DE$RegulationStart[epiMT_DEGs_DE$geneID %in% epiMT_DEGs_DE_Upwithin_60$geneID] <- "Induced 42-60hrs"
epiMT_DEGs_DE$RegulationStart[epiMT_DEGs_DE$geneID %in% epiMT_DEGs_DE_Downwithin_60$geneID] <- "Repressed 42-60hrs"

table(epiMT_DEGs_DE$RegulationStart)


#
#### heatmap/GO plotting, can skip for stats ####

mat <- epiMT_DEGs_DE[,17:28] #4721 genes - similarly dynamic model than SVSMC (which is ~4345)
rownames(mat) <- epiMT_DEGs_DE[,1]
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
mat <- t(apply(mat, 1, cal_z_score))

myColor <- colorRampPalette(c("steelblue", "white", "red"))(50)
myBreaks <- c(seq(min(mat), 0, 
                  length.out=ceiling(50/2)), 
              seq(max(mat)/50, 
                  max(mat), 
                  length.out=floor(50/2)))


condition<-c(rep("0hr",3),rep("4hr",3), rep("8hr",3),rep("16hr",3))
patient<- c(rep(c("1","2","3"),4))
data_cols<-data.frame(condition=condition,patient=patient)

data_colsHeat <- data.frame("Hours" = data_cols[,1], stringsAsFactors = T)
rownames(data_colsHeat) <- colnames(mat)
data_colsHeat$Hours <- factor(data_colsHeat$Hours, levels(data_colsHeat$Hours)[c(1,3,4,2)])

mati <- mat[rownames(mat) %in% epiMT_DEGs_DE_Upwithin_4$geneID,]
mati <- mat[rownames(mat) %in% epiMT_DEGs_DE_Downwithin_4$geneID,]
mati <- mat[rownames(mat) %in% epiMT_DEGs_DE_Upwithin_8$geneID,]
mati <- mat[rownames(mat) %in% epiMT_DEGs_DE_Downwithin_8$geneID,]
mati <- mat[rownames(mat) %in% epiMT_DEGs_DE_Upwithin_24$geneID,]
mati <- mat[rownames(mat) %in% epiMT_DEGs_DE_Downwithin_24$geneID,]

library(pheatmap)
p <- pheatmap(mati,clustering_method = "complete",annotation_legend = F,
              annotation_col = data_colsHeat,
              show_colnames = F, 
              show_rownames = F, 
              cluster_cols = F,
              cluster_rows = T,
              #cutree_rows = 3,
              treeheight_col = 0, 
              treeheight_row = 45,
              legend = F,
              color = myColor, 
              breaks = myBreaks,
              border_color = NA)


#build support for clustering - biological theming
epiMT_DEGs_DE$Ens_ID_merge <- sapply(sapply(epiMT_DEGs_DE$geneID, strsplit, "\\."), "[[", 1)

#background:
epiMT_DEGs_CPM_exprs <- filter(epiMT_DEGs_CPM, rowMax(as.matrix(epiMT_DEGs_CPM[,47:55]))>1)
epiMT_DEGs_CPM_exprs$Ens_ID_merge <- sapply(sapply(epiMT_DEGs_CPM_exprs$geneID, strsplit, "\\."), "[[", 1)

trial <- split(epiMT_DEGs_DE, epiMT_DEGs_DE$RegulationStart)

library(clusterProfiler)
library(org.Hs.eg.db)

#only really necessary for manuscript to show induction of immune response in an early phase
trial <- lapply(trial[c(1,4,5,9,12,13)], function(x){
  enrichGO(gene          = unique(x$Ens_ID_merge),
           universe      = unique(epiMT_DEGs_CPM_exprs$Ens_ID_merge),
           keyType       = "ENSEMBL",
           OrgDb         = org.Hs.eg.db,
           ont           = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05,
           readable      = TRUE)
})

sapply(trial, dim)

saveRDS(trial, "GO_list_epiMT_2026_firstphases.rds")

edox2 <- enrichplot::pairwise_termsim(trial[[1]])
enrichplot::treeplot(edox2, showCategory = 30, #fontSize =2, #extend = 0.1, #hilight =F #nWords = 3,
                     cluster.params = list(n = 7), 
                     #label_format_tiplab = function(x) stringr::str_wrap(x, width=40),
                     label_format = function(x) stringr::str_wrap(x, width=25))
View(data.frame(trial[[1]]))
View(data.frame(trial[[2]]))
View(data.frame(trial[[3]]))
View(data.frame(trial[[4]]))
View(data.frame(trial[[5]]))
View(data.frame(trial[[6]]))

dotplot(simplify(trial[[1]]), showCategory = 10)
dotplot(simplify(trial[[2]]), showCategory = 10)


#### temporal biases of PCGs and subclasses ####

#size of waves in general:
trial <- as.data.frame(table(epiMT_DEGs_DE$RegulationStart))
trial$Freq/dim(epiMT_DEGs_DE)[1]*100
#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var1), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var1), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-42hrs", "42-60hrs")
trial$UpDown <- sapply(sapply(as.character(trial$Var1), strsplit, " "),"[[" , 1)
trial$PercAllLncs <- trial$Freq/dim(epiMT_DEGs_DE)[1]*100

trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,5,7,8,2,3,6)])

AllGWaveBias <- trial

ggplot(AllGWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(AllGWaveBias, grepl("Induced", Var1)), 
           aes(y = PercAllLncs, fill = UpDown)) +
  geom_label(data = filter(AllGWaveBias, grepl("Induced", Var1)), 
             aes(y = PercAllLncs+1, label = Freq), size = 5) +
  geom_col(data = filter(AllGWaveBias, grepl("Repressed", Var1)), 
           aes(y = -PercAllLncs, fill = UpDown)) +
  geom_label(data = filter(AllGWaveBias, grepl("Repressed", Var1)), 
             aes(y = -PercAllLncs-1, label = Freq), size = 5) +
  ylab("% all DEGs") +
  xlab("") +
  #scale_y_continuous(limits = c(-25,25), breaks = seq(-20,20, by = 20),
  #                   labels = (c(seq(20, 0, by = -20), seq(20,20,by=20)))) +
  theme_minimal() +
  theme(text = element_text(size=15))


#plot genes per cluster
epiMT_DEGs_DE <- merge(GeneBiotypes[,c(1,2,4,9)], epiMT_DEGs_DE, by.x = "CAT_geneID", by.y = "geneID")

Cluster_biotype <- as.data.frame(table(epiMT_DEGs_DE$CAT_geneClassII, epiMT_DEGs_DE$RegulationStart))

table(epiMT_DEGs_DE$CAT_geneClassII)

#bias of PCGs:
Selected <- "coding_mRNA"
table(epiMT_DEGs_DE$CAT_geneClassII)[Selected]
table(epiMT_DEGs_DE$CAT_geneClassII)/dim(epiMT_DEGs_DE)[1]*100 

trial <- filter(Cluster_biotype, grepl(Selected, Var1))
trial$Freq/table(epiMT_DEGs_DE$RegulationStart)*100
#looks like possible biases, less than expected in 2x later points

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(epiMT_DEGs_DE$RegulationStart)

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

ClusterNames <- trial$Var2

PCGEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(epiMT_DEGs_DE$CAT_geneClassII)[Selected]
  d <- dim(epiMT_DEGs_DE)[1]
  
  PCGEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("PCG" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est,
                                       fisher.test(data.frame("PCG" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p)
}

names(PCGEnrich_cluster) <- ClusterNames
triali <- bind_rows(PCGEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")
#Induced first timepoint, rarely at second
#2x significant biases - opposite of lncRNAs

#percentage plots
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-42hrs", "42-60hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,5,7,8,2,3,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(epiMT_DEGs_DE$CAT_geneClassII)[Selected]*100
trial$PercBackground <- trial$selection/dim(epiMT_DEGs_DE)[1]*100

PCGWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(PCGWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(PCGWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(PCGWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(PCGWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 5) +
  geom_col(data = filter(PCGWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(PCGWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(PCGWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory, label = Freq), size = 5) +
  ylab("% DE PCGs") +
  xlab("") +
  #scale_y_continuous(limits = c(-30,35),breaks = seq(-30,30, by = 10),
  #                   labels = (c(seq(30, 0, by = -10), seq(10,30,by=10)))) +
  theme_minimal()+
  theme(text = element_text(size=15))


#CCs
CC_Freeman <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/FreemanCellCycle.csv", header = T)[,1:35]

epiMT_DEGs_DE$CC <- epiMT_DEGs_DE$GeneClassUpdate
epiMT_DEGs_DE$CC[gsub("\\.[0-9]*", "", epiMT_DEGs_DE$CAT_geneID) %in% CC_Freeman$Ensembl_ID] <- "CC"

Selected <- "CC"
table(epiMT_DEGs_DE$CC)[Selected]
table(epiMT_DEGs_DE$CC)/dim(epiMT_DEGs_DE)[1]*100

Cluster_biotype <- as.data.frame(table(epiMT_DEGs_DE$CC, epiMT_DEGs_DE$RegulationStart))
trial <- filter(Cluster_biotype, grepl(Selected, Var1))
trial$Freq/table(epiMT_DEGs_DE$RegulationStart)*100
#looks like definite biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(epiMT_DEGs_DE$RegulationStart)

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

ClusterNames <- trial$Var2

CCEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(epiMT_DEGs_DE$CC)[Selected]
  d <- dim(epiMT_DEGs_DE)[1]
  
  CCEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                             "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est,
                                      fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                             "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p)
}

names(CCEnrich_cluster) <- ClusterNames
triali <- bind_rows(CCEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")
#all significant biases

#percentage plots
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-42hrs", "42-60hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,5,7,8,2,3,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(epiMT_DEGs_DE$CC)[Selected]*100
trial$PercBackground <- trial$selection/dim(epiMT_DEGs_DE)[1]*100

CCWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(CCWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(CCWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(CCWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CCWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 5) +
  geom_col(data = filter(CCWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(CCWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CCWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory, label = Freq), size = 5) +
  ylab("% DE S/G2M") +
  xlab("") +
  #scale_y_continuous(limits = c(-30,35),breaks = seq(-30,30, by = 10),
  #                   labels = (c(seq(30, 0, by = -10), seq(10,30,by=10)))) +
  theme_minimal() +
  theme(text = element_text(size = 15)) + Seurat::RotatedAxis()


#IEGs
#inserting new gene categories:
IEGs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Gene lists/arner_2015_table_S5_IEGs_lit.csv")
IEGs_hs <- filter(IEGs, !Hs_symbol %in% c("", "-"))
IEGs_hs <- filter(IEGs_hs, !grepl("[a-z]", IEGs_hs$Total_count))

epiMT_DEGs_DE$IEG <- epiMT_DEGs_DE$GeneClassUpdate
epiMT_DEGs_DE$IEG[epiMT_DEGs_DE$CAT_geneName %in% IEGs_hs$Hs_symbol] <- "IEG"

Selected <- "IEG"
table(epiMT_DEGs_DE$IEG)[Selected]
table(epiMT_DEGs_DE$IEG)/dim(epiMT_DEGs_DE)[1]*100

Cluster_biotype <- as.data.frame(table(epiMT_DEGs_DE$IEG, epiMT_DEGs_DE$RegulationStart))
trial <- filter(Cluster_biotype, grepl(Selected, Var1))
trial$Freq/table(epiMT_DEGs_DE$RegulationStart)*100
#looks like definite biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(epiMT_DEGs_DE$RegulationStart)

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

ClusterNames <- trial$Var2

IEGEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(epiMT_DEGs_DE$IEG)[Selected]
  d <- dim(epiMT_DEGs_DE)[1]
  
  IEGEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p)
}

names(IEGEnrich_cluster) <- ClusterNames
triali <- bind_rows(IEGEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")
#all significant biases

#percentage plots
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-42hrs", "42-60hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,5,7,8,2,3,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(epiMT_DEGs_DE$IEG)[Selected]*100
trial$PercBackground <- trial$selection/dim(epiMT_DEGs_DE)[1]*100

IEGWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(IEGWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(IEGWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(IEGWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(IEGWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory+2, label = Freq), size = 4.2) +
  geom_col(data = filter(IEGWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(IEGWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(IEGWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory-2, label = Freq), size = 4.2) +
  ylab("% DE IEGs") +
  xlab("") +
  #scale_y_continuous(limits = c(-30,35),breaks = seq(-30,30, by = 10),
  #                   labels = (c(seq(30, 0, by = -10), seq(10,30,by=10)))) +
  theme_minimal() +
  theme(text = element_text(size = 15)) + Seurat::RotatedAxis()


#### (to review) combined figure ####

#combined figure:
AllTypesWaveBias <- rbind(IEGWaveBias, CCWaveBias)

#correct timeframe:
#AllTypesWaveBias$FirstRegulation <- gsub("24hr", "16hr", AllTypesWaveBias$FirstRegulation)

AllTypesWaveBias$OR_corrected <- AllTypesWaveBias$OR
AllTypesWaveBias$OR_corrected[AllTypesWaveBias$p_adj >0.1] <- NA
AllTypesWaveBias$OR_corrected <- AllTypesWaveBias$OR_corrected + 0.001
AllTypesWaveBias$`Log2(Odds Ratio)` <- log2(AllTypesWaveBias$OR_corrected)

AllTypesWaveBias$`Log2(Odds Ratio)`[AllTypesWaveBias$`Log2(Odds Ratio)` < -5] <- -5
AllTypesWaveBias$`Log2(Odds Ratio)`[AllTypesWaveBias$`Log2(Odds Ratio)` >5] <- 5

AllTypesWaveBias$padj_corrected <- AllTypesWaveBias$p_adj
AllTypesWaveBias$padj_corrected[AllTypesWaveBias$p_adj >0.1] <- NA

AllTypesWaveBias$padj_simple <- NA
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.1] <- "p<0.1"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.05] <- "p<0.05"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.01] <- "p<0.01"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.001] <- "p<0.001"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.0001] <- "p<0.0001"

AllTypesWaveBias$UpDownType <- paste(AllTypesWaveBias$Var1, AllTypesWaveBias$UpDown, sep = "-")

#insert spacers for plotting:
AllTypesWaveBias <- AllTypesWaveBias[,c(15,4,12,14)]

appendSpacers <- data.frame(UpDownType = c(rep("Spacer1",3),
                                           rep("Spacer2",3)#,
                                           #rep("Spacer3",3)#,
                                           #rep("Spacer4",3)#,
                                           #rep("Spacer5",3),
                                           #rep("Spacer6",3),
                                           #rep("Spacer7",3),
                                           #rep("Spacer8",3)
), 
"FR" = rep(c("Within \n4hrs",
             "Within \n8hrs",
             "Within \n16hrs"), 2), 
"OR" = NA,  
`padj_simple` =NA)

colnames(appendSpacers) <- colnames(AllTypesWaveBias)

AllTypesWaveBias <- rbind(AllTypesWaveBias, appendSpacers)

AllTypesWaveBias$UpDownType <- factor(AllTypesWaveBias$UpDownType)
AllTypesWaveBias$UpDownType <- factor(AllTypesWaveBias$UpDownType, 
                                      levels = levels(AllTypesWaveBias$UpDownType)[c(2,1,
                                                                                     4,3)])

AllTypesWaveBias$FirstRegulation <- factor(AllTypesWaveBias$FirstRegulation)
#AllTypesWaveBias$FirstRegulation <- factor(AllTypesWaveBias$FirstRegulation, 
#                                      levels = levels(AllTypesWaveBias$FirstRegulation)[c(2,3,1)])


myColor <- colorRampPalette(c("steelblue", "white", "red"))(50)
myBreaks <- c(seq(min(AllTypesWaveBias$`Log2(Odds Ratio)`, na.rm = T), 0, 
                  length.out=ceiling(50/2)), 
              seq(max(AllTypesWaveBias$`Log2(Odds Ratio)`, na.rm = T)/50, 
                  max(AllTypesWaveBias$`Log2(Odds Ratio)`, na.rm = T), 
                  length.out=floor(50/2)))

ggplot(AllTypesWaveBias) + aes(x = FirstRegulation, y = UpDownType, size = padj_simple, fill = `Log2(Odds Ratio)`) +
  geom_point(color = "grey30", shape = 21) +
  xlab("") +
  #ylab("") +
  scale_size_discrete(range = c(5,12), limits = c(#"p<0.05", 
                                                  "p<0.01", 
                                                  #"p<0.001", 
                                                  "p<0.0001")) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "red") +
  theme_minimal() +
  theme(legend.key.size = unit(1.4, "line"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        axis.text.x = element_text(size=22),
        #axis.text.y = element_blank()
  ) + Seurat::RotatedAxis()



#### biases of lncRNAs ####

table(epiMT_DEGs_DE$CAT_geneClassII)
table(epiMT_DEGs_DE$CAT_geneClassII)/dim(epiMT_DEGs_DE)[1]*100

Cluster_biotype <- as.data.frame(table(epiMT_DEGs_DE$CAT_geneClassII, epiMT_DEGs_DE$RegulationStart))
trial <- filter(Cluster_biotype, Var1 == "lncRNA")
trial$Freq/table(epiMT_DEGs_DE$RegulationStart)*100
#looks like weak biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(epiMT_DEGs_DE$RegulationStart)

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

ClusterNames <- trial$Var2

LncEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(epiMT_DEGs_DE$CAT_geneClassII)["lncRNA"]
  d <- dim(epiMT_DEGs_DE)[1]
  
  LncEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p)
}

names(LncEnrich_cluster) <- ClusterNames
triali <- bind_rows(LncEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")


#percentage plots
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-42hrs", "42-60hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,5,7,8,2,3,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(epiMT_DEGs_DE$CAT_geneClassII)["lncRNA"]*100
trial$PercBackground <- trial$selection/dim(epiMT_DEGs_DE)[1]*100

LncWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(LncWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(LncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(LncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(LncWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 5) +
  geom_col(data = filter(LncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(LncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(LncWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory, label = Freq), size = 5) +
  ylab("% DE LncRNAs") +
  xlab("") +
  #scale_y_continuous(limits = c(-30,35),breaks = seq(-30,30, by = 10),
  #                   labels = (c(seq(30, 0, by = -10), seq(10,30,by=10)))) +
  theme_minimal() +
  theme(text = element_text(size=15))


#
#### Define candidate cis-acting lncRNAs (+PCGs) ####

#include HiC? unlikely to have for all cell types in these timecourses...
#CASMC HiC info from Zhao showed unlikely to add in a huge no. of interactions beyond taking anything 250kbp either side of lnc TSS anyway
#do without for now - can update later if looks good

#need co-ordinates:
#FCAT_robust <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/FANTOM_CAT.lv3_robust.bed", header = F)
#FCAT_robust <- FCAT_robust[,c(1:4,6)]
#FCAT_robust$gene <- sapply(sapply(FCAT_robust$V4, strsplit, "\\|"), "[[", 1)

#use this to subset big table via eddie power
epiMT_libs <- c(unique(unlist(strsplit(Supp18_epiMT$ref_sample_ID, ","))),
                unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "02hr00min")$qry_sample_ID, ","))),
                unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "04hr00min")$qry_sample_ID, ","))),
                unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "06hr00min")$qry_sample_ID, ",")))[-2],
                unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "08hr00min")$qry_sample_ID, ","))),
                unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "12hr00min")$qry_sample_ID, ","))),
                unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "16hr00min")$qry_sample_ID, ","))),
                #unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "24hr00min")$qry_sample_ID, ","))),
                unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "42hr00min")$qry_sample_ID, ","))),
                unique(unlist(strsplit(filter(Supp18_epiMT, qry_replicate_ID == "60hr00min")$qry_sample_ID, ","))))
#write.table(epiMT_libs, quote = F, "epiMT_libs", row.names = F)

#CAGEs for the 24hr epiMT genes, take this list to FANTOM:
epiMT_exprs <- unique(filter(epiMT_CPM, Biobase::rowMax(as.matrix(epiMT_CPM[,26:34])) >1)$geneID)
#write.table(epiMT_exprs, quote = F, "epiMT_exprs", row.names = F)

#now use "Subset_FANTOM5_epiMT.R" in eddie

#import, counts for selected libraries + genes:
epiMT_CAGEs_selectTime <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GeneralEffect_FANTOM_ER/epiMT_exprs_CAGEs_2026.csv", header = T)
dim(epiMT_CAGEs_selectTime)
colnames(epiMT_CAGEs_selectTime)
length(epiMT_libs)
#matches

#double check remove zeroes has worked via unix:
trial <- epiMT_CAGEs_selectTime
length(trial[,1]) == sum(!apply(trial[,-1], 1, function(row) all(row == 0)))

#find best expressed CAGEs per gene, this mapping object created alongside in eddie code:
#likely will exclude lncRNAs if use a >1 cut-off?
epiMT_CAGEs_mapping <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GeneralEffect_FANTOM_ER/epiMT_CAGEs_mapping_2026.csv", header = T)
length(unique(epiMT_CAGEs_mapping$geneID))#22268
length(unique(epiMT_exprs))#match

#all genes have a CAGE site
#dont need transcripts:
epiMT_CAGEs_mapping <- unique(epiMT_CAGEs_mapping[,c(1,3)])

#there will be loads of crap TSS in these samples
#take all that provide >20% of a gene's output:
epiMT_CAGEs_selectTime <- unique(merge(epiMT_CAGEs_selectTime, epiMT_CAGEs_mapping, by.x = "prmtrID", by.y = "CAGEClusterID"))

trial <- split(epiMT_CAGEs_selectTime, epiMT_CAGEs_selectTime$geneID)

triali <- lapply(trial, function(x){
  x$prmtrID[rowMeans(x[,2:26])/sum(rowMeans(x[,2:26]))>0.2]
})

perc20 <- filter(epiMT_CAGEs_selectTime, prmtrID %in% unlist(triali) )

length(unique(perc20$geneID))#missing a tiny minority of genes (18 only)

#for remaining go to 10%:
perc10 <- filter(epiMT_CAGEs_selectTime, !geneID %in% perc20$geneID)

trial <- split(perc10, perc10$geneID)

triali <- lapply(trial, function(x){
  x$prmtrID[rowMeans(x[,2:26])/sum(rowMeans(x[,2:26]))>0.1]
})

perc10 <- filter(perc10, prmtrID %in% unlist(triali) )

perc20_10 <- rbind(perc20, perc10)

length(unique(perc20_10$geneID))
#22267 genes
#no need for remaining 1 - too diffuse and weak expression

#liftOver to hg38 and filter:
perc20_10$TSS_start <- as.numeric(sapply(strsplit(sapply(strsplit(perc20_10$prmtrID, "\\:"), "[[", 2), "\\.\\."), "[[", 1))
perc20_10$TSS_stop <- as.numeric(gsub(",[+-]", "", sapply(strsplit(sapply(strsplit(perc20_10$prmtrID, "\\:"), "[[", 2), "\\.\\."), "[[", 2)))
perc20_10$chr <- sapply(strsplit(perc20_10$prmtrID, "\\:"), "[[", 1)
head(perc20_10)

library(rtracklayer)
chain <- import.chain("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/hg19ToHg38.over.chain")

trial <- makeGRangesFromDataFrame(perc20_10[,c(1,28:30)], start.field = "TSS_start", end.field = "TSS_stop", 
                                  seqnames.field = "chr", keep.extra.columns = T)

x <- liftOver(trial, chain)

trial <- data.frame(unlist(x))

#midpoint of each CAGE:
trial$mid <- trial$start + trial$width/2
perc20_10 <- merge(perc20_10, trial[,6:7], by = "prmtrID")
length(unique(perc20_10$geneID))
#some lost in liftOver

#biotypes:
trial <- merge(perc20_10, GeneBiotypes[,c(1,9)], by.x = "geneID", by.y = "CAT_geneID")
table(trial$CAT_geneClassII)
length(unique(trial$geneID))

#all G
TSSperG <- unique(trial[,c(1,2,28:32)])

length(unique(TSSperG$geneID))

#relevant biotypes only:
table(TSSperG$CAT_geneClassII)
TSSperG <- filter(TSSperG, CAT_geneClassII %in% c("coding_mRNA", "lncRNA"))

#for each gene
trial <- split(TSSperG, TSSperG$geneID)

#find all genes starting within 250kbp
triali <- lapply(trial, function(x){
  filter(TSSperG, chr == unique(x$chr), mid > min(x$mid) - 250000, mid < max(x$mid) + 250000)
})
trialii <- bind_rows(triali, .id = "AnchorGene")

#isolate self matches, re-insert genes with no neighbours
selfMatch <- filter(trialii, AnchorGene == geneID)
trialii <- filter(trialii, !AnchorGene == geneID)
selfMatchOnly <- filter(selfMatch, !AnchorGene %in% trialii$AnchorGene)
selfMatchOnly$geneID <- NA
#add together:
trialii <- rbind(trialii, selfMatchOnly)

trialii <- unique(trialii)

#identify PCG-PCG co-reg pairs (for later use)
AllPCG_AllPCG_1 <- filter(trialii, AnchorGene %in% filter(GeneBiotypes, CAT_geneClassII == "coding_mRNA")$CAT_geneID)
AllPCG_AllPCG_1 <- filter(AllPCG_AllPCG_1, geneID %in% filter(GeneBiotypes, CAT_geneClassII == "coding_mRNA")$CAT_geneID)
AllPCG_AllPCG_1$pairs <- paste(AllPCG_AllPCG_1$AnchorGene, AllPCG_AllPCG_1$geneID, sep = "-")

#remove reverse matches - NO, makes it a pain later on:
#AllPCG_AllPCG_1$pairs_ordered <- apply(AllPCG_AllPCG_1[, c(1,2)], 1, function(x) paste(sort(x), collapse = "-"))

#mostly pairs, some higher
#table(table(AllPCG_AllPCG_1$pairs_ordered))

#for each sorted pair label, remove any after first row which are duplicates:
#trial <- split(AllPCG_AllPCG_1, AllPCG_AllPCG_1$pairs_ordered)

#triali <- lapply(trial, function(x){
#  x[1,]
#})

#triali <- bind_rows(triali)
#AllPCG_AllPCG_1 <- unique(triali[,c(1:2,8,10)])

table(table(AllPCG_AllPCG_1$AnchorGene))

length(unique(AllPCG_AllPCG_1$AnchorGene))#12310 anchor genes (bit meaningless if removing duplicate pairs)
length(unique(c(AllPCG_AllPCG_1$AnchorGene, AllPCG_AllPCG_1$geneID)))#12310 PCGs expressed near another PCG

length(unique(AllPCG_AllPCG_1$pairs))#86356 pairs, many will be duplicates in terms of PCGx -> PCGy and PCGy -> PCGx

triali <- unique(merge(AllPCG_AllPCG_1[,-c(3:7)], TSSperG[,c(1,6)], by.x = "AnchorGene", by.y = "geneID"))
triali <- unique(merge(triali, TSSperG[,c(1,6)], by.x = "geneID", by.y = "geneID"))
triali$AbsDistLnc_PCG <- abs(triali$mid.x - triali$mid.y)/1000
triali$DistLnc_PCG <- (triali$mid.x - triali$mid.y)/1000
triali$pairs <- paste(triali$AnchorGene, triali$geneID, sep = "-")

#shortest distance per pair:
trialii <- split(triali, triali$pairs)

trialii[[101]]
trialiii <- lapply(trialii, function(x){
  unique(x[x$AbsDistLnc_PCG == min((x$AbsDistLnc_PCG)),])
})

trial <- unique(bind_rows(trialiii))

AllPCG_AllPCG_1 <- trial
#write.csv(AllLNC_AllPCG_1, "AllLNC_AllPCG_1_epiMT.csv", row.names = F)


#identify lncRNA-PCG co-reg pairs
AllLNC_AllPCG_1 <- filter(trialii, AnchorGene %in% filter(GeneBiotypes, CAT_geneClassII == "lncRNA")$CAT_geneID)
AllLNC_AllPCG_1 <- filter(AllLNC_AllPCG_1, geneID %in% filter(GeneBiotypes, CAT_geneClassII == "coding_mRNA")$CAT_geneID)

length(unique(AllLNC_AllPCG_1$AnchorGene))#5313 lncRNAs near a PCG
length(unique(AllLNC_AllPCG_1$geneID))#11265 PCGs near a lncRNA

#plot distance between
#n.b. multiple TSS per gene, take min distance between genes
#add the mid (middle of CAGE cluster) for lnc and neighbour
triali <- unique(merge(AllLNC_AllPCG_1[,-c(3:7)], TSSperG[,c(1,6)], by.x = "AnchorGene", by.y = "geneID"))
triali <- unique(merge(triali, TSSperG[,c(1,6)], by.x = "geneID", by.y = "geneID"))
triali$AbsDistLnc_PCG <- abs(triali$mid.x - triali$mid.y)/1000
triali$DistLnc_PCG <- (triali$mid.x - triali$mid.y)/1000
triali$pairs <- paste(triali$AnchorGene, triali$geneID, sep = "-")

#shortest distance per pair:
trialii <- split(triali, triali$pairs)

trialii[[101]]
trialiii <- lapply(trialii, function(x){
  unique(x[x$AbsDistLnc_PCG == min((x$AbsDistLnc_PCG)),])
})

trial <- bind_rows(trialiii)

AllLNC_AllPCG_1 <- trial
#write.csv(AllLNC_AllPCG_1, "AllLNC_AllPCG_1_epiMT_2026.csv", row.names = F)

#option for a closest neighbours analysis:
closestNeighbour_Lncs_epiMT <- split(AllLNC_AllPCG_1, AllLNC_AllPCG_1$AnchorGene)

#some variations:
#closest neighbour
#closestNeighbour <- lapply(closestNeighbour, function(x){
#  filter(x, AbsDistLnc_PCG == min(AbsDistLnc_PCG ))
#  }
#  )

#closest x neighbours within 250kbp etc:
#closestNeighbour <- lapply(closestNeighbour, function(x){
#    x[order(x$DisLnc_PCG, decreasing = F),][1:5,]}
#  )

#closest "surrounding" neighbours (consider up and down-stream):
#this selected for best "naive" results
closestNeighbour_Lncs_epiMT <- lapply(closestNeighbour_Lncs_epiMT, function(x){
  upstream <- filter(x, DistLnc_PCG < 0)
  downstream <- filter(x, DistLnc_PCG > 0)
  
  return(rbind(upstream[order(upstream$DistLnc_PCG, decreasing = T),][1,],
               downstream[order(downstream$DistLnc_PCG, decreasing = F),][1,]))
}
)

closestNeighbour_Lncs_epiMT <- bind_rows(closestNeighbour_Lncs_epiMT)
closestNeighbour_Lncs_epiMT <- filter(closestNeighbour_Lncs_epiMT, !is.na(pairs))


#now get the cis candidates:
#same timeframe only for simplicity
CoRegPairs_sameTimeframe_epiMT <- filter(closestNeighbour_Lncs_epiMT, 
                                       (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_2$geneID, epiMT_DEGs_DE_Downwithin_2$geneID) & 
                                          geneID %in% c(epiMT_DEGs_DE_Upwithin_2$geneID, epiMT_DEGs_DE_Downwithin_2$geneID)) |
                                         
                                         (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_4$geneID, epiMT_DEGs_DE_Downwithin_4$geneID) & 
                                            geneID %in% c(epiMT_DEGs_DE_Upwithin_4$geneID, epiMT_DEGs_DE_Downwithin_4$geneID)) |
                                         
                                         (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_6$geneID, epiMT_DEGs_DE_Downwithin_6$geneID) & 
                                            geneID %in% c(epiMT_DEGs_DE_Upwithin_6$geneID, epiMT_DEGs_DE_Downwithin_6$geneID)) |
                                         
                                         (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_8$geneID, epiMT_DEGs_DE_Downwithin_8$geneID) & 
                                            geneID %in% c(epiMT_DEGs_DE_Upwithin_8$geneID, epiMT_DEGs_DE_Downwithin_8$geneID)) |
                                         
                                         (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_12$geneID, epiMT_DEGs_DE_Downwithin_12$geneID) & 
                                            geneID %in% c(epiMT_DEGs_DE_Upwithin_12$geneID, epiMT_DEGs_DE_Downwithin_12$geneID)) |
                                         
                                         (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_16$geneID, epiMT_DEGs_DE_Downwithin_16$geneID) & 
                                            geneID %in% c(epiMT_DEGs_DE_Upwithin_16$geneID, epiMT_DEGs_DE_Downwithin_16$geneID)) |
                                         
                                         (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_42$geneID, epiMT_DEGs_DE_Downwithin_42$geneID) & 
                                            geneID %in% c(epiMT_DEGs_DE_Upwithin_42$geneID, epiMT_DEGs_DE_Downwithin_42$geneID)) |
                                         
                                         (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_60$geneID, epiMT_DEGs_DE_Downwithin_60$geneID) & 
                                            geneID %in% c(epiMT_DEGs_DE_Upwithin_60$geneID, epiMT_DEGs_DE_Downwithin_60$geneID))
                                       )
#250 CCLnc-target pairs

length(unique(CoRegPairs_sameTimeframe_epiMT$AnchorGene))#237 CClncRNAs
length(unique(CoRegPairs_sameTimeframe_epiMT$geneID))#236 potential targets


#### (needs revisit) Compare to control lncs, and SVSMC pairs ####

#compare to control cis lncs
ControlCisLncs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/ControlCisLncs.csv")

expressed_epiMT <- filter(epiMT_CPM, rowMax(as.matrix(epiMT_CPM[,c(13:16)])) >1)
sum(gsub("\\.[0-9]*", "", expressed_epiMT$geneID) %in% ControlCisLncs$x)#29 expressed
sum(unique(gsub("\\.[0-9]*", "", epiMT_DEGs_DE$CAT_geneID)) %in% ControlCisLncs$x)#5 are DE
sum(unique(gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended_epiMT$AnchorGene)) %in% ControlCisLncs$x)#2 in CClncs


#shared cis candidates between datasets:
CoRegPairs_04_48_24_extended <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CoRegPairs_04_48_24_extendedIII.csv")

#insert FANTOM ID where possible:
Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime.csv", header = T)
CoRegPairs_04_48_24_extended <- unique(merge(Enhancer_lociII[,c(1,14)], CoRegPairs_04_48_24_extended, by = "EnsID", all.y = T))

CoRegPairs_04_48_24_extended_epiMT$pairsII <- paste(
  gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended_epiMT$AnchorGene),
  gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended_epiMT$geneID),
  sep = "-"
)

CoRegPairs_04_48_24_extended$pairsII <- paste(
  gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended$FANTOM_ID),
  gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended$EnsID.y),
  sep = "-"
)


#11 x pairs found in both
sum(CoRegPairs_04_48_24_extended_epiMT$pairsII %in% CoRegPairs_04_48_24_extended$pairsII)

filter(CoRegPairs_04_48_24_extended, pairsII %in% CoRegPairs_04_48_24_extended_epiMT$pairsII)[,c(32,11)]
#INKILN is notable, NR2F2-AS1 and PITPNA-AS1 key loci, latter hits SERPINF1 which is not classed as influential (or a cis target) in SMC (locus is highlighted by other approach)
11/248 #4% of SVSMC pairs are common
11/834 #1% of epiMT pairs are common


#12x cclncRNAs predicted in both
sum(unique(gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended_epiMT$AnchorGene)) %in% gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended$FANTOM_ID))
length(unique(CoRegPairs_04_48_24_extended$FANTOM_ID))#106, but one is NA
12/106 #11% of cclncRNAs are found in both
length(unique(CoRegPairs_04_48_24_extended_epiMT$AnchorGene))#324
12/324 #4% of cclncRNAs are fond in both


#25X cclncRNA targets in both
sum(unique(gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended_epiMT$geneID)) %in% gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended$EnsID.y))
length(unique(CoRegPairs_04_48_24_extended$EnsID.y))#106, but one is NA
25/236 #11% of cclncRNAs are found in both
length(unique(CoRegPairs_04_48_24_extended_epiMT$geneID))#532
25/532 #5% of cclncRNAs are fond in both

#consider overlapping the GO terms too - handle this later on through all datasets


#### Recreate the SVSMC analysis - import necessary table here ####

#Temporal bias vs. other DEGs
#Elevation in FC for lnc targets vs. other PCGs in early phase
#Increased chance of DE PCG neighbour for DE lnc targets vs other lncs* and DE PCGs in early phase
#GO terms
# *required HiC/eQTL to work

#need DEGs:
#write.csv(epiMT_DEGs_DE, "epiMT_DEGs_DE_2026.csv", row.names = F)

#lnc-PCG neighbour table:
#write.csv(AllLNC_AllPCG_1, "AllLNC_AllPCG_1_epiMT_2026.csv", row.names = F)

#PCG-PCG neighbour table:
#write.csv(AllPCG_AllPCG_1, "AllPCG_AllPCG_1_epiMT_2026.csv", row.names = F)

#expressed genes table:
#write.csv(epiMT_CPM[,c(35,26:34)], "epiMT_CPM_2026.csv", row.names = F)

#import (needs changing on revisit):
epiMT_DEGs_DE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GeneralEffect_FANTOM_ER/epiMT_DEGs_DE_2026.csv")
epiMT_CPM <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GeneralEffect_FANTOM_ER/epiMT_CPM_2026.csv")
AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GeneralEffect_FANTOM_ER/AllLNC_AllPCG_1_epiMT_2026.csv")
AllPCG_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GeneralEffect_FANTOM_ER/AllPCG_AllPCG_1_epiMT_2026.csv")
GeneBiotypes <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/nature21374-s2/table11.csv")
#simplify lncRNA class:
GeneBiotypes$CAT_geneClassII <- GeneBiotypes$CAT_geneClass
GeneBiotypes$CAT_geneClassII[
  grepl("lncRNA", GeneBiotypes$CAT_geneClassII)] <- "lncRNA"

#n.b. (START FROM HERE ON RETURN)
epiMT_DEGs_DE_Upwithin_2 <- filter(epiMT_DEGs_DE, RegulationStart == "Induced <2hrs")
epiMT_DEGs_DE_Downwithin_2 <- filter(epiMT_DEGs_DE, RegulationStart == "Repressed <2hrs")

epiMT_DEGs_DE_Upwithin_4 <- filter(epiMT_DEGs_DE, RegulationStart == "Induced 2-4hrs")
epiMT_DEGs_DE_Downwithin_4 <- filter(epiMT_DEGs_DE, RegulationStart == "Repressed 2-4hrs")

epiMT_DEGs_DE_Upwithin_6 <- filter(epiMT_DEGs_DE, RegulationStart == "Induced 4-6hrs")
epiMT_DEGs_DE_Downwithin_6 <- filter(epiMT_DEGs_DE, RegulationStart == "Repressed 4-6hrs")

epiMT_DEGs_DE_Upwithin_8 <- filter(epiMT_DEGs_DE, RegulationStart == "Induced 6-8hrs")
epiMT_DEGs_DE_Downwithin_8 <- filter(epiMT_DEGs_DE, RegulationStart == "Repressed 6-8hrs")

epiMT_DEGs_DE_Upwithin_12 <- filter(epiMT_DEGs_DE, RegulationStart == "Induced 8-12hrs")
epiMT_DEGs_DE_Downwithin_12 <- filter(epiMT_DEGs_DE, RegulationStart == "Repressed 8-12hrs")

epiMT_DEGs_DE_Upwithin_16 <- filter(epiMT_DEGs_DE, RegulationStart == "Induced 12-16hrs")
epiMT_DEGs_DE_Downwithin_16 <- filter(epiMT_DEGs_DE, RegulationStart == "Repressed 12-16hrs")

epiMT_DEGs_DE_Upwithin_42 <- filter(epiMT_DEGs_DE, RegulationStart == "Induced 16-42hrs")
epiMT_DEGs_DE_Downwithin_42 <- filter(epiMT_DEGs_DE, RegulationStart == "Repressed 16-42hrs")

epiMT_DEGs_DE_Upwithin_60 <- filter(epiMT_DEGs_DE, RegulationStart == "Induced 42-60hrs")
epiMT_DEGs_DE_Downwithin_60 <- filter(epiMT_DEGs_DE, RegulationStart == "Repressed 42-60hrs")


closestNeighbour_Lncs_epiMT <- split(AllLNC_AllPCG_1, AllLNC_AllPCG_1$AnchorGene)

#some variations:
#closest neighbour
#closestNeighbour <- lapply(closestNeighbour, function(x){
#  filter(x, AbsDistLnc_PCG == min(AbsDistLnc_PCG ))
#  }
#  )

#closest x neighbours within 250kbp etc:
#closestNeighbour <- lapply(closestNeighbour, function(x){
#    x[order(x$DisLnc_PCG, decreasing = F),][1:5,]}
#  )

#closest "surrounding" neighbours (consider up and down-stream):
#this selected for best "naive" results
closestNeighbour_Lncs_epiMT <- lapply(closestNeighbour_Lncs_epiMT, function(x){
  upstream <- filter(x, DistLnc_PCG < 0)
  downstream <- filter(x, DistLnc_PCG > 0)
  
  return(rbind(upstream[order(upstream$DistLnc_PCG, decreasing = T),][1,],
               downstream[order(downstream$DistLnc_PCG, decreasing = F),][1,]))
}
)

closestNeighbour_Lncs_epiMT <- bind_rows(closestNeighbour_Lncs_epiMT)
closestNeighbour_Lncs_epiMT <- filter(closestNeighbour_Lncs_epiMT, !is.na(pairs))


#now get the cis candidates:
#same timeframe only for simplicity
CoRegPairs_sameTimeframe_epiMT <- filter(closestNeighbour_Lncs_epiMT, 
                                         (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_2$CAT_geneID, epiMT_DEGs_DE_Downwithin_2$CAT_geneID) & 
                                            geneID %in% c(epiMT_DEGs_DE_Upwithin_2$CAT_geneID, epiMT_DEGs_DE_Downwithin_2$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_4$CAT_geneID, epiMT_DEGs_DE_Downwithin_4$CAT_geneID) & 
                                              geneID %in% c(epiMT_DEGs_DE_Upwithin_4$CAT_geneID, epiMT_DEGs_DE_Downwithin_4$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_6$CAT_geneID, epiMT_DEGs_DE_Downwithin_6$CAT_geneID) & 
                                              geneID %in% c(epiMT_DEGs_DE_Upwithin_6$CAT_geneID, epiMT_DEGs_DE_Downwithin_6$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_8$CAT_geneID, epiMT_DEGs_DE_Downwithin_8$CAT_geneID) & 
                                              geneID %in% c(epiMT_DEGs_DE_Upwithin_8$CAT_geneID, epiMT_DEGs_DE_Downwithin_8$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_12$CAT_geneID, epiMT_DEGs_DE_Downwithin_12$CAT_geneID) & 
                                              geneID %in% c(epiMT_DEGs_DE_Upwithin_12$CAT_geneID, epiMT_DEGs_DE_Downwithin_12$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_16$CAT_geneID, epiMT_DEGs_DE_Downwithin_16$CAT_geneID) & 
                                              geneID %in% c(epiMT_DEGs_DE_Upwithin_16$CAT_geneID, epiMT_DEGs_DE_Downwithin_16$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_42$CAT_geneID, epiMT_DEGs_DE_Downwithin_42$CAT_geneID) & 
                                              geneID %in% c(epiMT_DEGs_DE_Upwithin_42$CAT_geneID, epiMT_DEGs_DE_Downwithin_42$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_60$CAT_geneID, epiMT_DEGs_DE_Downwithin_60$CAT_geneID) & 
                                              geneID %in% c(epiMT_DEGs_DE_Upwithin_60$CAT_geneID, epiMT_DEGs_DE_Downwithin_60$CAT_geneID))
)
#250 CCLnc-target pairs

length(unique(CoRegPairs_sameTimeframe_epiMT$AnchorGene))#237 CClncRNAs
length(unique(CoRegPairs_sameTimeframe_epiMT$geneID))#236 potential targets

#expressed lncs:
length(unique(filter(GeneBiotypes, 
              CAT_geneClassII == "lncRNA", 
              CAT_geneID %in% epiMT_CPM$geneID[Biobase::rowMax(as.matrix(epiMT_CPM[,2:10])) >1]
              )$CAT_geneID))

epiMT_CPM$geneID[Biobase::rowMax(as.matrix(epiMT_CPM[,2:10])) >1]


#
#### bias of IEGs ####

#IEGs
#inserting new gene categories:
IEGs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Gene lists/arner_2015_table_S5_IEGs_lit.csv")
IEGs_hs <- filter(IEGs, !Hs_symbol %in% c("", "-"))
IEGs_hs <- filter(IEGs_hs, !grepl("[a-z]", IEGs_hs$Total_count))

epiMT_DEGs_DE$IEG <- epiMT_DEGs_DE$GeneClassUpdate
epiMT_DEGs_DE$IEG[epiMT_DEGs_DE$CAT_geneName %in% IEGs_hs$Hs_symbol] <- "IEG"

Selected <- "IEG"
table(epiMT_DEGs_DE$IEG)[Selected]
table(epiMT_DEGs_DE$IEG)/dim(epiMT_DEGs_DE)[1]*100

Cluster_biotype <- as.data.frame(table(epiMT_DEGs_DE$IEG, epiMT_DEGs_DE$RegulationStart))
trial <- filter(Cluster_biotype, grepl(Selected, Var1))
trial$Freq/table(epiMT_DEGs_DE$RegulationStart)*100
#looks like definite biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(epiMT_DEGs_DE$RegulationStart)

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

ClusterNames <- trial$Var2

IEGEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(epiMT_DEGs_DE$IEG)[Selected]
  d <- dim(epiMT_DEGs_DE)[1]
  
  IEGEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p)
}

names(IEGEnrich_cluster) <- ClusterNames
triali <- bind_rows(IEGEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")
#all significant biases

#percentage plots
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-42hrs", "42-60hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,5,7,8,2,3,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(epiMT_DEGs_DE$IEG)[Selected]*100
trial$PercBackground <- trial$selection/dim(epiMT_DEGs_DE)[1]*100

IEGWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(IEGWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(IEGWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(IEGWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(IEGWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 4.2) +
  geom_col(data = filter(IEGWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(IEGWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(IEGWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory, label = Freq), size = 4.2) +
  ylab("% DE IEGs") +
  xlab("") +
  #scale_y_continuous(limits = c(-30,35),breaks = seq(-30,30, by = 10),
  #                   labels = (c(seq(30, 0, by = -10), seq(10,30,by=10)))) +
  theme_minimal() +
  theme(text = element_text(size = 15))


#
#### biases of lncRNAs ####

table(epiMT_DEGs_DE$CAT_geneClassII)
table(epiMT_DEGs_DE$CAT_geneClassII)/dim(epiMT_DEGs_DE)[1]*100

Cluster_biotype <- as.data.frame(table(epiMT_DEGs_DE$CAT_geneClassII, epiMT_DEGs_DE$RegulationStart))
trial <- filter(Cluster_biotype, Var1 == "lncRNA")
trial$Freq/table(epiMT_DEGs_DE$RegulationStart)*100
#looks like weak biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(epiMT_DEGs_DE$RegulationStart)

ClusterNames <- trial$Var2

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

LncEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(epiMT_DEGs_DE$CAT_geneClassII)["lncRNA"]
  d <- dim(epiMT_DEGs_DE)[1]
  
  LncEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p)
}

names(LncEnrich_cluster) <- ClusterNames
triali <- bind_rows(LncEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")

#percentage plots
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-42hrs", "42-60hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,5,7,8,2,3,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(epiMT_DEGs_DE$CAT_geneClassII)["lncRNA"]*100
trial$PercBackground <- trial$selection/dim(epiMT_DEGs_DE)[1]*100

LncWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(LncWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(LncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(LncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(LncWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 5) +
  geom_col(data = filter(LncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(LncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(LncWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory, label = Freq), size = 5) +
  ylab("% DE LncRNAs") +
  xlab("") +
  #scale_y_continuous(limits = c(-30,35),breaks = seq(-30,30, by = 10),
  #                   labels = (c(seq(30, 0, by = -10), seq(10,30,by=10)))) +
  theme_minimal() +
  theme(text = element_text(size=15))


#
#### cis lnc biases ####

#plot enrichment of cis lnc, non cis lnc, e cis lnc, non-e cis lnc:
#plot genes per cluster

#cis lnc
epiMT_DEGs_DE$CisLnc <- epiMT_DEGs_DE$CAT_geneClassII
epiMT_DEGs_DE$CisLnc[epiMT_DEGs_DE$CAT_geneID %in% CoRegPairs_sameTimeframe_epiMT$AnchorGene] <- "CisLnc"
table(epiMT_DEGs_DE$CisLnc)

Cluster_biotype <- as.data.frame(table(epiMT_DEGs_DE$CisLnc, epiMT_DEGs_DE$RegulationStart))

#bias of category:
table(epiMT_DEGs_DE$CisLnc)["CisLnc"]
table(epiMT_DEGs_DE$CisLnc)/dim(epiMT_DEGs_DE)[1]*100

trial <- filter(Cluster_biotype, Var1 == "CisLnc")
trial$Freq/table(epiMT_DEGs_DE$RegulationStart)*100
#looks like strong biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(epiMT_DEGs_DE$RegulationStart)

ClusterNames <- trial$Var2

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

LncEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(epiMT_DEGs_DE$CisLnc)["CisLnc"]
  d <- dim(epiMT_DEGs_DE)[1]
  
  LncEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater"
                                       )$p)
}

names(LncEnrich_cluster) <- ClusterNames
triali <- bind_rows(LncEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")
#Induced first timepoint, rarely at second
#2x significant biases

#percentage plots
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-42hrs", "42-60hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,5,7,8,2,3,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(epiMT_DEGs_DE$CisLnc)["CisLnc"]*100
trial$PercBackground <- trial$selection/dim(epiMT_DEGs_DE)[1]*100

CisLncWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(CisLncWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(CisLncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(CisLncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CisLncWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 5) +
  geom_col(data = filter(CisLncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(CisLncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CisLncWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory, label = Freq), size = 5) +
  ylab("% DE CisLncRNAs") +
  xlab("") +
  #scale_y_continuous(limits = c(-30,60),breaks = seq(-20,60, by = 20),
  #                   labels = (c(seq(20, 0, by = -20), seq(20,60,by=20)))) +
  theme_minimal() +
  theme(text = element_text(size = 15))

#cis lncRNAs again biased to early activation, may have their primary role early on

#compare to all lncs:
sum(LncWaveBias$Freq)
sum(CisLncWaveBias$Freq)

a <- 80
b <- 237
c <- 171
d <- 1199

fisher.test(data.frame("LncRNA" = c(a,b-a),
                       "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")

trial <- data.frame(c/d*100, a/b*100)

colnames(trial) <- c("All DE lncRNAs", "CClncRNAs")

#melt(trial)

trial <- melt(trial)

trial$variable <- factor(trial$variable)
trial$variable <- factor(trial$variable, levels = levels(trial$variable)[c(1,3,2)])

ggplot(trial) + aes(x = value, y = variable) +
  geom_bar(stat = "identity", fill = "olivedrab3", color = "grey60") +
  xlab("% 0-2hr induced") +
  ylab("") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0,40,20), limits = c(0,40)) +
  theme(text = element_text(size =24))


#
#### elncRNAs ####

#alternate (less specific) indicator of cis-acting function would be transcription from an annotated enhancer

epiMT_DEGs_DE$CAT_geneClassIII <- epiMT_DEGs_DE$CAT_geneClassII
epiMT_DEGs_DE$CAT_geneClassIII[epiMT_DEGs_DE$CAT_geneID %in% filter(GeneBiotypes, CAT_geneCategory == "e_lncRNA")$CAT_geneID] <- "e_lncRNA"

table(epiMT_DEGs_DE$CAT_geneClassIII)

Cluster_biotype <- as.data.frame(table(epiMT_DEGs_DE$CAT_geneClassIII, epiMT_DEGs_DE$RegulationStart))

table(epiMT_DEGs_DE$CAT_geneClassIII)

#bias of category:
table(epiMT_DEGs_DE$CAT_geneClassIII)["e_lncRNA"]
table(epiMT_DEGs_DE$CAT_geneClassIII)/dim(epiMT_DEGs_DE)[1]*100

trial <- filter(Cluster_biotype, Var1 == "e_lncRNA")
trial$Freq/table(epiMT_DEGs_DE$RegulationStart)*100

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(epiMT_DEGs_DE$RegulationStart)

ClusterNames <- trial$Var2

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

LncEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(epiMT_DEGs_DE$CAT_geneClassIII)["e_lncRNA"]
  d <- dim(epiMT_DEGs_DE)[1]
  
  LncEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p)
}

names(LncEnrich_cluster) <- ClusterNames
triali <- bind_rows(LncEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")
#Induced first timepoint, rarely at second
#2x significant biases

#percentage plots
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-42hrs", "42-60hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,5,7,8,2,3,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(epiMT_DEGs_DE$CAT_geneClassIII)["e_lncRNA"]*100
trial$PercBackground <- trial$selection/dim(epiMT_DEGs_DE)[1]*100

ELncWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(ELncWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(ELncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(ELncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(ELncWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 5) +
  geom_col(data = filter(ELncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(ELncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(ELncWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory, label = Freq), size = 5) +
  ylab("% DE Enhancer\nlncRNAs") +
  xlab("") +
  #scale_y_continuous(limits = c(-30,60),breaks = seq(-20,60, by = 20),
  #                   labels = (c(seq(20, 0, by = -20), seq(20,60,by=20)))) +
  theme_minimal() +
  theme(text = element_text(size = 15))

#early bias

#percentage of cclncRNAs amongst DE elncRNAs vs other DE lncRNAs:
sum(ELncWaveBias$Freq)
sum(LncWaveBias$Freq)
sum(CisLncWaveBias$Freq)

sum(filter(epiMT_DEGs_DE, CAT_geneClassIII == "e_lncRNA")$CAT_geneID %in% filter(epiMT_DEGs_DE, CisLnc == "CisLnc")$CAT_geneID)

64/237 #27% of cclncs are elncs
277/1199 #23% of all DELs are elncs

aii <- 64
bii <- 237    
cii <- 277   
dii <- 1199   

fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

trial <- data.frame(cii/dii*100, aii/bii*100)

colnames(trial) <- c("All DE lncRNAs", "CClncRNAs (same)")

#melt(trial)

trial <- melt(trial)

trial$variable <- factor(trial$variable)
trial$variable <- factor(trial$variable, levels = levels(trial$variable)[c(1,3,2)])

ggplot(trial) + aes(x = value, y = variable) +
  geom_bar(stat = "identity", fill = "olivedrab3", color = "grey60") +
  xlab("% ElncRNAs") +
  ylab("") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0,40,20), limits = c(0,45)) +
  theme(text = element_text(size =20))



#
#### all co-regulated (time-agnostic) ####

#plot enrichment of cis lnc, non cis lnc, e cis lnc, non-e cis lnc:
#plot genes per cluster

CoReg <- filter(AllLNC_AllPCG_1, AnchorGene %in% epiMT_DEGs_DE$CAT_geneID, geneID %in% epiMT_DEGs_DE$CAT_geneID)

#cis lnc
epiMT_DEGs_DE$CoRegLnc <- epiMT_DEGs_DE$CAT_geneClassII
epiMT_DEGs_DE$CoRegLnc[epiMT_DEGs_DE$CAT_geneID %in% CoReg$AnchorGene] <- "CoRegLnc"
table(epiMT_DEGs_DE$CoRegLnc)

Cluster_biotype <- as.data.frame(table(epiMT_DEGs_DE$CoRegLnc, epiMT_DEGs_DE$RegulationStart))

#bias of category:
table(epiMT_DEGs_DE$CoRegLnc)["CoRegLnc"]
table(epiMT_DEGs_DE$CoRegLnc)/dim(epiMT_DEGs_DE)[1]*100

trial <- filter(Cluster_biotype, Var1 == "CoRegLnc")
trial$Freq/table(epiMT_DEGs_DE$RegulationStart)*100
#looks like strong biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(epiMT_DEGs_DE$RegulationStart)

ClusterNames <- trial$Var2

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

LncEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(epiMT_DEGs_DE$CoRegLnc)["CoRegLnc"]
  d <- dim(epiMT_DEGs_DE)[1]
  
  LncEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), 
                                                              row.names = c("Cluster", "other")), alternative = "greater")$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), 
                                                              row.names = c("Cluster", "other")), alternative = "greater")$p)
}
names(LncEnrich_cluster) <- ClusterNames
triali <- bind_rows(LncEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")
#Induced first timepoint, rarely at second
#2x significant biases

#percentage plots
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-42hrs", "42-60hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,5,7,8,2,3,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(epiMT_DEGs_DE$CoRegLnc)["CoRegLnc"]*100
trial$PercBackground <- trial$selection/dim(epiMT_DEGs_DE)[1]*100

CoRegLncWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(CoRegLncWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(CoRegLncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(CoRegLncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CoRegLncWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 5) +
  geom_col(data = filter(CoRegLncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(CoRegLncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CoRegLncWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory, label = Freq), size = 5) +
  ylab("% DE CoRegLncRNAs") +
  xlab("") +
  #scale_y_continuous(limits = c(-30,60),breaks = seq(-20,60, by = 20),
  #                   labels = (c(seq(20, 0, by = -20), seq(20,60,by=20)))) +
  theme_minimal() +
  theme(text = element_text(size = 15))

#basically same as lncs generally


#
#### combined figure (and table save) ####

#combined figure:
AllTypesWaveBias <- rbind(#IEGWaveBias,
                          LncWaveBias, 
                          ELncWaveBias, #NonELncWaveBias, 
                          #CoRegLncWaveBias, #NonCoRegLncWaveBias,
                          CisLncWaveBias#, CisLncSameWaveBias #NonCisLncWaveBias
                          )

AllTypesWaveBias$FirstRegulation <- as.character(AllTypesWaveBias$FirstRegulation)
#AllTypesWaveBias$FirstRegulation[grepl("24hrs", AllTypesWaveBias$FirstRegulation)] <- "Within \n16hrs"

AllTypesWaveBias$OR_corrected <- AllTypesWaveBias$OR
AllTypesWaveBias$OR_corrected[AllTypesWaveBias$p_adj >0.1] <- NA
AllTypesWaveBias$OR_corrected <- AllTypesWaveBias$OR_corrected + 0.001
AllTypesWaveBias$`Log2(Odds Ratio)` <- log2(AllTypesWaveBias$OR_corrected)

AllTypesWaveBias$`Log2(Odds Ratio)`[AllTypesWaveBias$`Log2(Odds Ratio)` < -5] <- -5
AllTypesWaveBias$`Log2(Odds Ratio)`[AllTypesWaveBias$`Log2(Odds Ratio)` >5] <- 5

AllTypesWaveBias$padj_corrected <- AllTypesWaveBias$p_adj
AllTypesWaveBias$padj_corrected[AllTypesWaveBias$p_adj >0.1] <- NA

AllTypesWaveBias$padj_simple <- NA
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.1] <- "p<0.1"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.05] <- "p<0.05"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.01] <- "p<0.01"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.001] <- "p<0.001"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.0001] <- "p<0.0001"

AllTypesWaveBias$UpDownType <- paste(AllTypesWaveBias$Var1, AllTypesWaveBias$UpDown, sep = "-")

#save a copy:
#write.csv(AllTypesWaveBias, "epiMT_AllTypesWaveBias_2026.csv", row.names = F)


#insert spacers for plotting:
AllTypesWaveBias <- AllTypesWaveBias[,c(16,6,13,15)]

#appendSpacers <- data.frame(UpDownType = c(rep("Spacer1",3),
#                                           rep("Spacer2",3),
#                                           rep("Spacer3",3),
#                                           rep("Spacer4",3)#,
                                           #rep("Spacer5",3),
                                           #rep("Spacer6",3),
                                           #rep("Spacer7",3),
                                           #rep("Spacer8",3)
#), 
#"FR" = rep(c("Within \n4hrs",
#             "Within \n8hrs",
#             "Within \n16hrs"), 4), 
#"OR" = NA,  
#`padj_simple` =NA)

#colnames(appendSpacers) <- colnames(AllTypesWaveBias)

#AllTypesWaveBias <- rbind(AllTypesWaveBias, appendSpacers)

AllTypesWaveBias$UpDownType <- factor(AllTypesWaveBias$UpDownType)
AllTypesWaveBias$UpDownType <- factor(AllTypesWaveBias$UpDownType, 
                                      levels = levels(AllTypesWaveBias$UpDownType)[c(2,1,4,3,8,7,6,5)])

AllTypesWaveBias$FirstRegulation <- factor(AllTypesWaveBias$FirstRegulation)
AllTypesWaveBias$FirstRegulation <- factor(AllTypesWaveBias$FirstRegulation, 
                                      levels = levels(AllTypesWaveBias$FirstRegulation)[c(1,4,5,7,8,2,3,6)])

#8 space version: c(16,15,17,26,25,18,10,9,19,2,1,20,12,11,21,4,3,22,14,13,23,6,5,24,8,7)

myColor <- colorRampPalette(c("steelblue", "white", "red"))(50)
myBreaks <- c(seq(min(AllTypesWaveBias$`Log2(Odds Ratio)`, na.rm = T), 0, 
                  length.out=ceiling(50/2)), 
              seq(max(AllTypesWaveBias$`Log2(Odds Ratio)`, na.rm = T)/50, 
                  max(AllTypesWaveBias$`Log2(Odds Ratio)`, na.rm = T), 
                  length.out=floor(50/2)))

ggplot(AllTypesWaveBias) + aes(x = FirstRegulation, y = UpDownType, size = padj_simple, fill = `Log2(Odds Ratio)`) +
  geom_point(color = "grey30", shape = 21) +
  xlab("") +
  #ylab("") +
  scale_size_discrete(range = c(3,12), limits = c(#"p<0.1", 
                                                  "p<0.05", "p<0.01", "p<0.001", "p<0.0001")) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "red") +
  theme_minimal() +
  theme(legend.key.size = unit(1.4, "line"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        axis.text.x = element_text(size=22),
        #axis.text.y = element_blank()
  ) + Seurat::RotatedAxis()


#
#### bias of CCPCGs ####

closestNeighbour_PCGs_epiMT <- split(AllPCG_AllPCG_1, AllPCG_AllPCG_1$AnchorGene)

#some variations:
#closest neighbour
#closestNeighbour <- lapply(closestNeighbour, function(x){
#  filter(x, AbsDistLnc_PCG == min(AbsDistLnc_PCG ))
#  }
#  )

#closest x neighbours within 250kbp etc:
#closestNeighbour <- lapply(closestNeighbour, function(x){
#    x[order(x$DisLnc_PCG, decreasing = F),][1:5,]}
#  )

#closest "surrounding" neighbours (consider up and down-stream):
#this selected for best "naive" results
closestNeighbour_PCGs_epiMT <- lapply(closestNeighbour_PCGs_epiMT, function(x){
  upstream <- filter(x, DistLnc_PCG < 0)
  downstream <- filter(x, DistLnc_PCG > 0)
  
  return(rbind(upstream[order(upstream$DistLnc_PCG, decreasing = T),][1,],
               downstream[order(downstream$DistLnc_PCG, decreasing = F),][1,]))
}
)

closestNeighbour_PCGs_epiMT <- bind_rows(closestNeighbour_PCGs_epiMT)
closestNeighbour_PCGs_epiMT <- filter(closestNeighbour_PCGs_epiMT, !is.na(pairs))


#now get the cis candidates:
#same timeframe only for simplicity
CoRegPairs_sameTimeframe_epiMT_PCG <- filter(closestNeighbour_PCGs_epiMT, 
                                         (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_2$CAT_geneID, epiMT_DEGs_DE_Downwithin_2$CAT_geneID) & 
                                            geneID %in% c(epiMT_DEGs_DE_Upwithin_2$CAT_geneID, epiMT_DEGs_DE_Downwithin_2$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_4$CAT_geneID, epiMT_DEGs_DE_Downwithin_4$CAT_geneID) & 
                                              geneID %in% c(epiMT_DEGs_DE_Upwithin_4$CAT_geneID, epiMT_DEGs_DE_Downwithin_4$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_6$CAT_geneID, epiMT_DEGs_DE_Downwithin_6$CAT_geneID) & 
                                              geneID %in% c(epiMT_DEGs_DE_Upwithin_6$CAT_geneID, epiMT_DEGs_DE_Downwithin_6$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_8$CAT_geneID, epiMT_DEGs_DE_Downwithin_8$CAT_geneID) & 
                                              geneID %in% c(epiMT_DEGs_DE_Upwithin_8$CAT_geneID, epiMT_DEGs_DE_Downwithin_8$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_12$CAT_geneID, epiMT_DEGs_DE_Downwithin_12$CAT_geneID) & 
                                              geneID %in% c(epiMT_DEGs_DE_Upwithin_12$CAT_geneID, epiMT_DEGs_DE_Downwithin_12$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_16$CAT_geneID, epiMT_DEGs_DE_Downwithin_16$CAT_geneID) & 
                                              geneID %in% c(epiMT_DEGs_DE_Upwithin_16$CAT_geneID, epiMT_DEGs_DE_Downwithin_16$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_42$CAT_geneID, epiMT_DEGs_DE_Downwithin_42$CAT_geneID) & 
                                              geneID %in% c(epiMT_DEGs_DE_Upwithin_42$CAT_geneID, epiMT_DEGs_DE_Downwithin_42$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_60$CAT_geneID, epiMT_DEGs_DE_Downwithin_60$CAT_geneID) & 
                                              geneID %in% c(epiMT_DEGs_DE_Upwithin_60$CAT_geneID, epiMT_DEGs_DE_Downwithin_60$CAT_geneID))
)
#913 CCLnc-target pairs

length(unique(CoRegPairs_sameTimeframe_epiMT_PCG$AnchorGene))#863 CClncRNAs
length(unique(CoRegPairs_sameTimeframe_epiMT_PCG$geneID))#843 potential targets

#cis lnc
epiMT_DEGs_DE$CisPCG <- epiMT_DEGs_DE$CAT_geneClassII
epiMT_DEGs_DE$CisPCG[epiMT_DEGs_DE$CAT_geneID %in% c(CoRegPairs_sameTimeframe_epiMT_PCG$AnchorGene)] <- "CisPCG"
table(epiMT_DEGs_DE$CisPCG)

Cluster_biotype <- as.data.frame(table(epiMT_DEGs_DE$CisPCG, epiMT_DEGs_DE$RegulationStart))

#bias of category:
table(epiMT_DEGs_DE$CisPCG)["CisPCG"]
table(epiMT_DEGs_DE$CisPCG)/dim(epiMT_DEGs_DE)[1]*100

trial <- filter(Cluster_biotype, Var1 == "CisPCG")
trial$Freq/table(epiMT_DEGs_DE$RegulationStart)*100
#looks like strong biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(epiMT_DEGs_DE$RegulationStart)

ClusterNames <- trial$Var2

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

PCGEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(epiMT_DEGs_DE$CisPCG)["CisPCG"]
  d <- dim(epiMT_DEGs_DE)[1]
  
  LncEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p)
}
names(LncEnrich_cluster) <- ClusterNames
triali <- bind_rows(LncEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")
#Induced first timepoint, rarely at second
#2x significant biases

#percentage plots
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-42hrs", "42-60hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,5,7,8,2,3,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(epiMT_DEGs_DE$CisPCG)["CisPCG"]*100
trial$PercBackground <- trial$selection/dim(epiMT_DEGs_DE)[1]*100

CisPCGWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(CisPCGWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(CisPCGWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(CisPCGWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CisPCGWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 5) +
  geom_col(data = filter(CisPCGWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(CisPCGWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CisPCGWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory, label = Freq), size = 5) +
  ylab("% DE CisPCG") +
  xlab("") +
  #scale_y_continuous(limits = c(-30,60),breaks = seq(-20,60, by = 20),
  #                   labels = (c(seq(20, 0, by = -20), seq(20,60,by=20)))) +
  theme_minimal() +
  theme(text = element_text(size = 15))

#cis PCGs somewhat early

#compare to DE lncs with co-reg neighbours:
sum(CisLncWaveBias$Freq)
sum(CisPCGWaveBias$Freq)

a <- 80 #hit: 0-2 induced
b <- 237 #selection - all CC lncs
c <- 144
d <- 863 #background - all CC PCGs
fisher.test(data.frame("cisLnc" = c(a, b-a),
                       "Not"   = c(c, d-c)), alternative = "greater")

#plot
trial <- data.frame("variable" = c("CClncRNAs", "DE PCG with\nco-regulated\nneighbour"),
                         "value" = c(a/b*100, c/d*100))

trial$variable <- as.factor(trial$variable)
trial$variable <- factor(trial$variable, levels(trial$variable)[c(2,1)])

ggplot(trial) + aes(y = variable, x = value, fill = variable) +
  geom_bar(stat = "identity", color = "grey60") +
  scale_fill_manual(values = c(`CClncRNAs` = "olivedrab3", `DE PCG with\nco-regulated\nneighbour` = "mediumorchid")) +
  xlab("% 0-2hr induced") +
  ylab(" ") +
  theme_minimal() + 
  theme(text = element_text(size=24)) + 
  theme(legend.position = "none") #+ Seurat::RotatedAxis()


#
#### same timeframe and PCG comparison ####

#combined figure:
AllTypesWaveBias <- rbind(LncWaveBias, #PCGWaveBias, #the PCG result is weird, unbalanced table so OR is way off...in future remove conditional ORs
                          CisLncWaveBias, CisLncSameWaveBias, 
                          CisPCGWaveBias, CisPCGSameWaveBias)

AllTypesWaveBias$OR_corrected <- AllTypesWaveBias$OR
AllTypesWaveBias$OR_corrected[AllTypesWaveBias$p_adj >0.1] <- NA
AllTypesWaveBias$OR_corrected <- AllTypesWaveBias$OR_corrected + 0.001
AllTypesWaveBias$`Log2(Odds Ratio)` <- log2(AllTypesWaveBias$OR_corrected)

AllTypesWaveBias$`Log2(Odds Ratio)`[AllTypesWaveBias$`Log2(Odds Ratio)` < -5] <- -5
AllTypesWaveBias$`Log2(Odds Ratio)`[AllTypesWaveBias$`Log2(Odds Ratio)` >5] <- 5

AllTypesWaveBias$padj_corrected <- AllTypesWaveBias$p_adj
AllTypesWaveBias$padj_corrected[AllTypesWaveBias$p_adj >0.1] <- NA

AllTypesWaveBias$padj_simple <- NA
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.1] <- "p<0.1"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.05] <- "p<0.05"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.01] <- "p<0.01"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.001] <- "p<0.001"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.0001] <- "p<0.0001"

AllTypesWaveBias$UpDownType <- paste(AllTypesWaveBias$Var1, AllTypesWaveBias$UpDown, sep = "-")

#might be clearer to do a comparison at time point of interest, 4 hours (induced) only:
PCG_PCG_4hr_timing <- filter(AllTypesWaveBias, FirstRegulation == "Within \n4hrs", UpDown == "Induced")

PCG_PCG_4hr_timing$Var1 <- gsub("Same", "\n(same timeframe)", PCG_PCG_4hr_timing$Var1)

PCG_PCG_4hr_timing$Var1 <- factor(PCG_PCG_4hr_timing$Var1)
PCG_PCG_4hr_timing$Var1 <- factor(PCG_PCG_4hr_timing$Var1, levels = levels(PCG_PCG_4hr_timing$Var1)[c(5,1,2,3,4)])

ggplot(PCG_PCG_4hr_timing) +
  #scale_y_continuous(limits = c(0,15)) +
  geom_col(aes(x = Var1, y = PercCategory), fill = "steelblue") +
  geom_col(aes(x = Var1, y = PercBackground), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  #geom_label(aes(x = Var1, y = DEP_Pairs, label = NoDEP_Pairs), size = 3) +  
  ylab("% in 0-4hr induced DEGs") +
  xlab("\nGene biotype\nvs.\nall Gene biotypes") +
  theme_minimal() + Seurat::RotatedAxis()

#or a direct comparison for CClncs:
#numbers from:
length(unique(c(CoRegPairs_04_48_24_extended_epiMT$AnchorGene)))#324 CClncs
length(unique(c(CoRegPairs_04_48_24_extended_epiMT_PCG$AnchorGene)))#2132 PCG
PCG_PCG_4hr_timing

a <- 95
b <- 324 #selection - all CC lncs
c <- 448
d <- 2132 #background - all CC PCGs
fisher.test(data.frame("cisLnc" = c(a, b-a),
                       "Not"   = c(c, d-c)), alternative = "greater")#***
statmod::power.fisher.test(a/b, c/d, c, d)


#or a direct comparison for CClncs in same timeframe:
length(unique(c(CoRegPairs_04_48_24_epiMT_same$AnchorGene)))#259 sametimeframe CClncs
length(unique(c(CoRegPairs_04_48_24_epiMT_samePCG$AnchorGene)))#1567 PCG
PCG_PCG_4hr_timing

ai <- 68
bi <- 259 #selection all cclncs, same timeframe
ci <- 240
di <- 1567  #background - all CC lncs, same timeframe
fisher.test(data.frame("cisLnc" = c(ai, bi-ai),
                       "Not"   = c(ci, di-ci)), alternative = "greater")#****
statmod::power.fisher.test(ai/bi, ci/di, ci, di)#so ~20-30% chance to find a stat difference

CClnc_CCPCG_direct <- data.frame("hr4_induced" = c(a/(b)*100, ai/(bi)*100), 
                                 "all_CClncs_CCPCGs" = c(c/(d)*100, ci/(di)*100))

CClnc_CCPCG_direct$PairType <- as.factor(c("All types", "Same\ntimeframe only"))
#CClnc_CCPCG_direct$PairType <- factor(CClnc_CCPCG_direct$PairType, levels = levels(CClnc_CCPCG_direct$PairType)[c(2,1)])

#CClnc_CCPCG_direct$NoDEP_Pairs <- c(a,ai)

ggplot(CClnc_CCPCG_direct) +
  #scale_y_continuous(limits = c(0,12)) +
  geom_col(aes(x = PairType, y = hr4_induced), fill = "steelblue") +
  geom_col(aes(x = PairType, y = all_CClncs_CCPCGs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  #geom_label(aes(x = PairType, y = DEP_Pairs, label = NoDEP_Pairs), size = 3) +  
  ylab("% in 0-4hr induced DEGs") +
  xlab("\nCCLncRNAs\nvs.\nCCPCGs") +
  theme_minimal()


fisher.test(data.frame("cisLnc" = c(a, b-a),
                       "Not"   = c(c, d-c)), alternative = "greater")$p *2
fisher.test(data.frame("cisLnc" = c(ai, bi-ai),
                       "Not"   = c(ci, di-ci)), alternative = "greater")$p * 2


#repressed
PCG_PCG_4hr_timing <- filter(AllTypesWaveBias, FirstRegulation == "Within \n4hrs", UpDown == "Repressed")

PCG_PCG_4hr_timing$Var1 <- gsub("Same", "\n(same timeframe)", PCG_PCG_4hr_timing$Var1)

PCG_PCG_4hr_timing$Var1 <- factor(PCG_PCG_4hr_timing$Var1)
PCG_PCG_4hr_timing$Var1 <- factor(PCG_PCG_4hr_timing$Var1, levels = levels(PCG_PCG_4hr_timing$Var1)[c(5,1,2,3,4)])

ggplot(PCG_PCG_4hr_timing) +
  #scale_y_continuous(limits = c(0,15)) +
  geom_col(aes(x = Var1, y = PercCategory), fill = "steelblue") +
  geom_col(aes(x = Var1, y = PercBackground), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  #geom_label(aes(x = Var1, y = DEP_Pairs, label = NoDEP_Pairs), size = 3) +  
  ylab("% in 0-4hr repressed DEGs") +
  xlab("Gene class\nvs.\nall Gene classes") +
  theme_minimal() + Seurat::RotatedAxis()

#weak enrichment amongst lncs as seen in the dotplot


#### enrichment of gene induction near induced lncs vs near expressed lncs ####

#background of totally non-DE lncs overall whole timecourse, expressed in given timepoint:
epiMT_CPM_02 <- filter(epiMT_CPM, (Hour0_meanCPM>1 | Hour2_meanCPM>1), 
                       !geneID %in% filter(epiMT_DEGs_DE, !grepl("<2hrs", RegulationStart))$CAT_geneID)

epiMT_CPM_24 <- filter(epiMT_CPM, (Hour2_meanCPM>1 | Hour4_meanCPM>1), 
                       !geneID %in% filter(epiMT_DEGs_DE, !grepl("2-4hrs", RegulationStart))$CAT_geneID)

epiMT_CPM_46 <- filter(epiMT_CPM, (Hour4_meanCPM>1 | Hour6_meanCPM>1), 
                       !geneID %in% filter(epiMT_DEGs_DE, !grepl("4-6hrs", RegulationStart))$CAT_geneID)

epiMT_CPM_68 <- filter(epiMT_CPM, (Hour6_meanCPM>1 | Hour8_meanCPM>1), 
                       !geneID %in% filter(epiMT_DEGs_DE, !grepl("6-8hrs", RegulationStart))$CAT_geneID)

epiMT_CPM_812 <- filter(epiMT_CPM, (Hour8_meanCPM>1 | Hour12_meanCPM>1), 
                        !geneID %in% filter(epiMT_DEGs_DE, !grepl("8-12hrs", RegulationStart))$CAT_geneID)

epiMT_CPM_1216 <- filter(epiMT_CPM, (Hour12_meanCPM>1 | Hour16_meanCPM>1), 
                         !geneID %in% filter(epiMT_DEGs_DE, !grepl("12-16hrs", RegulationStart))$CAT_geneID)

epiMT_CPM_1642 <- filter(epiMT_CPM, (Hour16_meanCPM>1 | Hour42_meanCPM>1), 
                         !geneID %in% filter(epiMT_DEGs_DE, !grepl("16-42hrs", RegulationStart))$CAT_geneID)

epiMT_CPM_4260 <- filter(epiMT_CPM, (Hour42_meanCPM>1 | Hour60_meanCPM>1), 
                         !geneID %in% filter(epiMT_DEGs_DE, !grepl("42-60hrs", RegulationStart))$CAT_geneID)

#(n.b. alternate for same-timeframe is to just take expressed and not DE in previous timeframe)
backgrounds_list <- list(epiMT_CPM_02, epiMT_CPM_24, epiMT_CPM_46, epiMT_CPM_68,
                         epiMT_CPM_812, epiMT_CPM_1216, epiMT_CPM_1642, epiMT_CPM_4260)

#just up DELs:
selection_list <- list(epiMT_DEGs_DE_Upwithin_2, epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_6, 
                       epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_12, epiMT_DEGs_DE_Upwithin_16,
                       epiMT_DEGs_DE_Upwithin_42, epiMT_DEGs_DE_Upwithin_60)
#with up targets:
hits_list <- list(epiMT_DEGs_DE_Upwithin_2, epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_6, 
                  epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_12, epiMT_DEGs_DE_Upwithin_16,
                  epiMT_DEGs_DE_Upwithin_42, epiMT_DEGs_DE_Upwithin_60)

results_list <- list()

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$CAT_geneID)
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], 
                            geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            #!geneID %in% nonSelection_list[[i]]$CAT_geneID
  )$geneID
  ))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], 
                            geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            #!geneID %in% nonSelection_list[[i]]$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID
  ))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], 
                            geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], 
                            geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, 
                                           results_list[[2]][1]/results_list[[2]][2]*100, 
                                           results_list[[3]][1]/results_list[[3]][2]*100,
                                           results_list[[4]][1]/results_list[[4]][2]*100,
                                           results_list[[5]][1]/results_list[[5]][2]*100,
                                           results_list[[6]][1]/results_list[[6]][2]*100,
                                           results_list[[7]][1]/results_list[[7]][2]*100,
                                           results_list[[8]][1]/results_list[[8]][2]*100),
                           
                           "EL_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, 
                                          results_list[[2]][3]/results_list[[2]][4]*100, 
                                          results_list[[3]][3]/results_list[[3]][4]*100,
                                          results_list[[4]][3]/results_list[[4]][4]*100,
                                          results_list[[5]][3]/results_list[[5]][4]*100,
                                          results_list[[6]][3]/results_list[[6]][4]*100,
                                          results_list[[7]][3]/results_list[[7]][4]*100,
                                          results_list[[8]][3]/results_list[[8]][4]*100))

DEL_PCG_type$RegulationBegins <- as.factor(c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-42hrs", "42-60hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1], results_list[[4]][1],
                        results_list[[5]][1], results_list[[6]][1], results_list[[7]][1], results_list[[8]][1])

g2bi <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = EL_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17,23,29,35,41,47)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18,24,30,36,42,48)]
#write.csv(DEL_PCG_type, "epiMT_inducedDEL_DEP_DEneighboursame_conc.csv", row.names = F)

p.adjust(DEL_PCG_type$p, method = "bonferroni")

DEL_PCG_type_plot <- unique(melt(DEL_PCG_type[,1:3]))

DEL_PCG_type_plot$variable <- as.character(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable[grepl("^DEL", DEL_PCG_type_plot$variable)] <- "Induced lncRNAs"
DEL_PCG_type_plot$variable[grepl("^EL", DEL_PCG_type_plot$variable)] <- "Non-DE lncRNAs"

DEL_PCG_type_plot$variable <- as.factor(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable <- factor(DEL_PCG_type_plot$variable, levels = levels(DEL_PCG_type_plot$variable)[c(1,2)])

DEL_PCG_type_plot$RegulationBegins <- as.factor(DEL_PCG_type_plot$RegulationBegins)
DEL_PCG_type_plot$RegulationBegins <- factor(DEL_PCG_type_plot$RegulationBegins, 
                                             levels = levels(DEL_PCG_type_plot$RegulationBegins)[c(1,4,5,7,8,2,3,6)])

ggplot(DEL_PCG_type_plot) + aes(x = RegulationBegins, y = value, fill = variable) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.7), color = "grey60") +
  theme_minimal() +
  theme(text = element_text(size = 24)) + #, legend.position = "none") +
  scale_fill_manual(values = c(`Induced lncRNAs` = "#D6604D", `Non-DE lncRNAs` = "grey30")) +
  xlab("") +
  ylab("% With induced\nneighbour") + Seurat::RotatedAxis()


#
#### enrichment of gene repression near induced lncs vs near expressed lncs ####

#background of totally non-DE lncs overall whole timecourse, expressed in given timepoint:
epiMT_CPM_02 <- filter(epiMT_CPM, (Hour0_meanCPM>1 | Hour2_meanCPM>1), 
                       !geneID %in% filter(epiMT_DEGs_DE, !grepl("<2hrs", RegulationStart))$CAT_geneID)

epiMT_CPM_24 <- filter(epiMT_CPM, (Hour2_meanCPM>1 | Hour4_meanCPM>1), 
                       !geneID %in% filter(epiMT_DEGs_DE, !grepl("2-4hrs", RegulationStart))$CAT_geneID)

epiMT_CPM_46 <- filter(epiMT_CPM, (Hour4_meanCPM>1 | Hour6_meanCPM>1), 
                       !geneID %in% filter(epiMT_DEGs_DE, !grepl("4-6hrs", RegulationStart))$CAT_geneID)

epiMT_CPM_68 <- filter(epiMT_CPM, (Hour6_meanCPM>1 | Hour8_meanCPM>1), 
                       !geneID %in% filter(epiMT_DEGs_DE, !grepl("6-8hrs", RegulationStart))$CAT_geneID)

epiMT_CPM_812 <- filter(epiMT_CPM, (Hour8_meanCPM>1 | Hour12_meanCPM>1), 
                        !geneID %in% filter(epiMT_DEGs_DE, !grepl("8-12hrs", RegulationStart))$CAT_geneID)

epiMT_CPM_1216 <- filter(epiMT_CPM, (Hour12_meanCPM>1 | Hour16_meanCPM>1), 
                         !geneID %in% filter(epiMT_DEGs_DE, !grepl("12-16hrs", RegulationStart))$CAT_geneID)

epiMT_CPM_1642 <- filter(epiMT_CPM, (Hour16_meanCPM>1 | Hour42_meanCPM>1), 
                         !geneID %in% filter(epiMT_DEGs_DE, !grepl("16-42hrs", RegulationStart))$CAT_geneID)

epiMT_CPM_4260 <- filter(epiMT_CPM, (Hour42_meanCPM>1 | Hour60_meanCPM>1), 
                         !geneID %in% filter(epiMT_DEGs_DE, !grepl("42-60hrs", RegulationStart))$CAT_geneID)

#(n.b. alternate for same-timeframe is to just take expressed and not DE in previous timeframe)
backgrounds_list <- list(epiMT_CPM_02, epiMT_CPM_24, epiMT_CPM_46, epiMT_CPM_68,
                         epiMT_CPM_812, epiMT_CPM_1216, epiMT_CPM_1642, epiMT_CPM_4260)

#just up DELs:
selection_list <- list(epiMT_DEGs_DE_Upwithin_2, epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_6, 
                       epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_12, epiMT_DEGs_DE_Upwithin_16,
                       epiMT_DEGs_DE_Upwithin_42, epiMT_DEGs_DE_Upwithin_60)
#with up targets:
hits_list <- list(epiMT_DEGs_DE_Downwithin_2, epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_6, 
                  epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_12, epiMT_DEGs_DE_Downwithin_16,
                  epiMT_DEGs_DE_Downwithin_42, epiMT_DEGs_DE_Downwithin_60)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$CAT_geneID)
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], 
                            geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            #!geneID %in% nonSelection_list[[i]]$CAT_geneID
  )$geneID
  ))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], 
                            geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            #!geneID %in% nonSelection_list[[i]]$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID
  ))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], 
                            geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], 
                            geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, 
                                           results_list[[2]][1]/results_list[[2]][2]*100, 
                                           results_list[[3]][1]/results_list[[3]][2]*100,
                                           results_list[[4]][1]/results_list[[4]][2]*100,
                                           results_list[[5]][1]/results_list[[5]][2]*100,
                                           results_list[[6]][1]/results_list[[6]][2]*100,
                                           results_list[[7]][1]/results_list[[7]][2]*100,
                                           results_list[[8]][1]/results_list[[8]][2]*100),
                           
                           "EL_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, 
                                          results_list[[2]][3]/results_list[[2]][4]*100, 
                                          results_list[[3]][3]/results_list[[3]][4]*100,
                                          results_list[[4]][3]/results_list[[4]][4]*100,
                                          results_list[[5]][3]/results_list[[5]][4]*100,
                                          results_list[[6]][3]/results_list[[6]][4]*100,
                                          results_list[[7]][3]/results_list[[7]][4]*100,
                                          results_list[[8]][3]/results_list[[8]][4]*100))

DEL_PCG_type$RegulationBegins <- as.factor(c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-42hrs", "42-60hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1], results_list[[4]][1],
                        results_list[[5]][1], results_list[[6]][1], results_list[[7]][1], results_list[[8]][1])

g2bi <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = EL_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17,23,29,35,41,47)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18,24,30,36,42,48)]
#write.csv(DEL_PCG_type, "epiMT_inducedDEL_DEP_DEneighboursame_conc.csv", row.names = F)

p.adjust(DEL_PCG_type$p, method = "bonferroni")

DEL_PCG_type_plot <- unique(melt(DEL_PCG_type[,1:3]))

DEL_PCG_type_plot$variable <- as.character(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable[grepl("^DEL", DEL_PCG_type_plot$variable)] <- "Induced lncRNAs"
DEL_PCG_type_plot$variable[grepl("^EL", DEL_PCG_type_plot$variable)] <- "Non-DE lncRNAs"

DEL_PCG_type_plot$variable <- as.factor(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable <- factor(DEL_PCG_type_plot$variable, levels = levels(DEL_PCG_type_plot$variable)[c(1,2)])

DEL_PCG_type_plot$RegulationBegins <- as.factor(DEL_PCG_type_plot$RegulationBegins)
DEL_PCG_type_plot$RegulationBegins <- factor(DEL_PCG_type_plot$RegulationBegins, 
                                             levels = levels(DEL_PCG_type_plot$RegulationBegins)[c(1,4,5,7,8,2,3,6)])

ggplot(DEL_PCG_type_plot) + aes(x = RegulationBegins, y = value, fill = variable) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.7), color = "grey60") +
  theme_minimal() +
  theme(text = element_text(size = 24)) + #, legend.position = "none") +
  scale_fill_manual(values = c(`Induced lncRNAs` = "#D6604D", `Non-DE lncRNAs` = "grey30")) +
  xlab("") +
  ylab("% With repressed\nneighbour") + Seurat::RotatedAxis()


#
#### enrichment of gene induction near induced lncs vs near induced PCGs ####

#just up DELs:
backgrounds_list <- list(epiMT_DEGs_DE_Upwithin_2, epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_6, 
                         epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_12, epiMT_DEGs_DE_Upwithin_16,
                         epiMT_DEGs_DE_Upwithin_42, epiMT_DEGs_DE_Upwithin_60)
#with up targets:
hits_list <- list(epiMT_DEGs_DE_Upwithin_2, epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_6, 
                  epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_12, epiMT_DEGs_DE_Upwithin_16,
                  epiMT_DEGs_DE_Upwithin_42, epiMT_DEGs_DE_Upwithin_60)

results_list <- list()

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$CAT_geneID)
  
  #PCG genes
  d <- length(unique(filter(backgrounds_list[[i]], CAT_geneID %in% filter(GeneBiotypes, grepl("coding_mRNA", CAT_geneClass))$CAT_geneID)$CAT_geneID))
  #all hits in PCGs
  c <- length(unique(filter(backgrounds_list[[i]], CAT_geneID %in% filter(GeneBiotypes, grepl("coding_mRNA", CAT_geneClass))$CAT_geneID,
                            CAT_geneID %in% filter(AllPCG_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$CAT_geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], CAT_geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$CAT_geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], CAT_geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            CAT_geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$CAT_geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "PCG" = c(c,d-c), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "PCG" = c(c,d-c), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, 
                                           results_list[[2]][1]/results_list[[2]][2]*100, 
                                           results_list[[3]][1]/results_list[[3]][2]*100,
                                           results_list[[4]][1]/results_list[[4]][2]*100,
                                           results_list[[5]][1]/results_list[[5]][2]*100,
                                           results_list[[6]][1]/results_list[[6]][2]*100,
                                           results_list[[7]][1]/results_list[[7]][2]*100,
                                           results_list[[8]][1]/results_list[[8]][2]*100),
                           
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, 
                                           results_list[[2]][3]/results_list[[2]][4]*100, 
                                           results_list[[3]][3]/results_list[[3]][4]*100,
                                           results_list[[4]][3]/results_list[[4]][4]*100,
                                           results_list[[5]][3]/results_list[[5]][4]*100,
                                           results_list[[6]][3]/results_list[[6]][4]*100,
                                           results_list[[7]][3]/results_list[[7]][4]*100,
                                           results_list[[8]][3]/results_list[[8]][4]*100))

DEL_PCG_type$RegulationBegins <- as.factor(c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-42hrs", "42-60hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1], results_list[[4]][1],
                        results_list[[5]][1], results_list[[6]][1], results_list[[7]][1], results_list[[8]][1])

g2bi <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17,23,29,35,41,47)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18,24,30,36,42,48)]
#write.csv(DEL_PCG_type, "epiMT_inducedDEL_DEP_DEneighboursame_conc.csv", row.names = F)

p.adjust(DEL_PCG_type$p, method = "bonferroni")

DEL_PCG_type_plot <- unique(melt(DEL_PCG_type[,1:3]))

DEL_PCG_type_plot$variable <- as.character(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable[grepl("^DEL", DEL_PCG_type_plot$variable)] <- "Induced lncRNAs"
DEL_PCG_type_plot$variable[grepl("DEP", DEL_PCG_type_plot$variable)] <- "Induced PCGs"

DEL_PCG_type_plot$variable <- as.factor(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable <- factor(DEL_PCG_type_plot$variable, levels = levels(DEL_PCG_type_plot$variable)[c(1,2)])

DEL_PCG_type_plot$RegulationBegins <- as.factor(DEL_PCG_type_plot$RegulationBegins)
DEL_PCG_type_plot$RegulationBegins <- factor(DEL_PCG_type_plot$RegulationBegins, 
                                             levels = levels(DEL_PCG_type_plot$RegulationBegins)[c(1,4,5,7,8,2,3,6)])

ggplot(DEL_PCG_type_plot) + aes(x = RegulationBegins, y = value, fill = variable) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.7), color = "grey60") +
  theme_minimal() +
  theme(text = element_text(size = 24)) + #, legend.position = "none") +
  scale_fill_manual(values = c(`Induced lncRNAs` = "#D6604D", `Induced PCGs` = "grey80")) +
  xlab("") +
  ylab("% With induced\nneighbour") + Seurat::RotatedAxis()


#
#### GO/KEGG/etc terms early-induced targets ####

library(clusterProfiler)
library(org.Hs.eg.db)

epiMT_DE_PCGs_DE_Upwithin_2 <- filter(epiMT_DEGs_DE_Upwithin_2, grepl("coding_mRNA", CAT_geneClassII))

EarlyInduced_CClncRNAs_CoUp_Targets_epiMT <- filter(CoRegPairs_sameTimeframe_epiMT, AnchorGene %in% epiMT_DEGs_DE_Upwithin_2$CAT_geneID,
                                              geneID %in% epiMT_DEGs_DE_Upwithin_2$CAT_geneID)

#targets of early lncs vs other DEGs in that timeframe
CoUp_DE <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", EarlyInduced_CClncRNAs_CoUp_Targets_epiMT$geneID)),
                    universe      = gsub("\\.[0-9]*", "", epiMT_DE_PCGs_DE_Upwithin_2$CAT_geneID),
                    keyType       = "ENSEMBL",
                    OrgDb         = org.Hs.eg.db,
                    ont           = "all",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
CoUp_DE_df <- as.data.frame(CoUp_DE)
#0x terms

#loosen the background, expressed PCGs:
epiMT_02_exprs_PCG <- filter(epiMT_CPM, geneID %in% filter(GeneBiotypes, CAT_geneClassII == "coding_mRNA")$CAT_geneID, 
                             (Hour0_meanCPM>1 | Hour2_meanCPM>1))

CoUp_DE2 <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", EarlyInduced_CClncRNAs_CoUp_Targets_epiMT$geneID)),
                     universe      = gsub("\\.[0-9]*", "", epiMT_02_exprs_PCG$geneID),
                     keyType       = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
CoUp_DE2_df <- as.data.frame(CoUp_DE2)
#much better
#endothelial/epithelial activation terms, immune activation terms

dotplot(simplify(CoUp_DE2), showCategory = 10)


#
#### specific target lists early-induced targets ####

#IEGs (inc. lots of immune)
IEGs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Gene lists/arner_2015_table_S5_IEGs_lit.csv")
IEGs_hs <- filter(IEGs, !Hs_symbol %in% c("", "-"))
IEGs_hs <- filter(IEGs_hs, !grepl("[a-z]", IEGs_hs$Total_count))
IEGs_hs$Hs_symbol[IEGs_hs$Hs_symbol == "IL8"] <- "CXCL8"

#Epithelial-biased
FANT_S10 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/nature21374-s2/table 10.csv")
FANT_S10_Epith <- filter(FANT_S10, grepl("epithelial", sample_ontology_term))
FANT_S10_Epith_G <- gsub("\\.[0-9]*", "", unlist(strsplit(FANT_S10_Epith$associated_geneID, ",")))

#Epigenetic modifiers (LISA has some definitions but not well described, CRdb seems recent and fine)
CRdb_dataB <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CRdb Data Browse.csv")

TF_Lambert2018 <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/TF_Lambert2018.txt", header = T, stringsAsFactors = F)

#add HMGA2 and FOXL1 to this list, unclear why not in there...

GeneLists <- list("TFs"= unique(filter(epiMT_DEGs_DE, gsub("\\.[0-9]*", "", CAT_geneID) %in% TF_Lambert2018$ID)$CAT_geneID),
                  "CRs" = unique(filter(epiMT_DEGs_DE, CAT_geneName %in% c(CRdb_dataB$CR, "FOXL1", "HMGA2"))$CAT_geneID),
                  "IEGs"= unique(filter(epiMT_DEGs_DE, CAT_geneName %in% c(IEGs_hs[,4]))$CAT_geneID),
                  "Epith"= unique(filter(epiMT_DEGs_DE, gsub("\\.[0-9]*", "", CAT_geneID) %in% FANT_S10_Epith_G)$CAT_geneID))

results_list <- list()

#same timeframe
epiMT_DE_PCGs_DE_Upwithin_2 <- filter(epiMT_DEGs_DE_Upwithin_2, grepl("coding_mRNA", CAT_geneClassII))

EarlyInduced_CClncRNAs_CoUp_Targets_epiMT <- filter(CoRegPairs_sameTimeframe_epiMT, AnchorGene %in% epiMT_DEGs_DE_Upwithin_2$CAT_geneID,
                                                    geneID %in% epiMT_DEGs_DE_Upwithin_2$CAT_geneID)

for (i in 1:length(GeneLists)){
  
  a <- length(unique(EarlyInduced_CClncRNAs_CoUp_Targets_epiMT$geneID[ EarlyInduced_CClncRNAs_CoUp_Targets_epiMT$geneID %in% GeneLists[[i]] ]))
  b <- length(unique(EarlyInduced_CClncRNAs_CoUp_Targets_epiMT$geneID))
  c <- length(unique(epiMT_DE_PCGs_DE_Upwithin_2$CAT_geneID[ epiMT_DE_PCGs_DE_Upwithin_2$CAT_geneID %in% GeneLists[[i]]]))
  d <- length(unique(epiMT_DE_PCGs_DE_Upwithin_2$CAT_geneID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate,
                                  paste(unique(EarlyInduced_CClncRNAs_CoUp_Targets_epiMT$geneID[ 
                                    EarlyInduced_CClncRNAs_CoUp_Targets_epiMT$geneID %in% GeneLists[[i]] 
                                    ]), collapse = "/"))
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:8] <- c("pval", "OR", "Genes")
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")
#some enrichment of IEG, others not bothered
GeneLists_enrichedSameT <- GeneLists_enriched

filter(GeneBiotypes, CAT_geneID %in% unique(EarlyInduced_CClncRNAs_CoUp_Targets_epiMT$geneID[ 
  EarlyInduced_CClncRNAs_CoUp_Targets_epiMT$geneID %in% GeneLists[["IEGs"]] 
  ]))[,1:3]



#
#### examine early lncs/targets ####

#add gene names
trial <- merge(CoRegPairs_sameTimeframe_epiMT, GeneBiotypes[,1:3], by.x = "AnchorGene", by.y = "CAT_geneID")
trial <- merge(trial, GeneBiotypes[,1:3], by.x = "geneID", by.y = "CAT_geneID")

CoRegPairs_sameTimeframe_epiMT_early <- filter(trial, AnchorGene %in% epiMT_DEGs_DE_Upwithin_2$CAT_geneID,
                                               geneID %in% epiMT_DEGs_DE_Upwithin_2$CAT_geneID)

#notable targets: HMGA2, LIF, CXCL3
#notable early lncs: RP11-221N13.3

#control cis lncs:
ControlCisLncs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Write-ups/supplement/ControlCisLncs.csv", header = T)

#2x identified in the spike of cis-acting lnc evidence
CoRegPairs_sameTimeframe_epiMT_early_controlCisLncs <- filter(CoRegPairs_sameTimeframe_epiMT_early, 
                                                              gsub("\\.[0-9]*", "", AnchorGene) %in% ControlCisLncs$Ens_ID)

#up to 4hrs:
CoRegPairs_sameTimeframe_epiMT_early2 <- filter(trial, AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_2$CAT_geneID, epiMT_DEGs_DE_Upwithin_4$CAT_geneID),
                                               geneID %in% c(epiMT_DEGs_DE_Upwithin_2$CAT_geneID, epiMT_DEGs_DE_Upwithin_4$CAT_geneID))

#nothing extra is obvious

#
#######################
#### previous code ####
#######################
#### group1a same/delayed cclncs - Up/down targets ####

#up/down DELs

#loop across all timeframes:
#backgrounds, all expressed lncs in sequential timeframes, without earlier/later DEGs
backgrounds_list <- list(epiMT_CPM_04, epiMT_CPM_48, epiMT_CPM_824)
#selections, all DELs in sequential timeframes
selection_list <- list(epiMT_DEGs_DE_within_4, epiMT_DEGs_DE_within_8, epiMT_DEGs_DE_within_24)
#hits, target type
hits_list <- list(epiMT_DEGs_DE_within_4, epiMT_DEGs_DE_within_8, epiMT_DEGs_DE_within_24)
#results:
results_list <- list()

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- as.character(unique(unlist(sapply(hits_list[i:length(hits_list)], 
                                                           function(x){x[,1]})
  )))
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g1ai <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_DEL_EL_DEneighbour.csv", row.names = F)


#just up DELs:
selection_list <- list(epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_24)
#remove down DELs from background:
nonSelection_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- as.character(unique(unlist(sapply(hits_list[i:length(hits_list)], 
                                                           function(x){x[,1]})
  )))
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            #!geneID %in% nonSelection_list[[i]]$CAT_geneID
                            )$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            #!geneID %in% nonSelection_list[[i]]$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g1aii <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_inducedDEL_EL_DEneighbour.csv", row.names = F)


#just down DELs:
selection_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- as.character(unique(unlist(sapply(hits_list[i:length(hits_list)], 
                                                           function(x){x[,1]})
  )))
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g1aiii <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_repressedDEL_EL_DEneighbour.csv", row.names = F)


grid.arrange(g1ai, g1aii, g1aiii, ncol = 3,
             left = textGrob("\n% With a candidate\ncis target", gp = gpar(fontface = 'bold', fontsize = 12), rot = 90))

#summary: some sig for early if grouping, some sig for repression if splitting

#### group1b same cclncs - Up/down targets ####

#up/down DELs

#loop across all timeframes:
#backgrounds, all expressed lncs in sequential timeframes, without earlier/later DEGs
backgrounds_list <- list(epiMT_CPM_04, epiMT_CPM_48, epiMT_CPM_824)
#selections, all DELs in sequential timeframes
selection_list <- list(epiMT_DEGs_DE_within_4, epiMT_DEGs_DE_within_8, epiMT_DEGs_DE_within_24)
#hits, target type
hits_list <- list(epiMT_DEGs_DE_within_4, epiMT_DEGs_DE_within_8, epiMT_DEGs_DE_within_24)
#results:
results_list <- list()

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$CAT_geneID)
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g1bi <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_DEL_EL_DEneighboursame.csv", row.names = F)


#just up DELs:
selection_list <- list(epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$CAT_geneID)
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g1bii <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_inducedDEL_EL_DEneighboursame.csv", row.names = F)


#just down DELs:
selection_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$CAT_geneID)
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g1biii <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_repressedDEL_EL_DEneighboursame.csv", row.names = F)


grid.arrange(g1bi, g1bii, g1biii, ncol = 3,
             left = textGrob("\n% With a candidate\ncis target", gp = gpar(fontface = 'bold', fontsize = 12), rot = 90))

#summary: some sig for early if grouping, some sig for induced and repressed if splitting


#### group1c delayed cclncs - Up/down targets ####

#up.down DELs

#just up DELs:
selection_list <- list(epiMT_DEGs_DE_within_4, epiMT_DEGs_DE_within_8, epiMT_DEGs_DE_within_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- as.character(unique(unlist(sapply(hits_list[c(1:length(hits_list)) >i], 
                                                           function(x){x[,1]})
  )))
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g1ci <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_DEL_EL_DEneighbourdelay.csv", row.names = F)


#just up DELs:
selection_list <- list(epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- as.character(unique(unlist(sapply(hits_list[c(1:length(hits_list)) >i], 
                                                           function(x){x[,1]})
  )))
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g1cii <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_inducedDEL_EL_DEneighbourdelay.csv", row.names = F)


#just down DELs:
selection_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- as.character(unique(unlist(sapply(hits_list[c(1:length(hits_list)) >i], 
                                                           function(x){x[,1]})
  )))
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g1ciii <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_repressedDEL_EL_DEneighbourdelay.csv", row.names = F)


#strong for down
grid.arrange(g1ci, g1cii, g1ciii, ncol = 3,
             left = textGrob("\n% With a candidate\ncis target", gp = gpar(fontface = 'bold', fontsize = 12), rot = 90))


#### group2a same/delayed cclncs - concordant targets ####

#just up DELs:
selection_list <- list(epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_24)
#with up targets:
hits_list <- list(epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- as.character(unique(unlist(sapply(hits_list[i:length(hits_list)], 
                                                           function(x){x[,1]})
  )))
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g2ai <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_inducedDEL_EL_DEneighbour_conc.csv", row.names = F)


#just down DELs:
selection_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)
#with down targets
hits_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- as.character(unique(unlist(sapply(hits_list[i:length(hits_list)], 
                                                           function(x){x[,1]})
  )))
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g2aii <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_repressedDEL_EL_DEneighbour_conc.csv", row.names = F)


grid.arrange(g2ai, g2aii, ncol = 3,
             left = textGrob("\n% With a candidate\ncis target", gp = gpar(fontface = 'bold', fontsize = 12), rot = 90))

#strong sig size for induced + repressed

#### group2b same cclncs - concordant targets ####

#just up DELs:
selection_list <- list(epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_24)
#with up targets:
hits_list <- list(epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_24)
#optional: remove down DELs from background:
#nonSelection_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$CAT_geneID)
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], 
                            geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            #!geneID %in% nonSelection_list[[i]]$CAT_geneID
                            )$geneID
                     ))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], 
                            geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            #!geneID %in% nonSelection_list[[i]]$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
                            )$geneID
                     ))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], 
                            geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], 
                            geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g2bi <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
#write.csv(DEL_PCG_type, "epiMT_inducedDEL_EL_DEneighboursame_conc.csv", row.names = F)

p.adjust(DEL_PCG_type$p, method = "bonferroni")


#just down DELs:
selection_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)
#with down targets:
hits_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$CAT_geneID)
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g2bii <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_repressedDEL_EL_DEneighboursame_conc.csv", row.names = F)

p.adjust(DEL_PCG_type$p, method = "bonferroni")


#strong for both up and down
grid.arrange(g2bi, g2bii, ncol = 3,
             left = textGrob("\n% With a candidate\ncis target", gp = gpar(fontface = 'bold', fontsize = 12), rot = 90))


#### group2c delayed cclncs - concordant targets ####

#just up DELs:
selection_list <- list(epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_24)
#with up targets:
hits_list <- list(epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- as.character(unique(unlist(sapply(hits_list[c(1:length(hits_list)) >i], 
                                                           function(x){x[,1]})
  )))
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g2ci <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_inducedDEL_EL_DEneighbourdelay_conc.csv", row.names = F)


#just down DELs:
selection_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)
#with down targets:
hits_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- as.character(unique(unlist(sapply(hits_list[c(1:length(hits_list)) >i], 
                                                           function(x){x[,1]})
  )))
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g2cii <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_repressedDEL_EL_DEneighbourdelay_conc.csv", row.names = F)


#strong for down
grid.arrange(g2ci, g2cii, ncol = 3,
             left = textGrob("\n% With a candidate\ncis target", gp = gpar(fontface = 'bold', fontsize = 12), rot = 90))



#### group3a same/delayed cclncs - disconcordant targets ####

#just up DELs:
selection_list <- list(epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_24)
#with down targets:
hits_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- as.character(unique(unlist(sapply(hits_list[i:length(hits_list)], 
                                                           function(x){x[,1]})
  )))
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g3ai <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_inducedDEL_EL_DEneighbour_disc.csv", row.names = F)


#just down DELs:
selection_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)
#with down targets
hits_list <- list(epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- as.character(unique(unlist(sapply(hits_list[i:length(hits_list)], 
                                                           function(x){x[,1]})
  )))
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g3aii <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_repressedDEL_EL_DEneighbour_disc.csv", row.names = F)

#ns across board
grid.arrange(g3ai, g3aii, ncol = 3,
             left = textGrob("\n% With a candidate\ncis target", gp = gpar(fontface = 'bold', fontsize = 12), rot = 90))


#### group3b same cclncs - disconcordant targets ####

#just up DELs:
selection_list <- list(epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_24)
#with down targets:
hits_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$CAT_geneID)
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g3bi <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_inducedDEL_EL_DEneighboursame_disc.csv", row.names = F)


#just down DELs:
selection_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)
#with up targets:
hits_list <- list(epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$CAT_geneID)
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g3bii <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_repressedDEL_EL_DEneighboursame_disc.csv", row.names = F)


#ns across board
grid.arrange(g3bi, g3bii, ncol = 3,
             left = textGrob("\n% With a candidate\ncis target", gp = gpar(fontface = 'bold', fontsize = 12), rot = 90))


#### group3c delayed cclncs - disconcordant targets ####

#just up DELs:
selection_list <- list(epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_24)
#with down targets:
hits_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- as.character(unique(unlist(sapply(hits_list[c(1:length(hits_list)) >i], 
                                                           function(x){x[,1]})
  )))
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g3ci <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_inducedDEL_EL_DEneighbourdelay_disc.csv", row.names = F)


#just down DELs:
selection_list <- list(epiMT_DEGs_DE_Downwithin_4, epiMT_DEGs_DE_Downwithin_8, epiMT_DEGs_DE_Downwithin_24)
#with up targets:
hits_list <- list(epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Upwithin_8, epiMT_DEGs_DE_Upwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- as.character(unique(unlist(sapply(hits_list[c(1:length(hits_list)) >i], 
                                                           function(x){x[,1]})
  )))
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID)$geneID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% selection_list[[i]]$CAT_geneID
  )$geneID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], geneID %in% filter(GeneBiotypes, grepl("lncRNA", CAT_geneClass))$CAT_geneID,
                            geneID %in% filter(AllLNC_AllPCG_1, AnchorGene %in% selection_list[[i]]$CAT_geneID, geneID %in% cisPotentialTargets)$AnchorGene
  )$geneID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-16hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g3cii <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18)]
write.csv(DEL_PCG_type, "epiMT_repressedDEL_EL_DEneighbourdelay_disc.csv", row.names = F)

#strong sig for down
grid.arrange(g3ci, g3cii, ncol = 3,
             left = textGrob("\n% With a candidate\ncis target", gp = gpar(fontface = 'bold', fontsize = 12), rot = 90))


#### exploring bigger FCs for lnc-targeted genes ####

PCGDE_Upwithin_4 <- filter(epiMT_DEGs_DE_Upwithin_4, grepl("coding_mRNA", CAT_geneClassII))

#now only co-induced pairs
PCGDE_Upwithin_4$DEneighbourTypeIII <- "None"

#identifies all PCG neighbour pairs with a change in 4hrs
CoInducedPCG_4hr <- filter(AllPCG_AllPCG_1, AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_4$CAT_geneID),
                           geneID %in% c(epiMT_DEGs_DE_Upwithin_4$CAT_geneID))

CoInducedLncRNA_4hr <- filter(CoRegPairs_04_48_24_extended_epiMT, AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_4$CAT_geneID), 
                              geneID %in% c(epiMT_DEGs_DE_Upwithin_4$CAT_geneID))

PCGDE_Upwithin_4$DEneighbourTypeIII[PCGDE_Upwithin_4$CAT_geneID %in% 
                                      c(CoInducedPCG_4hr$AnchorGene)] <- "PCG_only"
PCGDE_Upwithin_4$DEneighbourTypeIII[PCGDE_Upwithin_4$CAT_geneID %in% CoInducedLncRNA_4hr$geneID] <- "CCLncRNA"

table(PCGDE_Upwithin_4$DEneighbourTypeIII)#76 PCG coinduced with a CCLnc in 4hrs, 169 are coinduced with a PCG

PCGDE_Upwithin_4$`Target concordant with:` <- PCGDE_Upwithin_4$DEneighbourTypeIII
library(ggbeeswarm)
ggplot(PCGDE_Upwithin_4) + aes(x = `Target concordant with:`, y = abs(LFC_0_4), color = `Target concordant with:`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Log2FC (0-4hrs)")

#enough n for anova first
summary(aov(abs(LFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Upwithin_4))# ***
TukeyHSD(aov(abs(LFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Upwithin_4))# */ns for lncs, no change to PCG

#better results with correlation previously? T2 improved this in SVSMC


#repressed 
PCGDE_Downwithin_4 <- filter(epiMT_DEGs_DE_Downwithin_4, grepl("coding_mRNA", CAT_geneClassII))

#now only co-Repressed pairs
PCGDE_Downwithin_4$DEneighbourTypeIII <- "None"

#identifies all PCG neighbour pairs with a change in 4hrs
CoRepressedPCG_4hr <- filter(AllPCG_AllPCG_1, AnchorGene %in% c(epiMT_DEGs_DE_Downwithin_4$CAT_geneID),
                             geneID %in% c(epiMT_DEGs_DE_Downwithin_4$CAT_geneID))

CoRepressedLncRNA_4hr <- filter(CoRegPairs_04_48_24_extended_epiMT, AnchorGene %in% c(epiMT_DEGs_DE_Downwithin_4$CAT_geneID), 
                                geneID %in% c(epiMT_DEGs_DE_Downwithin_4$CAT_geneID))

PCGDE_Downwithin_4$DEneighbourTypeIII[PCGDE_Downwithin_4$CAT_geneID %in% 
                                        c(CoRepressedPCG_4hr$AnchorGene)] <- "PCG_only"
PCGDE_Downwithin_4$DEneighbourTypeIII[PCGDE_Downwithin_4$CAT_geneID %in% CoRepressedLncRNA_4hr$geneID] <- "CCLncRNA"

table(PCGDE_Downwithin_4$DEneighbourTypeIII)#16 PCG coRepressed with a CCLnc in 4hrs, 94 are coRepressed with a PCG

PCGDE_Downwithin_4$`Target concordant with:` <- PCGDE_Downwithin_4$DEneighbourTypeIII

ggplot(PCGDE_Downwithin_4) + aes(x = `Target concordant with:`, y = abs(LFC_0_4), color = `Target concordant with:`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Absolute Log2FC (0-4hrs)")

#KW low n
kruskal.test(abs(LFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Downwithin_4)# ns


#### exploring bigger FCs for lnc-targeted genes (4-8hr pairs) ####

PCGDE_Upwithin_8 <- filter(epiMT_DEGs_DE_Upwithin_8, grepl("coding_mRNA", CAT_geneClassII))

#now only co-induced pairs
PCGDE_Upwithin_8$DEneighbourTypeIII <- "None"

#identifies all PCG neighbour pairs with a change in 4hrs
CoInducedPCG_8hr <- filter(AllPCG_AllPCG_1, AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_8$CAT_geneID),
                           geneID %in% c(epiMT_DEGs_DE_Upwithin_8$CAT_geneID))

CoInducedLncRNA_8hr <- filter(CoRegPairs_04_48_24_extended_epiMT, AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_8$CAT_geneID), 
                              geneID %in% c(epiMT_DEGs_DE_Upwithin_8$CAT_geneID))

PCGDE_Upwithin_8$DEneighbourTypeIII[PCGDE_Upwithin_8$CAT_geneID %in% 
                                      c(CoInducedPCG_8hr$AnchorGene)] <- "PCG_only"
PCGDE_Upwithin_8$DEneighbourTypeIII[PCGDE_Upwithin_8$CAT_geneID %in% CoInducedLncRNA_8hr$geneID] <- "CCLncRNA"

table(PCGDE_Upwithin_8$DEneighbourTypeIII)#33 PCG coinduced with a CCLnc in 4hrs, 163 are coinduced with a PCG

PCGDE_Upwithin_8$`Target concordant with:` <- PCGDE_Upwithin_8$DEneighbourTypeIII

PCGDE_Upwithin_8_assess <- filter(PCGDE_Upwithin_8, LFC_0_8 >log2(1.5), pa_0_8 < 0.05)

library(ggbeeswarm)
ggplot(PCGDE_Upwithin_8_assess) + aes(x = `Target concordant with:`, y = abs(LFC_0_8), color = `Target concordant with:`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Log2FC (0-8hrs)")

#enough n for anova first
summary(aov(abs(LFC_0_8) ~ DEneighbourTypeIII, data = PCGDE_Upwithin_8_assess))# ns


#repressed 
PCGDE_Downwithin_8 <- filter(epiMT_DEGs_DE_Downwithin_8, grepl("coding_mRNA", CAT_geneClassII))

#now only co-Repressed pairs
PCGDE_Downwithin_8$DEneighbourTypeIII <- "None"

#identifies all PCG neighbour pairs with a change in 4hrs
CoRepressedPCG_8hr <- filter(AllPCG_AllPCG_1, AnchorGene %in% c(epiMT_DEGs_DE_Downwithin_8$CAT_geneID),
                             geneID %in% c(epiMT_DEGs_DE_Downwithin_8$CAT_geneID))

CoRepressedLncRNA_8hr <- filter(CoRegPairs_04_48_24_extended_epiMT, AnchorGene %in% c(epiMT_DEGs_DE_Downwithin_8$CAT_geneID), 
                                geneID %in% c(epiMT_DEGs_DE_Downwithin_8$CAT_geneID))

PCGDE_Downwithin_8$DEneighbourTypeIII[PCGDE_Downwithin_8$CAT_geneID %in% 
                                        c(CoRepressedPCG_8hr$AnchorGene)] <- "PCG_only"
PCGDE_Downwithin_8$DEneighbourTypeIII[PCGDE_Downwithin_8$CAT_geneID %in% CoRepressedLncRNA_8hr$geneID] <- "CCLncRNA"

table(PCGDE_Downwithin_8$DEneighbourTypeIII)#33 PCG coRepressed with a CCLnc in 4hrs, 210 are coRepressed with a PCG

PCGDE_Downwithin_8$`Target concordant with:` <- PCGDE_Downwithin_8$DEneighbourTypeIII

PCGDE_Downwithin_8_assess <- filter(PCGDE_Downwithin_8, LFC_0_8 < -log2(1.5), pa_0_8 < 0.05)

ggplot(PCGDE_Downwithin_8_assess) + aes(x = `Target concordant with:`, y = abs(LFC_0_8), color = `Target concordant with:`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Absolute Log2FC (0-8hrs)")

summary(aov(abs(LFC_0_8) ~ DEneighbourTypeIII, data = PCGDE_Downwithin_8_assess))# ns



#### exploring bigger FCs for lnc-targeted genes (8-24hr pairs) ####

PCGDE_Upwithin_24 <- filter(epiMT_DEGs_DE_Upwithin_24, grepl("coding_mRNA", CAT_geneClassII))

#now only co-induced pairs
PCGDE_Upwithin_24$DEneighbourTypeIII <- "None"

#identifies all PCG neighbour pairs with a change in 4hrs
CoInducedPCG_24hr <- filter(AllPCG_AllPCG_1, AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_24$CAT_geneID),
                            geneID %in% c(epiMT_DEGs_DE_Upwithin_24$CAT_geneID))

CoInducedLncRNA_24hr <- filter(CoRegPairs_04_48_24_extended_epiMT, AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_24$CAT_geneID), 
                               geneID %in% c(epiMT_DEGs_DE_Upwithin_24$CAT_geneID))

PCGDE_Upwithin_24$DEneighbourTypeIII[PCGDE_Upwithin_24$CAT_geneID %in% 
                                       c(CoInducedPCG_24hr$AnchorGene)] <- "PCG_only"
PCGDE_Upwithin_24$DEneighbourTypeIII[PCGDE_Upwithin_24$CAT_geneID %in% CoInducedLncRNA_24hr$geneID] <- "CCLncRNA"

table(PCGDE_Upwithin_24$DEneighbourTypeIII)#42 PCG coinduced with a CCLnc in 4hrs, 173 are coinduced with a PCG

PCGDE_Upwithin_24$`Target concordant with:` <- PCGDE_Upwithin_24$DEneighbourTypeIII

PCGDE_Upwithin_24_assess <- filter(PCGDE_Upwithin_24, LFC_0_24 >log2(1.5), pa_0_24 < 0.05)

library(ggbeeswarm)
ggplot(PCGDE_Upwithin_24_assess) + aes(x = `Target concordant with:`, y = abs(LFC_0_24), color = `Target concordant with:`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Log2FC (0-24hrs)")

#enough n for anova first
summary(aov(abs(LFC_0_24) ~ DEneighbourTypeIII, data = PCGDE_Upwithin_24_assess))# ns


#repressed 
PCGDE_Downwithin_24 <- filter(epiMT_DEGs_DE_Downwithin_24, grepl("coding_mRNA", CAT_geneClassII))

#now only co-Repressed pairs
PCGDE_Downwithin_24$DEneighbourTypeIII <- "None"

#identifies all PCG neighbour pairs with a change in 4hrs
CoRepressedPCG_24hr <- filter(AllPCG_AllPCG_1, AnchorGene %in% c(epiMT_DEGs_DE_Downwithin_24$CAT_geneID),
                              geneID %in% c(epiMT_DEGs_DE_Downwithin_24$CAT_geneID))

CoRepressedLncRNA_24hr <- filter(CoRegPairs_04_48_24_extended_epiMT, AnchorGene %in% c(epiMT_DEGs_DE_Downwithin_24$CAT_geneID), 
                                 geneID %in% c(epiMT_DEGs_DE_Downwithin_24$CAT_geneID))

PCGDE_Downwithin_24$DEneighbourTypeIII[PCGDE_Downwithin_24$CAT_geneID %in% 
                                         c(CoRepressedPCG_24hr$AnchorGene)] <- "PCG_only"
PCGDE_Downwithin_24$DEneighbourTypeIII[PCGDE_Downwithin_24$CAT_geneID %in% CoRepressedLncRNA_24hr$geneID] <- "CCLncRNA"

table(PCGDE_Downwithin_24$DEneighbourTypeIII)#49 PCG coRepressed with a CCLnc in 4hrs, 261 are coRepressed with a PCG

PCGDE_Downwithin_24$`Target concordant with:` <- PCGDE_Downwithin_24$DEneighbourTypeIII

PCGDE_Downwithin_24_assess <- filter(PCGDE_Downwithin_24, LFC_0_24 < -log2(1.5), pa_0_24 < 0.05)

ggplot(PCGDE_Downwithin_24_assess) + aes(x = `Target concordant with:`, y = abs(LFC_0_24), color = `Target concordant with:`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Absolute Log2FC (0-24hrs)")

summary(aov(abs(LFC_0_24) ~ DEneighbourTypeIII, data = PCGDE_Downwithin_24_assess))# ns


#### combined plotting ####

PCGDE_Upwithin_4$time <- "4hr Induced"
PCGDE_Upwithin_8_assess$time <- "8hr Induced"
PCGDE_Upwithin_24_assess$time <- "16hr Induced"

PCGDE_Upwithin_4$LFC_2_compare <- PCGDE_Upwithin_4$LFC_0_4
PCGDE_Upwithin_8_assess$LFC_2_compare <- PCGDE_Upwithin_8_assess$LFC_0_8
PCGDE_Upwithin_24_assess$LFC_2_compare <- PCGDE_Upwithin_24_assess$LFC_0_24

trial <- rbind(PCGDE_Upwithin_4,
               PCGDE_Upwithin_8_assess,
               PCGDE_Upwithin_24_assess)

trial$time <- factor(trial$time)
trial$time <- factor(trial$time, levels = levels(trial$time)[c(2,3,1)])

trial$DEneighbourTypeIII[trial$DEneighbourTypeIII == "None"] <- "No"
trial$DEneighbourTypeIII[trial$DEneighbourTypeIII == "PCG_only"] <- "No"
trial$DEneighbourTypeIII[trial$DEneighbourTypeIII == "CCLncRNA"] <- "Yes"

ggplot(trial) + aes(x = DEneighbourTypeIII, y = LFC_2_compare, color = DEneighbourTypeIII) +
  geom_quasirandom(alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha =0.4) +
  facet_wrap(~time) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "mediumorchid", "No" = "olivedrab")) +
  xlab("\nInduced LncRNA Neighbour")+
  theme_minimal() +
  ggtitle("EpiMT Induced PCGs") +
  theme(#axis.text.x = element_blank(),
    strip.background = element_rect(),
    legend.position = "none") +
  ylab("Log2FC \n(max. time vs. 0hr)")

trial$timechar <- as.character(trial$time)
wilcox.test(filter(trial, time == "0-4hrs Induced",  DEneighbourTypeIII == "Yes")$LFC_0_4,
            filter(trial, time == "0-4hrs Induced",  DEneighbourTypeIII == "No")$LFC_0_4, paired = F)
wilcox.test(filter(trial, time == "4-8hrs Induced",  DEneighbourTypeIII == "Yes")$LFC_0_8,
            filter(trial, time == "4-8hrs Induced",  DEneighbourTypeIII == "No")$LFC_0_8, paired = F)
wilcox.test(filter(trial, time == "8-16hrs Induced",  DEneighbourTypeIII == "Yes")$LFC_0_24,
            filter(trial, time == "8-16hrs Induced",  DEneighbourTypeIII == "No")$LFC_0_24, paired = F)
#sufficient for a t:
table(trial$DEneighbourTypeIII, trial$time)
t.test(filter(trial, time == "4hr Induced",  DEneighbourTypeIII == "Yes")$LFC_0_4,
       filter(trial, time == "4hr Induced",  DEneighbourTypeIII == "No")$LFC_0_4, var.equal = T)
t.test(filter(trial, time == "8hr Induced",  DEneighbourTypeIII == "Yes")$LFC_0_8,
       filter(trial, time == "8hr Induced",  DEneighbourTypeIII == "No")$LFC_0_8, var.equal = T)
t.test(filter(trial, time == "16hr Induced",  DEneighbourTypeIII == "Yes")$LFC_0_24,
       filter(trial, time == "16hr Induced",  DEneighbourTypeIII == "No")$LFC_0_24, var.equal = T)
#ns across board but early closest

#keep distinct y axes
trial$DEneighbourTypeIII[trial$DEneighbourTypeIII == "Yes"] <- "Co-induced\nwith\nLncRNA\nNeighbour"
trial$DEneighbourTypeIII[trial$DEneighbourTypeIII == "No"] <- "Other\nPCGs"

ggplot(filter(trial, time == "4hr Induced")) + aes(x = DEneighbourTypeIII, 
                                                   y = LFC_2_compare ,
                                                   color = DEneighbourTypeIII) +
  geom_quasirandom(alpha = 0.35, size = 4) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha =0.7) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Co-induced\nwith\nLncRNA\nNeighbour" = "mediumorchid", "Other\nPCGs" = "grey60")) +
  xlab("")+
  theme_minimal() +
  theme(text = element_text(size = 28),
        legend.position = "none") +
  ylab("Log2FC 0-4hr")

ggplot(filter(trial, time == "8hr Induced")) + aes(x = DEneighbourTypeIII, 
                                                   y = LFC_2_compare ,
                                                   color = DEneighbourTypeIII) +
  geom_quasirandom(alpha = 0.35, size = 4) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha =0.7) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Co-induced\nwith\nLncRNA\nNeighbour" = "mediumorchid", "Other\nPCGs" = "grey60")) +
  xlab("")+
  theme_minimal() +
  theme(text = element_text(size = 28),
        legend.position = "none") +
  ylab("Log2FC 0-8hr")

ggplot(filter(trial, time == "16hr Induced")) + aes(x = DEneighbourTypeIII, 
                                                    y = LFC_2_compare ,
                                                    color = DEneighbourTypeIII) +
  geom_quasirandom(alpha = 0.35, size = 4) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha =0.7) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Co-induced\nwith\nLncRNA\nNeighbour" = "mediumorchid", "Other\nPCGs" = "grey60")) +
  xlab("")+
  theme_minimal() +
  theme(text = element_text(size = 28),
        legend.position = "none") +
  ylab("Log2FC 0-16hr")


PCGDE_Downwithin_4$time <- "4hr Repressed"
PCGDE_Downwithin_8_assess$time <- "8hr Repressed"
PCGDE_Downwithin_24_assess$time <- "16hr Repressed"

PCGDE_Downwithin_4$abs(LFC_2_compare) <- PCGDE_Downwithin_4$LFC_0_4
PCGDE_Downwithin_8_assess$abs(LFC_2_compare) <- PCGDE_Downwithin_8_assess$LFC_0_8
PCGDE_Downwithin_24_assess$abs(LFC_2_compare) <- PCGDE_Downwithin_24_assess$LFC_0_24

trial <- rbind(PCGDE_Downwithin_4,
               PCGDE_Downwithin_8_assess,
               PCGDE_Downwithin_24_assess)

trial$time <- factor(trial$time)
trial$time <- factor(trial$time, levels = levels(trial$time)[c(2,3,1)])

trial$DEneighbourTypeIII[trial$DEneighbourTypeIII == "None"] <- "No"
trial$DEneighbourTypeIII[trial$DEneighbourTypeIII == "PCG_only"] <- "No"
trial$DEneighbourTypeIII[trial$DEneighbourTypeIII == "CCLncRNA"] <- "Yes"

ggplot(trial) + aes(x = DEneighbourTypeIII, y = abs(abs(LFC_2_compare)), color = DEneighbourTypeIII) +
  geom_quasirandom(alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha =0.4) +
  facet_wrap(~time) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "mediumorchid", "No" = "olivedrab")) +
  xlab("\nRepressed LncRNA Neighbour")+
  theme_minimal() +
  ggtitle("EpiMT Repressed PCGs") +
  theme(#axis.text.x = element_blank(),
    strip.background = element_rect(),
    legend.position = "none") +
  ylab("Absolute LFC \n(max. time vs. 0hr)")

trial$timechar <- as.character(trial$time)
wilcox.test(filter(trial, time == "4hr Repressed",  DEneighbourTypeIII == "Yes")$LFC_0_4,
            filter(trial, time == "4hr Repressed",  DEneighbourTypeIII == "No")$LFC_0_4, paired = F)
wilcox.test(filter(trial, time == "8hr Repressed",  DEneighbourTypeIII == "Yes")$LFC_0_8,
            filter(trial, time == "8hr Repressed",  DEneighbourTypeIII == "No")$LFC_0_8, paired = F)
wilcox.test(filter(trial, time == "16hr Repressed",  DEneighbourTypeIII == "Yes")$LFC_0_24,
            filter(trial, time == "16hr Repressed",  DEneighbourTypeIII == "No")$LFC_0_24, paired = F)
#sufficient for a t:
table(trial$DEneighbourTypeIII, trial$time)
t.test(filter(trial, time == "4hr Repressed",  DEneighbourTypeIII == "Yes")$LFC_0_4,
       filter(trial, time == "4hr Repressed",  DEneighbourTypeIII == "No")$LFC_0_4, var.equal = T)
t.test(filter(trial, time == "8hr Repressed",  DEneighbourTypeIII == "Yes")$LFC_0_8,
       filter(trial, time == "8hr Repressed",  DEneighbourTypeIII == "No")$LFC_0_8, var.equal = T)
t.test(filter(trial, time == "16hr Repressed",  DEneighbourTypeIII == "Yes")$LFC_0_24,
       filter(trial, time == "16hr Repressed",  DEneighbourTypeIII == "No")$LFC_0_24, var.equal = T)
#ns across board but early closest

#keep distinct y axes
trial$DEneighbourTypeIII[trial$DEneighbourTypeIII == "Yes"] <- "Co-Repressed\nwith\nLncRNA\nNeighbour"
trial$DEneighbourTypeIII[trial$DEneighbourTypeIII == "No"] <- "Other\nPCGs"

ggplot(filter(trial, time == "4hr Repressed")) + aes(x = DEneighbourTypeIII, 
                                                     y = abs(LFC_2_compare) ,
                                                     color = DEneighbourTypeIII) +
  geom_quasirandom(alpha = 0.35, size = 4) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha =0.7) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Co-Repressed\nwith\nLncRNA\nNeighbour" = "mediumorchid", "Other\nPCGs" = "grey60")) +
  xlab("")+
  theme_minimal() +
  theme(text = element_text(size = 28),
        legend.position = "none") +
  ylab("Absolute LFC 0-4hr")

ggplot(filter(trial, time == "8hr Repressed")) + aes(x = DEneighbourTypeIII, 
                                                     y = abs(LFC_2_compare) ,
                                                     color = DEneighbourTypeIII) +
  geom_quasirandom(alpha = 0.35, size = 4) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha =0.7) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Co-Repressed\nwith\nLncRNA\nNeighbour" = "mediumorchid", "Other\nPCGs" = "grey60")) +
  xlab("")+
  theme_minimal() +
  theme(text = element_text(size = 28),
        legend.position = "none") +
  ylab("Absolute LFC 0-8hr")

ggplot(filter(trial, time == "16hr Repressed")) + aes(x = DEneighbourTypeIII, 
                                                      y = abs(LFC_2_compare) ,
                                                      color = DEneighbourTypeIII) +
  geom_quasirandom(alpha = 0.35, size = 4) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha =0.7) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Co-Repressed\nwith\nLncRNA\nNeighbour" = "mediumorchid", "Other\nPCGs" = "grey60")) +
  xlab("")+
  theme_minimal() +
  theme(text = element_text(size = 28),
        legend.position = "none") +
  ylab("Absolute LFC 0-16hr")



#### enrichment of target DEGs vs. other DEGs ####
library(clusterProfiler)
library(org.Hs.eg.db)
TF_Lambert2018 <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/TF_Lambert2018.txt", header = T, stringsAsFactors = F)

#as in SMC, lnc targets vs. all DE
CoReg_DE <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended_epiMT$geneID)),
                     universe      = unique(gsub("\\.[0-9]*", "", epiMT_DEGs_DE$CAT_geneID)),
                     keyType       = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "all",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

CoReg_DE_df <- as.data.frame(CoReg_DE)
#no terms returned

#0-4hr targets vs. all DEGs:
CoReg_DE <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", filter(CoRegPairs_04_48_24_extended_epiMT, 
                                                                        AnchorGene %in% epiMT_DEGs_DE_within_4$CAT_geneID)$geneID)),
                     universe      = unique(gsub("\\.[0-9]*", "", epiMT_DEGs_DE$CAT_geneID)),
                     keyType       = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "all",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
CoReg_DE_df <- as.data.frame(CoReg_DE)
#1 term...


#0-4hr targets vs. 0-4hr DEGs:
CoReg_DE <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", filter(CoRegPairs_04_48_24_extended_epiMT, 
                                                                        AnchorGene %in% epiMT_DEGs_DE_within_4$CAT_geneID)$geneID)),
                     universe      = unique(gsub("\\.[0-9]*", "", epiMT_DEGs_DE_within_4$CAT_geneID)),
                     keyType       = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "all",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
CoReg_DE_df <- as.data.frame(CoReg_DE)
#1 term...


#focus on lambert TFs - no "CAT" IDs so remove these
length(unique(CoRegPairs_04_48_24_extended_epiMT$geneID))#532 targets
length(unique(filter(CoRegPairs_04_48_24_extended_epiMT, grepl("ENSG", geneID))$geneID))#528 with ENS IDs
length(unique(filter(CoRegPairs_04_48_24_extended_epiMT, gsub("\\.[0-9]*", "", geneID) %in% TF_Lambert2018$ID)$geneID))

43/528*100 #8.1% of lnc targets

epiMT_PCGs_DE <- filter(epiMT_DEGs_DE, CAT_geneClassII == "coding_mRNA")

length(unique(filter(epiMT_PCGs_DE, grepl("ENSG", CAT_geneID))$CAT_geneID))
length(unique(filter(epiMT_PCGs_DE, gsub("\\.[0-9]*", "", CAT_geneID) %in% TF_Lambert2018$ID)$CAT_geneID))
279/3805*100 #7.3% of all DEGs

a <- 43
b <- 528
c <- 279
d <- 3805

fisher.test(data.frame("lncT" = c(a,b-a),
                       "ElseW" = c(c-a,d-c-(b-a)), row.names = c("TF", "Not")), alternative = "greater")
#ns


#0-4hrs:

epiMT_PCGs_DE_within_4 <- filter(epiMT_DEGs_DE_within_4, CAT_geneClassII == "coding_mRNA")

length(unique(filter(epiMT_PCGs_DE_within_4, grepl("ENSG", CAT_geneID))$CAT_geneID))
length(unique(filter(epiMT_PCGs_DE_within_4, gsub("\\.[0-9]*", "", CAT_geneID) %in% TF_Lambert2018$ID)$CAT_geneID))
113/1081*100 #10.5% of 0-4hr DEGs are TFs

length(unique(filter(CoRegPairs_04_48_24_extended_epiMT, 
                     AnchorGene %in% epiMT_DEGs_DE_within_4$CAT_geneID,
                     geneID %in% epiMT_DEGs_DE_within_4$CAT_geneID,
                     grepl("ENSG", geneID))$geneID))#265 with ENS IDs
length(unique(filter(CoRegPairs_04_48_24_extended_epiMT, 
                     AnchorGene %in% epiMT_DEGs_DE_within_4$CAT_geneID,
                     geneID %in% epiMT_DEGs_DE_within_4$CAT_geneID,
                     gsub("\\.[0-9]*", "", geneID) %in% TF_Lambert2018$ID)$geneID))
14/106*100 #13.2% of 0-4hr CClnc targets are TFs

a <- 14
b <- 106
c <- 113
d <- 1081

fisher.test(data.frame("lncT" = c(a,b-a),
                       "ElseW" = c(c-a,d-c-(b-a)), row.names = c("TF", "Not")), alternative = "greater")
#ns but closeish


#0-4hrs induced cclnc targets:
length(unique(filter(CoRegPairs_04_48_24_extended_epiMT, 
                     AnchorGene %in% epiMT_DEGs_DE_Upwithin_4$CAT_geneID,
                     geneID %in% epiMT_DEGs_DE_within_4$CAT_geneID,
                     grepl("ENSG", geneID))$geneID))#207 with ENS IDs
length(unique(filter(CoRegPairs_04_48_24_extended_epiMT, 
                     AnchorGene %in% epiMT_DEGs_DE_Upwithin_4$CAT_geneID,
                     geneID %in% epiMT_DEGs_DE_within_4$CAT_geneID,
                     gsub("\\.[0-9]*", "", geneID) %in% TF_Lambert2018$ID)$geneID))
9/85*100 #10.6

length(unique(filter(epiMT_PCGs_DE_within_4, grepl("ENSG", CAT_geneID))$CAT_geneID))
length(unique(filter(epiMT_PCGs_DE_within_4, gsub("\\.[0-9]*", "", CAT_geneID) %in% TF_Lambert2018$ID)$CAT_geneID))
113/1081*100 #9.3%

a <- 9
b <- 85
c <- 113
d <- 1081

fisher.test(data.frame("lncT" = c(a,b-a),
                       "ElseW" = c(c-a,d-c-(b-a)), row.names = c("TF", "Not")), alternative = "greater")



#IEGs:
IEGs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Gene lists/arner_2015_table_S5_IEGs_lit.csv")
IEGs_hs <- filter(IEGs, !Hs_symbol %in% c("", "-"))
IEGs_hs <- filter(IEGs_hs, !grepl("[a-z]", IEGs_hs$Total_count))

#remove those without a gene symbol for fairness

CoRegPairs_04_48_24_extended_epiMTi <- merge(CoRegPairs_04_48_24_extended_epiMT, epiMT_DEGs_DE[,1:2], by.x = "geneID", by.y = "CAT_geneID")
CoRegPairs_04_48_24_extended_epiMTi <- filter(CoRegPairs_04_48_24_extended_epiMTi, !grepl("CATG0", CAT_geneName))

epiMT_PCGs_DE_within_4 <- filter(epiMT_DEGs_DE_within_4, CAT_geneClassII == "coding_mRNA")

length(unique(filter(epiMT_PCGs_DE_within_4, !grepl("CATG0", CAT_geneName))$CAT_geneID))
length(unique(filter(epiMT_PCGs_DE_within_4, CAT_geneName %in% IEGs_hs$Hs_symbol)$CAT_geneID))
90/1081*100 #8.3% of 0-4hr DEGs are IEGs


length(unique(filter(CoRegPairs_04_48_24_extended_epiMTi, 
                     AnchorGene %in% epiMT_DEGs_DE_within_4$CAT_geneID,
                     geneID %in% epiMT_DEGs_DE_within_4$CAT_geneID)$geneID))#106 with ENS IDs
length(unique(filter(CoRegPairs_04_48_24_extended_epiMTi, 
                     AnchorGene %in% epiMT_DEGs_DE_within_4$CAT_geneID,
                     geneID %in% epiMT_DEGs_DE_within_4$CAT_geneID,
                     CAT_geneName %in% IEGs_hs$Hs_symbol)$geneID))
12/106*100 #11.3% of 0-4hr CClncRNA targets are IEGs

a <- 12
b <- 106
c <- 90
d <- 1081

fisher.test(data.frame("lncT" = c(a,b-a),
                       "ElseW" = c(c-a,d-c-(b-a)), row.names = c("IEG", "Not")), alternative = "greater")


#0-4hrs induced cclnc targets:
length(unique(filter(CoRegPairs_04_48_24_extended_epiMTi, 
                     AnchorGene %in% epiMT_DEGs_DE_Upwithin_4$CAT_geneID,
                     geneID %in% epiMT_DEGs_DE_within_4$CAT_geneID)$geneID))#265 with ENS IDs
length(unique(filter(CoRegPairs_04_48_24_extended_epiMTi, 
                     AnchorGene %in% epiMT_DEGs_DE_Upwithin_4$CAT_geneID,
                     geneID %in% epiMT_DEGs_DE_within_4$CAT_geneID,
                     CAT_geneName %in% IEGs_hs$Hs_symbol)$geneID))
9/85*100 #10.6%

length(unique(filter(epiMT_PCGs_DE_within_4, !grepl("CATG0", CAT_geneName))$CAT_geneID))
length(unique(filter(epiMT_PCGs_DE_within_4, CAT_geneName %in% IEGs_hs$Hs_symbol)$CAT_geneID))
90/1081*100 #8.3% of 0-4hr DEGs

a <- 9
b <- 85
c <- 90
d <- 1081

fisher.test(data.frame("lncT" = c(a,b-a),
                       "ElseW" = c(c-a,d-c-(b-a)), row.names = c("IEG", "Not")), alternative = "greater")
#ns


#### attempts with more robustly expressed/bigger FC PCG ####

filtered_DEGs <- filter(epiMT_DEGs_DE_Upwithin_4, pa_0_4 <0.01, (Hour4_meanCPM>5 | Hour0_meanCPM >5))

PCGDE_Upwithin_4 <- filter(filtered_DEGs, grepl("coding_mRNA", CAT_geneClassII))

#now only co-induced pairs
PCGDE_Upwithin_4$DEneighbourTypeIII <- "None"

#identifies all PCG neighbour pairs with a change in 4hrs
CoInducedPCG_4hr <- filter(AllPCG_AllPCG_1, AnchorGene %in% c(filtered_DEGs$CAT_geneID),
                           geneID %in% c(filtered_DEGs$CAT_geneID))

CoInducedLncRNA_4hr <- filter(CoRegPairs_04_48_24_extended_epiMT, AnchorGene %in% c(filtered_DEGs$CAT_geneID), 
                              geneID %in% c(filtered_DEGs$CAT_geneID))

PCGDE_Upwithin_4$DEneighbourTypeIII[PCGDE_Upwithin_4$CAT_geneID %in% 
                                      c(CoInducedPCG_4hr$AnchorGene, CoInducedPCG_4hr$geneID)] <- "PCG_only"
PCGDE_Upwithin_4$DEneighbourTypeIII[PCGDE_Upwithin_4$CAT_geneID %in% CoInducedLncRNA_4hr$geneID] <- "CCLncRNA"

table(PCGDE_Upwithin_4$DEneighbourTypeIII)#76 PCG coinduced with a CCLnc in 4hrs, 169 are coinduced with a PCG

PCGDE_Upwithin_4$`Target concordant with:` <- PCGDE_Upwithin_4$DEneighbourTypeIII
library(ggbeeswarm)
ggplot(PCGDE_Upwithin_4) + aes(x = `Target concordant with:`, y = abs(LFC_0_4), color = `Target concordant with:`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Log2FC (0-4hrs)")

#anova first
summary(aov(abs(LFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Upwithin_4))# ***
summary(aov(abs(2^LFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Upwithin_4))# *, but non-normalised
kruskal.test(abs(LFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Upwithin_4)#** but close

pairwise.t.test(PCGDE_Upwithin_4$LFC_0_4, PCGDE_Upwithin_4$`Target concordant with:`, p.adjust.method = "bonferroni")
#Dunnets better for a multiple to control comparison? but worth knowing about the pcg-none comparison too...

#some signs...

#repressed 
filtered_DEGs <- filter(epiMT_DEGs_DE_Downwithin_4, pa_0_4 <0.01, (Hour4_meanCPM>5 | Hour0_meanCPM >5))

PCGDE_Downwithin_4 <- filter(filtered_DEGs, grepl("coding_mRNA", CAT_geneClassII))

#now only co-Repressed pairs
PCGDE_Downwithin_4$DEneighbourTypeIII <- "None"

#identifies all PCG neighbour pairs with a change in 4hrs
CoRepressedPCG_4hr <- filter(AllPCG_AllPCG_1, AnchorGene %in% c(filtered_DEGs$CAT_geneID),
                             geneID %in% c(filtered_DEGs$CAT_geneID))

CoRepressedLncRNA_4hr <- filter(CoRegPairs_04_48_24_extended_epiMT, AnchorGene %in% c(filtered_DEGs$CAT_geneID), 
                                geneID %in% c(filtered_DEGs$CAT_geneID))

PCGDE_Downwithin_4$DEneighbourTypeIII[PCGDE_Downwithin_4$CAT_geneID %in% 
                                        c(CoRepressedPCG_4hr$AnchorGene, CoRepressedPCG_4hr$geneID)] <- "PCG_only"
PCGDE_Downwithin_4$DEneighbourTypeIII[PCGDE_Downwithin_4$CAT_geneID %in% CoRepressedLncRNA_4hr$geneID] <- "CCLncRNA"

table(PCGDE_Downwithin_4$DEneighbourTypeIII)#6 PCG coRepressed with a CCLnc in 4hrs, 41 are coRepressed with a PCG

PCGDE_Downwithin_4$`Target concordant with:` <- PCGDE_Downwithin_4$DEneighbourTypeIII

ggplot(PCGDE_Downwithin_4) + aes(x = `Target concordant with:`, y = abs(LFC_0_4), color = `Target concordant with:`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Absolute Log2FC (0-4hrs)")

#anova first
summary(aov(abs(LFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Downwithin_4))# ns
summary(aov(abs(2^LFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Downwithin_4))# ns
kruskal.test(abs(LFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Downwithin_4)# ns


#### attempts with correlation ####

#Correlation:
#Spearman's on FPKM
CoRegPairs_04_48_24_epiMT_Corr <- merge(CoRegPairs_04_48_24_extended_epiMT[,c(1:2)], epiMT_CPM[,c(1:13)], by.x = "AnchorGene", by.y = "geneID")
CoRegPairs_04_48_24_epiMT_Corr <- merge(CoRegPairs_04_48_24_epiMT_Corr, epiMT_CPM[,c(1:13)], by.x = "geneID", by.y = "geneID")
CoRegPairs_04_48_24_epiMT_Corr <- unique(CoRegPairs_04_48_24_epiMT_Corr)
colnames(CoRegPairs_04_48_24_epiMT_Corr)
CoRegPairs_04_48_24_epiMT_Corr$pairs <- paste(CoRegPairs_04_48_24_epiMT_Corr$AnchorGene, CoRegPairs_04_48_24_epiMT_Corr$geneID, sep ="-")

trial <- split(CoRegPairs_04_48_24_epiMT_Corr, CoRegPairs_04_48_24_epiMT_Corr$pairs)
triali <- sapply(trial, function(x){
  cor.test(as.numeric(x[3:14]), as.numeric(x[15:26]), method = "spearman")$estimate
})
trialii <- sapply(trial, function(x){
  cor.test(as.numeric(x[3:14]), as.numeric(x[15:26]), method = "spearman")$p.value
})

trialiii <- data.frame("pairs" = names(trial), "spear_rho" = triali, "spear_p" = trialii, "spear_p_adj" = p.adjust(trialii, method = "BH"))

CoRegPairs_04_48_24_epiMT_Corr <- merge(trialiii, CoRegPairs_04_48_24_epiMT_Corr, by = "pairs")
dim(filter(CoRegPairs_04_48_24_epiMT_Corr, spear_p_adj < 0.05))#280, (used to be 90) from Q3 2024updates

CoRegPairs_04_48_24_epiMT_Corr$corSig <- "No"
CoRegPairs_04_48_24_epiMT_Corr$corSig[abs(CoRegPairs_04_48_24_epiMT_Corr$spear_rho) >0.2 & 
                                        CoRegPairs_04_48_24_epiMT_Corr$spear_p_adj <0.05] <- "Yes"
table(CoRegPairs_04_48_24_epiMT_Corr$corSig)#still 86 (with or without the rho filter)
CoRegPairs_04_48_24_epiMT_Corr <- filter(CoRegPairs_04_48_24_epiMT_Corr, corSig == "Yes")


#PCG-PCG:
CoRegPairs_04_48_24_extended_epiMT_PCG <- filter(AllPCG_AllPCG_1, 
                                                 (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_4$CAT_geneID, epiMT_DEGs_DE_Downwithin_4$CAT_geneID) & 
                                                    geneID %in% epiMT_DEGs_DE$CAT_geneID
                                                 ) |
                                                   (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_8$CAT_geneID, epiMT_DEGs_DE_Downwithin_8$CAT_geneID) & 
                                                      geneID %in% c(epiMT_DEGs_DE_Upwithin_8$CAT_geneID, epiMT_DEGs_DE_Downwithin_8$CAT_geneID
                                                                    , epiMT_DEGs_DE_Upwithin_24$CAT_geneID, epiMT_DEGs_DE_Downwithin_24$CAT_geneID
                                                      )) |
                                                   (AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_24$CAT_geneID, epiMT_DEGs_DE_Downwithin_24$CAT_geneID) & 
                                                      geneID %in% c(epiMT_DEGs_DE_Upwithin_24$CAT_geneID, epiMT_DEGs_DE_Downwithin_24$CAT_geneID
                                                      ))
)

CoRegPCGs_epiMT_Corr <- merge(CoRegPairs_04_48_24_extended_epiMT_PCG[,c(1:2)], epiMT_CPM[,c(1:13)], by.x = "AnchorGene", by.y = "geneID")
CoRegPCGs_epiMT_Corr <- merge(CoRegPCGs_epiMT_Corr, epiMT_CPM[,c(1:13)], by.x = "geneID", by.y = "geneID")
CoRegPCGs_epiMT_Corr <- unique(CoRegPCGs_epiMT_Corr)
colnames(CoRegPCGs_epiMT_Corr)
CoRegPCGs_epiMT_Corr$pairs <- paste(CoRegPCGs_epiMT_Corr$AnchorGene, CoRegPCGs_epiMT_Corr$geneID, sep ="-")

trial <- split(CoRegPCGs_epiMT_Corr, CoRegPCGs_epiMT_Corr$pairs)
triali <- sapply(trial, function(x){
  cor.test(as.numeric(x[3:14]), as.numeric(x[15:26]), method = "spearman")$estimate
})
trialii <- sapply(trial, function(x){
  cor.test(as.numeric(x[3:14]), as.numeric(x[15:26]), method = "spearman")$p.value
})

trialiii <- data.frame("pairs" = names(trial), "spear_rho" = triali, "spear_p" = trialii, "spear_p_adj" = p.adjust(trialii, method = "BH"))

CoRegPCGs_epiMT_Corr <- merge(trialiii, CoRegPCGs_epiMT_Corr, by = "pairs")
dim(filter(CoRegPCGs_epiMT_Corr, spear_p_adj < 0.05))#2269

CoRegPCGs_epiMT_Corr$corSig <- "No"
CoRegPCGs_epiMT_Corr$corSig[abs(CoRegPCGs_epiMT_Corr$spear_rho) >0.2 & 
                              CoRegPCGs_epiMT_Corr$spear_p_adj <0.05] <- "Yes"
table(CoRegPCGs_epiMT_Corr$corSig)
CoRegPCGs_epiMT_Corr <- filter(CoRegPCGs_epiMT_Corr, corSig == "Yes")


#
#### exploring bigger FCs for correlated cis-targeted genes ####

PCGDE_Upwithin_4 <- filter(epiMT_DEGs_DE_Upwithin_4, grepl("coding_mRNA", CAT_geneClassII))

#now only co-induced pairs
PCGDE_Upwithin_4$DEneighbourTypeIII <- "None"

#identifies all correlated PCG neighbour targets with a change in 4hrs
CoInducedPCG_4hr <- filter(CoRegPCGs_epiMT_Corr, AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_4$CAT_geneID),
                           geneID %in% c(epiMT_DEGs_DE_Upwithin_4$CAT_geneID))

CoInducedLncRNA_4hr <- filter(CoRegPairs_04_48_24_epiMT_Corr, AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_4$CAT_geneID), 
                              geneID %in% c(epiMT_DEGs_DE_Upwithin_4$CAT_geneID))

PCGDE_Upwithin_4$DEneighbourTypeIII[PCGDE_Upwithin_4$CAT_geneID %in% 
                                      c(CoInducedPCG_4hr$AnchorGene, CoInducedPCG_4hr$geneID)] <- "PCG_only"
PCGDE_Upwithin_4$DEneighbourTypeIII[PCGDE_Upwithin_4$CAT_geneID %in% CoInducedLncRNA_4hr$geneID] <- "CCLncRNA"

table(PCGDE_Upwithin_4$DEneighbourTypeIII)#41 PCG coinduced with a CCLnc in 4hrs, 89 are coinduced with a PCG

PCGDE_Upwithin_4$`Target concordant with:` <- PCGDE_Upwithin_4$DEneighbourTypeIII
library(ggbeeswarm)
ggplot(PCGDE_Upwithin_4) + aes(x = `Target concordant with:`, y = abs(LFC_0_4), color = `Target concordant with:`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Log2FC (0-4hrs)")

#anova first
summary(aov(abs(LFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Upwithin_4))# ***
TukeyHSD(aov(abs(LFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Upwithin_4))# ***

#better than none but not PCG 


#repressed 
PCGDE_Downwithin_4 <- filter(epiMT_DEGs_DE_Downwithin_4, grepl("coding_mRNA", CAT_geneClassII))

#now only co-Repressed pairs
PCGDE_Downwithin_4$DEneighbourTypeIII <- "None"

#identifies all PCG neighbour pairs with a change in 4hrs
CoRepressedPCG_4hr <- filter(CoRegPCGs_epiMT_Corr, AnchorGene %in% c(epiMT_DEGs_DE_Downwithin_4$CAT_geneID),
                             geneID %in% c(epiMT_DEGs_DE_Downwithin_4$CAT_geneID))

CoRepressedLncRNA_4hr <- filter(CoRegPairs_04_48_24_epiMT_Corr, AnchorGene %in% c(epiMT_DEGs_DE_Downwithin_4$CAT_geneID), 
                                geneID %in% c(epiMT_DEGs_DE_Downwithin_4$CAT_geneID))

PCGDE_Downwithin_4$DEneighbourTypeIII[PCGDE_Downwithin_4$CAT_geneID %in% 
                                        c(CoRepressedPCG_4hr$AnchorGene, CoRepressedPCG_4hr$geneID)] <- "PCG_only"
PCGDE_Downwithin_4$DEneighbourTypeIII[PCGDE_Downwithin_4$CAT_geneID %in% CoRepressedLncRNA_4hr$geneID] <- "CCLncRNA"

table(PCGDE_Downwithin_4$DEneighbourTypeIII)#10 PCG coRepressed with a CCLnc in 4hrs, 54 are coRepressed with a PCG

PCGDE_Downwithin_4$`Target concordant with:` <- PCGDE_Downwithin_4$DEneighbourTypeIII

ggplot(PCGDE_Downwithin_4) + aes(x = `Target concordant with:`, y = abs(LFC_0_4), color = `Target concordant with:`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Absolute Log2FC (0-4hrs)")

#KW, close
kruskal.test(abs(LFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Downwithin_4)# 0.06
FSA::dunnTest(PCGDE_Downwithin_4$LFC_0_4, PCGDE_Downwithin_4$`Target concordant with:`, method = "bonferroni")


#previously combined both to get sig? not a good approach but curious where the test has changed since last time
#repressed 
epiMT_DEGs_DE_within_4 <- rbind(epiMT_DEGs_DE_Upwithin_4, epiMT_DEGs_DE_Downwithin_4)

PCGDE_within_4 <- filter(epiMT_DEGs_DE_within_4, grepl("coding_mRNA", CAT_geneClassII))

#now only co-Repressed pairs
PCGDE_within_4$DEneighbourTypeIII <- "None"

#identifies all PCG neighbour pairs with a change in 4hrs
PCGDE_within_4$DEneighbourTypeIII[PCGDE_within_4$CAT_geneID %in% 
                                        c(CoInducedPCG_4hr$AnchorGene, CoInducedPCG_4hr$geneID,
                                          CoRepressedPCG_4hr$AnchorGene, CoRepressedPCG_4hr$geneID)] <- "PCG_only"
PCGDE_within_4$DEneighbourTypeIII[PCGDE_within_4$CAT_geneID %in% c(CoInducedLncRNA_4hr$geneID, 
                                                                       CoRepressedLncRNA_4hr$geneID)] <- "CCLncRNA"

table(PCGDE_within_4$DEneighbourTypeIII)

PCGDE_within_4$`Target concordant with:` <- PCGDE_within_4$DEneighbourTypeIII

ggplot(PCGDE_within_4) + aes(x = `Target concordant with:`, y = abs(LFC_0_4), color = `Target concordant with:`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Absolute Log2FC (0-4hrs)")

#anova first
summary(aov(abs(LFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_within_4))# ***
TukeyHSD(aov(abs(LFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_within_4))# still no sig vs. PCG

#KW?
kruskal.test(abs(LFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_within_4)#***
FSA::dunnTest(PCGDE_within_4$LFC_0_4, PCGDE_within_4$`Target concordant with:`, method = "bh") # close if using BH


#previous test worked due to a) old method b) incorrect approach of combining induced/repressed genes



#### exploring bigger FCs for lnc-targeted genes (prior) ####

epiMT_DEGs_DE_Upwithin_4_PCG <- filter(epiMT_DEGs_DE_Upwithin_4, 
                                       geneID %in% filter(epiMT_DEGs_DE, CAT_geneClassII == "coding_mRNA")$CAT_geneID)

epiMT_DEGs_DE_Upwithin_4_PCG$Partner <- "Not_targeted"

#targets of 4hr induced PCGs
CisCandidateLoci_epiMT_PCG4 <- filter(CisCandidateLoci_epiMT_PCG, AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_4$geneID, 
                                                                                    epiMT_DEGs_DE_Downwithin_4$geneID))

epiMT_DEGs_DE_Upwithin_4_PCG$Partner[epiMT_DEGs_DE_Upwithin_4_PCG$geneID %in%
                                       filter(CisCandidateLoci_epiMT_PCG4, 
                                              AnchorGene %in% epiMT_DEGs_DE_Upwithin_4$geneID)$gene] <- "PCG_coinduced"

#targets of 4hr induced lncRNAs
epiMT_DEGs_DE_Upwithin_4_PCG$Partner[epiMT_DEGs_DE_Upwithin_4_PCG$geneID %in%
                                       filter(CisCandidateLoci_epiMT_4, 
                                              AnchorGene %in% epiMT_DEGs_DE_Upwithin_4$geneID)$gene] <- "Other_lncRNA_coinduced"

#targets of 4hr induced elncRNAs
epiMT_DEGs_DE_Upwithin_4_PCG$Partner[epiMT_DEGs_DE_Upwithin_4_PCG$geneID %in%
                                       filter(CisCandidateLoci_epiMT_4, 
                                              AnchorGene %in% epiMT_DEGs_DE_Upwithin_4$geneID,
                                              AnchorGene %in% filter(GeneBiotypes, CAT_geneCategory == "e_lncRNA")$CAT_geneID)$gene] <- "ElncRNA_coinduced"
table(epiMT_DEGs_DE_Upwithin_4_PCG$Partner)

library(ggbeeswarm)
ggplot(epiMT_DEGs_DE_Upwithin_4_PCG) + aes(x = Partner, y = LFC_0_4, color = Partner) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Absolute Log2FC")

kruskal.test(LFC_0_4 ~ Partner, data = epiMT_DEGs_DE_Upwithin_4_PCG)

pairwise.wilcox.test(epiMT_DEGs_DE_Upwithin_4_PCG$LFC_0_4, epiMT_DEGs_DE_Upwithin_4_PCG$Partner, p.adjust.method = "BH")



epiMT_DEGs_DE_Downwithin_4_PCG <- filter(epiMT_DEGs_DE_Downwithin_4, 
                                       geneID %in% filter(epiMT_DEGs_DE, CAT_geneClassII == "coding_mRNA")$CAT_geneID)

epiMT_DEGs_DE_Downwithin_4_PCG$Partner <- "Not_targeted"

#targets of 4hr repressed PCGs
CisCandidateLoci_epiMT_PCG4 <- filter(CisCandidateLoci_epiMT_PCG, AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_4$geneID, 
                                                                                    epiMT_DEGs_DE_Downwithin_4$geneID))

epiMT_DEGs_DE_Downwithin_4_PCG$Partner[epiMT_DEGs_DE_Downwithin_4_PCG$geneID %in%
                                       filter(CisCandidateLoci_epiMT_PCG4, 
                                              AnchorGene %in% epiMT_DEGs_DE_Downwithin_4$geneID)$gene] <- "PCG_corepressed"

epiMT_DEGs_DE_Downwithin_4_PCG$Partner[epiMT_DEGs_DE_Downwithin_4_PCG$geneID %in%
                                       filter(CisCandidateLoci_epiMT_4, 
                                              AnchorGene %in% epiMT_DEGs_DE_Downwithin_4$geneID)$gene] <- "Other_lncRNA_corepressed"

epiMT_DEGs_DE_Downwithin_4_PCG$Partner[epiMT_DEGs_DE_Downwithin_4_PCG$geneID %in%
                                       filter(CisCandidateLoci_epiMT_4, 
                                              AnchorGene %in% epiMT_DEGs_DE_Downwithin_4$geneID,
                                              AnchorGene %in% filter(GeneBiotypes, CAT_geneCategory == "e_lncRNA")$CAT_geneID)$gene] <- "ElncRNA_corepressed"
table(epiMT_DEGs_DE_Downwithin_4_PCG$Partner)

ggplot(epiMT_DEGs_DE_Downwithin_4_PCG) + aes(x = Partner, y = abs(LFC_0_4), color = Partner) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Absolute Log2FC")

kruskal.test(LFC_0_4 ~ Partner, data = epiMT_DEGs_DE_Downwithin_4_PCG)

pairwise.wilcox.test(epiMT_DEGs_DE_Downwithin_4_PCG$LFC_0_4, epiMT_DEGs_DE_Downwithin_4_PCG$Partner, p.adjust.method = "BH")

#elncRNAs and PCG targeted genes elevated FCs

#lncRNA-targeted genes are elevated vs. untargted genes but not PCG-targeted genes

#no sign of elevated FCs at lncRNA loci


#could combine up and down genes?
epiMT_DEGs_DE_4_PCG <- filter(rbind(epiMT_DEGs_DE_Upwithin_4, 
                                    epiMT_DEGs_DE_Downwithin_4), 
                              geneID %in% filter(epiMT_DEGs_DE, CAT_geneClassII == "coding_mRNA")$CAT_geneID)
epiMT_DEGs_DE_4_PCG$Partner <- "Not_targeted"

#targets of 4hr repressed PCGs
CisCandidateLoci_epiMT_PCG4 <- filter(CisCandidateLoci_epiMT_PCG, AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_4$geneID, 
                                                                                    epiMT_DEGs_DE_Downwithin_4$geneID))

epiMT_DEGs_DE_4_PCG$Partner[epiMT_DEGs_DE_4_PCG$geneID %in%
                                         CisCandidateLoci_epiMT_PCG4$gene] <- "PCG_targeted"

CisCandidateLoci_epiMT_4 <- filter(CisCandidateLoci_epiMT, AnchorGene %in% c(epiMT_DEGs_DE_Upwithin_4$geneID, 
                                                                             epiMT_DEGs_DE_Downwithin_4$geneID))

epiMT_DEGs_DE_4_PCG$Partner[epiMT_DEGs_DE_4_PCG$geneID %in%
                                         CisCandidateLoci_epiMT_4$gene] <- "LncRNA_targeted"

#epiMT_DEGs_DE_4_PCG$Partner[epiMT_DEGs_DE_4_PCG$geneID %in%
#                                         filter(CisCandidateLoci_epiMT_4,
#                                                AnchorGene %in% filter(GeneBiotypes, CAT_geneCategory == "e_lncRNA")$CAT_geneID)$gene] <- "ElncRNA_targeted"
table(epiMT_DEGs_DE_4_PCG$Partner)

ggplot(epiMT_DEGs_DE_4_PCG) + aes(x = Partner, y = abs(LFC_0_4), color = Partner) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Absolute Log2FC")

kruskal.test(abs(LFC_0_4) ~ Partner, data = epiMT_DEGs_DE_4_PCG)
pairwise.wilcox.test(abs(epiMT_DEGs_DE_4_PCG$LFC_0_4), epiMT_DEGs_DE_4_PCG$Partner, p.adjust.method = "BH")
#lncRNA-targeted genes have elevated FCs compared to PCG-targeted or not-targeted, elncRNAs in particular
#stats look better if not splitting enhancer/non-enhancer - t test does not like the lnc-PCG comparison:
pairwise.t.test(abs(epiMT_DEGs_DE_4_PCG$LFC_0_4), epiMT_DEGs_DE_4_PCG$Partner, p.adjust.method = "bonferroni")
#seems overly stringent...
t.test(abs(filter(epiMT_DEGs_DE_4_PCG, Partner == "LncRNA_targeted")$LFC_0_4), 
       abs(filter(epiMT_DEGs_DE_4_PCG, Partner == "PCG_targeted")$LFC_0_4))

#### missing pairs of interest ####

#expecting more RP11 and HMGA2? ENSG00000256268.1 ENSG00000149948.9
#NR2F2 ENSG00000185551.8 and AS ENSG00000247809.3
exampleCounts <- plotCounts(dds, "ENSG00000247809.3", intgroup = c("rep", "time"), returnData = T)

ggplot(exampleCounts,
             aes(x = time, y = count, color = rep, group = rep)) + 
  geom_jitter(width=0.05) + stat_summary(fun=mean, geom="line") #+
  #scale_y_log10()

dds_assay <- assay(rld)
exampleCounts$rld_counts <- dds_assay["ENSG00000247809.3",]

ggplot(exampleCounts,
       aes(x = time, y = rld_counts, color = rep, group = rep)) + 
  geom_jitter(width=0.05) + stat_summary(fun=mean, geom="line") #+
#scale_y_log10()

DEGs_4_0["ENSG00000247809.3",]

#DESeq2 is filtering ENSG00000256268.1 for some reason, it does pass edgeR in their original analysis though...