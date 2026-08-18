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
Supp18_aoSMCil1b <- filter(Supp18_ER, grepl("IL1b", series_name))

#include some finer increments, 2hrs, can subset later if needed
aoSMCil1b_counts <- ER_counts[,c("geneID", unique(unlist(strsplit(Supp18_aoSMCil1b$ref_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_aoSMCil1b, qry_replicate_ID == "01hr")$qry_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_aoSMCil1b, qry_replicate_ID == "03hr")$qry_sample_ID, ","))),
                             unique(unlist(strsplit(filter(Supp18_aoSMCil1b, qry_replicate_ID == "05hr")$qry_sample_ID, ",")))
                             )]
dim(aoSMCil1b_counts)
#n.b. 6 per timepoint...timepoints have been merged, v unhelpful...

#possiblydue to QC, beware
aoSMCil1b_meta <- filter(Timecourse_FANTOM_samples, CAGE_lib_ID %in% colnames(aoSMCil1b_counts))
colnames(aoSMCil1b_counts)[-1] == aoSMCil1b_meta$CAGE_lib_ID
colnames(aoSMCil1b_counts)[-1] %in% aoSMCil1b_meta$CAGE_lib_ID

#sort the ordering, donor and rep
aoSMCil1b_meta <- aoSMCil1b_meta[order(aoSMCil1b_meta$time, aoSMCil1b_meta$rep),]
aoSMCil1b_counts <- aoSMCil1b_counts[,c(1, match(aoSMCil1b_meta$CAGE_lib_ID, colnames(aoSMCil1b_counts)[-1])+1)]
colnames(aoSMCil1b_counts)[-1] == aoSMCil1b_meta$CAGE_lib_ID

#create sample name
aoSMCil1b_meta$time <- gsub(" ","", aoSMCil1b_meta$time)
aoSMCil1b_meta$sample_nameII <- paste(aoSMCil1b_meta$time, aoSMCil1b_meta$rep, sep = "-")
aoSMCil1b_meta$sample_nameII <- gsub("00min-biol", "", aoSMCil1b_meta$sample_nameII)
rownames(aoSMCil1b_meta) <- aoSMCil1b_meta$sample_nameII
colnames(aoSMCil1b_counts)[2:(1+length(rownames(aoSMCil1b_meta)))] <- aoSMCil1b_meta$sample_nameII

aoSMCil1b_counts <- as.matrix(aoSMCil1b_counts[,2:(1+length(rownames(aoSMCil1b_meta)))])
rownames(aoSMCil1b_counts) <- genes$geneID

#essential basic QC step not to be overlooked, library size:
aoSMCil1b_libSizes <- data.frame("sample" = colnames(aoSMCil1b_counts), "libSize" = colSums(aoSMCil1b_counts))
ggplot(aoSMCil1b_libSizes) + aes(x = sample, y = libSize) +
  geom_bar(stat = "identity") + Seurat::RotatedAxis()
#some low, none terrible

summary(aoSMCil1b_libSizes$libSize/1000000)
#median 1.2mill counts

#
#### for aoSMCil1b, find DEGs ####

#obtain table of 0,4,8,24:
Supp18_aoSMCil1b <- filter(Supp18_ER, grepl("IL1b", series_name))

#ideally should obtain table of 0,4,8,24 BUT

#depth issue in 24hrs, and 1 of 16hr
#there are only 2 12hr samples
aoSMCil1b_counts <- ER_counts[,c("geneID", unique(unlist(strsplit(Supp18_aoSMCil1b$ref_sample_ID, ","))),
                                 unique(unlist(strsplit(filter(Supp18_aoSMCil1b, qry_replicate_ID == "01hr")$qry_sample_ID, ","))),
                                 unique(unlist(strsplit(filter(Supp18_aoSMCil1b, qry_replicate_ID == "03hr")$qry_sample_ID, ","))),
                                 unique(unlist(strsplit(filter(Supp18_aoSMCil1b, qry_replicate_ID == "05hr")$qry_sample_ID, ",")))
                                 )]

aoSMCil1b_meta <- filter(Timecourse_FANTOM_samples, CAGE_lib_ID %in% colnames(aoSMCil1b_counts))
#sort the ordering, donor and rep
aoSMCil1b_meta <- aoSMCil1b_meta[order(aoSMCil1b_meta$time, aoSMCil1b_meta$rep),]
aoSMCil1b_counts <- aoSMCil1b_counts[,c(1, match(aoSMCil1b_meta$CAGE_lib_ID, colnames(aoSMCil1b_counts)[-1])+1)]
colnames(aoSMCil1b_counts)[-1] == aoSMCil1b_meta$CAGE_lib_ID
#create sample name
aoSMCil1b_meta$time <- gsub(" ","", aoSMCil1b_meta$time)
aoSMCil1b_meta$sample_nameII <- paste(aoSMCil1b_meta$time, aoSMCil1b_meta$rep, sep = "-")
aoSMCil1b_meta$sample_nameII <- gsub("00min-biol", "", aoSMCil1b_meta$sample_nameII)
rownames(aoSMCil1b_meta) <- aoSMCil1b_meta$sample_nameII
colnames(aoSMCil1b_counts)[2:(1+length(rownames(aoSMCil1b_meta)))] <- aoSMCil1b_meta$sample_nameII

aoSMCil1b_counts <- as.matrix(aoSMCil1b_counts[,2:(1+length(rownames(aoSMCil1b_meta)))])
rownames(aoSMCil1b_counts) <- genes$geneID

#cut to 0,2,4,6
colnames(aoSMCil1b_counts)
colnames(aoSMCil1b_counts)[-c(4:9,13:15,19:21)]
rownames(aoSMCil1b_meta)[-c(4:9,13:15,19:21)]

aoSMCil1b_counts <- aoSMCil1b_counts[,-c(4:9,13:15,19:21)]
aoSMCil1b_meta <- aoSMCil1b_meta[-c(4:9,13:15,19:21),]

library(DESeq2)
dds <- DESeqDataSetFromMatrix(aoSMCil1b_counts[,], 
                              aoSMCil1b_meta[,], design=~rep#+batch
                              +time)
dds <- dds[rowSums(counts(dds))>10, ]
dds <- DESeq(dds)

rld <- rlog(dds, blind=FALSE)

pcaData <- plotPCA(rld, intgroup = c("rep", "time"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#colnames(pcaData)[c(4,5)] <- c("rep", "time")

ggplot(pcaData, aes(x = PC1, y = PC2, color = time, shape = rep)) +
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
#pretty convincing gradual tx'omic shift + donor pairing (bit of a pull back after 4hr)

#can cut to 0, 2, 4, 6 for simplicity (done in lines above)

#with this many timepoints, stick to a 0hr comparison
DEGs_2_0 <- results(dds, contrast=c("time","02hr","00hr00min"))
DEGs_4_0 <- results(dds, contrast=c("time","04hr","00hr00min"))
DEGs_6_0 <- results(dds, contrast=c("time","06hr","00hr00min"))

#identify DEGs WT v mut in either treatment:
aoSMCil1b_DEGs <- data.frame(geneID = rownames(DEGs_4_0), 
                         LFC_0_2 = DEGs_2_0$log2FoldChange, pa_0_2 = DEGs_2_0$padj,
                         LFC_0_4 = DEGs_4_0$log2FoldChange, pa_0_4 = DEGs_4_0$padj,
                         LFC_0_6 = DEGs_6_0$log2FoldChange, pa_0_6 = DEGs_6_0$padj)

#get a CPM table alongside:
trial <- lapply(as.data.frame(aoSMCil1b_counts), function(x){
  x/sum(x)*1000000
})
aoSMCil1b_CPM <- bind_cols(trial)
colnames(aoSMCil1b_CPM)
#careful with col indices, not all n=3:
aoSMCil1b_CPM$Hour0_meanCPM <- rowMeans(aoSMCil1b_CPM[,1:3])
aoSMCil1b_CPM$Hour2_meanCPM <- rowMeans(aoSMCil1b_CPM[,4:6])
aoSMCil1b_CPM$Hour4_meanCPM <- rowMeans(aoSMCil1b_CPM[,7:9])
aoSMCil1b_CPM$Hour6_meanCPM <- rowMeans(aoSMCil1b_CPM[,10:12])

aoSMCil1b_CPM$geneID <- genes$geneID
colnames(aoSMCil1b_CPM)
aoSMCil1b_DEGs_CPM <- merge(aoSMCil1b_DEGs, aoSMCil1b_CPM[,], by = "geneID")

#count how many samples reach >1 (too far, no diff to DEGs)
#aoSMCil1b_CPM$Hour0_CPM1_no <-  apply(aoSMCil1b_CPM[,1:3], 1, function(x){sum(x >1)})
#aoSMCil1b_CPM$Hour4_CPM1_no <-  apply(aoSMCil1b_CPM[,4:6], 1, function(x){sum(x >1)})
#aoSMCil1b_CPM$Hour8_CPM1_no <-  apply(aoSMCil1b_CPM[,7:9], 1, function(x){sum(x >1)})
#aoSMCil1b_CPM$Hour24_CPM1_no <-  apply(aoSMCil1b_CPM[,1:3], 1, function(x){sum(x >1)})

aoSMCil1b_DEGs_DE <- filter(aoSMCil1b_DEGs_CPM, 
                        #to find DEGs in WT sham vs. cryo:
                        ((Hour0_meanCPM >1 | Hour2_meanCPM >1) & abs(LFC_0_2) > log2(1.5) & pa_0_2 <0.05) |
                          ((Hour0_meanCPM >1 | Hour4_meanCPM >1) & abs(LFC_0_4) > log2(1.5) & pa_0_4 <0.05) |
                          ((Hour0_meanCPM >1 | Hour6_meanCPM >1) & abs(LFC_0_6) > log2(1.5) & pa_0_6 <0.05))
dim(aoSMCil1b_DEGs_DE)[1] #no. DEGs total
sum(Biobase::rowMax(as.matrix(aoSMCil1b_CPM[,13:16])) >1)
169/18331 #9%

#pretty weak IL1b response in these aoSMCs

#biotypes:
trial <- merge(aoSMCil1b_CPM, GeneBiotypes, by.x = "geneID", by.y = "CAT_geneID")
table(filter(trial, Biobase::rowMax(as.matrix(trial[,14:17])) >1)$CAT_geneClassII)

#temporal ordering:
#regulated within 2 hours:
aoSMCil1b_DEGs_DE_Upwithin_2 <- filter(aoSMCil1b_DEGs_DE, (Hour0_meanCPM >=1 | Hour2_meanCPM >=1) &
                                     (LFC_0_2 >= log2(1.5) & pa_0_2 <0.05))
aoSMCil1b_DEGs_DE_Downwithin_2 <- filter(aoSMCil1b_DEGs_DE, (Hour0_meanCPM >=1 | Hour2_meanCPM >=1) &
                                       (LFC_0_2 < -log2(1.5) & pa_0_2 <0.05))

#of remainder, regulated within 4 hours:
aoSMCil1b_DEGs_DE_Upwithin_4 <- filter(aoSMCil1b_DEGs_DE, !geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_2$geneID, aoSMCil1b_DEGs_DE_Downwithin_2$geneID),
                                   (Hour0_meanCPM >=1 | Hour4_meanCPM >=1) &
                                     (LFC_0_4 >= log2(1.5) & pa_0_4 <0.05))
aoSMCil1b_DEGs_DE_Downwithin_4 <- filter(aoSMCil1b_DEGs_DE, !geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_2$geneID, aoSMCil1b_DEGs_DE_Downwithin_2$geneID),
                                     (Hour0_meanCPM >=1 | Hour4_meanCPM >=1) &
                                       (LFC_0_4 < -log2(1.5) & pa_0_4 <0.05))

#of remainder, regulated within 6 hours:
aoSMCil1b_DEGs_DE_Upwithin_6 <- filter(aoSMCil1b_DEGs_DE, !geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_2$geneID, aoSMCil1b_DEGs_DE_Downwithin_2$geneID,
                                                                 aoSMCil1b_DEGs_DE_Upwithin_4$geneID, aoSMCil1b_DEGs_DE_Downwithin_4$geneID),
                                   (Hour0_meanCPM >=1 | Hour6_meanCPM >=1) &
                                     (LFC_0_6 >= log2(1.5) & pa_0_6 <0.05))
aoSMCil1b_DEGs_DE_Downwithin_6 <- filter(aoSMCil1b_DEGs_DE, !geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_2$geneID, aoSMCil1b_DEGs_DE_Downwithin_2$geneID,
                                                                   aoSMCil1b_DEGs_DE_Upwithin_4$geneID, aoSMCil1b_DEGs_DE_Downwithin_4$geneID),
                                     (Hour0_meanCPM >=1 | Hour6_meanCPM >=1) &
                                       (LFC_0_6 < -log2(1.5) & pa_0_6 <0.05))


#for figures, assign each gene to wave, should have all the genes seperated by the above:
aoSMCil1b_DEGs_DE$RegulationStart[aoSMCil1b_DEGs_DE$geneID %in% aoSMCil1b_DEGs_DE_Upwithin_2$geneID] <- "Induced <2hrs"
aoSMCil1b_DEGs_DE$RegulationStart[aoSMCil1b_DEGs_DE$geneID %in% aoSMCil1b_DEGs_DE_Downwithin_2$geneID] <- "Repressed <2hrs"

aoSMCil1b_DEGs_DE$RegulationStart[aoSMCil1b_DEGs_DE$geneID %in% aoSMCil1b_DEGs_DE_Upwithin_4$geneID] <- "Induced 2-4hrs"
aoSMCil1b_DEGs_DE$RegulationStart[aoSMCil1b_DEGs_DE$geneID %in% aoSMCil1b_DEGs_DE_Downwithin_4$geneID] <- "Repressed 2-4hrs"

aoSMCil1b_DEGs_DE$RegulationStart[aoSMCil1b_DEGs_DE$geneID %in% aoSMCil1b_DEGs_DE_Upwithin_6$geneID] <- "Induced 4-6hrs"
aoSMCil1b_DEGs_DE$RegulationStart[aoSMCil1b_DEGs_DE$geneID %in% aoSMCil1b_DEGs_DE_Downwithin_6$geneID] <- "Repressed 4-6hrs"

table(aoSMCil1b_DEGs_DE$RegulationStart)


#
#### heatmap, can skip for stats ####

mat <- aoSMCil1b_DEGs_DE[,17:28] #4721 genes - similarly dynamic model than SVSMC (which is ~4345)
rownames(mat) <- aoSMCil1b_DEGs_DE[,1]
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

mati <- mat[rownames(mat) %in% aoSMCil1b_DEGs_DE_Upwithin_4$geneID,]
mati <- mat[rownames(mat) %in% aoSMCil1b_DEGs_DE_Downwithin_4$geneID,]
mati <- mat[rownames(mat) %in% aoSMCil1b_DEGs_DE_Upwithin_8$geneID,]
mati <- mat[rownames(mat) %in% aoSMCil1b_DEGs_DE_Downwithin_8$geneID,]
mati <- mat[rownames(mat) %in% aoSMCil1b_DEGs_DE_Upwithin_24$geneID,]
mati <- mat[rownames(mat) %in% aoSMCil1b_DEGs_DE_Downwithin_24$geneID,]

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
aoSMCil1b_DEGs_DE$Ens_ID_merge <- sapply(sapply(aoSMCil1b_DEGs_DE$geneID, strsplit, "\\."), "[[", 1)

#background:
aoSMCil1b_DEGs_CPM_exprs <- filter(aoSMCil1b_DEGs_CPM, rowMax(as.matrix(aoSMCil1b_DEGs_CPM[,26:29]))>1)

aoSMCil1b_DEGs_CPM_exprs$Ens_ID_merge <- sapply(sapply(aoSMCil1b_DEGs_CPM_exprs$geneID, strsplit, "\\."), "[[", 1)

trial <- split(aoSMCil1b_DEGs_DE, aoSMCil1b_DEGs_DE$RegulationStart)

library(clusterProfiler)
library(org.Hs.eg.db)
trial <- lapply(trial, function(x){
  enrichGO(gene          = unique(x$Ens_ID_merge),
           universe      = unique(aoSMCil1b_DEGs_CPM_exprs$Ens_ID_merge),
           keyType       = "ENSEMBL",
           OrgDb         = org.Hs.eg.db,
           ont           = "ALL",
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05,
           readable      = TRUE)
})

sapply(trial, dim)#only lost within 8hrs is weakly themed

#saveRDS(trial, "GO_list_aoSMCil1b.rds")

edox2 <- enrichplot::pairwise_termsim(trial[[6]])
enrichplot::treeplot(edox2, showCategory = 30, #fontSize =2, #extend = 0.1, #hilight =F #nWords = 3,
                     cluster.params = list(n = 7), 
                     #label_format_tiplab = function(x) stringr::str_wrap(x, width=40),
                     label_format = function(x) stringr::str_wrap(x, width=25))
View(data.frame(trial[[6]]))


#### temporal biases of PCGs and subclasses ####

#size of waves in general:
trial <- as.data.frame(table(aoSMCil1b_DEGs_DE$RegulationStart))
trial$Freq/dim(aoSMCil1b_DEGs_DE)[1]*100

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var1), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var1), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-24hrs", "24-36hrs", "36-48hrs")
trial$UpDown <- sapply(sapply(as.character(trial$Var1), strsplit, " "),"[[" , 1)
trial$PercAllLncs <- trial$Freq/dim(aoSMCil1b_DEGs_DE)[1]*100

trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,7,8,9,2,3,5,6)])

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
aoSMCil1b_DEGs_DE <- merge(GeneBiotypes[,c(1,2,4,9)], aoSMCil1b_DEGs_DE, by.x = "CAT_geneID", by.y = "geneID")

Cluster_biotype <- as.data.frame(table(aoSMCil1b_DEGs_DE$CAT_geneClassII, aoSMCil1b_DEGs_DE$RegulationStart))

table(aoSMCil1b_DEGs_DE$CAT_geneClassII)

#bias of PCGs:
Selected <- "coding_mRNA"
table(aoSMCil1b_DEGs_DE$CAT_geneClassII)[Selected]
table(aoSMCil1b_DEGs_DE$CAT_geneClassII)/dim(aoSMCil1b_DEGs_DE)[1]*100 

trial <- filter(Cluster_biotype, grepl(Selected, Var1))
trial$Freq/table(aoSMCil1b_DEGs_DE$RegulationStart)*100
#looks like possible biases, less than expected in 2x later points

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(aoSMCil1b_DEGs_DE$RegulationStart)

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

ClusterNames <- trial$Var2

PCGEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(aoSMCil1b_DEGs_DE$CAT_geneClassII)[Selected]
  d <- dim(aoSMCil1b_DEGs_DE)[1]
  
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
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-24hrs", "24-36hrs", "36-48hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,7,8,9,2,3,5,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(aoSMCil1b_DEGs_DE$CAT_geneClassII)[Selected]*100
trial$PercBackground <- trial$selection/dim(aoSMCil1b_DEGs_DE)[1]*100

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

aoSMCil1b_DEGs_DE$CC <- aoSMCil1b_DEGs_DE$GeneClassUpdate
aoSMCil1b_DEGs_DE$CC[gsub("\\.[0-9]*", "", aoSMCil1b_DEGs_DE$CAT_geneID) %in% CC_Freeman$Ensembl_ID] <- "CC"

Selected <- "CC"
table(aoSMCil1b_DEGs_DE$CC)[Selected]
table(aoSMCil1b_DEGs_DE$CC)/dim(aoSMCil1b_DEGs_DE)[1]*100

Cluster_biotype <- as.data.frame(table(aoSMCil1b_DEGs_DE$CC, aoSMCil1b_DEGs_DE$RegulationStart))
trial <- filter(Cluster_biotype, grepl(Selected, Var1))
trial$Freq/table(aoSMCil1b_DEGs_DE$RegulationStart)*100
#looks like definite biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(aoSMCil1b_DEGs_DE$RegulationStart)

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

ClusterNames <- trial$Var2

CCEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(aoSMCil1b_DEGs_DE$CC)[Selected]
  d <- dim(aoSMCil1b_DEGs_DE)[1]
  
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
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-24hrs", "24-36hrs", "36-48hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,7,8,9,2,3,5,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(aoSMCil1b_DEGs_DE$CC)[Selected]*100
trial$PercBackground <- trial$selection/dim(aoSMCil1b_DEGs_DE)[1]*100

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
  theme(text = element_text(size = 15))


#IEGs
#inserting new gene categories:
IEGs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Gene lists/arner_2015_table_S5_IEGs_lit.csv")
IEGs_hs <- filter(IEGs, !Hs_symbol %in% c("", "-"))
IEGs_hs <- filter(IEGs_hs, !grepl("[a-z]", IEGs_hs$Total_count))

aoSMCil1b_DEGs_DE$IEG <- aoSMCil1b_DEGs_DE$GeneClassUpdate
aoSMCil1b_DEGs_DE$IEG[aoSMCil1b_DEGs_DE$CAT_geneName %in% IEGs_hs$Hs_symbol] <- "IEG"

Selected <- "IEG"
table(aoSMCil1b_DEGs_DE$IEG)[Selected]
table(aoSMCil1b_DEGs_DE$IEG)/dim(aoSMCil1b_DEGs_DE)[1]*100

Cluster_biotype <- as.data.frame(table(aoSMCil1b_DEGs_DE$IEG, aoSMCil1b_DEGs_DE$RegulationStart))
trial <- filter(Cluster_biotype, grepl(Selected, Var1))
trial$Freq/table(aoSMCil1b_DEGs_DE$RegulationStart)*100
#looks like definite biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(aoSMCil1b_DEGs_DE$RegulationStart)

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

ClusterNames <- trial$Var2

IEGEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(aoSMCil1b_DEGs_DE$IEG)[Selected]
  d <- dim(aoSMCil1b_DEGs_DE)[1]
  
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
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-24hrs", "24-36hrs", "36-48hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,7,8,9,2,3,5,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(aoSMCil1b_DEGs_DE$IEG)[Selected]*100
trial$PercBackground <- trial$selection/dim(aoSMCil1b_DEGs_DE)[1]*100

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


#### (to review) combined figure ####

#combined figure:
AllTypesWaveBias <- rbind(TFWaveBias, IEGWaveBias, CCWaveBias)

#correct timeframe:
AllTypesWaveBias$FirstRegulation <- gsub("24hr", "16hr", AllTypesWaveBias$FirstRegulation)

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
                                      levels = levels(AllTypesWaveBias$UpDownType)[c(2,1,5,
                                                                                     8,7,6,
                                                                                     4,3)])

AllTypesWaveBias$FirstRegulation <- factor(AllTypesWaveBias$FirstRegulation)
AllTypesWaveBias$FirstRegulation <- factor(AllTypesWaveBias$FirstRegulation, 
                                           levels = levels(AllTypesWaveBias$FirstRegulation)[c(2,3,1)])


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
  scale_size_discrete(range = c(5,12), limits = c("p<0.05", "p<0.01", "p<0.001", "p<0.0001")) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "red") +
  theme_minimal() +
  theme(legend.key.size = unit(1.4, "line"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        axis.text.x = element_text(size=22),
        #axis.text.y = element_blank()
  )



#### biases of lncRNAs - hang on - 3 only! ####

table(aoSMCil1b_DEGs_DE$CAT_geneClassII) #3 lncRNAs only!
table(aoSMCil1b_DEGs_DE$CAT_geneClassII)/dim(aoSMCil1b_DEGs_DE)[1]*100

Cluster_biotype <- as.data.frame(table(aoSMCil1b_DEGs_DE$CAT_geneClassII, aoSMCil1b_DEGs_DE$RegulationStart))
trial <- filter(Cluster_biotype, Var1 == "lncRNA")
trial$Freq/table(aoSMCil1b_DEGs_DE$RegulationStart)*100
#looks like weak biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(aoSMCil1b_DEGs_DE$RegulationStart)

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

ClusterNames <- trial$Var2

LncEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(aoSMCil1b_DEGs_DE$CAT_geneClassII)["lncRNA"]
  d <- dim(aoSMCil1b_DEGs_DE)[1]
  
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
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-24hrs", "24-36hrs", "36-48hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,7,8,9,2,3,5,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(aoSMCil1b_DEGs_DE$CAT_geneClassII)["lncRNA"]*100
trial$PercBackground <- trial$selection/dim(aoSMCil1b_DEGs_DE)[1]*100

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
aoSMCil1b_libs <- c("geneID", aoSMCil1b_meta$CAGE_lib_ID)
#write.table(aoSMCil1b_libs, quote = F, "aoSMCil1b_libs_2026", row.names = F)

#CAGEs for the 24hr aoSMCil1b genes, take this list to FANTOM:
aoSMCil1b_exprs <- unique(filter(aoSMCil1b_CPM, Biobase::rowMax(as.matrix(aoSMCil1b_CPM[,13:16])) >1)$geneID)
#write.table(aoSMCil1b_exprs, quote = F, "aoSMCil1b_exprs_2026", row.names = F)

#now use "Subset_FANTOM5_aoSMCil1b.R" in eddie

#import, counts for selected libraries + genes:
aoSMCil1b_CAGEs_selectTime <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GeneralEffect_FANTOM_ER/aoSMCil1b_exprs_CAGEs_2026.csv", header = T)
dim(aoSMCil1b_CAGEs_selectTime)
colnames(aoSMCil1b_CAGEs_selectTime)
length(aoSMCil1b_libs)
#matches

#double check remove zeroes has worked via unix:
trial <- aoSMCil1b_CAGEs_selectTime
length(trial[,1]) == sum(!apply(trial[,-1], 1, function(row) all(row == 0)))

#find best expressed CAGEs per gene, this mapping object created alongside in eddie code:
#likely will exclude lncRNAs if use a >1 cut-off?
aoSMCil1b_CAGEs_mapping <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GeneralEffect_FANTOM_ER/aoSMCil1b_CAGEs_mapping_2026.csv", header = T)
length(unique(aoSMCil1b_CAGEs_mapping$geneID))#18331
length(unique(aoSMCil1b_exprs))#18331
#also means all genes have a CAGE site

#dont need transcripts:
aoSMCil1b_CAGEs_mapping <- unique(aoSMCil1b_CAGEs_mapping[,c(1,3)])

#there will be loads of crap TSS in these samples
#take all that provide >20% of a gene's output:
aoSMCil1b_CAGEs_selectTime <- unique(merge(aoSMCil1b_CAGEs_selectTime, aoSMCil1b_CAGEs_mapping, by.x = "prmtrID", by.y = "CAGEClusterID"))

trial <- split(aoSMCil1b_CAGEs_selectTime, aoSMCil1b_CAGEs_selectTime$geneID)

triali <- lapply(trial, function(x){
  x$prmtrID[rowMeans(x[,2:13])/sum(rowMeans(x[,2:13]))>0.2]
})

perc20 <- filter(aoSMCil1b_CAGEs_selectTime, prmtrID %in% unlist(triali) )

length(unique(perc20$geneID))#missing a tiny minority of genes (11 only)


#for remaining go to 10%:
perc10 <- filter(aoSMCil1b_CAGEs_selectTime, !geneID %in% perc20$geneID)

trial <- split(perc10, perc10$geneID)

triali <- lapply(trial, function(x){
  x$prmtrID[rowMeans(x[,2:13])/sum(rowMeans(x[,2:13]))>0.1]
})

perc10 <- filter(perc10, prmtrID %in% unlist(triali) )

perc20_10 <- rbind(perc20, perc10)

length(unique(perc20_10$geneID))
#18331 genes


#liftOver to hg38 and filter:
perc20_10$TSS_start <- as.numeric(sapply(strsplit(sapply(strsplit(perc20_10$prmtrID, "\\:"), "[[", 2), "\\.\\."), "[[", 1))
perc20_10$TSS_stop <- as.numeric(gsub(",[+-]", "", sapply(strsplit(sapply(strsplit(perc20_10$prmtrID, "\\:"), "[[", 2), "\\.\\."), "[[", 2)))
perc20_10$chr <- sapply(strsplit(perc20_10$prmtrID, "\\:"), "[[", 1)
head(perc20_10)

library(rtracklayer)
chain <- import.chain("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/hg19ToHg38.over.chain")

trial <- makeGRangesFromDataFrame(perc20_10[,c(1,15:17)], start.field = "TSS_start", end.field = "TSS_stop", 
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
TSSperG <- unique(trial[,c(1,2,15:19)])

length(unique(TSSperG$geneID))#18329 (2 lost in liftover)

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


#identify lncRNA-PCG co-reg pairs
AllLNC_AllPCG_1 <- filter(trialii, AnchorGene %in% filter(GeneBiotypes, CAT_geneClassII == "lncRNA")$CAT_geneID)
AllLNC_AllPCG_1 <- filter(AllLNC_AllPCG_1, geneID %in% filter(GeneBiotypes, CAT_geneClassII == "coding_mRNA")$CAT_geneID)
AllLNC_AllPCG_1$pairs <- paste(AllLNC_AllPCG_1$AnchorGene, AllLNC_AllPCG_1$geneID, sep = "-")

#identify PCG-PCG co-reg pairs (for later use)
AllPCG_AllPCG_1 <- filter(trialii, AnchorGene %in% filter(GeneBiotypes, CAT_geneClassII == "coding_mRNA")$CAT_geneID)
AllPCG_AllPCG_1 <- filter(AllPCG_AllPCG_1, geneID %in% filter(GeneBiotypes, CAT_geneClassII == "coding_mRNA")$CAT_geneID)
AllPCG_AllPCG_1$pairs <- paste(AllPCG_AllPCG_1$AnchorGene, AllPCG_AllPCG_1$geneID, sep = "-")

length(unique(AllLNC_AllPCG_1$AnchorGene))#3012 lncRNAs near a PCG
length(unique(AllLNC_AllPCG_1$geneID))#8611 PCGs near a lncRNA

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
#write.csv(AllLNC_AllPCG_1, "AllLNC_AllPCG_1_aoSMCil1b_2026.csv", row.names = F)


#PCG formating
length(unique(AllPCG_AllPCG_1$AnchorGene))#11906 anchor genes (bit meaningless if removing duplicate pairs)
length(unique(c(AllPCG_AllPCG_1$AnchorGene, AllPCG_AllPCG_1$geneID)))#11906 PCGs expressed near another PCG
length(unique(AllPCG_AllPCG_1$pairs))#82080 pairs, many will be duplicates in terms of PCGx -> PCGy and PCGy -> PCGx

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
#write.csv(AllPCG_AllPCG_1, "AllPCG_AllPCG_1_aoSMCil1b_2026.csv", row.names = F)


#
#### (needs revisit) Compare to control lncs, and SVSMC pairs ####

#compare to control cis lncs
ControlCisLncs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/ControlCisLncs.csv")

expressed_aoSMCil1b <- filter(aoSMCil1b_CPM, rowMax(as.matrix(aoSMCil1b_CPM[,c(13:16)])) >1)
sum(gsub("\\.[0-9]*", "", expressed_aoSMCil1b$geneID) %in% ControlCisLncs$x)#29 expressed
sum(unique(gsub("\\.[0-9]*", "", aoSMCil1b_DEGs_DE$CAT_geneID)) %in% ControlCisLncs$x)#5 are DE
sum(unique(gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended_aoSMCil1b$AnchorGene)) %in% ControlCisLncs$x)#2 in CClncs


#shared cis candidates between datasets:
CoRegPairs_04_48_24_extended <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CoRegPairs_04_48_24_extendedIII.csv")

#insert FANTOM ID where possible:
Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime.csv", header = T)
CoRegPairs_04_48_24_extended <- unique(merge(Enhancer_lociII[,c(1,14)], CoRegPairs_04_48_24_extended, by = "EnsID", all.y = T))

CoRegPairs_04_48_24_extended_aoSMCil1b$pairsII <- paste(
  gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended_aoSMCil1b$AnchorGene),
  gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended_aoSMCil1b$geneID),
  sep = "-"
)

CoRegPairs_04_48_24_extended$pairsII <- paste(
  gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended$FANTOM_ID),
  gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended$EnsID.y),
  sep = "-"
)


#11 x pairs found in both
sum(CoRegPairs_04_48_24_extended_aoSMCil1b$pairsII %in% CoRegPairs_04_48_24_extended$pairsII)

filter(CoRegPairs_04_48_24_extended, pairsII %in% CoRegPairs_04_48_24_extended_aoSMCil1b$pairsII)[,c(32,11)]
#INKILN is notable, NR2F2-AS1 and PITPNA-AS1 key loci, latter hits SERPINF1 which is not classed as influential (or a cis target) in SMC (locus is highlighted by other approach)
11/248 #4% of SVSMC pairs are common
11/834 #1% of aoSMCil1b pairs are common


#12x cclncRNAs predicted in both
sum(unique(gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended_aoSMCil1b$AnchorGene)) %in% gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended$FANTOM_ID))
length(unique(CoRegPairs_04_48_24_extended$FANTOM_ID))#106, but one is NA
12/106 #11% of cclncRNAs are found in both
length(unique(CoRegPairs_04_48_24_extended_aoSMCil1b$AnchorGene))#324
12/324 #4% of cclncRNAs are fond in both


#25X cclncRNA targets in both
sum(unique(gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended_aoSMCil1b$geneID)) %in% gsub("\\.[0-9]*", "", CoRegPairs_04_48_24_extended$EnsID.y))
length(unique(CoRegPairs_04_48_24_extended$EnsID.y))#106, but one is NA
25/236 #11% of cclncRNAs are found in both
length(unique(CoRegPairs_04_48_24_extended_aoSMCil1b$geneID))#532
25/532 #5% of cclncRNAs are fond in both

#consider overlapping the GO terms too - handle this later on through all datasets


#### Recreate the SVSMC analysis - import necessary table here ####

#Temporal bias vs. other DEGs
#Elevation in FC for lnc targets vs. other PCGs in early phase
#Increased chance of DE PCG neighbour for DE lnc targets vs other lncs* and DE PCGs in early phase
#GO terms
# *required HiC/eQTL to work

#need DEGs:
write.csv(aoSMCil1b_DEGs_DE, "aoSMCil1b_DEGs_DE_2026.csv", row.names = F)

#lnc-PCG neighbour table:
write.csv(AllLNC_AllPCG_1, "AllLNC_AllPCG_1_aoSMCil1b_2026.csv", row.names = F)

#PCG-PCG neighbour table:
write.csv(AllPCG_AllPCG_1, "AllPCG_AllPCG_1_aoSMCil1b_2026.csv", row.names = F)

#expressed genes table:
write.csv(aoSMCil1b_CPM[,c(17,13:16)], "aoSMCil1b_CPM_2026.csv", row.names = F)

#import (needs changing on revisit):
aoSMCil1b_DEGs_DE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GeneralEffect_FANTOM_ER/aoSMCil1b_DEGs_DE_2026.csv")
aoSMCil1b_CPM <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GeneralEffect_FANTOM_ER/aoSMCil1b_CPM_2026.csv")
AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GeneralEffect_FANTOM_ER/AllLNC_AllPCG_1_aoSMCil1b_2026.csv")
AllPCG_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GeneralEffect_FANTOM_ER/AllPCG_AllPCG_1_aoSMCil1b_2026.csv")
GeneBiotypes <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/nature21374-s2/table11.csv")

#simplify lncRNA class:
GeneBiotypes$CAT_geneClassII <- GeneBiotypes$CAT_geneClass
GeneBiotypes$CAT_geneClassII[
  grepl("lncRNA", GeneBiotypes$CAT_geneClassII)] <- "lncRNA"

#n.b. (START FROM HERE ON RETURN)
aoSMCil1b_DEGs_DE$CAT_geneID <- aoSMCil1b_DEGs_DE$geneID

aoSMCil1b_DEGs_DE_Upwithin_2 <- filter(aoSMCil1b_DEGs_DE, RegulationStart == "Induced <2hrs")
aoSMCil1b_DEGs_DE_Downwithin_2 <- filter(aoSMCil1b_DEGs_DE, RegulationStart == "Repressed <2hrs")

aoSMCil1b_DEGs_DE_Upwithin_4 <- filter(aoSMCil1b_DEGs_DE, RegulationStart == "Induced 2-4hrs")
aoSMCil1b_DEGs_DE_Downwithin_4 <- filter(aoSMCil1b_DEGs_DE, RegulationStart == "Repressed 2-4hrs")

aoSMCil1b_DEGs_DE_Upwithin_6 <- filter(aoSMCil1b_DEGs_DE, RegulationStart == "Induced 4-6hrs")
aoSMCil1b_DEGs_DE_Downwithin_6 <- filter(aoSMCil1b_DEGs_DE, RegulationStart == "Repressed 4-6hrs")

closestNeighbour_Lncs_aoSMCil1b <- split(AllLNC_AllPCG_1, AllLNC_AllPCG_1$AnchorGene)

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
closestNeighbour_Lncs_aoSMCil1b <- lapply(closestNeighbour_Lncs_aoSMCil1b, function(x){
  upstream <- filter(x, DistLnc_PCG < 0)
  downstream <- filter(x, DistLnc_PCG > 0)
  
  return(rbind(upstream[order(upstream$DistLnc_PCG, decreasing = T),][1,],
               downstream[order(downstream$DistLnc_PCG, decreasing = F),][1,]))
}
)

closestNeighbour_Lncs_aoSMCil1b <- bind_rows(closestNeighbour_Lncs_aoSMCil1b)
closestNeighbour_Lncs_aoSMCil1b <- filter(closestNeighbour_Lncs_aoSMCil1b, !is.na(pairs))


#now get the cis candidates:
#same timeframe only for simplicity
CoRegPairs_sameTimeframe_aoSMCil1b <- filter(closestNeighbour_Lncs_aoSMCil1b, 
                                         (AnchorGene %in% c(aoSMCil1b_DEGs_DE_Upwithin_2$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_2$CAT_geneID) & 
                                            geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_2$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_2$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(aoSMCil1b_DEGs_DE_Upwithin_4$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_4$CAT_geneID) & 
                                              geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_4$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_4$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(aoSMCil1b_DEGs_DE_Upwithin_6$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_6$CAT_geneID) & 
                                              geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_6$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_6$CAT_geneID)) 
                                         )
#0 CCLnc-target pairs

CoRegPairs_sameTimeframe_aoSMCil1b <- filter(closestNeighbour_Lncs_aoSMCil1b, 
                                             (AnchorGene %in% c(aoSMCil1b_DEGs_DE_Upwithin_2$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_2$CAT_geneID) & 
                                                geneID %in% aoSMCil1b_DEGs_DE$geneID))
#not even any with later change
CoRegPairs_sameTimeframe_aoSMCil1b <- filter(closestNeighbour_Lncs_aoSMCil1b, 
                                             (AnchorGene %in% aoSMCil1b_DEGs_DE$geneID & 
                                                geneID %in% aoSMCil1b_DEGs_DE$geneID))
#only 1 co-regulated pair found...


#expressed lncs:
length(unique(filter(GeneBiotypes, 
                     CAT_geneClassII == "lncRNA", 
                     CAT_geneID %in% aoSMCil1b_CPM$geneID[Biobase::rowMax(as.matrix(aoSMCil1b_CPM[,2:11])) >1]
)$CAT_geneID)) #3911

aoSMCil1b_CPM$geneID[Biobase::rowMax(as.matrix(aoSMCil1b_CPM[,2:11])) >1]


#
#### bias of IEGs ####

#IEGs
#inserting new gene categories:
IEGs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Gene lists/arner_2015_table_S5_IEGs_lit.csv")
IEGs_hs <- filter(IEGs, !Hs_symbol %in% c("", "-"))
IEGs_hs <- filter(IEGs_hs, !grepl("[a-z]", IEGs_hs$Total_count))

aoSMCil1b_DEGs_DE$IEG <- aoSMCil1b_DEGs_DE$GeneClassUpdate
aoSMCil1b_DEGs_DE$IEG[aoSMCil1b_DEGs_DE$CAT_geneName %in% IEGs_hs$Hs_symbol] <- "IEG"

Selected <- "IEG"
table(aoSMCil1b_DEGs_DE$IEG)[Selected]
table(aoSMCil1b_DEGs_DE$IEG)/dim(aoSMCil1b_DEGs_DE)[1]*100

Cluster_biotype <- as.data.frame(table(aoSMCil1b_DEGs_DE$IEG, aoSMCil1b_DEGs_DE$RegulationStart))
trial <- filter(Cluster_biotype, grepl(Selected, Var1))
trial$Freq/table(aoSMCil1b_DEGs_DE$RegulationStart)*100
#looks like definite biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(aoSMCil1b_DEGs_DE$RegulationStart)

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

ClusterNames <- trial$Var2

IEGEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(aoSMCil1b_DEGs_DE$IEG)[Selected]
  d <- dim(aoSMCil1b_DEGs_DE)[1]
  
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
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-24hrs", "24-36hrs", "36-48hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,7,8,9,2,3,5,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(aoSMCil1b_DEGs_DE$IEG)[Selected]*100
trial$PercBackground <- trial$selection/dim(aoSMCil1b_DEGs_DE)[1]*100

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

table(aoSMCil1b_DEGs_DE$CAT_geneClassII)
table(aoSMCil1b_DEGs_DE$CAT_geneClassII)/dim(aoSMCil1b_DEGs_DE)[1]*100

Cluster_biotype <- as.data.frame(table(aoSMCil1b_DEGs_DE$CAT_geneClassII, aoSMCil1b_DEGs_DE$RegulationStart))
trial <- filter(Cluster_biotype, Var1 == "lncRNA")
trial$Freq/table(aoSMCil1b_DEGs_DE$RegulationStart)*100
#looks like weak biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(aoSMCil1b_DEGs_DE$RegulationStart)

ClusterNames <- trial$Var2

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

LncEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(aoSMCil1b_DEGs_DE$CAT_geneClassII)["lncRNA"]
  d <- dim(aoSMCil1b_DEGs_DE)[1]
  
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
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-24hrs", "24-36hrs", "36-48hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,7,8,9,2,3,5,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(aoSMCil1b_DEGs_DE$CAT_geneClassII)["lncRNA"]*100
trial$PercBackground <- trial$selection/dim(aoSMCil1b_DEGs_DE)[1]*100

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

#plot enrichment of cis lnc
#plot genes per cluster

#cis lnc
aoSMCil1b_DEGs_DE$CisLnc <- aoSMCil1b_DEGs_DE$CAT_geneClassII
aoSMCil1b_DEGs_DE$CisLnc[aoSMCil1b_DEGs_DE$CAT_geneID %in% CoRegPairs_sameTimeframe_aoSMCil1b$AnchorGene] <- "CisLnc"
table(aoSMCil1b_DEGs_DE$CisLnc)

Cluster_biotype <- as.data.frame(table(aoSMCil1b_DEGs_DE$CisLnc, aoSMCil1b_DEGs_DE$RegulationStart))

#bias of category:
table(aoSMCil1b_DEGs_DE$CisLnc)["CisLnc"]
table(aoSMCil1b_DEGs_DE$CisLnc)/dim(aoSMCil1b_DEGs_DE)[1]*100

trial <- filter(Cluster_biotype, Var1 == "CisLnc")
trial$Freq/table(aoSMCil1b_DEGs_DE$RegulationStart)*100
#looks like strong biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(aoSMCil1b_DEGs_DE$RegulationStart)

ClusterNames <- trial$Var2

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

LncEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(aoSMCil1b_DEGs_DE$CisLnc)["CisLnc"]
  d <- dim(aoSMCil1b_DEGs_DE)[1]
  
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
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-24hrs", "24-36hrs", "36-48hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,7,8,9,2,3,5,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(aoSMCil1b_DEGs_DE$CisLnc)["CisLnc"]*100
trial$PercBackground <- trial$selection/dim(aoSMCil1b_DEGs_DE)[1]*100

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
  #scale_y_continuous(limits = c(-30,48),breaks = seq(-20,48, by = 20),
  #                   labels = (c(seq(20, 0, by = -20), seq(20,48,by=20)))) +
  theme_minimal() +
  theme(text = element_text(size = 15))

#cis lncRNAs again biased to early activation, may have their primary role early on

#compare to all lncs:
sum(LncWaveBias$Freq)
sum(CisLncWaveBias$Freq)

a <- 50
b <- 72
c <- 99
d <- 361

fisher.test(data.frame("LncRNA" = c(a,b-a),
                       "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")

trial <- data.frame(c/d*100, a/b*100)

colnames(trial) <- c("All DE lncRNAs", "CClncRNAs (same)")

#melt(trial)

trial <- melt(trial)

trial$variable <- factor(trial$variable)
trial$variable <- factor(trial$variable, levels = levels(trial$variable)[c(1,3,2)])

ggplot(trial) + aes(x = value, y = variable) +
  geom_bar(stat = "identity", fill = "olivedrab3", color = "grey48") +
  xlab("% 0-4hr induced") +
  ylab("") +
  theme_minimal() +
  #scale_x_continuous(breaks = seq(0,50,25), limits = c(0,50)) +
  theme(text = element_text(size =20))

#
#### elncRNAs ####

#alternate (less specific) indicator of cis-acting function would be transcription from an annotated enhancer

aoSMCil1b_DEGs_DE$CAT_geneClassIII <- aoSMCil1b_DEGs_DE$CAT_geneClassII
aoSMCil1b_DEGs_DE$CAT_geneClassIII[aoSMCil1b_DEGs_DE$CAT_geneID %in% filter(GeneBiotypes, CAT_geneCategory == "e_lncRNA")$CAT_geneID] <- "e_lncRNA"

table(aoSMCil1b_DEGs_DE$CAT_geneClassIII)

Cluster_biotype <- as.data.frame(table(aoSMCil1b_DEGs_DE$CAT_geneClassIII, aoSMCil1b_DEGs_DE$RegulationStart))

table(aoSMCil1b_DEGs_DE$CAT_geneClassIII)

#bias of category:
table(aoSMCil1b_DEGs_DE$CAT_geneClassIII)["e_lncRNA"]
table(aoSMCil1b_DEGs_DE$CAT_geneClassIII)/dim(aoSMCil1b_DEGs_DE)[1]*100

trial <- filter(Cluster_biotype, Var1 == "e_lncRNA")
trial$Freq/table(aoSMCil1b_DEGs_DE$RegulationStart)*100

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(aoSMCil1b_DEGs_DE$RegulationStart)

ClusterNames <- trial$Var2

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

LncEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(aoSMCil1b_DEGs_DE$CAT_geneClassIII)["e_lncRNA"]
  d <- dim(aoSMCil1b_DEGs_DE)[1]
  
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
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-24hrs", "24-36hrs", "36-48hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,7,8,9,2,3,5,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(aoSMCil1b_DEGs_DE$CAT_geneClassIII)["e_lncRNA"]*100
trial$PercBackground <- trial$selection/dim(aoSMCil1b_DEGs_DE)[1]*100

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
  #scale_y_continuous(limits = c(-30,48),breaks = seq(-20,48, by = 20),
  #                   labels = (c(seq(20, 0, by = -20), seq(20,48,by=20)))) +
  theme_minimal() +
  theme(text = element_text(size = 15))

#early bias

#percentage of cclncRNAs amongst DE elncRNAs vs other DE lncRNAs:
sum(ELncWaveBias$Freq)
sum(LncWaveBias$Freq)
sum(CisLncWaveBias$Freq)

sum(filter(aoSMCil1b_DEGs_DE, CAT_geneClassIII == "e_lncRNA")$CAT_geneID %in% filter(aoSMCil1b_DEGs_DE, CisLnc == "CisLnc")$CAT_geneID)

64/237 #27% of cclncs are elncs
277/1199 #23% of all DELs are elncs

aii <- 25
bii <- 72
cii <- 97   
dii <- 361   

fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

trial <- data.frame(cii/dii*100, aii/bii*100)

colnames(trial) <- c("All DE lncRNAs", "CClncRNAs (same)")

#melt(trial)

trial <- melt(trial)

trial$variable <- factor(trial$variable)
trial$variable <- factor(trial$variable, levels = levels(trial$variable)[c(1,3,2)])

ggplot(trial) + aes(x = value, y = variable) +
  geom_bar(stat = "identity", fill = "olivedrab3", color = "grey48") +
  xlab("% ElncRNAs") +
  ylab("") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0,40,20), limits = c(0,45)) +
  theme(text = element_text(size =20))


#
#### all co-regulated (time-agnostic) ####

#plot enrichment of cis lnc, non cis lnc, e cis lnc, non-e cis lnc:
#plot genes per cluster

CoReg <- filter(AllLNC_AllPCG_1, AnchorGene %in% aoSMCil1b_DEGs_DE$CAT_geneID, geneID %in% aoSMCil1b_DEGs_DE$CAT_geneID)

#cis lnc
aoSMCil1b_DEGs_DE$CoRegLnc <- aoSMCil1b_DEGs_DE$CAT_geneClassII
aoSMCil1b_DEGs_DE$CoRegLnc[aoSMCil1b_DEGs_DE$CAT_geneID %in% CoReg$AnchorGene] <- "CoRegLnc"
table(aoSMCil1b_DEGs_DE$CoRegLnc)

Cluster_biotype <- as.data.frame(table(aoSMCil1b_DEGs_DE$CoRegLnc, aoSMCil1b_DEGs_DE$RegulationStart))

#bias of category:
table(aoSMCil1b_DEGs_DE$CoRegLnc)["CoRegLnc"]
table(aoSMCil1b_DEGs_DE$CoRegLnc)/dim(aoSMCil1b_DEGs_DE)[1]*100

trial <- filter(Cluster_biotype, Var1 == "CoRegLnc")
trial$Freq/table(aoSMCil1b_DEGs_DE$RegulationStart)*100
#looks like strong biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(aoSMCil1b_DEGs_DE$RegulationStart)

ClusterNames <- trial$Var2

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

LncEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(aoSMCil1b_DEGs_DE$CoRegLnc)["CoRegLnc"]
  d <- dim(aoSMCil1b_DEGs_DE)[1]
  
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
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-24hrs", "24-36hrs", "36-48hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,7,8,9,2,3,5,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(aoSMCil1b_DEGs_DE$CoRegLnc)["CoRegLnc"]*100
trial$PercBackground <- trial$selection/dim(aoSMCil1b_DEGs_DE)[1]*100

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
  #scale_y_continuous(limits = c(-30,48),breaks = seq(-20,48, by = 20),
  #                   labels = (c(seq(20, 0, by = -20), seq(20,48,by=20)))) +
  theme_minimal() +
  theme(text = element_text(size = 15))

#basically same as lncs generally


#
#### combined figure (and table save) ####

#combined figure:
AllTypesWaveBias <- rbind(IEGWaveBias,
                          LncWaveBias, 
                          ELncWaveBias, #NonELncWaveBias, 
                          #CoRegLncWaveBias, #NonCoRegLncWaveBias,
                          CisLncWaveBias#, CisLncSameWaveBias #NonCisLncWaveBias
)

AllTypesWaveBias$FirstRegulation <- as.character(AllTypesWaveBias$FirstRegulation)

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
write.csv(AllTypesWaveBias, "aoSMCil1b_AllTypesWaveBias.csv", row.names = F)


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
                                           levels = levels(AllTypesWaveBias$FirstRegulation)[c(1,4,7,8,9,2,3,5,6)])

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
  scale_size_discrete(range = c(3,12), limits = c("p<0.1", "p<0.05", "p<0.01", "p<0.001", 
                                                  "p<0.0001")) +
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

closestNeighbour_PCGs_aoSMCil1b <- split(AllPCG_AllPCG_1, AllPCG_AllPCG_1$AnchorGene)

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
closestNeighbour_PCGs_aoSMCil1b <- lapply(closestNeighbour_PCGs_aoSMCil1b, function(x){
  upstream <- filter(x, DistLnc_PCG < 0)
  downstream <- filter(x, DistLnc_PCG > 0)
  
  return(rbind(upstream[order(upstream$DistLnc_PCG, decreasing = T),][1,],
               downstream[order(downstream$DistLnc_PCG, decreasing = F),][1,]))
}
)

closestNeighbour_PCGs_aoSMCil1b <- bind_rows(closestNeighbour_PCGs_aoSMCil1b)
closestNeighbour_PCGs_aoSMCil1b <- filter(closestNeighbour_PCGs_aoSMCil1b, !is.na(pairs))


#now get the cis candidates:
#same timeframe only for simplicity
CoRegPairs_sameTimeframe_aoSMCil1b_PCG <- filter(closestNeighbour_PCGs_aoSMCil1b, 
                                             (AnchorGene %in% c(aoSMCil1b_DEGs_DE_Upwithin_2$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_2$CAT_geneID) & 
                                                geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_2$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_2$CAT_geneID)) |
                                               
                                               (AnchorGene %in% c(aoSMCil1b_DEGs_DE_Upwithin_4$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_4$CAT_geneID) & 
                                                  geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_4$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_4$CAT_geneID)) |
                                               
                                               (AnchorGene %in% c(aoSMCil1b_DEGs_DE_Upwithin_6$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_6$CAT_geneID) & 
                                                  geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_6$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_6$CAT_geneID)) |
                                               
                                               (AnchorGene %in% c(aoSMCil1b_DEGs_DE_Upwithin_8$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_8$CAT_geneID) & 
                                                  geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_8$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_8$CAT_geneID)) |
                                               
                                               (AnchorGene %in% c(aoSMCil1b_DEGs_DE_Upwithin_12$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_12$CAT_geneID) & 
                                                  geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_12$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_12$CAT_geneID)) |
                                               
                                               (AnchorGene %in% c(aoSMCil1b_DEGs_DE_Upwithin_16$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_16$CAT_geneID) & 
                                                  geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_16$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_16$CAT_geneID)) |
                                               
                                               (AnchorGene %in% c(aoSMCil1b_DEGs_DE_Upwithin_24$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_24$CAT_geneID) & 
                                                  geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_24$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_24$CAT_geneID)) |
                                               
                                               (AnchorGene %in% c(aoSMCil1b_DEGs_DE_Upwithin_36$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_36$CAT_geneID) & 
                                                  geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_36$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_36$CAT_geneID)) |
                                               
                                               (AnchorGene %in% c(aoSMCil1b_DEGs_DE_Upwithin_48$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_48$CAT_geneID) & 
                                                  geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_48$CAT_geneID, aoSMCil1b_DEGs_DE_Downwithin_48$CAT_geneID))
)
#597 CCLnc-target pairs

length(unique(CoRegPairs_sameTimeframe_aoSMCil1b_PCG$AnchorGene))#566 CClncRNAs
length(unique(CoRegPairs_sameTimeframe_aoSMCil1b_PCG$geneID))#556 potential targets


#cis lnc
aoSMCil1b_DEGs_DE$CisPCG <- aoSMCil1b_DEGs_DE$CAT_geneClassII
aoSMCil1b_DEGs_DE$CisPCG[aoSMCil1b_DEGs_DE$CAT_geneID %in% c(CoRegPairs_sameTimeframe_aoSMCil1b_PCG$AnchorGene)] <- "CisPCG"
table(aoSMCil1b_DEGs_DE$CisPCG)

Cluster_biotype <- as.data.frame(table(aoSMCil1b_DEGs_DE$CisPCG, aoSMCil1b_DEGs_DE$RegulationStart))

#bias of category:
table(aoSMCil1b_DEGs_DE$CisPCG)["CisPCG"]
table(aoSMCil1b_DEGs_DE$CisPCG)/dim(aoSMCil1b_DEGs_DE)[1]*100

trial <- filter(Cluster_biotype, Var1 == "CisPCG")
trial$Freq/table(aoSMCil1b_DEGs_DE$RegulationStart)*100
#looks like strong biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(aoSMCil1b_DEGs_DE$RegulationStart)

ClusterNames <- trial$Var2

#order by time:
trial$timeQuery <- as.numeric(gsub("hrs", "", sapply(strsplit(as.character(trial$Var2), "[<-]"), "[[",2)))
trial$UpDown <- sapply(strsplit(as.character(trial$Var2), " "), "[[",1)
trial <- trial[order(trial$UpDown, trial$timeQuery),]

PCGEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(aoSMCil1b_DEGs_DE$CisPCG)["CisPCG"]
  d <- dim(aoSMCil1b_DEGs_DE)[1]
  
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
trial$FirstRegulation <- c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-24hrs", "24-36hrs", "36-48hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(1,4,7,8,9,2,3,5,6)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(aoSMCil1b_DEGs_DE$CisPCG)["CisPCG"]*100
trial$PercBackground <- trial$selection/dim(aoSMCil1b_DEGs_DE)[1]*100

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
  #scale_y_continuous(limits = c(-30,48),breaks = seq(-20,48, by = 20),
  #                   labels = (c(seq(20, 0, by = -20), seq(20,48,by=20)))) +
  theme_minimal() +
  theme(text = element_text(size = 15))

#cis PCGs somewhat early

#compare to DE lncs with co-reg neighbours:
sum(CisLncWaveBias$Freq)
sum(CisPCGWaveBias$Freq)

a <- 50 #hit: 0-4 induced
b <- 72 #selection - all CC lncs
c <- 88
d <- 566 #background - all CC PCGs
fisher.test(data.frame("cisLnc" = c(a, b-a),
                       "Not"   = c(c, d-c)), alternative = "greater")

#plot
trial <- data.frame("variable" = c("DE LncRNA w/ DE PCG\nneighbour(same)", "DE PCG w/ DE PCG\nneighbour(same)"),
                    "value" = c(a/b*100, c/d*100))

trial$variable <- as.factor(trial$variable)
trial$variable <- factor(trial$variable, levels(trial$variable)[c(2,1)])

ggplot(trial) + aes(y = variable, x = value, fill = variable) +
  geom_bar(stat = "identity", color = "grey48") +
  scale_fill_manual(values = c(`DE LncRNA w/ DE PCG\nneighbour(same)` = "olivedrab3", `DE PCG w/ DE PCG\nneighbour(same)` = "mediumorchid")) +
  xlab("% 0-4hr induced") +
  ylab(" ") +
  theme_minimal() + 
  theme(text = element_text(size=24)) + 
  theme(legend.position = "none") #+ Seurat::RotatedAxis()


#
#### (skip) same timeframe and PCG comparison

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
length(unique(c(CoRegPairs_04_48_24_extended_aoSMCil1b$AnchorGene)))#324 CClncs
length(unique(c(CoRegPairs_04_48_24_extended_aoSMCil1b_PCG$AnchorGene)))#2132 PCG
PCG_PCG_4hr_timing

a <- 95
b <- 324 #selection - all CC lncs
c <- 448
d <- 2132 #background - all CC PCGs
fisher.test(data.frame("cisLnc" = c(a, b-a),
                       "Not"   = c(c, d-c)), alternative = "greater")#***
statmod::power.fisher.test(a/b, c/d, c, d)


#or a direct comparison for CClncs in same timeframe:
length(unique(c(CoRegPairs_04_48_24_aoSMCil1b_same$AnchorGene)))#259 sametimeframe CClncs
length(unique(c(CoRegPairs_04_48_24_aoSMCil1b_samePCG$AnchorGene)))#1567 PCG
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
aoSMCil1b_CPM_02 <- filter(aoSMCil1b_CPM, (Hour0_meanCPM>1 | Hour2_meanCPM>1), 
                       !geneID %in% filter(aoSMCil1b_DEGs_DE, !grepl("<2hrs", RegulationStart))$CAT_geneID)

aoSMCil1b_CPM_24 <- filter(aoSMCil1b_CPM, (Hour2_meanCPM>1 | Hour4_meanCPM>1), 
                       !geneID %in% filter(aoSMCil1b_DEGs_DE, !grepl("2-4hrs", RegulationStart))$CAT_geneID)

aoSMCil1b_CPM_46 <- filter(aoSMCil1b_CPM, (Hour4_meanCPM>1 | Hour6_meanCPM>1), 
                       !geneID %in% filter(aoSMCil1b_DEGs_DE, !grepl("4-6hrs", RegulationStart))$CAT_geneID)

aoSMCil1b_CPM_68 <- filter(aoSMCil1b_CPM, (Hour6_meanCPM>1 | Hour8_meanCPM>1), 
                       !geneID %in% filter(aoSMCil1b_DEGs_DE, !grepl("6-8hrs", RegulationStart))$CAT_geneID)

aoSMCil1b_CPM_812 <- filter(aoSMCil1b_CPM, (Hour8_meanCPM>1 | Hour12_meanCPM>1), 
                        !geneID %in% filter(aoSMCil1b_DEGs_DE, !grepl("8-12hrs", RegulationStart))$CAT_geneID)

aoSMCil1b_CPM_1216 <- filter(aoSMCil1b_CPM, (Hour12_meanCPM>1 | Hour16_meanCPM>1), 
                         !geneID %in% filter(aoSMCil1b_DEGs_DE, !grepl("12-16hrs", RegulationStart))$CAT_geneID)

aoSMCil1b_CPM_1624 <- filter(aoSMCil1b_CPM, (Hour16_meanCPM>1 | Hour24_meanCPM>1), 
                         !geneID %in% filter(aoSMCil1b_DEGs_DE, !grepl("16-24hrs", RegulationStart))$CAT_geneID)

aoSMCil1b_CPM_2436 <- filter(aoSMCil1b_CPM, (Hour24_meanCPM>1 | Hour36_meanCPM>1), 
                         !geneID %in% filter(aoSMCil1b_DEGs_DE, !grepl("24-36hrs", RegulationStart))$CAT_geneID)

aoSMCil1b_CPM_3648 <- filter(aoSMCil1b_CPM, (Hour36_meanCPM>1 | Hour48_meanCPM>1), 
                         !geneID %in% filter(aoSMCil1b_DEGs_DE, !grepl("36-48hrs", RegulationStart))$CAT_geneID)

#(n.b. alternate for same-timeframe is to just take expressed and not DE in previous timeframe)
backgrounds_list <- list(aoSMCil1b_CPM_02, aoSMCil1b_CPM_24, aoSMCil1b_CPM_46, aoSMCil1b_CPM_68,
                         aoSMCil1b_CPM_812, aoSMCil1b_CPM_1216, aoSMCil1b_CPM_1624, aoSMCil1b_CPM_2436, aoSMCil1b_CPM_3648)

#just up DELs:
selection_list <- list(aoSMCil1b_DEGs_DE_Upwithin_2, aoSMCil1b_DEGs_DE_Upwithin_4, aoSMCil1b_DEGs_DE_Upwithin_6, 
                       aoSMCil1b_DEGs_DE_Upwithin_8, aoSMCil1b_DEGs_DE_Upwithin_12, aoSMCil1b_DEGs_DE_Upwithin_16,
                       aoSMCil1b_DEGs_DE_Upwithin_24, aoSMCil1b_DEGs_DE_Upwithin_36, aoSMCil1b_DEGs_DE_Upwithin_48)
#with up targets:
hits_list <- list(aoSMCil1b_DEGs_DE_Upwithin_2, aoSMCil1b_DEGs_DE_Upwithin_4, aoSMCil1b_DEGs_DE_Upwithin_6, 
                  aoSMCil1b_DEGs_DE_Upwithin_8, aoSMCil1b_DEGs_DE_Upwithin_12, aoSMCil1b_DEGs_DE_Upwithin_16,
                  aoSMCil1b_DEGs_DE_Upwithin_24, aoSMCil1b_DEGs_DE_Upwithin_36, aoSMCil1b_DEGs_DE_Upwithin_48)

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
                                           results_list[[8]][1]/results_list[[8]][2]*100,
                                           results_list[[9]][1]/results_list[[9]][2]*100),
                           
                           "EL_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, 
                                          results_list[[2]][3]/results_list[[2]][4]*100, 
                                          results_list[[3]][3]/results_list[[3]][4]*100,
                                          results_list[[4]][3]/results_list[[4]][4]*100,
                                          results_list[[5]][3]/results_list[[5]][4]*100,
                                          results_list[[6]][3]/results_list[[6]][4]*100,
                                          results_list[[7]][3]/results_list[[7]][4]*100,
                                          results_list[[8]][3]/results_list[[8]][4]*100,
                                          results_list[[9]][1]/results_list[[9]][2]*100))

DEL_PCG_type$RegulationBegins <- as.factor(c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-24hrs", "24-36hrs", "36-48hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1], results_list[[4]][1],
                        results_list[[5]][1], results_list[[6]][1], results_list[[7]][1], results_list[[8]][1],
                        results_list[[9]][1])

g2bi <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = EL_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17,23,29,35,41,47,53)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18,24,30,36,36,48,54)]
#write.csv(DEL_PCG_type, "aoSMCil1b_inducedDEL_DEP_DEneighboursame_conc.csv", row.names = F)

p.adjust(DEL_PCG_type$p, method = "bonferroni")

DEL_PCG_type_plot <- unique(melt(DEL_PCG_type[,1:3]))

DEL_PCG_type_plot$variable <- as.character(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable[grepl("^DEL", DEL_PCG_type_plot$variable)] <- "Induced lncRNAs"
DEL_PCG_type_plot$variable[grepl("^EL", DEL_PCG_type_plot$variable)] <- "Non-DE lncRNAs"

DEL_PCG_type_plot$variable <- as.factor(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable <- factor(DEL_PCG_type_plot$variable, levels = levels(DEL_PCG_type_plot$variable)[c(1,2)])

DEL_PCG_type_plot$RegulationBegins <- as.factor(DEL_PCG_type_plot$RegulationBegins)
DEL_PCG_type_plot$RegulationBegins <- factor(DEL_PCG_type_plot$RegulationBegins, 
                                             levels = levels(DEL_PCG_type_plot$RegulationBegins)[c(1,4,7,8,9,2,3,5,6)])

ggplot(DEL_PCG_type_plot) + aes(x = RegulationBegins, y = value, fill = variable) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.7), color = "grey48") +
  theme_minimal() +
  theme(text = element_text(size = 24)) + #, legend.position = "none") +
  scale_fill_manual(values = c(`Induced lncRNAs` = "#D6484D", `Non-DE lncRNAs` = "grey30")) +
  xlab("") +
  ylab("% With induced\nneighbour") + Seurat::RotatedAxis()


#
#### enrichment of gene repression near induced lncs vs near expressed lncs ####

#(n.b. alternate for same-timeframe is to just take expressed and not DE in previous timeframe)
backgrounds_list <- list(aoSMCil1b_CPM_02, aoSMCil1b_CPM_24, aoSMCil1b_CPM_46, aoSMCil1b_CPM_68,
                         aoSMCil1b_CPM_812, aoSMCil1b_CPM_1216, aoSMCil1b_CPM_1624, aoSMCil1b_CPM_2436, aoSMCil1b_CPM_3648)

#just up DELs:
selection_list <- list(aoSMCil1b_DEGs_DE_Upwithin_2, aoSMCil1b_DEGs_DE_Upwithin_4, aoSMCil1b_DEGs_DE_Upwithin_6, 
                       aoSMCil1b_DEGs_DE_Upwithin_8, aoSMCil1b_DEGs_DE_Upwithin_12, aoSMCil1b_DEGs_DE_Upwithin_16,
                       aoSMCil1b_DEGs_DE_Upwithin_24, aoSMCil1b_DEGs_DE_Upwithin_36, aoSMCil1b_DEGs_DE_Upwithin_48)
#with up targets:
hits_list <- list(aoSMCil1b_DEGs_DE_Downwithin_2, aoSMCil1b_DEGs_DE_Downwithin_4, aoSMCil1b_DEGs_DE_Downwithin_6, 
                  aoSMCil1b_DEGs_DE_Downwithin_8, aoSMCil1b_DEGs_DE_Downwithin_12, aoSMCil1b_DEGs_DE_Downwithin_16,
                  aoSMCil1b_DEGs_DE_Upwithin_24, aoSMCil1b_DEGs_DE_Downwithin_36, aoSMCil1b_DEGs_DE_Downwithin_48)

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
                                           results_list[[8]][1]/results_list[[8]][2]*100,
                                           results_list[[9]][1]/results_list[[9]][2]*100),
                           
                           "EL_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, 
                                          results_list[[2]][3]/results_list[[2]][4]*100, 
                                          results_list[[3]][3]/results_list[[3]][4]*100,
                                          results_list[[4]][3]/results_list[[4]][4]*100,
                                          results_list[[5]][3]/results_list[[5]][4]*100,
                                          results_list[[6]][3]/results_list[[6]][4]*100,
                                          results_list[[7]][3]/results_list[[7]][4]*100,
                                          results_list[[8]][3]/results_list[[8]][4]*100,
                                          results_list[[9]][1]/results_list[[9]][2]*100))

DEL_PCG_type$RegulationBegins <- as.factor(c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-24hrs", "24-36hrs", "36-48hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1], results_list[[4]][1],
                        results_list[[5]][1], results_list[[6]][1], results_list[[7]][1], results_list[[8]][1],
                        results_list[[9]][1])

g2bi <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = EL_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17,23,29,35,41,47,53)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18,24,30,36,36,48,54)]
#write.csv(DEL_PCG_type, "aoSMCil1b_inducedDEL_DEP_DEneighboursame_conc.csv", row.names = F)

p.adjust(DEL_PCG_type$p, method = "bonferroni")

DEL_PCG_type_plot <- unique(melt(DEL_PCG_type[,1:3]))

DEL_PCG_type_plot$variable <- as.character(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable[grepl("^DEL", DEL_PCG_type_plot$variable)] <- "Induced lncRNAs"
DEL_PCG_type_plot$variable[grepl("^EL", DEL_PCG_type_plot$variable)] <- "Non-DE lncRNAs"

DEL_PCG_type_plot$variable <- as.factor(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable <- factor(DEL_PCG_type_plot$variable, levels = levels(DEL_PCG_type_plot$variable)[c(1,2)])

DEL_PCG_type_plot$RegulationBegins <- as.factor(DEL_PCG_type_plot$RegulationBegins)
DEL_PCG_type_plot$RegulationBegins <- factor(DEL_PCG_type_plot$RegulationBegins, 
                                             levels = levels(DEL_PCG_type_plot$RegulationBegins)[c(1,4,7,8,9,2,3,5,6)])

ggplot(DEL_PCG_type_plot) + aes(x = RegulationBegins, y = value, fill = variable) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.7), color = "grey48") +
  theme_minimal() +
  theme(text = element_text(size = 24)) + #, legend.position = "none") +
  scale_fill_manual(values = c(`Induced lncRNAs` = "#D6484D", `Non-DE lncRNAs` = "grey30")) +
  xlab("") +
  ylab("% With repressed\nneighbour") + Seurat::RotatedAxis()


#
#### enrichment of gene induction near induced lncs vs near induced PCGs ####

#just up DELs:
backgrounds_list <- list(aoSMCil1b_DEGs_DE_Upwithin_2, aoSMCil1b_DEGs_DE_Upwithin_4, aoSMCil1b_DEGs_DE_Upwithin_6, 
                         aoSMCil1b_DEGs_DE_Upwithin_8, aoSMCil1b_DEGs_DE_Upwithin_12, aoSMCil1b_DEGs_DE_Upwithin_16,
                         aoSMCil1b_DEGs_DE_Upwithin_24, aoSMCil1b_DEGs_DE_Upwithin_36, aoSMCil1b_DEGs_DE_Upwithin_48)
#with up targets:
hits_list <- list(aoSMCil1b_DEGs_DE_Upwithin_2, aoSMCil1b_DEGs_DE_Upwithin_4, aoSMCil1b_DEGs_DE_Upwithin_6, 
                  aoSMCil1b_DEGs_DE_Upwithin_8, aoSMCil1b_DEGs_DE_Upwithin_12, aoSMCil1b_DEGs_DE_Upwithin_16,
                  aoSMCil1b_DEGs_DE_Upwithin_24, aoSMCil1b_DEGs_DE_Upwithin_36, aoSMCil1b_DEGs_DE_Upwithin_48)

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
                                           results_list[[8]][1]/results_list[[8]][2]*100,
                                           results_list[[9]][1]/results_list[[9]][2]*100),
                           
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, 
                                           results_list[[2]][3]/results_list[[2]][4]*100, 
                                           results_list[[3]][3]/results_list[[3]][4]*100,
                                           results_list[[4]][3]/results_list[[4]][4]*100,
                                           results_list[[5]][3]/results_list[[5]][4]*100,
                                           results_list[[6]][3]/results_list[[6]][4]*100,
                                           results_list[[7]][3]/results_list[[7]][4]*100,
                                           results_list[[8]][3]/results_list[[8]][4]*100,
                                           results_list[[9]][1]/results_list[[9]][2]*100))

DEL_PCG_type$RegulationBegins <- as.factor(c("0-2hrs", "2-4hrs", "4-6hrs", "6-8hrs", "8-12hrs", "12-16hrs", "16-24hrs", "24-36hrs", "36-48hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1], results_list[[4]][1],
                        results_list[[5]][1], results_list[[6]][1], results_list[[7]][1], results_list[[8]][1], results_list[[9]][1])

g2bi <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "steelblue") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  ylab("") +   xlab("") +
  
  theme_minimal()

DEL_PCG_type$p <- unlist(results_list)[c(5,11,17,23,29,35,41,47,53)]
DEL_PCG_type$OR <- unlist(results_list)[c(6,12,18,24,30,36,36,48,54)]
#write.csv(DEL_PCG_type, "aoSMCil1b_inducedDEL_DEP_DEneighboursame_conc.csv", row.names = F)

p.adjust(DEL_PCG_type$p, method = "bonferroni")


DEL_PCG_type_plot <- unique(melt(DEL_PCG_type[,1:3]))

DEL_PCG_type_plot$variable <- as.character(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable[grepl("^DEL", DEL_PCG_type_plot$variable)] <- "Induced lncRNAs"
DEL_PCG_type_plot$variable[grepl("DEP", DEL_PCG_type_plot$variable)] <- "Induced PCGs"

DEL_PCG_type_plot$variable <- as.factor(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable <- factor(DEL_PCG_type_plot$variable, levels = levels(DEL_PCG_type_plot$variable)[c(1,2)])

DEL_PCG_type_plot$RegulationBegins <- as.factor(DEL_PCG_type_plot$RegulationBegins)
DEL_PCG_type_plot$RegulationBegins <- factor(DEL_PCG_type_plot$RegulationBegins, 
                                             levels = levels(DEL_PCG_type_plot$RegulationBegins)[c(1,4,7,8,9,2,3,5,6)])

ggplot(DEL_PCG_type_plot) + aes(x = RegulationBegins, y = value, fill = variable) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.7), color = "grey48") +
  theme_minimal() +
  theme(text = element_text(size = 24)) + #, legend.position = "none") +
  scale_fill_manual(values = c(`Induced lncRNAs` = "#D6484D", `Induced PCGs` = "grey80")) +
  xlab("") +
  ylab("% With induced\nneighbour") + Seurat::RotatedAxis()


#
#### GO/KEGG/etc terms early-induced targets ####

library(clusterProfiler)
library(org.Hs.eg.db)

aoSMCil1b_DE_PCGs_DE_Upwithin_2 <- filter(aoSMCil1b_DEGs_DE_Upwithin_2, grepl("coding_mRNA", CAT_geneClassII))

EarlyInduced_CClncRNAs_CoUp_Targets_aoSMCil1b <- filter(CoRegPairs_sameTimeframe_aoSMCil1b, AnchorGene %in% aoSMCil1b_DEGs_DE_Upwithin_2$CAT_geneID,
                                                    geneID %in% aoSMCil1b_DEGs_DE_Upwithin_2$CAT_geneID)

#targets of early lncs vs other DEGs in that timeframe
CoUp_DE <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", EarlyInduced_CClncRNAs_CoUp_Targets_aoSMCil1b$geneID)),
                    universe      = gsub("\\.[0-9]*", "", aoSMCil1b_DE_PCGs_DE_Upwithin_2$CAT_geneID),
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
aoSMCil1b_02_exprs_PCG <- filter(aoSMCil1b_CPM, geneID %in% filter(GeneBiotypes, CAT_geneClassII == "coding_mRNA")$CAT_geneID, 
                             (Hour0_meanCPM>1 | Hour2_meanCPM>1))

CoUp_DE2 <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", EarlyInduced_CClncRNAs_CoUp_Targets_aoSMCil1b$geneID)),
                     universe      = gsub("\\.[0-9]*", "", aoSMCil1b_02_exprs_PCG$geneID),
                     keyType       = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "all",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
CoUp_DE2_df <- as.data.frame(CoUp_DE2)
#much better
#endothelial/epithelial activation terms, immune activation terms

dotplot(simplify(CoUp_DE2))


#### specific target lists early-induced targets ####

#IEGs (inc. lots of immune)
IEGs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Gene lists/arner_2015_table_S5_IEGs_lit.csv")
IEGs_hs <- filter(IEGs, !Hs_symbol %in% c("", "-"))
IEGs_hs <- filter(IEGs_hs, !grepl("[a-z]", IEGs_hs$Total_count))
IEGs_hs$Hs_symbol[IEGs_hs$Hs_symbol == "IL8"] <- "CXCL8"

#-biased
FANT_S10 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/nature21374-s2/table 10.csv")
FANT_S10_macroMono <- filter(FANT_S10, grepl("macrop|monocy", sample_ontology_term))
FANT_S10_macroMono_G <- gsub("\\.[0-9]*", "", unlist(strsplit(FANT_S10_macroMono$associated_geneID, ",")))

#Epigenetic modifiers (LISA has some definitions but not well described, CRdb seems recent and fine)
CRdb_dataB <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CRdb Data Browse.csv")

TF_Lambert2018 <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/TF_Lambert2018.txt", header = T, stringsAsFactors = F)

#add HMGA2 and FOXL1 to this list, unclear why not in there...

GeneLists <- list("TFs"= unique(filter(aoSMCil1b_DEGs_DE, gsub("\\.[0-9]*", "", CAT_geneID) %in% TF_Lambert2018$ID)$CAT_geneID),
                  "CRs" = unique(filter(aoSMCil1b_DEGs_DE, CAT_geneName %in% c(CRdb_dataB$CR, "FOXL1", "HMGA2"))$CAT_geneID),
                  "IEGs"= unique(filter(aoSMCil1b_DEGs_DE, CAT_geneName %in% c(IEGs_hs[,4]))$CAT_geneID),
                  "MacMon"= unique(filter(aoSMCil1b_DEGs_DE, gsub("\\.[0-9]*", "", CAT_geneID) %in% FANT_S10_macroMono_G)$CAT_geneID))

results_list <- list()

#same timeframe
aoSMCil1b_DE_PCGs_DE_Upwithin_2 <- filter(aoSMCil1b_DEGs_DE_Upwithin_2, grepl("coding_mRNA", CAT_geneClassII))

EarlyInduced_CClncRNAs_CoUp_Targets_aoSMCil1b <- filter(CoRegPairs_sameTimeframe_aoSMCil1b, AnchorGene %in% aoSMCil1b_DEGs_DE_Upwithin_2$CAT_geneID,
                                                    geneID %in% aoSMCil1b_DEGs_DE_Upwithin_2$CAT_geneID)

for (i in 1:length(GeneLists)){
  
  a <- length(unique(EarlyInduced_CClncRNAs_CoUp_Targets_aoSMCil1b$geneID[ EarlyInduced_CClncRNAs_CoUp_Targets_aoSMCil1b$geneID %in% GeneLists[[i]] ]))
  b <- length(unique(EarlyInduced_CClncRNAs_CoUp_Targets_aoSMCil1b$geneID))
  c <- length(unique(aoSMCil1b_DE_PCGs_DE_Upwithin_2$CAT_geneID[ aoSMCil1b_DE_PCGs_DE_Upwithin_2$CAT_geneID %in% GeneLists[[i]]]))
  d <- length(unique(aoSMCil1b_DE_PCGs_DE_Upwithin_2$CAT_geneID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate,
                                  paste(unique(EarlyInduced_CClncRNAs_CoUp_Targets_aoSMCil1b$geneID[ 
                                    EarlyInduced_CClncRNAs_CoUp_Targets_aoSMCil1b$geneID %in% GeneLists[[i]] 
                                  ]), collapse = "/"))
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:8] <- c("pval", "OR", "Genes")
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")
#some enrichment of IEG, others not bothered
GeneLists_enrichedSameT <- GeneLists_enriched

filter(GeneBiotypes, CAT_geneID %in% unique(EarlyInduced_CClncRNAs_CoUp_Targets_aoSMCil1b$geneID[ 
  EarlyInduced_CClncRNAs_CoUp_Targets_aoSMCil1b$geneID %in% GeneLists[["IEGs"]] 
]))[,1:3]


#
#### examine early lncs/targets ####

#add gene names
trial <- merge(CoRegPairs_sameTimeframe_aoSMCil1b, GeneBiotypes[,1:3], by.x = "AnchorGene", by.y = "CAT_geneID")
trial <- merge(trial, GeneBiotypes[,1:3], by.x = "geneID", by.y = "CAT_geneID")

CoRegPairs_sameTimeframe_aoSMCil1b_early <- filter(trial, AnchorGene %in% aoSMCil1b_DEGs_DE_Upwithin_2$CAT_geneID,
                                               geneID %in% aoSMCil1b_DEGs_DE_Upwithin_2$CAT_geneID)

#notable targets: few/unclear
#notable early lncs: few/unclear

#control cis lncs:
ControlCisLncs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Write-ups/supplement/ControlCisLncs.csv", header = T)

#0x identified in the spike of cis-acting lnc evidence
CoRegPairs_sameTimeframe_aoSMCil1b_early_controlCisLncs <- filter(CoRegPairs_sameTimeframe_aoSMCil1b_early, 
                                                              gsub("\\.[0-9]*", "", AnchorGene) %in% ControlCisLncs$Ens_ID)

#up to 4hrs:
CoRegPairs_sameTimeframe_aoSMCil1b_early2 <- filter(trial, AnchorGene %in% c(aoSMCil1b_DEGs_DE_Upwithin_2$CAT_geneID, aoSMCil1b_DEGs_DE_Upwithin_4$CAT_geneID),
                                                geneID %in% c(aoSMCil1b_DEGs_DE_Upwithin_2$CAT_geneID, aoSMCil1b_DEGs_DE_Upwithin_4$CAT_geneID))

#nothing extra is obvious, maybe TNFAIP3 showing up all the time... otherwise v distinct gene set

