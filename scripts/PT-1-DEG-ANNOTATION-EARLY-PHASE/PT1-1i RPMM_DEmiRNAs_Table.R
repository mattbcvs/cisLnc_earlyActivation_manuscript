# shortstack general outputs:
library(DESeq2)
library(dplyr)
library(ggplot2)

Counts <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker_group-bioinformatics/sRNAseq_HSVMC_IP_TimeCourse_Effie_analysis/shortstack_miRNAs_HSVSMC_sRNAseq/Counts.txt", 
                     sep = "\t", header = T)

#results has info on miR calling, main sequence etc
Results <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker_group-bioinformatics/sRNAseq_HSVMC_IP_TimeCourse_Effie_analysis/shortstack_miRNAs_HSVSMC_sRNAseq/Results.txt", 
                      sep = "\t", header = T)
# Counts txt file
countdata <- read.csv( "\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker_group-bioinformatics/sRNAseq_HSVMC_IP_TimeCourse_Effie_analysis/All_Patients.txt", header=TRUE )
countdataframe <- DataFrame(countdata)

#  Metadata (metadata.txt file attached)
metaData <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker_group-bioinformatics/sRNAseq_HSVMC_IP_TimeCourse_Effie_analysis/metadata.txt", header = TRUE, stringsAsFactors = T)
metadataframe <- matrix(metaData)

# Setting samples labeled as 'control' in the metadata file as the reference level:
metaData$condition <- relevel(metaData$condition, "control") 
metadataframe <- DataFrame(metaData)

#metadataframe

#LRT including patient

# DE expression 
dds_m <- DESeqDataSetFromMatrix(countData=countdataframe, colData=metaData, design=~patient+condition, tidy = TRUE)###added patient into model
dds_m <-dds_m[rowSums(counts(dds_m))>10, ]
dds <- DESeq(dds_m, test="LRT", full=~patient+condition, reduced=~patient) ####reduced model contains just patient

res <- results(dds)

res.table <- as.data.frame(res)#733 miRs with a p value (with the model before there was 1416)
res.tableDE <- filter(res.table, padj <0.05) #30 miRs with a significant p value instead of 6, get alot more significance cos takes into account patient variability


rld_data <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(rld_data, intgroup = c("condition", "patient"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = patient)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) #+coord_fixed()

colnames(pcaData)[c(4,5)] <- c("Timepoint", "Patient")

ggplot(pcaData, aes(x = PC1, y = PC2, shape = Patient, color = Timepoint)) +
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

assay(rld_data) <- limma::removeBatchEffect(assay(rld_data), rld_data$patient)

pcaData <- plotPCA(rld_data, intgroup = c("patient", "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

colnames(pcaData)[c(4,5)] <- c("Patient", "Timepoint")

pcaData$Timepoint <- rep(c("0hr", "4hr", "8hr", "24hr"), 4)

pcaData$Timepoint <- factor(pcaData$Timepoint)
pcaData$Timepoint <- factor(pcaData$Timepoint, levels = levels(pcaData$Timepoint)[c(1,3,4,2)])

ggplot(pcaData, aes(x = PC1, y = PC2, color = Timepoint, shape = patient)) +
  geom_point(size =3) +
  xlab(paste0("\nPC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance\n")) +
  #coord_fixed() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size =16.5),
        axis.title.x = element_text(size = 16.5),
        legend.text = element_text(size =15),
        legend.title = element_text(size = 16.5))


#### Results from DESeq2 ####

dds_pairs <- DESeq(dds_m) ####reduced model contains just patient

res0h4h_contrast <- results(dds_pairs, contrast = c("condition", "group4h", "control"))
res0h4h.table <- as.data.frame(res0h4h_contrast)

res4h8h_contrast <- results(dds_pairs, contrast = c("condition", "group8h", "group4h"))
res4h8h.table <- as.data.frame(res4h8h_contrast)

res8h24h_contrast <- results(dds_pairs, contrast = c("condition", "group24h", "group8h"))
res8h24h.table <- as.data.frame(res8h24h_contrast)

res0h8h_contrast <- results(dds_pairs, contrast = c("condition", "group8h", "control"))
res0h8h.table <- as.data.frame(res0h8h_contrast)

res0h24h_contrast <- results(dds_pairs, contrast = c("condition", "group24h", "control"))
res0h24h.table <- as.data.frame(res0h24h_contrast)

res4h24h_contrast <- results(dds_pairs, contrast = c("condition", "group24h", "group4h"))
res4h24h.table <- as.data.frame(res4h24h_contrast)

resAllPairs <- as.data.frame(rbind(res0h4h.table, res4h8h.table, res8h24h.table, res0h8h.table, res0h24h.table, res4h24h.table))

resAllPairs$preadj <- p.adjust(resAllPairs$pvalue, method = "BH")

#recapitulate NA values from padj:
resAllPairs$preadj[is.na(resAllPairs$padj)] <- NA

#assemble table
res.table$LogFC_0_4 <- res0h4h.table$log2FoldChange
res.table$preadj_0_4 <- resAllPairs$preadj[1:1107]

res.table$LogFC_4_8 <- res4h8h.table$log2FoldChange
res.table$preadj_4_8 <- resAllPairs$preadj[1108:(1107*2)]

res.table$LogFC_8_24 <- res8h24h.table$log2FoldChange
res.table$preadj_8_24 <- resAllPairs$preadj[2215:(1107*3)]

res.table$LogFC_0_8 <- res0h8h.table$log2FoldChange
res.table$preadj_0_8 <- resAllPairs$preadj[3322:(1107*4)]

res.table$LogFC_0_24 <- res0h24h.table$log2FoldChange
res.table$preadj_0_24 <- resAllPairs$preadj[4429:(1107*5)]

res.table$LogFC_4_24 <- res4h24h.table$log2FoldChange
res.table$preadj_4_24 <- (resAllPairs$preadj[5536:(1107*6)])

miRTimecourse_DE <- data.frame("ENSEMBL" = rownames(res.table), res.table[,6:18])



#results from shortstack
shortstacknames <- list.files( "\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker_group-bioinformatics/Effie/shortstack_miRNAs_HSVSMC_sRNAseq/shortstack_conditions/")
shortstackfiles <- paste(list.files( "\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker_group-bioinformatics/Effie/shortstack_miRNAs_HSVSMC_sRNAseq/shortstack_conditions/", full.names = T),
                         "/Results.txt", sep = "")

shortstackcounts <- paste(list.files( "\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker_group-bioinformatics/Effie/shortstack_miRNAs_HSVSMC_sRNAseq/shortstack_conditions/", full.names = T),
                          "/Counts.txt", sep = "")
#counts table to calculate rpm using mapped reads rather than all reads (as default in shortstack)
#should be same as if calculating "rpmm" (there is an option for this in shortstack)
miRs<-read.table(shortstackcounts[1],header=TRUE,sep="\t",stringsAsFactors = FALSE)[,1:2]
counts <-do.call(cbind,lapply(shortstackcounts,function(fn)read.table(fn,header=TRUE,sep="\t",stringsAsFactors = FALSE)[,4]))
#function to make RPMM:
rpm <- sapply(as.data.frame(counts), function(x){x/sum(x)*1000000})
rpm <- data.frame(miRs,rpm,stringsAsFactors = FALSE)
colnames(rpm)<-c("locus","name", shortstacknames)
rpm <- rpm[,c(1:2, (c(1,5,9,13,3,7,11,15,4,8,12,16,2,6,10,14)+2))]

rpm_list_ctrl <- as.list(as.data.frame(t(rpm[,3:6])))
rpm_list_pd <- as.list(as.data.frame(t(rpm[,7:10])))
rpm_list_il <- as.list(as.data.frame(t(rpm[,11:14])))
rpm_list_bo <- as.list(as.data.frame(t(rpm[,15:18])))

rpm_mean_ctrl <- sapply(rpm_list_ctrl, mean)
rpm_mean_pd <- sapply(rpm_list_pd, mean)
rpm_mean_il <- sapply(rpm_list_il, mean)
rpm_mean_bo <- sapply(rpm_list_bo, mean)

std <- function(x) sd(x)/sqrt(length(x))

rpm_se_basal <- sapply(rpm_list_ctrl, std)
rpm_se_chol <- sapply(rpm_list_pd, std)
rpm_se_mig <- sapply(rpm_list_il, std)
rpm_se_pdgf <- sapply(rpm_list_bo, std)

rpm_mean_treatment <- data.frame("Hour0_meanrpm" = rpm_mean_ctrl, "Hour0_serpm" = rpm_se_basal,
                                 "Hour4_meanrpm" = rpm_mean_pd, "Hour4_serpm" = rpm_se_chol,
                                 "Hour8_meanrpm" = rpm_mean_il, "Hour8_serpm" = rpm_se_mig,
                                 "Hour24_meanrpm" = rpm_mean_bo, "Hour24_serpm" = rpm_se_pdgf)

trial <- as.list(as.data.frame(t(rpm_mean_treatment)))
rpm_max_treatment <- as.numeric(sapply(trial, max))

rpm <- cbind(rpm, rpm_mean_treatment, rpm_max_treatment)

rpm$mergeCol <- paste(rpm$locus, rpm$name, sep = "+")

miRTimecourse_DEnonDE_rpm <- merge(rpm, miRTimecourse_DE, by.x = "mergeCol", by.y = "ENSEMBL", all.x = T)


#### miRBase annotation########
#correct strands, add in a primary transcript column, remove/deal with affected miRs
miRBaseGFF3 <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/hsa.gff3", sep = "\t", skip = 13, head = F)#v21 downloaded 3/3/25

miRBaseGFF3pri <- filter(miRBaseGFF3, !V3 == "miRNA")
miRBaseGFF3 <- filter(miRBaseGFF3, V3 == "miRNA")
miRBaseGFF3$miRname <- sapply(sapply(sapply(sapply(miRBaseGFF3$V9, strsplit, "\\;"), "[[", 3), strsplit, "="), "[[", 2)
miRBaseGFF3$primiRID <- sapply(sapply(sapply(sapply(miRBaseGFF3$V9, strsplit, "\\;"), "[[", 4), strsplit, "="), "[[", 2)
miRBaseGFF3pri$primiRID <- sapply(sapply(sapply(sapply(miRBaseGFF3pri$V9, strsplit, "\\;"), "[[", 2), strsplit, "="), "[[", 2)
miRBaseGFF3pri$primiRname <- sapply(sapply(sapply(sapply(miRBaseGFF3pri$V9, strsplit, "\\;"), "[[", 3), strsplit, "="), "[[", 2)

library(miRBaseConverter)

#some example functions, find version:
checkMiRNAVersion(miRBaseGFF3pri$primiRname, verbose = T)#v22 downloaded

#match format to shortstack annotation used
miRBase_annotation <- merge(miRBaseGFF3pri[,c(1,4,5,7,9:11)], miRBaseGFF3[,c(1,4,5,7,9:11)], by = "primiRID")
miRBase_annotation$mergeCol <- paste(miRBase_annotation$V1.y, miRBase_annotation$V4.y, sep = ":")
miRBase_annotation$mergeCol <- paste(miRBase_annotation$mergeCol, miRBase_annotation$V5.y, sep = "-")
miRBase_annotation$mergeCol <- paste(miRBase_annotation$mergeCol, miRBase_annotation$miRname, sep = "+")

#match up
rpm_allmiRs <- filter(miRTimecourse_DEnonDE_rpm, rpm_max_treatment >0)
trial <- merge(miRBase_annotation[,c(1,7,11,13,14)], rpm_allmiRs, by = "mergeCol", all.y = T)
#0 entries are not found in exact same way in miRBase format
unmatched <- filter(trial, is.na(primiRID))

rpm_allmiRs_annotated <- trial


#### annotating guide/pass ####

#calculate mature and *
trial <- split(rpm_allmiRs_annotated, rpm_allmiRs_annotated$primiRID)
table(sapply(trial, function(x){length(x[,1])})) #550 hairpins only 1 arm detected, 601 miRs both arms detected
#select miRs with both arms expressed:
triali <- trial[sapply(trial, function(x){length(x[,1])}) > 1]
#average diff in each sample:
trialii <- lapply(triali, function(x){data.frame(x[which.min(x$rpm_max_treatment),c(1:3,5)],
                                                 "StartvMature" = mean(as.numeric((x[which.min(x$rpm_max_treatment),8:23]/x[which.max(x$rpm_max_treatment),8:23])))*100)})
trialiii <- bind_rows(trialii)
#most pairs have one arm which is on average 8-40% of RPMM:
boxplot(trialiii$StartvMature)

#compare with just rpm_max
trialiimax <- lapply(triali, function(x){data.frame(x[which.min(x$rpm_max_treatment),c(1:3,5)],
                                                    "StartvMatureMax" = (x[which.min(x$rpm_max_treatment),32]/x[which.max(x$rpm_max_treatment),32])*100)})
trialiiimax <- bind_rows(trialiimax)

compareStarArmMeasure <- merge(trialiii, trialiiimax[,c(1,5)], by = "mergeCol")
#check for DE miRs, e.g. miR-222-5p

#as miR-222-5p is of interest, shows not to discount any * arms
rpm_allmiRs_annotated$Stars[rpm_allmiRs_annotated$mergeCol %in% filter(compareStarArmMeasure, StartvMature >= 20)$mergeCol] <- ">20% RPMM"
rpm_allmiRs_annotated$Stars[rpm_allmiRs_annotated$mergeCol %in% filter(compareStarArmMeasure, StartvMature < 20)$mergeCol] <- "<20% RPMM"
rpm_allmiRs_annotated$Stars[rpm_allmiRs_annotated$mergeCol %in% filter(compareStarArmMeasure, StartvMature < 10)$mergeCol] <- "<10% RPMM"
rpm_allmiRs_annotated$Stars[rpm_allmiRs_annotated$mergeCol %in% filter(compareStarArmMeasure, StartvMature < 1)$mergeCol] <- "<1% RPMM"
rpm_allmiRs_annotated$Stars[!rpm_allmiRs_annotated$mergeCol %in% compareStarArmMeasure$mergeCol] <- "Single miR only"
table(rpm_allmiRs_annotated$Stars)
#name with stars will be more informative ID than mergeCol from now on:
rpm_allmiRs_annotated$nameStars <- rpm_allmiRs_annotated$name
rpm_allmiRs_annotated$nameStars[grepl("<20", rpm_allmiRs_annotated$Stars)] <- 
  paste(rpm_allmiRs_annotated$name[grepl("<20", rpm_allmiRs_annotated$Stars)], "*", sep = "")

rpm_allmiRs_annotated$nameStars[grepl("<10", rpm_allmiRs_annotated$Stars)] <- 
  paste(rpm_allmiRs_annotated$name[grepl("<10", rpm_allmiRs_annotated$Stars)], "**", sep = "")

rpm_allmiRs_annotated$nameStars[grepl("<1%", rpm_allmiRs_annotated$Stars)] <- 
  paste(rpm_allmiRs_annotated$name[grepl("<1%", rpm_allmiRs_annotated$Stars)], "***", sep = "")


#### DE miRs ####

colnames(rpm_allmiRs_annotated)
rpm_allmiRs_lrt <- unique(filter(rpm_allmiRs_annotated, padj < 0.05)[,1])

rpm_allmiRs_pairs <-  unique(filter(rpm_allmiRs_annotated, (preadj_0_4 <0.05 & (LogFC_0_4 > log2(1.5) | LogFC_0_4 < -log2(1.5))) |
                                      (preadj_0_8 <0.05 & (LogFC_0_8 > log2(1.5) | LogFC_0_8 < -log2(1.5))) |
                                      (preadj_0_24 <0.05 & (LogFC_0_24 > log2(1.5) | LogFC_0_24 < -log2(1.5))) |
                                      (preadj_4_8 <0.05 & (LogFC_4_8 > log2(1.5) | LogFC_4_8 < -log2(1.5))) |
                                      (preadj_4_24 <0.05 & (LogFC_4_24 > log2(1.5) | LogFC_4_24 < -log2(1.5))) |
                                      (preadj_8_24 <0.05 & (LogFC_8_24 > log2(1.5) | LogFC_8_24 < -log2(1.5))))[,1])

rpm_allmiRs_lrt %in% rpm_allmiRs_pairs
rpm_allmiRs_annotated$DE_LRT <- NA
rpm_allmiRs_annotated$DE_pairs <- NA

rpm_allmiRs_annotated$DE_LRT[rpm_allmiRs_annotated$mergeCol %in% rpm_allmiRs_lrt] <- "DE"
rpm_allmiRs_annotated$DE_pairs[rpm_allmiRs_annotated$mergeCol %in% rpm_allmiRs_pairs] <- "DE"

rpm_allmiRs_annotated$DE_consensus <- "Stable"
rpm_allmiRs_annotated$DE_consensus[!is.na(rpm_allmiRs_annotated$DE_LRT) & !is.na(rpm_allmiRs_annotated$DE_pairs)] <- "DE"

rpm_allmiRs_annotated_DE <- filter(rpm_allmiRs_annotated, DE_consensus == "DE")
#27 in absence of expression filter


#overall filter - per 10 RPMM miRs
rpm_express_miR <- filter(rpm_allmiRs_annotated, rpm_max_treatment >10) #358 over 10 rpm
rpm_express_miR_DE <- filter(rpm_express_miR, DE_consensus == "DE")

22/358 #only 6%

#per 100 RPMM miRs (i.e. are we being too incusive of low expressed miRs?)
rpm_express_miRstrict <- filter(rpm_allmiRs_annotated, rpm_max_treatment >100) #358 over 10 rpm
rpm_express_miRstrict_DE <- filter(rpm_express_miR, DE_consensus == "DE")
7/177 #4%, even less

#good conclusion, even if being strict on expression, v. few expressed miRs are DE


#### final table #### 

#write.csv(rpm_allmiRs_annotated, "rpm_allmiRs_annotated_2025.csv", row.names = F)

rpm_allmiRs_annotated <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/rpm_allmiRs_annotated_2025.csv")


#### expected miRs ####

#list from https://academic.oup.com/cardiovascres/article/114/4/611/4781705 - expected miRs in SMC

SMC_miRs <- c("miR-1-", "miR-133a-", "miR-21-", "miR-23b-", "miR-24-", "miR-27b-", "miR-26a-", "miR-29a-", "miR-29b-", "miR-29c-",
  "miR-34a-", "miR-130a-", "miR-138-", "miR-143-", "miR-145-", "miR-146a-", "miR-155-", "miR-195-", "miR-204-", "miR-205-",
  "miR-206-", "miR-210-", "miR-221-", "miR-222-", "miR-424-", "miR-322-", "miR-663-")

trial <- lapply(SMC_miRs, function(x){
  rpm_express_miR[(grepl(x, rpm_express_miR$mergeCol)),]
                  })

triali <- bind_rows(trial)
length(unique(triali$primiRID))
23/27 #85% of these loci have a mature miR detected >10 RPMM

trial <- lapply(SMC_miRs, function(x){
  rpm_express_miR_DE[(grepl(x, rpm_express_miR_DE$mergeCol)),]
})
triali <- bind_rows(trial)
length(unique(triali$primiRID))
7/27 #26% have a mature miR detected >10 RPMM and DE

unique(triali$primiRname)
unique(triali$nameStars)
#miR-222 only one with 2x arms

#2x other passenger strands in there too

known_SMC_miRs <- triali


#### heatmap ####

batched_rld <- assay(rld_data)
dim(batched_rld)

batched_rld_miRs <- batched_rld[rownames(batched_rld) %in% rpm_express_miR_DE$mergeCol,]

batched_rld_miRs <- as.data.frame(batched_rld_miRs)

batched_rld_miRs <- batched_rld_miRs[,c(1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16)]

mat <- batched_rld_miRs[,1:16]

MatchRowNamesTable <- unique(rpm_express_miR_DE[,c(1,50)])
MatchRowNamesTable <- MatchRowNamesTable[match(rownames(mat), MatchRowNamesTable$mergeCol),]

rownames(mat) <- MatchRowNamesTable[,2]

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


condition<-c(rep("0hr",4),rep("4hr",4),rep("8hr",4),rep("24hr",4))
patient<-rep(c("1","2","3","4"),4)
data_cols<-data.frame(condition=condition,patient=patient)

data_colsHeat <- data.frame("Hours" = data_cols[,1], stringsAsFactors = T)
rownames(data_colsHeat) <- colnames(mat)
data_colsHeat$Hours <- factor(data_colsHeat$Hours, levels(data_colsHeat$Hours)[c(1,3,4,2)])

library(pheatmap)
p <- pheatmap(mat,clustering_method = "complete",annotation_legend = F,
              annotation_col = data_colsHeat,
              show_colnames = F, 
              show_rownames = T, 
              cluster_cols = F,
              cluster_rows = T, gaps_col = c(4,8,12),
              #cutree_rows = 2,
              treeheight_col = 0, 
              treeheight_row = 0,
              legend = F,
              color = myColor, 
              breaks = myBreaks,
              border_color = NA)


known_SMC_miRs$SMC_pheno_marker <- c("Differentiation", "Activation", "Activation", "Activation", NA, "Activation", "Activation", "Activation")

data_rowsHeat <- as.data.frame(known_SMC_miRs[c(1:4,6:8),c(51)])
rownames(data_rowsHeat) <- known_SMC_miRs$nameStars[c(1:4,6:8)]
colnames(data_rowsHeat) <- "Marker in SMC for" 


p <- pheatmap(mat,clustering_method = "complete",annotation_legend = T,
              #annotation_col = data_colsHeat,
              annotation_row = data_rowsHeat,
              show_colnames = F, 
              show_rownames = T, 
              cluster_cols = F,
              cluster_rows = T, gaps_col = c(4,8,12),
              #cutree_rows = 2,
              treeheight_col = 0, 
              treeheight_row = 0,
              legend = T,
              color = myColor, 
              breaks = myBreaks,
              border_color = NA)

#can just plot activating markers..

p <- pheatmap(mat[rownames(mat) %in% rownames(data_rowsHeat)[-1],],
              clustering_method = "complete",annotation_legend = F,
              #annotation_col = data_colsHeat,
              #annotation_row = data_rowsHeat,
              show_colnames = F, 
              show_rownames = T, 
              cluster_cols = F,
              cluster_rows = T, gaps_col = c(4,8,12),
              #cutree_rows = 2,
              treeheight_col = 0, 
              treeheight_row = 0,
              legend = F,
              color = myColor, 
              breaks = myBreaks,
              border_color = NA)
