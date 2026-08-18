#obtain CAGE sites for timecourse lncRNAs for
#a) confirming annotations
#b) enhancer annotation
#c) novel lncRNAs in FANTOM if needed
#d) comparing to their timecourse data
#e) more accurate definitions of neighbours (TSS-TSS distance)

#libraries
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(reshape2)

#starting data, fpkm_allG with artefacts removed:
fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026filt.csv", header = T, stringsAsFactors = F)
fpkm_PLARG <- filter(fpkm_allG, V55 == "Bona fide lncRNA")
length(unique(fpkm_PLARG$EnsID))#597 unique entries
length(unique(fpkm_PLARG$MSTRG_Tx_ID))#1575 unique entries

#finding CAGE clusters for 1575 transcripts linked to 597 genes (all newly assembled and expressed lncRNAs)

#isolate tx co-ordinates
novel_overlay <- fpkm_PLARG
trial <- novel_overlay[,c(47:48)]
trialii <- sapply(as.list(trial[,2]), strsplit, ":")
trial[,3] <- sapply(trialii, "[[", 1)
trial[,4] <- sapply(trialii, "[[", 2)
trialii <- sapply(as.list(trial[,4]), strsplit, "-")
trial[,5] <- sapply(trialii, "[[", 1)
trial[,6] <- sapply(trialii, "[[", 2)
trialii <- sapply(as.list(trial[,6]), strsplit, " ")
trial[,7] <- sapply(trialii, "[[", 1)
lnc_Coords <- data.frame("chr" = trial[,3], "start" = as.numeric(trial[,5]), "stop" = as.numeric(trial[,7]), "MSTRG_Tx_ID" = trial[,1], stringsAsFactors = F)

#retrieve genomic co-ordinates of areas 1kbp around 5' range, for use in determining surrounding CAGE clusters
lnc5range_hSVSMC <- unique(lnc_Coords)
lnc5range_hSVSMC[,5] <- 1000 #arbitrary number to mean table is .bed format
lnc5range_hSVSMC <- merge(lnc5range_hSVSMC, fpkm_PLARG[,c(8,47)], by = "MSTRG_Tx_ID")#insert strand info
lnc5range_hSVSMC <- lnc5range_hSVSMC[,c(2:4,1,5,6)] #proper .bed format order now
lnc5range_hSVSMC[lnc5range_hSVSMC$str == "+",2] <- lnc5range_hSVSMC[lnc5range_hSVSMC$str == "+",2] -1000
lnc5range_hSVSMC[lnc5range_hSVSMC$str == "+",3] <- lnc5range_hSVSMC[lnc5range_hSVSMC$str == "+",2] +2000
lnc5range_hSVSMC[lnc5range_hSVSMC$str == "-",2] <- lnc5range_hSVSMC[lnc5range_hSVSMC$str == "-",3] -1000
lnc5range_hSVSMC[lnc5range_hSVSMC$str == "-",3] <- lnc5range_hSVSMC[lnc5range_hSVSMC$str == "-",3] +1000

#take this object to bed tools to find CAGE sites overlapping the area 1kbp around the 5' end (same strand overlaps only)
#write.table(lnc5range_hSVSMC, "\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/allexprslnc5range_hSVSMCtime_2026.bed", sep = "\t", quote = F, row.names = F)

FANTOM_CAGE2allLncs_hSVSMC <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/bedtools_error_nonMergeTime_Fantom", sep = "\t", stringsAsFactors = F, header = F)
length(unique(FANTOM_CAGE2allLncs_hSVSMC$V13))#1485 tx have a CAGE
length(unique(FANTOM_CAGE2allLncs_hSVSMC$V4))#9183 CAGEs

trial <- merge(fpkm_PLARG[,c(47,49,50)], FANTOM_CAGE2allLncs_hSVSMC[,c(1:8,11,13)], by.x = "MSTRG_Tx_ID", by.y = "V13", all.x = T)

FANTOM_CAGE2allLncs <- unique(trial)

colnames(FANTOM_CAGE2allLncs) <- c("MSTRG_Tx_ID", "Exon_count", "Spliced_length", "CTSSchr", "CTSSstart", "CTSSstop", "CTSSname", "CTSSexp", "CTSSstrand",
                                   "CTSSpeakstart", "CTSSpeakstop", "5pStringtie")
FANTOM_CAGEclusterCounts <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/FANTOM_CAT.source_CAGE_cluster.count.tsv", sep = "\t")

trial <- FANTOM_CAGE2allLncs
triali <- merge(trial, FANTOM_CAGEclusterCounts[,c(1,3,6)], by.x = "CTSSname", by.y = "DPIClstrID", all.x = T)
triali$CTSSexp <- triali$clstrCount_pm
FANTOM_CAGE2allLncs <- triali[,c(2:7,1,8:12)]

########determine TIEScore###########

#adjust to get 5' end of tx according to StringTie
FANTOM_CAGE2allLncs[,12] <- FANTOM_CAGE2allLncs[,12] +1000 #not strand sensitive

#calculate TIEScore based on FANTOM CAT supplement:
#weighted length score
FANTOM_CAGE2allLncs[,13] <- 0.1*(5.45*(log2(FANTOM_CAGE2allLncs$Exon_count * FANTOM_CAGE2allLncs$Spliced_length)) + 33.25)

#weighted distance score, n.b. not clear whether they used peak or nearest part of cluster
FANTOM_CAGE2allLncs[,14] <- FANTOM_CAGE2allLncs$`5pStringtie` - FANTOM_CAGE2allLncs$CTSSpeakstart #expecting positive or small neg numbers for + stranded tx
#liftOver alterations means that CTSS >1000bp away are returned, get rid of silly ones that creep in:
FANTOM_CAGE2allLncs <- filter(FANTOM_CAGE2allLncs, is.na(V14) | (V14 <1050 & V14 >-1050))
FANTOM_CAGE2allLncs[,15] <- sqrt(FANTOM_CAGE2allLncs[,14]*FANTOM_CAGE2allLncs[,14]) #CTSS peaks should be upstream or within exon but not downstream, otherwise implies exon is not validated
FANTOM_CAGE2allLncs[,15] <- 0.55*(-1.19*(sqrt(FANTOM_CAGE2allLncs[,15])) + 85.30)

#weighted expression score
FANTOM_CAGE2allLncs[,16] <- 0.35*(6.63*(log2(FANTOM_CAGE2allLncs$CTSSexp)) + 77.89)

FANTOM_CAGE2allLncs[,17] <- FANTOM_CAGE2allLncs[,13] + FANTOM_CAGE2allLncs[,15] + FANTOM_CAGE2allLncs[,16]

FANTOM_CAGE2allLncs[,18] <- as.factor(FANTOM_CAGE2allLncs$MSTRG_Tx_ID)
colnames(FANTOM_CAGE2allLncs)[13:18] <- c("E", "5p_CTSS_Dist", "D", "C", "TIEScore", "Tx_Factor")

#max TIEScore for each transcript
trial <- split(FANTOM_CAGE2allLncs, FANTOM_CAGE2allLncs$Tx_Factor)
trialii <- lapply(trial, function(x){x$CTSSname[which.max(x$TIEScore)]})
trialiii <- lapply(trial, function(x){max(x$TIEScore)})
trialiv <- lapply(trial, function(x){x$`5p_CTSS_Dist`[which.max(x$TIEScore)]})
trialv <- lapply(trial, function(x){x$CTSSstart[which.max(x$TIEScore)]})

BestCTSS <- unlist(trialii)
BestCTSS <- BestCTSS[!is.na(BestCTSS)]
BEstTIEScore <- unlist(trialiii)
BEstTIEScore <- BEstTIEScore[!is.na(BEstTIEScore)]
BEstDist <- unlist(trialiv)
BEstDist <- BEstDist[!is.na(BEstDist)]
BestStart <- unlist(trialv)
BestStart <- BestStart[!is.na(BestStart)]

Best <- data.frame("MSTRG_Tx_ID" = names(BestCTSS), BestCTSS, BEstTIEScore, BEstDist, BestStart)

allLncs_BestCAGE <- merge(unique(FANTOM_CAGE2allLncs[,c(1:3,9,12)]), Best, by.x = "MSTRG_Tx_ID", all.x = T)

summary(allLncs_BestCAGE$BEstTIEScore)

#numbers so far
length(unique(allLncs_BestCAGE$MSTRG_Tx_ID))#1574 lost 1 somewhere
length(unique(allLncs_BestCAGE$MSTRG_Tx_ID[is.na(allLncs_BestCAGE$BestCTSS)]))#90 of these, no CTSS in 1000bp
#easier ifget rid of NAs
allLncs_BestCAGE <- filter(allLncs_BestCAGE, !is.na(CTSSstrand))

#get rid of those with silly exon 1 changes:
#get first exon length per transcript
hSVSMC_annot <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/PLAR_4timepoints/Timecourse3PLAR.annot.bed", header=F, sep = "\t", stringsAsFactors = F)
trial <- sapply(hSVSMC_annot$V11[hSVSMC_annot$V6 == "+"], strsplit, ",")
trial <- sapply(trial, "[[", 1)
hSVSMC_annot[hSVSMC_annot$V6 == "+",13] <- unlist(trial)
trial <- sapply(hSVSMC_annot$V11[hSVSMC_annot$V6 == "-"], strsplit, ",")
trial <- sapply(trial, function(x){x[length(x)]})
hSVSMC_annot[hSVSMC_annot$V6 == "-",13] <- unlist(trial)

trial <- sapply(hSVSMC_annot$V4, strsplit, ":")
trial <- sapply(trial, "[[", 1)
trial <- sapply(trial, strsplit, "\\|")
trial <- sapply(trial, function(x){x[length(x)]})
hSVSMC_annot[,14] <- unlist(trial)

allLncs_BestCAGE <- unique(merge(allLncs_BestCAGE, hSVSMC_annot[,13:14], by.x = "MSTRG_Tx_ID", by.y = "V14", all.x = T))
allLncs_BestCAGE[,10] <- as.numeric(allLncs_BestCAGE[,10])
class(allLncs_BestCAGE[,10])

#first exons ~86-302bp IQR based on above annotations
#for pos strand tx, CTSS peak should be within 5p + exon 1
allLncs_BestCAGE[allLncs_BestCAGE$CTSSstrand == "+",11] <- ((as.numeric(allLncs_BestCAGE[allLncs_BestCAGE$CTSSstrand == "+",10]) + 
                                                               allLncs_BestCAGE[allLncs_BestCAGE$CTSSstrand == "+",5])
                                                            > allLncs_BestCAGE[allLncs_BestCAGE$CTSSstrand == "+",9])#if exon + 5p is larger than CTSS peak (TRUE) then makes sense as TSS is before end of first exon

allLncs_BestCAGE[allLncs_BestCAGE$CTSSstrand == "-",11] <- ((allLncs_BestCAGE[allLncs_BestCAGE$CTSSstrand == "-",5] - 
                                                               as.numeric(allLncs_BestCAGE[allLncs_BestCAGE$CTSSstrand == "-",10]))
                                                            < allLncs_BestCAGE[allLncs_BestCAGE$CTSSstrand == "-",9])#if 5p - exon 1 is smaller than CTSS peak (TRUE) then makes sense as TSS is after end of first exon                   
table(allLncs_BestCAGE$V11)
#4% false links, way lower:
boxplot(allLncs_BestCAGE$BEstTIEScore~allLncs_BestCAGE$V11)

#length of distance between CTSS peak and first splice site
allLncs_BestCAGE[allLncs_BestCAGE$CTSSstrand == "+",11] <- ((as.numeric(allLncs_BestCAGE[allLncs_BestCAGE$CTSSstrand == "+",10]) + 
                                                               allLncs_BestCAGE[allLncs_BestCAGE$CTSSstrand == "+",5])
                                                            - allLncs_BestCAGE[allLncs_BestCAGE$CTSSstrand == "+",9])# -ve number indicates CTSS starts after first splice site

allLncs_BestCAGE[allLncs_BestCAGE$CTSSstrand == "-",11] <- (-(allLncs_BestCAGE[allLncs_BestCAGE$CTSSstrand == "-",5] - 
                                                                as.numeric(allLncs_BestCAGE[allLncs_BestCAGE$CTSSstrand == "-",10])) +
                                                              allLncs_BestCAGE[allLncs_BestCAGE$CTSSstrand == "-",9])# -ve number indicates CTSS starts after first splice site

allLncs_BestCAGE[,12] <- allLncs_BestCAGE$V11/as.numeric(allLncs_BestCAGE$V13)
colnames(allLncs_BestCAGE)[11:12] <- c("Distance Upstream of SS1", "PercentExon1Retained")

allLncs_BestCAGE <- merge(fpkm_PLARG[,c(2,47)], allLncs_BestCAGE, by = "MSTRG_Tx_ID", all = T)

#shows vast majority of linked CAGE sites identified are pretty bang on 
#even without cutting anything with low TIEScore, retain close to 100% of exon (removes unmapped and definitely false CTSSs)
trial <- allLncs_BestCAGE
trial$type[grepl("MSTRG", trial$EnsID)] <- "Newly-assembled"
trial$type[!grepl("MSTRG", trial$EnsID)] <- "GENCODEv26"

summary(100-(100*trial$PercentExon1Retained))

ggplot(trial, aes(x = type, y = (100-(100*PercentExon1Retained)))) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.title.y = element_blank()) +
  ylab("First exon resize after CAGEseq (%)") +
  scale_y_continuous(limits= c(-200,200), breaks= c(-200,-100,0,100,200)) + 
  coord_flip()

ggplot(trial, aes(x = type, y = `Distance Upstream of SS1`)) +
  geom_boxplot() +
  theme_minimal() 

#######Filter out weaker TIEScores############

#FANTOM genes already tethered to many of these sites
#expression info available for all but RAMPAGE validation only for permissive/robust/stringent boundaries set by TIEScore
conversion <- read.table("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/FANTOM_CAT.lv1_raw.info_table.ID_mapping.tsv", header=T, sep = "\t", stringsAsFactors = F)
allFANTOMG <- read.table("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/allTx", header=T, sep = "\t", stringsAsFactors = F)
robFANTOMG <- read.table("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/robustTx", header=T, sep = "\t", stringsAsFactors = F)
strinFANTOMG <- read.table("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/stringentTx", header=T, sep = "\t", stringsAsFactors = F)

allLncs_BestCAGE <- unique(merge(allLncs_BestCAGE, conversion[,c(1,3)], by.x = "BestCTSS", by.y = "CAGEClusterID", all.x = T))
colnames(allLncs_BestCAGE)[14] <- "FANTOM_ID"
length(unique(allLncs_BestCAGE$BestCTSS))#866 CTSS for
length(unique(allLncs_BestCAGE$MSTRG_Tx_ID))#1575 transcripts (some have none here)

allLncs_BestCAGE[is.na(allLncs_BestCAGE$BestCTSS),15] <- "Absent"
allLncs_BestCAGE[is.na(allLncs_BestCAGE$V15),15] <- "NoFANTLinkedG"
allLncs_BestCAGE[allLncs_BestCAGE$FANTOM_ID %in% allFANTOMG$geneID,15] <- "Raw"
allLncs_BestCAGE[allLncs_BestCAGE$FANTOM_ID %in% robFANTOMG$geneID,15] <- "Robust"
allLncs_BestCAGE[allLncs_BestCAGE$FANTOM_ID %in% strinFANTOMG$geneID,15] <- "Stringent"

#write.csv(allLncs_BestCAGE, "allLncs_FANTOMCAGE_3PLARtimecourse_2026.csv", row.names = F)

allLncs_BestCAGE_NameFinder <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/FANTOMNameFinder.csv", stringsAsFactors = F,
                                          header = F, sep = "\t")
head(allLncs_BestCAGE_NameFinder)

trial <- sapply(allLncs_BestCAGE_NameFinder$V1, strsplit, ";")
allLncs_BestCAGE_NameFinder$V1 <- unlist(trial)
trial <- sapply(allLncs_BestCAGE_NameFinder$V2, strsplit, ";")
allLncs_BestCAGE_NameFinder$V2 <- unlist(trial)

allLncs_BestCAGE <- merge(allLncs_BestCAGE, allLncs_BestCAGE_NameFinder, by.x = "FANTOM_ID", by.y = "V1", all.x = T)
allLncs_BestCAGE <- allLncs_BestCAGE[,c(2:15,1,16)]

length(unique(allLncs_BestCAGE$EnsID))#597 lncRNAs genes
length(unique(allLncs_BestCAGE$MSTRG_Tx_ID))#1575 lncRNA tx
#as expected

table(unique(allLncs_BestCAGE[,c(2,14)])$V15)

#some tx with CAGEs with 2x classes:
#e.g:
trial <- unique(allLncs_BestCAGE[,c(2,14)])
trial$MSTRG_Tx_ID[duplicated(trial$MSTRG_Tx_ID)]

#select CAGE validated according to parameters set in FANTOM atlas and IJMS paper
allLncs_BestCAGE_Valid <- filter(allLncs_BestCAGE, BEstTIEScore >60)
length(unique(allLncs_BestCAGE_Valid$EnsID)) #529 genes stringent

#if the CAGE is within the first exon it is bit problematic
allLncs_BestCAGE_Valid <- filter(allLncs_BestCAGE_Valid, PercentExon1Retained >0.1)
length(unique(allLncs_BestCAGE_Valid$EnsID)) #526 genes

allLncs_BestCAGE$Valid <- "No CAGEseq Match"
allLncs_BestCAGE$Valid[allLncs_BestCAGE$EnsID %in% allLncs_BestCAGE_Valid$EnsID] <- "CAGEseq Match"
allLncs_BestCAGE$type <- "GENCODE"
allLncs_BestCAGE$type[grepl("MSTRG", allLncs_BestCAGE$EnsID)] <- "Newly-assembled"
allLncs_BestCAGE <- allLncs_BestCAGE[,-16]

#write.csv(allLncs_BestCAGE, "allLncs_BestCAGE_2026.csv", row.names = F)
#write.csv(allLncs_BestCAGE_Valid, "allLncs_BestCAGE_Valid_2026.csv", row.names = F)
