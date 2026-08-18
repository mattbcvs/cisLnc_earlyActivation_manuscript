library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(ggbeeswarm)

#### (skip on revisit) annotate extra biotype detail ####

fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026.csv", header = T)
table(unique(fpkm_allG[,c(2,55)])$V55)
length(unique(fpkm_allG$EnsID))#14235 expressed genes, 847 bona fide lncs, 11887 PCGs

#add in TFs from lambert:
TF_Lambert2018 <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/TF_Lambert2018.txt", header = T, stringsAsFactors = F)
fpkm_allG$TF[gsub("\\.[0-9]*", "", fpkm_allG$EnsID) %in% TF_Lambert2018$ID & !grepl("artefacts|Too low", fpkm_allG$V55)] <- "TF"
table(filter(fpkm_allG, !is.na(TF))$V55)
#all PCGs

#add in cell cycle genes:
CC_Freeman <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/FreemanCellCycle.csv", header = T)[,1:35]
fpkm_allG$CC[gsub("\\.[0-9]*", "", fpkm_allG$EnsID) %in% CC_Freeman$Ensembl_ID & !grepl("artefacts|Too low", fpkm_allG$V55)] <- "Cell cycle"
table(filter(fpkm_allG, !is.na(CC))$V55)
#some lncs in there

fpkm_allG$GeneClassUpdate <- fpkm_allG$V55
fpkm_allG$GeneClassUpdate[!is.na(fpkm_allG$TF)] <- "TF"
fpkm_allG$GeneClassUpdate[!is.na(fpkm_allG$CC) & !fpkm_allG$V55 %in% c("Bona fide lncRNA", "Putative lncRNA")] <- "CC" #excluding lncs from cell cycle def
fpkm_allG$GeneClassUpdate[!is.na(fpkm_allG$TF) & !is.na(fpkm_allG$CC) & !fpkm_allG$V55 == "Bona fide lncRNA"] <- "TF + CC"

table(unique(fpkm_allG[,c(1,58)])$GeneClassUpdate)


#### (skip on revisit) remove AS artefacts ####

#data is unstranded, lncRNAs particularly prone to being mis-quantified via this type of data
#if looking at cisLncs, overlapping pairs could be source of artefacts

fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026.csv", header = T)
table(unique(fpkm_allG[,c(2,55)])$V55)
length(unique(fpkm_allG$EnsID))#14235 expressed genes, 847 bona fide lncs, 11887 PCGs

#manual check of all overlapping lncRNA loci in IGV here:
FPKM_CQV_OVERLAP_fpkm <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SameoverlapsG_lncs.csv")
table(FPKM_CQV_OVERLAP_fpkm$IGV)

# remove artefacts from PLAR too:
fpkm_allG_filt <- filter(fpkm_allG, grepl("chr", chr), !grepl("artefacts|Too low", GeneClassUpdate), !is.na(GeneClassUpdate))

#remove manual found fails
fpkm_allG_filt_manual <- filter(fpkm_allG_filt, 
                                !EnsID %in% filter(FPKM_CQV_OVERLAP_fpkm, IGV == "fail")$EnsID.x, #remove manual fails
                                )

table(unique(fpkm_allG_filt[,c(2,58)])$GeneClassUpdate)
#847 lncs
table(unique(fpkm_allG_filt_manual[,c(2,58)])$GeneClassUpdate)
#597 after filtering

fpkm_allG <- fpkm_allG_filt_manual[,-c(5:8)]
#write.csv(fpkm_allG, "fpkm_allG_2026filt.csv", row.names = F)

#repeat DEGs, post QC of gene/lncRNA artefacts:

#start matrix to find DEGs
fpkm_allGDE_lrt <- unique(filter(fpkm_allG, padj < 0.05)[,-c(1,43:50)])
fpkm_allGDE_pairs <-  unique(filter(fpkm_allG, (preadj_0_4 <0.05 & (LogFC_0_4 > log2(1.5) | LogFC_0_4 < -log2(1.5))) |
                                      (preadj_0_8 <0.05 & (LogFC_0_8 > log2(1.5) | LogFC_0_8 < -log2(1.5))) |
                                      (preadj_0_24 <0.05 & (LogFC_0_24 > log2(1.5) | LogFC_0_24 < -log2(1.5))) |
                                      (preadj_4_8 <0.05 & (LogFC_4_8 > log2(1.5) | LogFC_4_8 < -log2(1.5))) |
                                      (preadj_4_24 <0.05 & (LogFC_4_24 > log2(1.5) | LogFC_4_24 < -log2(1.5))) |
                                      (preadj_8_24 <0.05 & (LogFC_8_24 > log2(1.5) | LogFC_8_24 < -log2(1.5))))[,-c(1,43:50)])

length(unique(filter(fpkm_allGDE_lrt, EnsID %in% fpkm_allGDE_pairs$EnsID)$EnsID)) #5082 overlap
length(unique(fpkm_allGDE_lrt$EnsID)) #5308 LRT DEGs (may contain some low FCs across borad)
length(unique(fpkm_allGDE_pairs$EnsID)) #5750 pairwise DEGs (may contains some inconsistent patient-patient genes)

#Overlap probs best, kicks out low fold changing genes (LRT only) or inconsistent changing genes (Pairwise only)
fpkm_allGDE <- filter(fpkm_allGDE_lrt, EnsID %in% fpkm_allGDE_pairs$EnsID)

#no duplicates
dupEnsID <- fpkm_allGDE$EnsID[duplicated(fpkm_allGDE[,1])]

