#define elnc status for timecourse lncs
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(GenomicRanges)
library(ggpubr)
library(FSA)
library(rcompanion)
library(DescTools)

#all lncRNAs + CAGE:
allLncs_BestCAGE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/allLncs_FANTOMCAGE_3PLARtimecourse_2026.csv",
                             header = T, stringsAsFactors = F)

FANTOM_OntEnrich <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/nature21374-s2/table11.csv", header = T, stringsAsFactors = F)

Enhancer_loci <- allLncs_BestCAGE
length(unique(Enhancer_loci$FANTOM_ID))
#591 valid FANTOM IDs checked
Enhancer_loci$CAGEvalidity[Enhancer_loci$BEstTIEScore >60 & Enhancer_loci$PercentExon1Retained > 0.1] <- "Valid CAGE"

#use CAGE TSS to get promoter region for these(strand-sensitive)
Enhancer_loci$Enh_Start[!is.na(Enhancer_loci$CAGEvalidity) & 
                          Enhancer_loci$CTSSstrand == "+"] <- Enhancer_loci$BestStart[!is.na(Enhancer_loci$CAGEvalidity) & 
                                                                                        Enhancer_loci$CTSSstrand == "+"] -2000
Enhancer_loci$Enh_Stop[!is.na(Enhancer_loci$CAGEvalidity) & 
                         Enhancer_loci$CTSSstrand == "+"] <- Enhancer_loci$BestStart[!is.na(Enhancer_loci$CAGEvalidity) & 
                                                                                       Enhancer_loci$CTSSstrand == "+"] +200

Enhancer_loci$Enh_Stop[!is.na(Enhancer_loci$CAGEvalidity) & 
                         Enhancer_loci$CTSSstrand == "-"] <- Enhancer_loci$BestStart[!is.na(Enhancer_loci$CAGEvalidity) & 
                                                                                       Enhancer_loci$CTSSstrand == "-"] +2000
Enhancer_loci$Enh_Start[!is.na(Enhancer_loci$CAGEvalidity) & 
                          Enhancer_loci$CTSSstrand == "-"] <- Enhancer_loci$BestStart[!is.na(Enhancer_loci$CAGEvalidity) & 
                                                                                        Enhancer_loci$CTSSstrand == "-"] -200

#otherwise just use the regular 5'
hSVSMC_5p <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/allexprslnc5range_hSVSMCtime_2026.bed", header = T, sep = "\t", stringsAsFactors = F)

trial <- unique(merge(Enhancer_loci, hSVSMC_5p[,c(1,2,4,6)], by = "MSTRG_Tx_ID"))
trial$start <- trial$start+1000 #real TSS, central point of range
trial$X5pStringtie <- trial$start
trial$CTSSstrand <- trial$str
#now no more na values for rows without a CTSS
Enhancer_lociII <- trial

Enhancer_lociII$Enh_Start[is.na(Enhancer_lociII$CAGEvalidity) & 
                            Enhancer_lociII$CTSSstrand == "+"] <- Enhancer_lociII$start[is.na(Enhancer_lociII$CAGEvalidity) & 
                                                                                          Enhancer_lociII$CTSSstrand == "+"] -2000
Enhancer_lociII$Enh_Stop[is.na(Enhancer_lociII$CAGEvalidity) & 
                           Enhancer_lociII$CTSSstrand == "+"] <- Enhancer_lociII$start[is.na(Enhancer_lociII$CAGEvalidity) & 
                                                                                         Enhancer_lociII$CTSSstrand == "+"] +200

Enhancer_lociII$Enh_Stop[is.na(Enhancer_lociII$CAGEvalidity) & 
                           Enhancer_lociII$CTSSstrand == "-"] <- Enhancer_lociII$start[is.na(Enhancer_lociII$CAGEvalidity) & 
                                                                                         Enhancer_lociII$CTSSstrand == "-"] +2000
Enhancer_lociII$Enh_Start[is.na(Enhancer_lociII$CAGEvalidity) & 
                            Enhancer_lociII$CTSSstrand == "-"] <- Enhancer_lociII$start[is.na(Enhancer_lociII$CAGEvalidity) & 
                                                                                          Enhancer_lociII$CTSSstrand == "-"] -200

#create .bed file to extract Genehancer annotations
Enhancer_lociIIbed <- unique(Enhancer_lociII[,c(19,17:18,1,21)])

#write.table(Enhancer_lociIIbed[,-5], "AllLncs_CAGE_NonCAGE_Promoterstime_2026.bed", quote = F, sep = "\t", row.names = F)

#Import from table browser UCSC:
AllLncs_Genehancer <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GHtimecourseLncs_2026.csv", header = T)

length(unique(AllLncs_Genehancer$name))#450 GH annotations (383 previously)
AllLncs_Genehancer <- unique(AllLncs_Genehancer)

AllLncs_Genehancer_Coords <- sort(makeGRangesFromDataFrame(unique(AllLncs_Genehancer[,c(1:5,10:11)]), 
                                                           seqnames.field = "chrom", 
                                                           start.field = "chromStart", end.field = "chromEnd",
                                                           keep.extra.columns = T))

Enhancer_lociII_Coords <- makeGRangesFromDataFrame(Enhancer_lociIIbed, keep.extra.columns = T, strand.field = "str")
Enhancer_lociII_Coords <- sort(Enhancer_lociII_Coords)
AllLncs2GenHan <- subsetByOverlaps(Enhancer_lociII_Coords, AllLncs_Genehancer_Coords) #1064 tx overlap a Genehancer annotation
GenHan2AllLncs <- subsetByOverlaps(AllLncs_Genehancer_Coords, Enhancer_lociII_Coords) #450 enhancers overlap 1064 tx

GenHan2AllLncs_Index <- findOverlaps(AllLncs_Genehancer_Coords, Enhancer_lociII_Coords) #1089 links
#query hits, which genehancer rows match subject hits, transcript rows

trial <- Enhancer_lociII_Coords$MSTRG_Tx_ID[subjectHits(GenHan2AllLncs_Index)]
trialii <- AllLncs_Genehancer_Coords[queryHits(GenHan2AllLncs_Index),]
AllLncs2GenHan_Info <- data.frame("MSTRG_Tx_ID" = trial, 
                                  "Genehancer_name" = trialii$name, 
                                  "Genehancer_type" = trialii$elementType, 
                                  "Genehancer_score" = trialii$score, 
                                  "Genehancer_source" = trialii$evidenceSources,
                                  stringsAsFactors = F)

#total number of available enhancer annotations
Enhancer_lociII <- merge(Enhancer_lociII, AllLncs2GenHan_Info, by = "MSTRG_Tx_ID", all.x = T)
Enhancer_lociII$FANTOMenhancer[Enhancer_lociII$FANTOM_ID %in% filter(FANTOM_OntEnrich, CAT_geneCategory == "e_lncRNA")$CAT_geneID] <- "FANTOM Enhancer"
Enhancer_lociII$FANTOMpromoter[Enhancer_lociII$FANTOM_ID %in% filter(FANTOM_OntEnrich, grepl("p_lncRNA", CAT_geneCategory))$CAT_geneID] <- "FANTOM Promoter"

#add usual gene info
fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026filt.csv", header = T, stringsAsFactors = F)

fpkm_PLARG <- filter(fpkm_allG, V55 == "Bona fide lncRNA")
length(unique(fpkm_PLARG$EnsID))#597 unique entries
length(unique(fpkm_PLARG$MSTRG_Tx_ID))#1575 unique entries

trial <- merge(Enhancer_lociII, fpkm_PLARG[,c(2,34:46)], by = "EnsID", all.x = T)
Enhancer_lociII <- trial

Enhancer_lociII$EnhancerVerdict <- NA
summary(unique(filter(Enhancer_lociII, Genehancer_type == "Enhancer")[c(22,24)])$Genehancer_score) #select 1st quartile for stringency or ignore (they are "double elite")
Enhancer_lociII$EnhancerVerdict[!is.na(Enhancer_lociII$FANTOMenhancer) | 
                                  (Enhancer_lociII$Genehancer_type == "Enhancer" & Enhancer_lociII$Genehancer_score > 215)] <- "Enhancer"

length(unique(Enhancer_lociII$EnsID))#597 genes total (some MSTRGs will overlap)
length(unique(Enhancer_lociII$EnsID[grepl("Enhancer", Enhancer_lociII$EnhancerVerdict)]))#98 enhancer lncRNAs expressed

#will require going through to get DE genes in Timeoucrse4timepoints.... etc to get this:
Enhancer_lociII$DiffExprs[Enhancer_lociII$EnsID %in% fpkm_allGDE$EnsID] <- "DiffExprs"

#write.csv(Enhancer_lociII, "Enhancer_lociIItime_2026.csv", row.names = F)
Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime.csv", header = T)

#OUT OF DATE NUMBERS
#597 lncs expressed
length(unique(Enhancer_lociII$EnsID))
#221 DE lncs over 24 hours:
length(unique(filter(Enhancer_lociII, !is.na(DiffExprs))$EnsID))
#98 enhancer lncs expressed
length(unique(filter(Enhancer_lociII, !is.na(EnhancerVerdict))$EnsID))
#54 enhancer lncs expressed and DE
length(unique(filter(Enhancer_lociII, !is.na(DiffExprs), !is.na(EnhancerVerdict))$EnsID))

221/597 #37% all lncs DE
54/98   #55% of elncs DE

#elncs are particularly IP responsive
