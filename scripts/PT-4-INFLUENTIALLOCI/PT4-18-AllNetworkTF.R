library(dplyr)
library(ggplot2)

#### influential genes import ####

ModuleGenes <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/moduleGenesMM_Feb2026.csv")
HubGenes <- filter(ModuleGenes, moduleOfInterest == "Yes", grepl("Top", Hubness))
table(HubGenes$ModuleSummary)
table(HubGenes$ModuleSummary, HubGenes$CCLnc_association)#7 in sustained immune, 14 in trans immune/muscle prolif
table(HubGenes$ModuleSummary, HubGenes$Hubness)

LisaGenes <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/LISA_testInfluence2026.csv")
LisaGenes <- filter(LisaGenes, fpkm_max_treatment >1)
rownames(LisaGenes) <- NULL
LisaGenes <- LisaGenes[order(LisaGenes$lowestP, decreasing = F),]
#valid predictions only, expressed in the system to some degree:
#top 30:
dim(LisaGenes)[1]*0.1
LisaGenes_top30 <- LisaGenes[1:200,]
LisaGenes_top20 <- LisaGenes[1:133,]
LisaGenes_top10 <- LisaGenes[1:67,]

#alt LISA, top 30 in any DEG set:
LisaGenes <- LisaGenes[order(LisaGenes$upP04, decreasing = F),]
LisaGenes$upP04_rank <- 1:length(LisaGenes$Transcription.Factor)
LisaGenes <- LisaGenes[order(LisaGenes$downP04, decreasing = F),]
LisaGenes$downP04_rank <- 1:length(LisaGenes$Transcription.Factor)

LisaGenes <- LisaGenes[order(LisaGenes$upP08, decreasing = F),]
LisaGenes$upP08_rank <- 1:length(LisaGenes$Transcription.Factor)
LisaGenes <- LisaGenes[order(LisaGenes$downP08, decreasing = F),]
LisaGenes$downP08_rank <- 1:length(LisaGenes$Transcription.Factor)

LisaGenes <- LisaGenes[order(LisaGenes$upP024, decreasing = F),]
LisaGenes$upP024_rank <- 1:length(LisaGenes$Transcription.Factor)
LisaGenes <- LisaGenes[order(LisaGenes$downP024, decreasing = F),]
LisaGenes$downP024_rank <- 1:length(LisaGenes$Transcription.Factor)

LisaGenes$bestRank <- Biobase::rowMin(as.matrix(LisaGenes[,25:30]))

LisaGenes_top30_perSet <- filter(LisaGenes, bestRank <= 200)
386/665 #58% of expressed TFs are highly influential at one time or other... doubtful

ISMARAGenes <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/TF_ISMARA_ranking2026.csv")
ISMARAGenes <- ISMARAGenes[order(ISMARAGenes$Zscore, decreasing = T),]
#valid predictions only, well expressed:
ISMARAGenes <- filter(ISMARAGenes, fpkm_max_treatment >1)
#top 30:
dim(ISMARAGenes)[1]*0.1
ISMARAGenes_top30 <- ISMARAGenes[1:111,]
ISMARAGenes_top20 <- ISMARAGenes[1:74,]
ISMARAGenes_top10 <- ISMARAGenes[1:37,]


#### enrichment of muscle TFs ####

#isolate Lambert TFs for each
ModuleGenes_TFs <- filter(ModuleGenes, grepl("TF", GeneClassUpdate))
HubGenes_TFs <- filter(HubGenes, grepl("TF", GeneClassUpdate))

LisaGenes_lambert <- filter(LisaGenes, Transcription.Factor %in% filter(fpkm_allG, grepl("TF", GeneClassUpdate))$EnsName)
LisaGenes_top30_lambert <- filter(LisaGenes_top30, Transcription.Factor %in% filter(fpkm_allG, grepl("TF", GeneClassUpdate))$EnsName)

ISMARAGenes_lambert <- filter(ISMARAGenes, Symbols %in% filter(fpkm_allG, grepl("TF", GeneClassUpdate))$EnsName)
ISMARAGenes_top30_lambert <- filter(ISMARAGenes_top30, Symbols %in% filter(fpkm_allG, grepl("TF", GeneClassUpdate))$EnsName)

#expected TFs driving SMC phenos from lit knowledge:
muscle_TFs <- c("YY1", "KLF4", "SRF", "FOS", "MYOCD", "TET2", "SMAD3", "TCF21", "TEAD1", #my knowledge, TFs involved in promoting SMC mat. or dediff.
                #Miller snATAC paper, heatmap in fig2, motifs enriched in SMC marker genes in control art/athero tissue
                "MEF2A", "MEF2B", "MEF2C", "MEF2D", "TEAD", "TEAD2", "TEAD4", 
                "EBF1", "EBF", "BATF", "FRA1")

a <- sum(HubGenes_TFs$EnsName %in% muscle_TFs)
b <- length(HubGenes_TFs$EnsName)
c <- sum(ModuleGenes_TFs$EnsName %in% muscle_TFs)
d <- length(ModuleGenes_TFs$EnsName)

fisher.test(data.frame("muscleTF" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")

a <- sum(LisaGenes_top30_lambert$Transcription.Factor %in% muscle_TFs)
b <- length(LisaGenes_top30_lambert$Transcription.Factor)
c <- sum(LisaGenes_lambert$Transcription.Factor %in% muscle_TFs)
d <- length(LisaGenes_lambert$Transcription.Factor)

fisher.test(data.frame("muscleTF" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")

a <- sum(ISMARAGenes_top30$Symbols %in% muscle_TFs)
b <- length(ISMARAGenes_top30$Symbols)
c <- sum(ISMARAGenes$Symbols %in% muscle_TFs)
d <- length(ISMARAGenes$Symbols)

fisher.test(data.frame("muscleTF" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#ISMARA yes + decent, others no, LISA close (find more muscle TFs? better list?)


#IEG TFs
IEGs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Gene lists/arner_2015_table_S5_IEGs_lit.csv")
IEGs_hs <- filter(IEGs, !Hs_symbol %in% c("", "-"))
IEGs_hs <- filter(IEGs_hs, !grepl("[a-z]", IEGs_hs$Total_count))

#filter to IEG TFs - for all to have a fair display
fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026filt.csv")
length(unique(fpkm_allG$EnsID))#12740

IEGs_TFs <- IEGs_hs$Hs_symbol[IEGs_hs$Hs_symbol %in% filter(fpkm_allG, grepl("TF", GeneClassUpdate))$EnsName]

a <- sum(HubGenes_TFs$EnsName %in% IEGs_TFs)
b <- length(HubGenes_TFs$EnsName)
c <- sum(ModuleGenes_TFs$EnsName %in% IEGs_TFs)
d <- length(ModuleGenes_TFs$EnsName)

fisher.test(data.frame("muscleTF" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")


ai <- sum(LisaGenes_top30_lambert$Transcription.Factor %in% IEGs_TFs)
bi <- length(LisaGenes_top30_lambert$Transcription.Factor)
ci <- sum(LisaGenes_lambert$Transcription.Factor %in% IEGs_TFs)
di <- length(LisaGenes_lambert$Transcription.Factor)

fisher.test(data.frame("muscleTF" = c(ai, bi-ai),
                       "Not"   = c(ci-ai, di-ci-(bi-ai))), alternative = "greater")

aii <- sum(ISMARAGenes_top30_lambert$Symbols %in% IEGs_TFs)
bii <- length(ISMARAGenes_top30_lambert$Symbols)
cii <- sum(ISMARAGenes_lambert$Symbols %in% IEGs_TFs)
dii <- length(ISMARAGenes_lambert$Symbols)

fisher.test(data.frame("muscleTF" = c(aii, bii-aii),
                       "Not"   = c(cii-aii, dii-cii-(bii-aii))), alternative = "greater")

#all strong on IEG TFs

#use for validation:

TF_IEG_type <- data.frame("Top 30%" = c(a/b, 
                                           ai/bi, 
                                           aii/bii), 
                           "All" = c(c/d, 
                                           ci/di, 
                                           cii/dii))

TF_IEG_type$Approach <- as.factor(c("WGCNA", "LISA", "ISMARA"))
TF_IEG_type$NoTF <- c(a, ai, aii)

ggplot(TF_IEG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = Approach, y = Top.30.), fill = "#D6604D") +
  geom_col(aes(x = Approach, y = All), fill = NA, color = "grey30", linetype = "dashed", linewidth = 1.2) +
  geom_label(aes(x = Approach, y = Top.30., label = NoTF), size = 3) +  
  xlab("") +   ylab("\n% IEG TFs") +
  theme_minimal()

colnames(TF_IEG_type)[1:2] <- c("Top 30% Ranked TFs", "All Ranked TFs")
trial <- reshape2::melt(TF_IEG_type[,1:3])

ggplot(trial) + aes(x = Approach, y = value, fill = variable) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c(`Top 30% Ranked TFs` = "grey10", `All Ranked TFs` = "grey70")) +
  theme_minimal() +
  ylab("% IEG TFs") +
  xlab("") +
  theme(text = element_text(size=24))


#### build a heatmap of evidence for influence, and evidence for cis regulation (as with SCClncRNA supp figure) ####

SCClncRNAs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SCClncRNAs.csv")
SCClncRNAs$lncNameCol <- SCClncRNAs$EnsID
SCClncRNAs$lncNameCol[!is.na(SCClncRNAs$EnsName.x)] <- SCClncRNAs$EnsName.x[!is.na(SCClncRNAs$EnsName.x)]
SCClncRNAs$pairsNice <- paste(SCClncRNAs$lncNameCol, SCClncRNAs$EnsName.y, sep = "-")

#add influence info:
SCClncRNAs_inf <- merge(SCClncRNAs, HubGenes[,c(2,6,12,14)], by.x = "EnsID", by.y = "EnsID", all.x = T)
SCClncRNAs_inf <- merge(SCClncRNAs_inf, HubGenes[,c(2,6,12,14)], by.x = "EnsID.y", by.y = "EnsID", all.x = T)
SCClncRNAs_inf <- merge(SCClncRNAs_inf, ISMARAGenes_top30[,c(2,4)], by.x = "EnsName.y", by.y = "Symbols", all.x = T)
SCClncRNAs_inf <- merge(SCClncRNAs_inf, LisaGenes_top30[,c(1:8)], by.x = "EnsName.y", by.y = "Transcription.Factor", all.x = T)

#% influential TFs amongst targeted:
length(unique(filter(SCClncRNAs_inf, grepl("TF", GeneClassUpdate.y) | 
                       EnsName.y %in% c(ISMARAGenes$Symbols, LisaGenes$Transcription.Factor))$EnsName.y)) #allow in anything tested by LISA/ISMARA outside Lambert
length(unique(filter(SCClncRNAs_inf, (!is.na(Zscore) | !is.na(lowestP)),
                     grepl("TF", GeneClassUpdate.y) | 
                       EnsName.y %in% c(ISMARAGenes$Symbols, LisaGenes$Transcription.Factor))$EnsName.y)) #allow in anything tested by LISA/ISMARA outside Lambert

compare_overlap <- data.frame("HiC" = as.numeric(!SCClncRNAs_inf$loopMethod == "Neither"), 
                              "Correlated" = as.numeric(SCClncRNAs_inf$spear_p_adj < 0.05), 
                              "eQTL" = as.numeric(SCClncRNAs_inf$eQTLvalidations >0 | SCClncRNAs_inf$FANTOM_eQTL_Hon %in% "Significant Hon et al. 2017"),
                              "GapmeR-based prediction" = as.numeric(!is.na(SCClncRNAs_inf$Sources)),
                              #color by WGCNA module?
                              "WGNCA_central_lncRNA" = as.numeric(!is.na(SCClncRNAs_inf$Module.x)),
                              "WGNCA_central_PCG" = as.numeric(!is.na(SCClncRNAs_inf$Module.y)),#as.numeric(as.factor(SCClncRNAs_inf$Module)),
                              "ISMARA" = as.numeric(!is.na(SCClncRNAs_inf$Zscore)),
                              "LISA" = as.numeric(!is.na(SCClncRNAs_inf$lowestP))
                              )
rownames(compare_overlap) <- SCClncRNAs_inf$pairsNice

compare_overlap$Correlated[is.na(compare_overlap$Correlated)] <- 0
compare_overlap$eQTL[is.na(compare_overlap$eQTL)] <- 0

#rank by evidence
trial <- as.matrix(compare_overlap[,1:( length(colnames(compare_overlap)))])
rownames(trial) <- rownames(compare_overlap)
compare_overlap <- trial

compare_overlap <- compare_overlap[order(compare_overlap[,1], compare_overlap[,2], 
                                         compare_overlap[,3], compare_overlap[,4],
                                         decreasing = T),]

pheatmap::pheatmap(t(compare_overlap), cluster_cols = F, cluster_rows = F, angle_col = "315", legend = F, fontsize_row = 10)


compare_overlap_early <- compare_overlap[rownames(compare_overlap) %in% filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs")$pairsNice,]

#for supplement:
pheatmap::pheatmap(t(compare_overlap_early[,c(5,6,7,1,2,3,4)]), gaps_row = 3,
                   cluster_cols = F, cluster_rows = F, angle_col = "315", legend = F, fontsize_row = 9)


#just "influential"

compare_overlap_inf <- compare_overlap[rowSums(compare_overlap[,5:7])>0,]
pheatmap::pheatmap(t(compare_overlap_inf[,c(5,6,7,1,2,3,4)]), gaps_row = 3, cluster_cols = F, cluster_rows = F, angle_col = "315", legend = F, fontsize_row = 10)

compare_overlap_earlyUp_inf <- compare_overlap_early[rowSums(compare_overlap_early[,5:7])>0,]
pheatmap::pheatmap(t(compare_overlap_earlyUp_inf[,c(5,6,7,1,2,3,4)]), gaps_row = 3,
                   cluster_cols = F, cluster_rows = F, angle_col = "315", legend = F, fontsize_row = 9)

compare_overlap_earlyDown <- compare_overlap[rownames(compare_overlap) %in% filter(SCClncRNAs, Lnc_Cluster == "Repressed <4hrs")$pairsNice,]
compare_overlap_earlyDown_inf <- compare_overlap_earlyDown[rowSums(compare_overlap_earlyDown[,5:7])>0,]
pheatmap::pheatmap(t(compare_overlap_earlyDown_inf[,c(5,6,7,1,2,3,4)]), gaps_row = 3,
                   cluster_cols = F, cluster_rows = F, angle_col = "315", legend = F, fontsize_row = 9)

compare_overlap_other <- compare_overlap[!rownames(compare_overlap) %in% filter(SCClncRNAs, grepl("<4hrs", Lnc_Cluster))$pairsNice,]
compare_overlap_other_inf <- compare_overlap_other[rowSums(compare_overlap_other[,5:7])>0,]
pheatmap::pheatmap(t(compare_overlap_other_inf[,c(5,6,7,1,2,3,4)]), gaps_row = 3,
                   cluster_cols = F, cluster_rows = F, angle_col = "315", legend = F, fontsize_row = 9)


#### selecting for the lab ####

#just influential, just early induced:

SCClncRNAs_inf_lab <- filter(SCClncRNAs_inf, Lnc_Cluster == "Induced <4hrs", 
                               !is.na(SCClncRNAs_inf$Module.x) | 
                               !is.na(SCClncRNAs_inf$Module.y) | !is.na(SCClncRNAs_inf$Zscore) | !is.na(SCClncRNAs_inf$lowestP))

#practical, well-expressed lncRNAs
fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026filt.csv")
length(unique(fpkm_allG$EnsID))#12740

SCClncRNAs_inf_lab <- merge(SCClncRNAs_inf_lab, unique(fpkm_allG[,c(2,25,27,29,31,33)]), by = "EnsID")

#no simple way to reconstruct the selection we did which was based on:

#>3 FPKM, but then added in IPL-IL6
#previous SCClncRNAs, which included 12913 (now lost)
#wider early responses (not mentioned yet)


