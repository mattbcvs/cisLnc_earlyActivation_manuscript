library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

#most likely CClncRNAs to actually be acting via a cis mechanism

#### The 20x HiC-linked pairs within 1mbp from 0-4hr CClncRNAs ####
AllLNC_AllPCG_2d3d <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_2d3d_2026.csv")
HiC_pairs <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <1000, !loopMethod == "Neither", Lnc_Cluster == "Induced <4hrs", PCG_Timeframe == "<4hrs")

#### The 67x eQTL-linked pairs within 250kbp ####
CoRegPairs_eQTL_supported <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CoRegPairs_eQTL_supported.csv")


#### 86x correlated co-regulated pairs within 250kbp: ####
#Spearman's on FPKM
AllLNC_AllPCG_2d3d <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_2d3d_2026.csv")
fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026filt.csv")
fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE_2026filt.csv", header = T)
length(unique(fpkm_allGDE$EnsID))

fpkm_allGDE_Upwithin_4 <- filter(fpkm_allGDE, RegulationStart == "Induced <4hrs")
fpkm_allGDE_Downwithin_4 <- filter(fpkm_allGDE, RegulationStart == "Repressed <4hrs")

fpkm_allGDE_Upwithin_8 <- filter(fpkm_allGDE, RegulationStart == "Induced 4-8hrs")
fpkm_allGDE_Downwithin_8 <- filter(fpkm_allGDE, RegulationStart == "Repressed 4-8hrs")

fpkm_allGDE_Upwithin_24 <- filter(fpkm_allGDE, RegulationStart == "Induced 8-24hrs")
fpkm_allGDE_Downwithin_24 <- filter(fpkm_allGDE, RegulationStart == "Repressed 8-24hrs")

CoRegPairs_04_48_24_extended <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <250,
                                             #AllLNC_AllPCG_1,
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                           fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                               (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                             fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                              fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                               (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                             fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
                                             )

AllLNC_AllPCG_FPKM <- merge(CoRegPairs_04_48_24_extended[,c(2,1,3)], fpkm_allG[,c(2,9:24)], by.x = "EnsID", by.y = "EnsID")
AllLNC_AllPCG_FPKM <- merge(AllLNC_AllPCG_FPKM, fpkm_allG[,c(2,9:24)], by.x = "EnsID.y", by.y = "EnsID")
AllLNC_AllPCG_FPKM <- unique(AllLNC_AllPCG_FPKM)
colnames(AllLNC_AllPCG_FPKM)
trial <- split(AllLNC_AllPCG_FPKM, AllLNC_AllPCG_FPKM$pair)
triali <- sapply(trial, function(x){
  cor.test(as.numeric(x[4:19]), as.numeric(x[20:35]), method = "spearman")$estimate
})
trialii <- sapply(trial, function(x){
  cor.test(as.numeric(x[4:19]), as.numeric(x[20:35]), method = "spearman")$p.value
})
trialiii <- data.frame("pairs" = names(trial), "spear_rho" = triali, "spear_p" = trialii, "spear_p_adj" = p.adjust(trialii, method = "BH"))

CoRegPairs_04_48_24_extendedcorr <- merge(trialiii, CoRegPairs_04_48_24_extended, by = "pairs")
dim(filter(CoRegPairs_04_48_24_extendedcorr, spear_p_adj < 0.05))#86 pairs

CoRegPairs_04_48_24_extendedcorr$corSig <- "No"
CoRegPairs_04_48_24_extendedcorr$corSig[abs(CoRegPairs_04_48_24_extendedcorr$spear_rho) >0.5 & 
                                          CoRegPairs_04_48_24_extendedcorr$spear_p_adj <0.05] <- "Yes"
table(CoRegPairs_04_48_24_extendedcorr$corSig)#still 86 (with or without the rho filter)

correlatedPairs <- filter(CoRegPairs_04_48_24_extendedcorr, corSig == "Yes")

correlatedPairs <- filter(AllLNC_AllPCG_2d3d, pairs %in% correlatedPairs$pairs)


#### 25x previously found as eQTL-linked from FANTOM (within 250kbp) ####

FANTOM_eQTL_pairs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/nature21374-s2/table16.csv")
FANTOM_eQTL_pairs_cis <- filter(FANTOM_eQTL_pairs, cis_correlated_candidate == "yes")
#5264 pairs as reported in Hon et al 2017

Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime.csv", header = T)

#add FANTOM ID to Alllnc:
trial <- unique(merge(AllLNC_AllPCG_2d3d, Enhancer_lociII[,c(1,14)], by = "EnsID"))

#pair up columns:
trial$pairs_noSuff <- paste(gsub("\\.[0-9]*", "", trial$EnsID), gsub("\\.[0-9]*", "", trial$EnsID.y), sep = "-")
trial$FANTpairs_noSuff <- paste(gsub("\\.[0-9]*", "", trial$FANTOM_ID), gsub("\\.[0-9]*", "", trial$EnsID.y), sep = "-")

FANTOM_eQTL_pairs_cis$pairs_noSuff <- paste(gsub("\\.[0-9]*", "", FANTOM_eQTL_pairs_cis$lncRNA_geneID), 
                                            gsub("\\.[0-9]*", "", FANTOM_eQTL_pairs_cis$mRNA_geneID), sep = "-")

#no. pairs found in dataset:
triali <- filter(trial, FANTpairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff |
                   pairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff)
length(unique(triali$pairs))#298 pairs both found expressed in SVSMC

#part of the co-reg set, allow any range:
CoRegPairs_04_48_24_extended <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <250,
                                       #AllLNC_AllPCG_1,
                                       (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                     fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                       fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                        fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                       fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
                                       )

triali <- filter(trial, pairs %in% CoRegPairs_04_48_24_extended$pairs, 
                 FANTpairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff |
                   pairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff)
length(unique(triali$pairs))#25 found both expressed in SVSMC and in a co-reg pairing

FANTOM_eQTL_pairs <- filter(CoRegPairs_04_48_24_extended, pairs %in% triali$pairs)


#### expected cis-acting lncRNAs (within 250kbp otherwise too permissive in terms of pairs)####

#any FANTOM match up here:
Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime.csv", header = T)

ControlCisLncs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Write-ups/supplement/ControlCisLncs.csv", header = T)

#73 genes in total
unique(ControlCisLncs$Ens_ID)
table(ControlCisLncs$source)

fpkm_allG$EnsID_merge <- gsub("\\.[0-9]*", "", fpkm_allG$EnsID)
ControlCisLncs <- merge(ControlCisLncs, fpkm_allG[,c(2,3,59)], by.x = "Ens_ID", by.y = "EnsID_merge", all.x = T)

Enhancer_lociII$FANTID_merge <- gsub("\\.[0-9]*", "", Enhancer_lociII$FANTOM_ID)
ControlCisLncs <- merge(ControlCisLncs, Enhancer_lociII[,c(1,43)], by.x = "Ens_ID", by.y = "FANTID_merge", all.x = T)

ControlCisLncs <- unique(ControlCisLncs)

#EnsID for 42 lncRNAs expressed in SVSMC with exp evidence suggesting cis acting abiity
ControlCisLncs_exprs <- unique(c(ControlCisLncs$EnsID.x, ControlCisLncs$EnsID.y))
ControlCisLncs_exprs <- ControlCisLncs_exprs[!is.na(ControlCisLncs_exprs)]

CoRegPairs_04_48_24_extended <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <250,
                                       #AllLNC_AllPCG_1,
                                       (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                     fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                       fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                        fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                       fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
                                       )

expectedCisActing_pairs <- filter(CoRegPairs_04_48_24_extended, EnsID %in% ControlCisLncs_exprs)

#previously all expected cis-acting lncs induced within 4hr were hic/eqtl/corr candidates:
length(unique(expectedCisActing_pairs$EnsID)) #14
length(unique(CoRegPairs_04_48_24_extended$EnsID)) #124

length(unique(filter(expectedCisActing_pairs, Lnc_Timeframe == "<4hrs")$EnsID)) #9

#3/4 HiC
sum(unique(filter(expectedCisActing_pairs, Lnc_Timeframe == "<4hrs")$EnsName.x) %in% HiC_pairs$EnsName.x)
#2/4 eQTL
unique(filter(expectedCisActing_pairs, Lnc_Timeframe == "<4hrs")$EnsName.x) %in% CoRegPairs_eQTL_supported$EnsName.x
#4/4 correlated
unique(filter(expectedCisActing_pairs, Lnc_Cluster == "Induced <4hrs")$EnsName.x) %in% correlatedPairs$EnsName.x


#### aggregate all info ####

SCClncRNAs <- unique(rbind(HiC_pairs,
                    CoRegPairs_eQTL_supported[,1:19],
                    FANTOM_eQTL_pairs,
                    correlatedPairs[,1:19],
                    expectedCisActing_pairs
                    ))
length(unique(SCClncRNAs$EnsID))

unique(filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs")$EnsID)#37
unique(filter(SCClncRNAs, Lnc_Cluster == "Repressed <4hrs")$EnsID)#21
unique(filter(SCClncRNAs, Lnc_Cluster == "Induced 4-8hrs")$EnsID)#7
unique(filter(SCClncRNAs, Lnc_Cluster == "Repressed 4-8hrs")$EnsID)#11
unique(filter(SCClncRNAs, Lnc_Cluster == "Induced 8-24hrs")$EnsID)#7
unique(filter(SCClncRNAs, Lnc_Cluster == "Repressed 8-24hrs")$EnsID)#4

#columns:
#HiC covered x
#eQTL, number of validations:
SCClncRNAs <- merge(SCClncRNAs, CoRegPairs_eQTL_supported[,c(3,61)], by = "pairs", all.x = T)
#FANTOM eQTL yes/no
SCClncRNAs$FANTOM_eQTL_Hon <- "Not significant"
SCClncRNAs$FANTOM_eQTL_Hon[SCClncRNAs$pairs %in% FANTOM_eQTL_pairs$pairs] <- "Significant Hon et al. 2017"
#correlation, p/rho
SCClncRNAs <- merge(SCClncRNAs, CoRegPairs_04_48_24_extendedcorr[,c(1:4)], by = "pairs", all.x = T)
#expected cis, author source
#sort lncs with multiple sources:
trial <- ControlCisLncs
trial$EnsID.x[is.na(trial$EnsID.x)] <- trial$EnsID.y[is.na(trial$EnsID.x)]
trial <- as.data.frame(sapply(split(trial, trial$EnsID.x), function(x){paste(x$source, collapse = "|")}))
trial$EnsID.x <- rownames(trial)
colnames(trial)[1] <- "Sources"
SCClncRNAs <- merge(SCClncRNAs, trial, by.x = "EnsID", by.y = "EnsID.x", all.x = T)
#class HOTAIR with the HOTAIR paper, Rinn 2007:
SCClncRNAs$Sources[SCClncRNAs$EnsName.x == "HOTAIR"] <- "Rinn"

#write.csv(SCClncRNAs, "SCClncRNAs.csv", row.names = F)

#
#### import the SCClncRNAs ####

SCClncRNAs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SCClncRNAs.csv")
SCClncRNAs$lncNameCol <- SCClncRNAs$EnsID
SCClncRNAs$lncNameCol[!is.na(SCClncRNAs$EnsName.x)] <- SCClncRNAs$EnsName.x[!is.na(SCClncRNAs$EnsName.x)]
SCClncRNAs$pairsNice <- paste(SCClncRNAs$lncNameCol, SCClncRNAs$EnsName.y, sep = "-")

#total no:
length(unique(SCClncRNAs$EnsID))

#venn diagram:
colnames(SCClncRNAs)

write.csv(unique(filter(SCClncRNAs, !loopMethod == "Neither")$EnsID), "HiC_SCC.csv")
write.csv(unique(filter(SCClncRNAs, eQTLvalidations >0 | FANTOM_eQTL_Hon != "Not significant")$EnsID), "eQTL_SCC.csv")
write.csv(unique(filter(SCClncRNAs, abs(spear_rho) > 0.5, spear_p_adj < 0.05)$EnsID), "corr_SCC.csv")
write.csv(unique(filter(SCClncRNAs, !is.na(Sources))$EnsID), "lit_SCC.csv")

filter(SCClncRNAs, !EnsID %in%
unique(c(
  filter(SCClncRNAs, !loopMethod == "Neither")$EnsID,
  filter(SCClncRNAs, eQTLvalidations >0 | FANTOM_eQTL_Hon != "Not significant")$EnsID,
  filter(SCClncRNAs, abs(spear_rho) > 0.5, spear_p_adj < 0.05)$EnsID,
  filter(SCClncRNAs, !is.na(Sources))$EnsID
)))

unique(filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs")$EnsID)
28/37

#
#### supplementary table 6 ####

#add back in the eQTL data:
colnames(CoRegPairs_eQTL_supported)
trial <- merge(SCClncRNAs, CoRegPairs_eQTL_supported[,c(3,22:60)], by = "pairs", all.x = T)

write.csv(trial, "SuppTable6_SCClncRNAsi.csv", row.names = F)


#
#### % and number per cluster ####

SCClncRNAs_time <- as.data.frame(table(unique(SCClncRNAs[,c(1,9)])$Lnc_Cluster))
SCClncRNAs_time$selection <- table(fpkm_allGDE$RegulationStart)

sum(SCClncRNAs_time$Freq)
37/87
1105/5081

LncEnrich_cluster <- list()
for (i in 1:length(SCClncRNAs_time$Var1)){
  a <- SCClncRNAs_time[i,2]
  b <- SCClncRNAs_time[i,3]
  c <- sum(SCClncRNAs_time$Freq)
  d <- 5081
  
  LncEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p)
}
names(LncEnrich_cluster) <- SCClncRNAs_time$Var1
triali <- bind_rows(LncEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")

#percentage plots
SCClncRNAs_time$FirstRegulation <- c("Within \n4hrs", "Within \n8hrs", "Within \n24hrs", "Within \n4hrs", "Within \n8hrs", "Within \n24hrs")
SCClncRNAs_time$FirstRegulation <- as.factor(SCClncRNAs_time$FirstRegulation)
SCClncRNAs_time$FirstRegulation <- factor(SCClncRNAs_time$FirstRegulation, levels = levels(SCClncRNAs_time$FirstRegulation)[c(2,3,1)])

SCClncRNAs_time$UpDown <- sapply(sapply(as.character(SCClncRNAs_time$Var1), strsplit, " "),"[[" , 1)

SCClncRNAs_time$PercCategory <- SCClncRNAs_time$Freq/sum(SCClncRNAs_time$Freq)*100
SCClncRNAs_time$PercBackground <- table(fpkm_allGDE$RegulationStart)/5081*100

#ELncWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(SCClncRNAs_time, aes(x = FirstRegulation)) +
  geom_col(data = filter(SCClncRNAs_time, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown), color = "black") +
  geom_col(data = filter(SCClncRNAs_time, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed", linewidth = 1) +
  geom_label(data = filter(SCClncRNAs_time, grepl("Induced", UpDown)), 
             aes(y = PercCategory-1, label = Freq), size = 3.2) +
  geom_col(data = filter(SCClncRNAs_time, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown), color = "black") +
  geom_col(data = filter(SCClncRNAs_time, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed", linewidth = 1) +
  geom_label(data = filter(SCClncRNAs_time, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory+1, label = Freq), size = 3.2) +
  ylab("% SCClncRNAs") +
  xlab("") +
  scale_fill_manual(values = c("Induced" = "#D6604D", "Repressed" = "#67A9CF")) +
  #scale_y_continuous(limits = c(-30,60),breaks = seq(-20,60, by = 20),
  #                   labels = (c(seq(20, 0, by = -20), seq(20,60,by=20)))) +
  theme_minimal() +
  theme(text = element_text(size = 15), legend.position = "none")


#### Heatmap of cis-acting evidence for 0-4hr induced lncRNAs ####

compare_overlap <- data.frame("HiC" = as.numeric(!SCClncRNAs$loopMethod == "Neither"), 
                              "Correlated" = as.numeric(SCClncRNAs$spear_p_adj < 0.05), 
                              "eQTL" = as.numeric(SCClncRNAs$eQTLvalidations >0 | SCClncRNAs$FANTOM_eQTL_Hon %in% "Significant Hon et al. 2017"),
                              "GapmeR-based prediction" = as.numeric(!is.na(SCClncRNAs$Sources)))
rownames(compare_overlap) <- SCClncRNAs$pairsNice

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

rownames(compare_overlap_early)[rowSums(compare_overlap_early) >1]
26/85 #26 pairs with more than one piece of evidence
unique(sapply(strsplit(rownames(compare_overlap_early)[rowSums(compare_overlap_early) >1], "-"), function(x){paste(x[1:(length(x)-1)], collapse = "-")}))
19/37

#for supplement:
pheatmap::pheatmap((compare_overlap_early), cluster_cols = F, cluster_rows = F, angle_col = "315", legend = F, fontsize_row = 9)


#### enrichment of expected cis in HiC/corr/eQTL lncs? ####

#only fair to use 250kbp
length(unique(filter(SCClncRNAs, AbsDistLnc_PCG <250,
                     (!loopMethod == "Neither" | 
                       eQTLvalidations >0 | FANTOM_eQTL_Hon %in% "Significant Hon et al. 2017" | 
                       spear_p_adj < 0.05))$EnsID))#83 lncs
length(unique(filter(SCClncRNAs, AbsDistLnc_PCG <250,
                     (!loopMethod == "Neither" | 
                       eQTLvalidations >0 | FANTOM_eQTL_Hon %in% "Significant Hon et al. 2017" | 
                       spear_p_adj < 0.05) &  !is.na(Sources))$EnsID))#10 lncs
10/81 #12.3%

14/124 #11.3%

#not a great deal really


#### enrichment of SCClncRNAs in 0-4hr induced lncRNAs vs. other DE lncRNAs ####

length(unique(SCClncRNAs$EnsID))#87
unique(filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs")$EnsID)#37
table(fpkm_allGDE_Upwithin_4$GeneClassUpdate)#65
table(fpkm_allGDE$GeneClassUpdate)#221

a <- 37
b <- 87
c <- 65
d <- 221

a/b #40% of SCClncRNAs are 0-4hr up
c/d #29% of DELs are 0-4hr up

#big enrichment - but already covered this maybe?
fisher.test(data.frame("cisLnc" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")


#### enrichment of enhancer lncRNAs amongst 0-4hr induced SCClncRNAs vs other 0-4hr induced lncs ####

#add in eLncs
Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime_2026.csv", header = T)
fpkm_allGDE$GeneClassUpdate[fpkm_allGDE$EnsID %in% filter(Enhancer_lociII, !is.na(EnhancerVerdict))$EnsID & fpkm_allGDE$V55 == "Bona fide lncRNA"] <- "ELnc"

unique(filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs")$EnsID)#37
unique(filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs", EnsID %in% filter(fpkm_allGDE, GeneClassUpdate == "ELnc")$EnsID)$EnsID)#15

filter(fpkm_allGDE_Upwithin_4, GeneClassUpdate == "Bona fide lncRNA", EnsID %in% filter(fpkm_allGDE, GeneClassUpdate == "ELnc")$EnsID)$EnsID#23
table(fpkm_allGDE_Upwithin_4$GeneClassUpdate)#65

a <- 15
b <- 37
c <- 23
d <- 65

a/b #40% of 0-4hr SCClncRNAs are eLncs
c/d #35% of 0-4hr DELs are elncs

#not big enough enrich for sig
fisher.test(data.frame("cisLnc" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")



#### GO/KEGG/REACTOME analysis same/later of 4hr lncs ####

#potential for stronger bio signal amongst higher confidence pairs
#will try same/later pairs first
#note: potential interference from and MSTRG.24277 (all the CXCLs)

#simple GO check, vs. all DEGs:

#background of vs. all DEGs may give specificity (but also reduce power)
fpkm_PCGDE <- filter(fpkm_allGDE, grepl("protein_coding", EnsType))

#filter the scclncs
EarlyInduced_SCClncRNAs_Targets <- filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs")

#targets of early lncs vs all DEGs
CoReg_DE <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", EarlyInduced_SCClncRNAs_Targets$EnsID.y)),
                     universe      = gsub("\\.[0-9]*", "", fpkm_PCGDE$EnsID),
                     keyType       = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "all",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
CoReg_DE_df <- as.data.frame(CoReg_DE)
#terms of interest, but purely driven by cxcl

#loosen the background, expressed PCGs:
fpkm_PCG <- filter(fpkm_allG, grepl("protein_coding", EnsType))

CoReg_DE2 <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", EarlyInduced_SCClncRNAs_Targets$EnsID.y)),
                     universe      = gsub("\\.[0-9]*", "", fpkm_PCG$EnsID),
                     keyType       = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "all",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
CoReg_DE2_df <- as.data.frame(CoReg_DE2)
#few more genes, 5x loci involved in cytokine receptor binding (but still CXCL driven 4x genes)

#could select one from each (highest abundance @4hrs)
#AllTargets_loci1 <- filter(AllTargets_T2, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName.y))
#AllTargets_loci <- filter(AllTargets_T2, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName.y))
#AllTargets_loci_T2 <- rbind(AllTargets_loci, AllTargets_loci1)

#fpkm_PCGDE_loci1 <- filter(fpkm_PCGDE, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName))
#fpkm_PCGDE_loci <- filter(fpkm_PCGDE, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName))
#fpkm_PCGDE_loci <- rbind(fpkm_PCGDE_loci, fpkm_PCGDE_loci1)

#but feels messy
#try KEGG/REACTOME, then move on

#KEGG
convertEnsEnt <- bitr(unique(fpkm_allG$EnsName), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

CoReg_DE2_K <- enrichKEGG(gene = filter(convertEnsEnt, SYMBOL %in% EarlyInduced_SCClncRNAs_Targets$EnsName.y)$ENTREZID,
                       universe = filter(convertEnsEnt, SYMBOL %in% fpkm_PCG$EnsName)$ENTREZID,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
CoReg_DE2_K_df <- data.frame(CoReg_DE2_K)
unique(CoReg_DE2_K_df$Description)
filter(convertEnsEnt, 
       ENTREZID %in% unlist(strsplit(filter(CoReg_DE2_K_df, Description == unique(CoReg_DE2_K_df$Description)[2])$geneID, 
                                     split = "/")))
#essentially identifies same genes pathways

#REACTOME
library(ReactomePA)
CoReg_DE2_R <- enrichPathway(gene          = unique(filter(convertEnsEnt,SYMBOL %in% EarlyInduced_SCClncRNAs_Targets$EnsName.y)$ENTREZID),
                          universe      = unique(filter(convertEnsEnt,SYMBOL %in% fpkm_PCG$EnsName)$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)
CoReg_DE2_R_df <- as.data.frame(CoReg_DE2_R)
#as previous


#### GO/KEGG/REACTOME analysis co-induced with 4hr lncs ####

#try just same (n.b. almost all are co-induction)
#note: potential interference from and MSTRG.24277 (all the CXCLs)

#simple GO check, vs. all DEGs:

#background of vs. all DEGs may give specificity (but also reduce power)
#fpkm_allGDE_within_4 <- rbind(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Downwithin_4)
fpkm_PCGDE_Upwithin_4 <- filter(fpkm_allGDE_Upwithin_4, grepl("protein_coding", EnsType))

#filter the scclncs
EarlyInduced_SCClncRNAs_CoUp_Targets <- filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs", PCG_Cluster == "Induced <4hrs")

#targets of early lncs vs all DEGs
CoUp_DE <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", EarlyInduced_SCClncRNAs_CoUp_Targets$EnsID.y)),
                     universe      = gsub("\\.[0-9]*", "", fpkm_PCGDE_Upwithin_4$EnsID),
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
fpkm_PCG <- filter(fpkm_allG, grepl("protein_coding", EnsType))

CoUp_DE2 <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", EarlyInduced_SCClncRNAs_CoUp_Targets$EnsID.y)),
                      universe      = gsub("\\.[0-9]*", "", fpkm_PCG$EnsID),
                      keyType       = "ENSEMBL",
                      OrgDb         = org.Hs.eg.db,
                      ont           = "all",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
CoUp_DE2_df <- as.data.frame(CoUp_DE2)
#much better, despite clear CXCL locus drive
#gets FOXL1 and HMGA2 involved too

dotplot(simplify(CoUp_DE2))


#KEGG
convertEnsEnt <- bitr(unique(fpkm_allG$EnsName), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

CoUp_DE2_K <- enrichKEGG(gene = filter(convertEnsEnt, SYMBOL %in% EarlyInduced_SCClncRNAs_CoUp_Targets$EnsName.y)$ENTREZID,
                          universe = filter(convertEnsEnt, SYMBOL %in% fpkm_PCG$EnsName)$ENTREZID,
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)
CoUp_DE2_K_df <- data.frame(CoUp_DE2_K)
unique(CoUp_DE2_K$Description)
filter(convertEnsEnt, 
       ENTREZID %in% unlist(strsplit(filter(CoUp_DE2_K_df, Description == unique(CoUp_DE2_K_df$Description)[7])$geneID, 
                                     split = "/")))
#essentially identifies same genes pathways

#REACTOME
library(ReactomePA)
CoUp_DE2_R <- enrichPathway(gene          = unique(filter(convertEnsEnt,SYMBOL %in% EarlyInduced_SCClncRNAs_CoUp_Targets$EnsName.y)$ENTREZID),
                             universe      = unique(filter(convertEnsEnt,SYMBOL %in% fpkm_PCG$EnsName)$ENTREZID),
                             organism = "human",
                             pvalueCutoff = 0.05,
                             qvalueCutoff  = 0.05,
                             readable      = TRUE)
CoUp_DE2_R_df <- as.data.frame(CoUp_DE2_R)
#as previous

#GO terms work well:
CoUp_DE2_dfi <- data.frame(simplify(CoUp_DE2))
CoUp_DE2_dfi$selectHits <- as.numeric(sapply(strsplit(CoUp_DE2_dfi$GeneRatio, "\\/"), "[[", 1))
CoUp_DE2_dfi$select <- as.numeric(sapply(strsplit(CoUp_DE2_dfi$GeneRatio, "\\/"), "[[", 2))
CoUp_DE2_dfi$geneRatio <- CoUp_DE2_dfi$selectHits/CoUp_DE2_dfi$select*100

CoUp_DE2_dfi_MF <- filter(CoUp_DE2_dfi, ONTOLOGY == "MF")

CoUp_DE2_dfi_MF <- CoUp_DE2_dfi_MF[order(CoUp_DE2_dfi_MF$Description),]

#CoUp_DE2_df_MF$DescriptionII <- stringr::str_wrap(CoUp_DE2_df$Description, width = 20)

CoUp_DE2_dfi_MF$Description <- factor(CoUp_DE2_dfi_MF$Description, labels = CoUp_DE2_dfi_MF$Description)
CoUp_DE2_dfi_MF$Description <- factor(CoUp_DE2_dfi_MF$Description,
                                           levels = levels(CoUp_DE2_dfi_MF$Description)[order(CoUp_DE2_dfi_MF$geneRatio, decreasing = F)])

CoUp_DE2_dfi_MF <- CoUp_DE2_dfi_MF[,c(3,10,16)]

colnames(CoUp_DE2_dfi_MF)[3] <- "% SCClncRNA targets with GO term"

ggplot(CoUp_DE2_dfi_MF) + aes(x =`% SCClncRNA targets with GO term`, fill = -log10(p.adjust), y = Description) +
  geom_bar(stat = "identity", color = "grey60")+
  theme_bw() +
  scale_x_continuous(breaks = seq(0,30,5)) +
  theme(text = element_text(size=24)) +
  xlab("% SCClncRNA targets\nwith GO term") +
  ylab("")

#
#### Blunter checks on TFs, IEGs, Cell cycle genes, SMC-biased genes ####

#GO/KEGG/REACTOME reveals little unexpected, though can provide some info, check more specific lists

#TFs and CCs already annotated, others of interest here:

#IEGs (inc. lots of immune)
IEGs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Gene lists/arner_2015_table_S5_IEGs_lit.csv")
IEGs_hs <- filter(IEGs, !Hs_symbol %in% c("", "-"))
IEGs_hs <- filter(IEGs_hs, !grepl("[a-z]", IEGs_hs$Total_count))
IEGs_hs$Hs_symbol[IEGs_hs$Hs_symbol == "IL8"] <- "CXCL8"

#SMC-biased
FANT_S10 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/nature21374-s2/table 10.csv")
FANT_S10_SMC <- filter(FANT_S10, grepl("smooth", sample_ontology_term))
FANT_S10_SMC_G <- gsub("\\.[0-9]*", "", unlist(strsplit(FANT_S10_SMC$associated_geneID, ",")))

FANT_S10_VSMC <- FANT_S10_SMC[-c(1,12:17),]
FANT_S10_VSMC_G <- gsub("\\.[0-9]*", "", unlist(strsplit(FANT_S10_VSMC$associated_geneID, ",")))

#Epigenetic modifiers (LISA has some definitions but not well described, CRdb seems recent and fine)
CRdb_dataB <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CRdb Data Browse.csv")

#add HMGA2 and FOXL1 to this list, unclear why not in there...

GeneLists <- list("TFs"= unique(filter(fpkm_allG, grepl("TF", GeneClassUpdate))$EnsID),
                  "CRs" = unique(filter(fpkm_allG, EnsName %in% c(CRdb_dataB$CR, "FOXL1", "HMGA2"))$EnsID),
                  "CC"= unique(filter(fpkm_allG, grepl("CC", GeneClassUpdate))$EnsID),
                  "IEGs"= unique(filter(fpkm_allG, EnsName %in% c(IEGs_hs[,4]))$EnsID),
                  "SMC"= unique(filter(fpkm_allG, gsub("\\.[0-9]*", "", EnsID) %in% FANT_S10_SMC_G)$EnsID),
                  "VSMC"= unique(filter(fpkm_allG, gsub("\\.[0-9]*", "", EnsID) %in% FANT_S10_VSMC_G)$EnsID)
                  )

results_list <- list()

#same/later timeframe:
EarlyInduced_SCClncRNAs_Targets <- filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs")

for (i in 1:length(GeneLists)){
  
  a <- length(unique(EarlyInduced_SCClncRNAs_Targets$EnsID.y[ EarlyInduced_SCClncRNAs_Targets$EnsID.y %in% GeneLists[[i]] ]))
  b <- length(unique(EarlyInduced_SCClncRNAs_Targets$EnsID.y))
  c <- length(unique(fpkm_PCGDE$EnsID[ fpkm_PCGDE$EnsID %in% GeneLists[[i]]]))
  d <- length(unique(fpkm_PCGDE$EnsID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate,
                                  paste(unique(EarlyInduced_SCClncRNAs_Targets$EnsName.y[ EarlyInduced_SCClncRNAs_Targets$EnsID.y %in% GeneLists[[i]] ]), collapse = "/"))
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:7] <- c("pval", "OR")
#cut non-vasc SMC
GeneLists_enriched <- GeneLists_enriched[-5,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")
#seems like essentially all these groups (aside from cell cycle) are enriched to quite strong degree 
GeneLists_enrichedSameLaterT_T2 <- GeneLists_enriched


#with non-T2 the result is much weaker right?
CoRegPairs_04_48_24_extended <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <250,
                                       #AllLNC_AllPCG_1,
                                       (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                     fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                       fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                        fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                       fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
)

EarlyInduced_CClncRNAs_Targets <- unique(filter(CoRegPairs_04_48_24_extended,
                                                EnsID %in% fpkm_allGDE_Upwithin_4$EnsID))

for (i in 1:length(GeneLists)){
  
  a <- length(unique(EarlyInduced_CClncRNAs_Targets$EnsID.y[ EarlyInduced_CClncRNAs_Targets$EnsID.y %in% GeneLists[[i]] ]))
  b <- length(unique(EarlyInduced_CClncRNAs_Targets$EnsID.y))
  c <- length(unique(fpkm_PCGDE$EnsID[ fpkm_PCGDE$EnsID %in% GeneLists[[i]]]))
  d <- length(unique(fpkm_PCGDE$EnsID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate,
                                  paste(unique(EarlyInduced_CClncRNAs_Targets$EnsName.y[ EarlyInduced_CClncRNAs_Targets$EnsID.y %in% GeneLists[[i]] ]), collapse = "/"))
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:7] <- c("pval", "OR")
#cut non-vasc SMC
GeneLists_enriched <- GeneLists_enriched[-5,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")
#weaker/non-sig results without T2
GeneLists_enrichedSameLaterT <- GeneLists_enriched


#not sure if need a same timeframe result too:
EarlyInduced_SCClncRNAs_Targets <- filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs", PCG_Cluster == "Induced <4hrs")

for (i in 1:length(GeneLists)){
  
  a <- length(unique(EarlyInduced_SCClncRNAs_Targets$EnsID.y[ EarlyInduced_SCClncRNAs_Targets$EnsID.y %in% GeneLists[[i]] ]))
  b <- length(unique(EarlyInduced_SCClncRNAs_Targets$EnsID.y))
  c <- length(unique(fpkm_PCGDE_Upwithin_4$EnsID[ fpkm_PCGDE_Upwithin_4$EnsID %in% GeneLists[[i]]]))
  d <- length(unique(fpkm_PCGDE_Upwithin_4$EnsID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate,
                                  paste(unique(EarlyInduced_SCClncRNAs_Targets$EnsName.y[ EarlyInduced_SCClncRNAs_Targets$EnsID.y %in% GeneLists[[i]] ]), collapse = "/"))
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:7] <- c("pval", "OR")
#cut non-vasc SMC
GeneLists_enriched <- GeneLists_enriched[-5,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")
#seems like essentially all these groups (aside from cell cycle) are enriched to quite strong degree 
GeneLists_enrichedSameT_T2 <- GeneLists_enriched


#using just CClncRNAs:
CoRegPairs_04_48_24_extended <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <250,
                                       #AllLNC_AllPCG_1,
                                       (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                     fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                       fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                        fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                       fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
)
EarlyInduced_CClncRNAs_Targets <- unique(filter(CoRegPairs_04_48_24_extended,
                                                EnsID %in% fpkm_allGDE_Upwithin_4$EnsID, EnsID.y %in% fpkm_allGDE_Upwithin_4$EnsID))

for (i in 1:length(GeneLists)){
  
  a <- length(unique(EarlyInduced_CClncRNAs_Targets$EnsID.y[ EarlyInduced_CClncRNAs_Targets$EnsID.y %in% GeneLists[[i]] ]))
  b <- length(unique(EarlyInduced_CClncRNAs_Targets$EnsID.y))
  c <- length(unique(fpkm_PCGDE_Upwithin_4$EnsID[ fpkm_PCGDE_Upwithin_4$EnsID %in% GeneLists[[i]]]))
  d <- length(unique(fpkm_PCGDE_Upwithin_4$EnsID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate,
                                  paste(unique(EarlyInduced_CClncRNAs_Targets$EnsName.y[ EarlyInduced_CClncRNAs_Targets$EnsID.y %in% GeneLists[[i]] ]), collapse = "/"))
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:7] <- c("pval", "OR")
#cut non-vasc SMC
GeneLists_enriched <- GeneLists_enriched[-5,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")
#seems like essentially all these groups (aside from cell cycle) are enriched to quite strong degree 
GeneLists_enrichedSameT <- GeneLists_enriched


#plot enrichment of interesting subtypes in same tier2:
GeneLists_enrichedSameT_T2$selectHit <- GeneLists_enrichedSameT_T2$a/GeneLists_enrichedSameT_T2$b
GeneLists_enrichedSameT_T2$backHit <- GeneLists_enrichedSameT_T2$c/GeneLists_enrichedSameT_T2$d

library(reshape2)
plot_Gclass <- melt(GeneLists_enrichedSameT_T2[,c(1,10,11)])

plot_Gclass$GeneList <- c("TFs", "Chromatin Regulators", "Core S/G2M", "IEGs", "VSMC-enriched", 
                          "TFs", "Chromatin Regulators", "Core S/G2M", "IEGs", "VSMC-enriched")
plot_Gclass$GeneList <- as.factor(plot_Gclass$GeneList)
plot_Gclass$GeneList <- factor(plot_Gclass$GeneList, levels = levels(plot_Gclass$GeneList)[c(2,1,3,5,4)])

plot_Gclass$variable <- c(rep("0-4hr co-induced with SCClncRNA", 5), rep("All 0-4hr induced PCGs", 5))
plot_Gclass$variable <- as.factor(plot_Gclass$variable)
plot_Gclass$variable <- factor(plot_Gclass$variable, levels = levels(plot_Gclass$variable)[c(2,1)])

plot_Gclass$value <- plot_Gclass$value*100

ggplot(plot_Gclass) + aes(y = GeneList, x = value, fill = variable) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c(`0-4hr co-induced with SCClncRNA` = "mediumorchid",`All 0-4hr induced PCGs` = "grey60")) +
  theme_minimal() +
  ylab("") +
  xlab("%") +
  theme(text = element_text(size=24))

#plot enrichment of interesting subtypes in same/later tier2:
GeneLists_enrichedSameT$selectHit <- GeneLists_enrichedSameT$a/GeneLists_enrichedSameT$b
GeneLists_enrichedSameT$backHit <- GeneLists_enrichedSameT$c/GeneLists_enrichedSameT$d

library(reshape2)
plot_Gclass <- melt(GeneLists_enrichedSameT[,c(1,10,11)])

plot_Gclass$GeneList <- c("TFs", "Chromatin Regulators", "Core S/G2M", "IEGs", "VSMC-enriched", 
                          "TFs", "Chromatin Regulators", "Core S/G2M", "IEGs", "VSMC-enriched")
plot_Gclass$GeneList <- as.factor(plot_Gclass$GeneList)
plot_Gclass$GeneList <- factor(plot_Gclass$GeneList, levels = levels(plot_Gclass$GeneList)[c(2,1,3,5,4)])

plot_Gclass$variable <- c(rep("0-4hr co-induced with lncRNA", 5), rep("All 0-4hr induced PCGs", 5))
plot_Gclass$variable <- as.factor(plot_Gclass$variable)
plot_Gclass$variable <- factor(plot_Gclass$variable, levels = levels(plot_Gclass$variable)[c(2,1)])

plot_Gclass$value <- plot_Gclass$value*100

ggplot(plot_Gclass) + aes(y = GeneList, x = value, fill = variable) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c(`0-4hr co-induced with lncRNA` = "mediumorchid",`All 0-4hr induced PCGs` = "grey60")) +
  theme_minimal() +
  ylab("") +
  xlab("%") +
  theme(text = element_text(size=24))

