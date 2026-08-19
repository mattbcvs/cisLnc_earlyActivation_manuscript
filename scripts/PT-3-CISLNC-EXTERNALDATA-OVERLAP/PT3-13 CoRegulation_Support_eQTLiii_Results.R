#### overview of all tests + selection of eQTL validated lncRNAs ####
library(dplyr)
library(ggplot2)

GTEX_SuppTable_df <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX_SuppTable_df250.csv")
Shu_SuppTable_df <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Shu_SuppTable_df_250_p025.csv")

eQTL_SuppTable <- rbind(GTEX_SuppTable_df, Shu_SuppTable_df)

#remove "all lnc" clusters - useless for now (too hetergenous in terms of representation across tissue is current hypo)
eQTL_SuppTable <- filter(eQTL_SuppTable, !grepl("_All", Run))

#save for supplement:
#write.csv(eQTL_SuppTable, "eQTL_SuppTable_2026.csv", row.names = F)


#
#### no. successful tests per lnc cluster (Fig.3E) ####

#split by lnc timing:
eQTL_SuppTable$LncCluster <- sapply(strsplit(eQTL_SuppTable$Run, "_"), "[[", 2)
eQTL_SuppTable$LncCluster <- gsub("Up", "Induced ", eQTL_SuppTable$LncCluster)
eQTL_SuppTable$LncCluster <- gsub("Down", "Repressed ", eQTL_SuppTable$LncCluster)
eQTL_SuppTable$LncCluster <- gsub(" 4", " 0-4hr", eQTL_SuppTable$LncCluster)
eQTL_SuppTable$LncCluster <- gsub(" 8", " 4-8hr", eQTL_SuppTable$LncCluster)
eQTL_SuppTable$LncCluster <- gsub(" 24", " 8-24hr", eQTL_SuppTable$LncCluster)
eQTL_SuppTable$LncCluster <- gsub("Del", "\\.", eQTL_SuppTable$LncCluster)
table(eQTL_SuppTable$LncCluster)
eQTL_SuppTable$LncCluster <- as.factor(eQTL_SuppTable$LncCluster)
eQTL_SuppTable$LncCluster <- factor(eQTL_SuppTable$LncCluster, levels(eQTL_SuppTable$LncCluster)[c(1,3,4,2,5,7,8,6)])

#split by overlap type:
eQTL_SuppTable$OverlapType <- as.factor(sapply(strsplit(eQTL_SuppTable$Run, "_"), "[[", 3))
eQTL_SuppTable$OverlapType <- factor(eQTL_SuppTable$OverlapType, levels(eQTL_SuppTable$OverlapType)[c(3,4,1,2,5,6)])

successTests <- unique(filter(eQTL_SuppTable, p < 0.05))#39 


#good for supplement:
ggplot(successTests) + aes(x = OverlapType, fill = OverlapType) +
  geom_bar() + 
  theme_bw() +
  theme(text = element_text(size = 16)) +
  facet_wrap(~LncCluster, drop=F, ncol = 4) +
  scale_x_discrete(drop=F) +
  Seurat::RotatedAxis() +
  xlab("") + 
  ylab("No. eQTL Datasets predictive\nof CClncRNA-target pairs")

#no tests in other lnc clusters, focus on 0-4 up
ggplot(filter(successTests, LncCluster == "Induced 0-4hr")) + aes(x = OverlapType, fill = OverlapType) +
  geom_bar() + 
  theme_bw() +
  theme(text = element_text(size = 18)) +
  #facet_wrap(~LncCluster, drop=F, ncol = 4) +
  scale_x_discrete(drop=F) +
  Seurat::RotatedAxis() +
  xlab("") + 
  ylab("No. eQTL Datasets predictive\nof CClncRNA-target pairs")

eQTL_SuppTable$ORii <- eQTL_SuppTable$OR
eQTL_SuppTable$ORii[eQTL_SuppTable$ORii == Inf] <- max(eQTL_SuppTable$ORii[!eQTL_SuppTable$ORii == Inf], na.rm = T)
eQTL_SuppTable$ORii[eQTL_SuppTable$OR >10] <- 10
eQTL_SuppTable$ORii[eQTL_SuppTable$p > 0.05] <- NA

ggplot(eQTL_SuppTable) + aes(x = OverlapType, y = -log10(p)+0.001, color = ORii) +
  geom_jitter(alpha = 0.6, width =0.06, size = 1.5) + 
  theme_bw() +
  facet_wrap(~LncCluster, drop=F, ncol = 4) +
  scale_color_gradient2(low = "steelblue", mid = "orange", high = "firebrick4", na.value = "grey") +
  scale_x_discrete(drop=F) +
  Seurat::RotatedAxis() +
  geom_hline(yintercept = -log10(0.05+0.001), linetype = "dashed", color = "grey60")+ 
  ylab("eQTL Enrichment in\nCClncRNA-targets (-log10p)") +
  xlab("eQTL-lncRNA Overlap Type")


#
#### confirmed pairs (requires codes PT2-11/PT2-12) ####

#reimport for ease
GTEX_SuppTable_df <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX_SuppTable_df250.csv")
Shu_SuppTable_df <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Shu_SuppTable_df_250_p025.csv")

eQTL_SuppTable <- rbind(GTEX_SuppTable_df, Shu_SuppTable_df)

#remove "all lnc" clusters - useless for now (too hetergenous in terms of representation across tissue is current hypo)
eQTL_SuppTable <- filter(eQTL_SuppTable, !grepl("_All", Run))

#split by overlap type:
eQTL_SuppTable$OverlapType <- as.character(sapply(strsplit(eQTL_SuppTable$Run, "_"), "[[", 3))
eQTL_SuppTable$OverlapType_grep <- as.character(eQTL_SuppTable$OverlapType)
eQTL_SuppTable$OverlapType_grep[grepl("Locus", eQTL_SuppTable$OverlapType)] <- "Promoter|Exon|Intron|TTS|Splice"

#split by lnc timing:
eQTL_SuppTable$LncCluster <- sapply(strsplit(eQTL_SuppTable$Run, "_"), "[[", 2)

#runs to use to select candidates
#optional: including a requirement for certain number of eQTL linked pairs to be found
#trialling without
#and a quite decent p considering lack of adjustment
#successTests_GTEX <- unique(filter(eQTL_SuppTable, p < 0.05, c >=5, !tissue == "SMC_biobank"))#24
#successTests_GTEX <- unique(filter(eQTL_SuppTable, p < 0.05, c >=4, !tissue == "SMC_biobank"))
successTests_GTEX <- unique(filter(eQTL_SuppTable, p < 0.05, !tissue == "SMC_biobank")) #35 

#re-confirm pairs for pairs amongst 1mbp:
AllLNC_AllPCG_2d3d <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_2d3d_2026.csv")
AllLNC_AllPCG_2d3d <- filter(AllLNC_AllPCG_2d3d, pair_range < 1000000)
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
260/8116 #numbers for within 1Mbp

#250kbp pairs used ultimately:
AllLNC_AllPCG_2d3d_250 <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <250)
dim(filter(AllLNC_AllPCG_2d3d_250, !loopMethod == "Neither"))
185/2563

#co-reg genes
GTEX_runs_selection <- list("Up4" = filter(AllLNC_AllPCG_2d3d_250,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                        fpkm_allGDE_Downwithin_4$EnsID))),
                            "Down4" = filter(AllLNC_AllPCG_2d3d_250,
                                             (EnsID %in% c(fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                            fpkm_allGDE_Downwithin_4$EnsID))),
                            "Up8" = filter(AllLNC_AllPCG_2d3d_250,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                        fpkm_allGDE_Downwithin_8$EnsID))),
                            "Down8" = filter(AllLNC_AllPCG_2d3d_250,
                                             (EnsID %in% c(fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                            fpkm_allGDE_Downwithin_8$EnsID))),
                            "Up24" = filter(AllLNC_AllPCG_2d3d_250,
                                            (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                          fpkm_allGDE_Downwithin_24$EnsID))),
                            "Down24" = filter(AllLNC_AllPCG_2d3d_250,
                                              (EnsID %in% c(fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                              fpkm_allGDE_Downwithin_24$EnsID))),
                            "Up4Del" = filter(AllLNC_AllPCG_2d3d_250,
                                              (EnsID %in% fpkm_allGDE_Upwithin_4$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                        fpkm_allGDE_within_24$EnsID))),
                            "Down4Del" = filter(AllLNC_AllPCG_2d3d_250,
                                                (EnsID %in% fpkm_allGDE_Downwithin_4$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                            fpkm_allGDE_within_24$EnsID))),
                            "All-same" = filter(AllLNC_AllPCG_2d3d_250,
                                                (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                              fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                               fpkm_allGDE_Downwithin_4$EnsID)) |
                                                  (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                                 fpkm_allGDE_Downwithin_8$EnsID)) |
                                                  (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                                  fpkm_allGDE_Downwithin_24$EnsID))),
                            "All-delayed"= filter(AllLNC_AllPCG_2d3d_250,
                                                  (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                                 fpkm_allGDE_within_24$EnsID)) |
                                                    (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                  fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_within_24$EnsID)))
                            )

eQTLconfirmed_pairs <- list()

for (i in 1:length(successTests_GTEX$tissue)){
  
  triali <- AllLNC_AllPCG_2d3d_250
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(successTests_GTEX$tissue[i], tissueType), 
                                                        grepl(successTests_GTEX$OverlapType_grep[i], OverlapType), 
                                                        pval_nominal < 0.05,
                                                        TotalOverlaps == 1)$pairs] <- "Yes"
  
  #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
  triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[successTests_GTEX$tissue[i]]], 
                                               pval_nominal < 0.05)$gene_id)
  
  #remove anything un-expressed >1 TPM in this tissue:
  findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", successTests_GTEX$tissue[i])), 
                                                      gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
  table(findExprs)
  triali <- filter(triali, EnsID.y %in% GTEX_exprs[findExprs, 1])
  
  triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% GTEX_runs_selection[[successTests_GTEX$LncCluster[i]]]$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  DELNC_DEPCG_1_eQTL$tissueName <- successTests_GTEX$tissue[i]
  DELNC_DEPCG_1_eQTL$overlapType <- successTests_GTEX$OverlapType[i]
  DELNC_DEPCG_1_eQTL$LncCluster <- successTests_GTEX$LncCluster[i]
  eQTLconfirmed_pairs[[i]] <- unique(DELNC_DEPCG_1_eQTL[,c(3,20:23)])
}


#### add pairs for Shu Biobank ####

successTests_Shu <- filter(eQTL_SuppTable, tissue == "SMC_biobank", p < 0.05)

AllLNC_AllPCG_2d3d_250$pairs_merge <- gsub("\\.[0-9]*$", "", AllLNC_AllPCG_2d3d_250$pairs)
AllLNC_AllPCG_2d3d_250$EnsID_merge <- gsub("\\.[0-9]*", "", AllLNC_AllPCG_2d3d_250$EnsID.y)
Shu_LncVar$pairs <- paste(Shu_LncVar$EnsID, Shu_LncVar$EnsID.y, sep = "-")
Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]

Shuconfirmed_pairs <- list()

for (i in 1:length(successTests_Shu$tissue)){
  
  triali <- filter(AllLNC_AllPCG_2d3d_250, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% filter(Shu_allVar_p, pvalue < 0.025)$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < 0.025,
                                                              grepl(successTests_Shu$OverlapType_grep[i], OverlapType), 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% GTEX_runs_selection[[successTests_Shu$LncCluster[i]]]$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  DELNC_DEPCG_1_eQTL$tissueName <- successTests_Shu$tissue[i]
  DELNC_DEPCG_1_eQTL$overlapType <- successTests_Shu$OverlapType[i]
  DELNC_DEPCG_1_eQTL$LncCluster <- successTests_Shu$LncCluster[i]
  Shuconfirmed_pairs[[i]] <- unique(DELNC_DEPCG_1_eQTL[,c(3,22:25)])
}


#### add columns to table, save all cclnc eQTL pairs #### 

CoRegPairs_eQTL_checked <- filter(AllLNC_AllPCG_2d3d_250,
                                     (EnsID %in% fpkm_allGDE_within_4$EnsID & EnsID.y %in% fpkm_allGDE$EnsID) |
                                       (EnsID %in% fpkm_allGDE_within_8$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                               fpkm_allGDE_within_24$EnsID)) |
                                       (EnsID %in% fpkm_allGDE_within_24$EnsID & EnsID.y %in% fpkm_allGDE_within_24$EnsID))

trial <- CoRegPairs_eQTL_checked

for (i in 1:length(eQTLconfirmed_pairs)){
  trial$newCol <- FALSE
  trial$newCol[trial$pairs %in% eQTLconfirmed_pairs[[i]]$pairs] <- TRUE
  colnames(trial)[length(colnames(trial))] <- paste(#eQTLconfirmed_pairs[[i]][1,2],
    eQTLconfirmed_pairs[[i]][1,3],
    eQTLconfirmed_pairs[[i]][1,4],
    eQTLconfirmed_pairs[[i]][1,5],
    sep = "_")
}

for (i in 1:length(Shuconfirmed_pairs)){
  trial$newCol <- FALSE
  trial$newCol[trial$pairs %in% Shuconfirmed_pairs[[i]]$pairs] <- TRUE
  colnames(trial)[length(colnames(trial))] <- paste(#eQTLconfirmed_pairs[[i]][1,2],
    Shuconfirmed_pairs[[i]][1,3],
    Shuconfirmed_pairs[[i]][1,4],
    Shuconfirmed_pairs[[i]][1,5],
    sep = "_")
}

CoRegPairs_eQTL_checked <- trial
colnames(CoRegPairs_eQTL_checked)
CoRegPairs_eQTL_checked$eQTLvalidations <- rowSums(
  CoRegPairs_eQTL_checked[,c(
    22:dim(CoRegPairs_eQTL_checked)[2]
  )])


CoRegPairs_eQTL_supported <- filter(CoRegPairs_eQTL_checked, eQTLvalidations > 0)

#write.csv(CoRegPairs_eQTL_supported, "CoRegPairs_eQTL_supported.csv",row.names = F)

# some stats
length(unique(CoRegPairs_eQTL_supported$pairs)) #67 eQTL-validated pairs
length(unique(CoRegPairs_eQTL_supported$EnsID)) #43 lncRNAs

# distance of eQTL-confirmed pairs:
summary(CoRegPairs_eQTL_supported$AbsDistLnc_PCG) #90kbp median
summary(CoRegPairs_eQTL_checked$AbsDistLnc_PCG) #120kbp median, bit longer


#### overlap identified eQTL-supported pairs with those from FANTOM analysis (Fig.3F) ####

CoRegPairs_eQTL_supported <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CoRegPairs_eQTL_supported.csv")

fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE_2026filt.csv", header = T)
length(unique(fpkm_allGDE$EnsID))

fpkm_allGDE_Upwithin_4 <- filter(fpkm_allGDE, RegulationStart == "Induced <4hrs")
fpkm_allGDE_Downwithin_4 <- filter(fpkm_allGDE, RegulationStart == "Repressed <4hrs")
fpkm_allGDE_within_4 <- rbind(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Downwithin_4)

fpkm_allGDE_Upwithin_8 <- filter(fpkm_allGDE, RegulationStart == "Induced 4-8hrs")
fpkm_allGDE_Downwithin_8 <- filter(fpkm_allGDE, RegulationStart == "Repressed 4-8hrs")
fpkm_allGDE_within_8 <- rbind(fpkm_allGDE_Upwithin_8, fpkm_allGDE_Downwithin_8)

fpkm_allGDE_Upwithin_24 <- filter(fpkm_allGDE, RegulationStart == "Induced 8-24hrs")
fpkm_allGDE_Downwithin_24 <- filter(fpkm_allGDE, RegulationStart == "Repressed 8-24hrs")
fpkm_allGDE_within_24 <- rbind(fpkm_allGDE_Upwithin_24, fpkm_allGDE_Downwithin_24)

CoRegPairs_eQTL_checked <- filter(AllLNC_AllPCG_2d3d_250,
                                  (EnsID %in% fpkm_allGDE_within_4$EnsID & EnsID.y %in% fpkm_allGDE$EnsID) |
                                    (EnsID %in% fpkm_allGDE_within_8$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                            fpkm_allGDE_within_24$EnsID)) |
                                    (EnsID %in% fpkm_allGDE_within_24$EnsID & EnsID.y %in% fpkm_allGDE_within_24$EnsID))
#compare to Hon 2017
#slight diff method - they looked for any eQTL-linked pairs which were significantly more co-expressed 
#than background shuffled pairs
#within similar distance, with same orientation
#correlation across all 1.8k samples - non-specific

#we compared to a set of eQTL linked lncRNA-mRNA from previous analysis
#which was a) more generic, not focused on confirming co-expression in one cell type
#b) looking only for linear correlation (not timing based synchronising like done here)

FANTOM_eQTL_pairs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/nature21374-s2/table16.csv")
FANTOM_eQTL_pairs_cis <- filter(FANTOM_eQTL_pairs, cis_correlated_candidate == "yes")
#5264 pairs - same as reported in Hon et al 2017

#pairs found as neighbours:
#add FANTOM ID to Alllnc:
trial <- unique(merge(AllLNC_AllPCG_2d3d_250, Enhancer_lociII[,c(1,14)], by = "EnsID"), all.x = T)

#pair up columns - check via EnsID for lnc
trial$pairs_noSuff <- paste(gsub("\\.[0-9]*", "", trial$EnsID), gsub("\\.[0-9]*", "", trial$EnsID.y), sep = "-")
#as well as FANT-ID for lnc - maybe subtle diffs
trial$FANTpairs_noSuff <- paste(gsub("\\.[0-9]*", "", trial$FANTOM_ID), gsub("\\.[0-9]*", "", trial$EnsID.y), sep = "-")

FANTOM_eQTL_pairs_cis$pairs_noSuff <- paste(gsub("\\.[0-9]*", "", FANTOM_eQTL_pairs_cis$lncRNA_geneID), 
                                            gsub("\\.[0-9]*", "", FANTOM_eQTL_pairs_cis$mRNA_geneID), sep = "-")

#no. pairs found in dataset:
triali <- filter(trial, FANTpairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff |
                   pairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff)
length(unique(triali$pairs))#183 pairs where both found expressed in SVSMC timecourse with TSS in 250kbp, whole locus in 1mbp

#no. pairs in the co-reg set checked for eQTLs:
triali <- filter(trial, pairs %in% CoRegPairs_eQTL_checked$pairs, 
                 FANTpairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff |
                   pairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff)
length(unique(triali$pairs))#15 found both expressed in SVSMC and in a co-reg pairing

#confirmed by this eQTL analysis? from any tissue at any lnc site:
triali <- filter(trial, pairs %in% CoRegPairs_eQTL_supported$pairs, 
                 FANTpairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff |
                   pairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff)
length(unique(triali$pairs))#7 recovered

#fantom pairs make up
7/64 # 10.4% of our eQTL-validated co regs
15/282 # 5.3% of the testable co-regs

#nice overlap

a <- 7
b <- 67
c <- 15
d <- 282

fisher.test(data.frame("cisLnc" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#p = 0.04, OR = 3


DEL_PCG_type <- rbind("eQTL-supported CClncRNA-\n-target pairs (this study)" = c(a,b, a/b),
                      "All CClncRNA-target pairs" = c(c,d, c/d))
colnames(DEL_PCG_type) <- c("eQTL-supported (FANTOM)", "Other", "%")

library(reshape2)
DEL_PCG_typei <- melt(DEL_PCG_type)[5:6,]

DEL_PCG_typei$value <- DEL_PCG_typei$value*100

DEL_PCG_typei$Var1 <- as.factor(DEL_PCG_typei$Var1)
DEL_PCG_typei$Var1 <- factor(DEL_PCG_typei$Var1, levels = levels(DEL_PCG_typei$Var1)[c(2,1)])

ggplot(DEL_PCG_typei) + aes(y = Var1, x = value) + 
  geom_bar(stat = "identity", color = "grey60") +
  xlab("% eQTL-supported\n(FANTOM, Hon et al. 2017)") +
  
  #xlab("") +
  ylab("") +
  theme_minimal() +
  scale_x_continuous(breaks = c(0,5,10)) +
  theme(text = element_text(size=24))# + Seurat::RotatedAxis()



#
#### pairs per tissue heatmap (Fig.S5) ####

#visualise which tissues have best predictive value - number of pairs recovered

#lncRNA names col:
CoRegPairs_eQTL_supported$lncNames <- CoRegPairs_eQTL_supported$EnsName.x
CoRegPairs_eQTL_supported$lncNames[grepl("MSTRG", CoRegPairs_eQTL_supported$EnsID)] <- 
  CoRegPairs_eQTL_supported$EnsID[grepl("MSTRG", CoRegPairs_eQTL_supported$EnsID)]

#nicely named pair col:
CoRegPairs_eQTL_supported$pairsNice <- paste(CoRegPairs_eQTL_supported$lncNames, 
                                                   CoRegPairs_eQTL_supported$EnsName.y, sep = "-")

#find cols for a pheatmap
colnames(CoRegPairs_eQTL_supported)

mat <- lapply(CoRegPairs_eQTL_supported[,22:60], as.numeric)
mat <- as.matrix(bind_rows(mat))
rownames(mat) <- CoRegPairs_eQTL_supported$pairsNice

pheatmap::pheatmap(mat, angle_col = 315, border_color = "grey60", fontsize_row = 6, fontsize_col = 12, cluster_rows = F, 
                   legend = F)

#certain pairs found in certain tissues
#some pairs widely validated

#aggregate tissue overlap - remove cluster not useful info:
cluster_2_agg <- unique(sapply(strsplit(colnames(mat), "_"), function(x){
  paste(x[1:(length(x)-1)], collapse = "_")
}))

cluster_2_agg_res <- list()

for(i in 1:length(cluster_2_agg)){
  cluster_2_agg_res[[i]] <- as.numeric(rowSums(as.matrix(mat[,grep(cluster_2_agg[i], colnames(mat))])) >0)
  }

names(cluster_2_agg_res) <- cluster_2_agg

mat_clust <- as.matrix(bind_cols(cluster_2_agg_res))
rownames(mat_clust) <- CoRegPairs_eQTL_supported$pairsNice

#put the cluster info instead into an annotation for the rows:
GTEX_runs_selection <- list("Up4" = filter(AllLNC_AllPCG_2d3d_250,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                        fpkm_allGDE_Downwithin_4$EnsID))),
                            "Down4" = filter(AllLNC_AllPCG_2d3d_250,
                                             (EnsID %in% c(fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                            fpkm_allGDE_Downwithin_4$EnsID))),
                            "Up8" = filter(AllLNC_AllPCG_2d3d_250,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                        fpkm_allGDE_Downwithin_8$EnsID))),
                            "Down8" = filter(AllLNC_AllPCG_2d3d_250,
                                             (EnsID %in% c(fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                            fpkm_allGDE_Downwithin_8$EnsID))),
                            "Up24" = filter(AllLNC_AllPCG_2d3d_250,
                                            (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                          fpkm_allGDE_Downwithin_24$EnsID))),
                            "Down24" = filter(AllLNC_AllPCG_2d3d_250,
                                              (EnsID %in% c(fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                              fpkm_allGDE_Downwithin_24$EnsID))),
                            "Up4Del" = filter(AllLNC_AllPCG_2d3d_250,
                                              (EnsID %in% fpkm_allGDE_Upwithin_4$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                        fpkm_allGDE_within_24$EnsID))),
                            "Down4Del" = filter(AllLNC_AllPCG_2d3d_250,
                                                (EnsID %in% fpkm_allGDE_Downwithin_4$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                            fpkm_allGDE_within_24$EnsID))),
                            "All-same" = filter(AllLNC_AllPCG_2d3d_250,
                                                (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                              fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                               fpkm_allGDE_Downwithin_4$EnsID)) |
                                                  (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                                 fpkm_allGDE_Downwithin_8$EnsID)) |
                                                  (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                                  fpkm_allGDE_Downwithin_24$EnsID))),
                            "All-delayed"= filter(AllLNC_AllPCG_2d3d_250,
                                                  (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                                 fpkm_allGDE_within_24$EnsID)) |
                                                    (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                  fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_within_24$EnsID)))
)

names(GTEX_runs_selection)

CoRegPairs_eQTL_supported$CClncRNA_Type[CoRegPairs_eQTL_supported$pairs %in% GTEX_runs_selection[[1]]$pairs] <- "Induced 0-4hr"
CoRegPairs_eQTL_supported$CClncRNA_Type[CoRegPairs_eQTL_supported$pairs %in% GTEX_runs_selection[[2]]$pairs] <- "Repressed 0-4hr"
CoRegPairs_eQTL_supported$CClncRNA_Type[CoRegPairs_eQTL_supported$pairs %in% GTEX_runs_selection[[3]]$pairs] <- "Induced 4-8hr"
CoRegPairs_eQTL_supported$CClncRNA_Type[CoRegPairs_eQTL_supported$pairs %in% GTEX_runs_selection[[4]]$pairs] <- "Repressed 4-8hr"
CoRegPairs_eQTL_supported$CClncRNA_Type[CoRegPairs_eQTL_supported$pairs %in% GTEX_runs_selection[[5]]$pairs] <- "Induced 8-24hr"
CoRegPairs_eQTL_supported$CClncRNA_Type[CoRegPairs_eQTL_supported$pairs %in% GTEX_runs_selection[[6]]$pairs] <- "Repressed 8-24hr"
CoRegPairs_eQTL_supported$CClncRNA_Type[CoRegPairs_eQTL_supported$pairs %in% GTEX_runs_selection[[7]]$pairs] <- "Induced 0-4hr"
CoRegPairs_eQTL_supported$CClncRNA_Type[CoRegPairs_eQTL_supported$pairs %in% GTEX_runs_selection[[8]]$pairs] <- "Repressed 0-4hr"
table(CoRegPairs_eQTL_supported$CClncRNA_Type)

CoRegPairs_eQTL_supported$PairType[CoRegPairs_eQTL_supported$pairs %in% GTEX_runs_selection[[1]]$pairs] <- "Same timeframe"
CoRegPairs_eQTL_supported$PairType[CoRegPairs_eQTL_supported$pairs %in% GTEX_runs_selection[[2]]$pairs] <- "Same timeframe"
CoRegPairs_eQTL_supported$PairType[CoRegPairs_eQTL_supported$pairs %in% GTEX_runs_selection[[3]]$pairs] <- "Same timeframe"
CoRegPairs_eQTL_supported$PairType[CoRegPairs_eQTL_supported$pairs %in% GTEX_runs_selection[[4]]$pairs] <- "Same timeframe"
CoRegPairs_eQTL_supported$PairType[CoRegPairs_eQTL_supported$pairs %in% GTEX_runs_selection[[5]]$pairs] <- "Same timeframe"
CoRegPairs_eQTL_supported$PairType[CoRegPairs_eQTL_supported$pairs %in% GTEX_runs_selection[[6]]$pairs] <- "Same timeframe"
CoRegPairs_eQTL_supported$PairType[CoRegPairs_eQTL_supported$pairs %in% GTEX_runs_selection[[7]]$pairs] <- "Later timeframe"
CoRegPairs_eQTL_supported$PairType[CoRegPairs_eQTL_supported$pairs %in% GTEX_runs_selection[[8]]$pairs] <- "Later timeframe"
table(CoRegPairs_eQTL_supported$PairType)

annotate_heatmap_rows <- CoRegPairs_eQTL_supported[,((dim(CoRegPairs_eQTL_supported)[2]-1):dim(CoRegPairs_eQTL_supported)[2])]
rownames(annotate_heatmap_rows) <- CoRegPairs_eQTL_supported$pairsNice

annotate_heatmap_rows

mat_clust_order <- mat_clust[order(annotate_heatmap_rows$PairType, annotate_heatmap_rows$CClncRNA_Type),]

pheatmap::pheatmap(mat_clust_order, annotation_row = annotate_heatmap_rows, cluster_rows = F,
                   angle_col = 315, border_color = "grey60", fontsize_row = 6, fontsize_col = 12, legend = F)


#aggregate overlap type, can maybe group the tissues/lncRNA pairs more easily
#tissues being profiled:
tiss_2_agg <- unique(sapply(strsplit(colnames(mat_clust), "_"), function(x){
  paste(x[1:length(x)-1], collapse = "_")
}))

tiss_2_agg_res <- list()

for(i in 1:length(tiss_2_agg)){
  tiss_2_agg_res[[i]] <- as.numeric(rowSums(as.matrix(mat_clust[,grep(tiss_2_agg[i], colnames(mat_clust))])) >0)
}

names(tiss_2_agg_res) <- tiss_2_agg

mat_tiss <- as.matrix(bind_cols(tiss_2_agg_res))
rownames(mat_tiss) <- CoRegPairs_eQTL_supported$pairsNice

pheatmap::pheatmap(mat_tiss, angle_col = 315, border_color = "grey60", fontsize_row = 7, cluster_rows = F)

mat_tiss_order <- mat_tiss[order(annotate_heatmap_rows$PairType, annotate_heatmap_rows$CClncRNA_Type),]

pheatmap::pheatmap(mat_tiss_order, annotation_row = annotate_heatmap_rows, cluster_rows = F, #cluster_cols = F,
                   angle_col = 315, border_color = "grey60", fontsize_row = 6, fontsize_col = 12, legend = F)

#order by tissues:
mat_tiss_order2 <- mat_tiss_order[,order(colSums(mat_tiss_order), decreasing = T)]

#put some breaks in for ease:
table(CoRegPairs_eQTL_supported$CClncRNA_Type)
table(CoRegPairs_eQTL_supported$PairType)

pheatmap::pheatmap(mat_tiss_order2, annotation_row = annotate_heatmap_rows, cluster_rows = F, cluster_cols = F, gaps_row = c(32,47,52,58,61),
                   angle_col = 315, border_color = "grey60", fontsize_row = 8, fontsize_col = 10.5, legend = F)


#aggregate overlap type
overlap_2_agg <- unique(sapply(strsplit(colnames(mat_clust), "_"), function(x){
  x[length(x)]
}))

overlap_2_agg_res <- list()

for(i in 1:length(overlap_2_agg)){
  overlap_2_agg_res[[i]] <- as.numeric(rowSums(as.matrix(mat_clust[,grep(overlap_2_agg[i], colnames(mat_clust))])) >0)
}

names(overlap_2_agg_res) <- overlap_2_agg

mat_overlap <- as.matrix(bind_cols(overlap_2_agg_res))
rownames(mat_overlap) <- CoRegPairs_eQTL_supported$pairsNice

mat_overlap_order <- mat_overlap[order(annotate_heatmap_rows$PairType, annotate_heatmap_rows$CClncRNA_Type),]

pheatmap::pheatmap(mat_overlap_order, annotation_row = annotate_heatmap_rows, angle_col = 315, border_color = "grey60", fontsize_row = 7, cluster_rows = F)

#order by overlaps:
mat_overlap_order2 <- mat_overlap_order[,order(colSums(mat_overlap_order), decreasing = T)]

pheatmap::pheatmap(mat_overlap_order2, annotation_row = annotate_heatmap_rows, angle_col = 315, border_color = "grey60", fontsize_row = 8, fontsize_col = 10.5, legend = F,
                   cluster_rows = F, cluster_cols = F, , gaps_row = c(32,47,52,58,61))


