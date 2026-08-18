library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)


#current tier2:
#HiC as well as eQTLs found at high enough threshold
CoRegPairs_04_48_24_extended <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CoRegPairs_04_48_24_extended_eQTL_Oct25.csv")
CoRegPairs_04_48_24_extended$lncNameCol <- CoRegPairs_04_48_24_extended$EnsID
CoRegPairs_04_48_24_extended$lncNameCol[!is.na(CoRegPairs_04_48_24_extended$EnsName.x)] <- CoRegPairs_04_48_24_extended$EnsName.x[!is.na(CoRegPairs_04_48_24_extended$EnsName.x)]

colnames(CoRegPairs_04_48_24_extended)
unique(CoRegPairs_04_48_24_extended$lncNameCol)

fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGClassUpdate.csv", header = T)

FPKM_CQV_OVERLAP_fpkm <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/FPKM_CQV_OVERLAP_fpkm.csv")
table(FPKM_CQV_OVERLAP_fpkm$IGV)#413 pass, 168 fail

#remove artefacts (push back to step 7?)
fpkm_allG_filt <- filter(fpkm_allG, grepl("chr", chr), !grepl("artefacts", GeneClassUpdate), !is.na(GeneClassUpdate))
#remove manual or Thresh4 fails
fpkm_allG_filt_manual <- filter(fpkm_allG_filt, 
                                !EnsID %in% filter(FPKM_CQV_OVERLAP_fpkm, IGV == "fail")$EnsID, #remove manual fails
)

fpkm_allG <- fpkm_allG_filt_manual
length(unique(fpkm_allG$EnsID))#10761

#
table(CoRegPairs_04_48_24_extended$eQTLvalidations)

#Correlation:
#Spearman's on FPKM
AllLNC_AllPCG_FPKM <- merge(CoRegPairs_04_48_24_extended[,c(2,1,3)], fpkm_allG[,c(2,11:26)], by.x = "EnsID", by.y = "EnsID")
AllLNC_AllPCG_FPKM <- merge(AllLNC_AllPCG_FPKM, fpkm_allG[,c(2,11:26)], by.x = "EnsID.y", by.y = "EnsID")
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

CoRegPairs_04_48_24_extended <- merge(trialiii, CoRegPairs_04_48_24_extended, by = "pairs")
dim(filter(CoRegPairs_04_48_24_extended, spear_p_adj < 0.05))#110

CoRegPairs_04_48_24_extended$corSig <- "No"
CoRegPairs_04_48_24_extended$corSig[abs(CoRegPairs_04_48_24_extended$spear_rho) >0.5 & 
                                      CoRegPairs_04_48_24_extended$spear_p_adj <0.05] <- "Yes"
table(CoRegPairs_04_48_24_extended$corSig)#still 110 (with or without the rho filter)

correlatedPairs <- filter(CoRegPairs_04_48_24_extended, corSig == "Yes")


#### additional SCClncRNAs from FANTOM ####

FANTOM_eQTL_pairs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/nature21374-s2/table16.csv")
FANTOM_eQTL_pairs_cis <- filter(FANTOM_eQTL_pairs, cis_correlated_candidate == "yes")
#5264 pairs as reported in Hon et al 2017

#pairs found as neighbours:
#pairs post-HiC
AllLNC_AllPCG_2d3d <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_2d3d_Aug2025.csv")
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
length(unique(triali$pairs))#143 pairs both found expressed in SVSMC

#part of the co-reg set:
triali <- filter(trial, pairs %in% CoRegPairs_04_48_24_extended$pairs, 
                 FANTpairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff |
                   pairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff)
length(unique(triali$pairs))#18 found both expressed in SVSMC and in a co-reg pairing

#any additional?
CoRegPairs_04_48_24_extended$FANTOM_eQTL <- "No"
CoRegPairs_04_48_24_extended$FANTOM_eQTL[CoRegPairs_04_48_24_extended$pairs %in% triali$pairs] <- "Yes"


#### known cis-acting lncRNAs ####

#any FANTOM match up here:
Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime.csv", header = T)

ControlCisLncs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Write-ups/supplement/ControlCisLncs.csv", header = T)

#73 genes in total
unique(ControlCisLncs$Ens_ID)
table(ControlCisLncs$source)

fpkm_allG$EnsID_merge <- gsub("\\.[0-9]*", "", fpkm_allG$EnsID)
ControlCisLncs <- merge(ControlCisLncs, fpkm_allG[,c(2,5,61)], by.x = "Ens_ID", by.y = "EnsID_merge", all.x = T)

Enhancer_lociII$FANTID_merge <- gsub("\\.[0-9]*", "", Enhancer_lociII$FANTOM_ID)
ControlCisLncs <- merge(ControlCisLncs, Enhancer_lociII[,c(1,43)], by.x = "Ens_ID", by.y = "FANTID_merge", all.x = T)

ControlCisLncs <- unique(ControlCisLncs)

#EnsID for 33 lncRNAs expressed in SVSMC with exp evidence suggesting cis acting abiity
ControlCisLncs_exprs <- unique(c(ControlCisLncs$EnsID.x, ControlCisLncs$EnsID.y))
ControlCisLncs_exprs <- ControlCisLncs_exprs[!is.na(ControlCisLncs_exprs)]

CoRegPairs_04_48_24_extended$ExpectedCis <- "No"
CoRegPairs_04_48_24_extended$ExpectedCis[CoRegPairs_04_48_24_extended$EnsID %in% ControlCisLncs_exprs] <- "Yes"

table(CoRegPairs_04_48_24_extended$ExpectedCis)


#### save all info to make SCClncRNAs ####

#write.csv(CoRegPairs_04_48_24_extended, "CoRegPairs_04_48_24_extended_SCClnc_Nov25.csv", row.names = F)
CoRegPairs_04_48_24_extended <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CoRegPairs_04_48_24_extended_SCClnc_Nov25.csv")


#### Heatmap of cis-acting evidence ####

#set per gene:
colnames(CoRegPairs_04_48_24_extended)
CoRegPairs_04_48_24_extended_G <- unique(filter(CoRegPairs_04_48_24_extended[,c(6,37,7,36,38,39)]))

eQTL_G <- unique(filter(CoRegPairs_04_48_24_extended_G, eQTLvalidations > 0 | FANTOM_eQTL == "Yes")$lncNameCol)
Corr_G <- unique(filter(CoRegPairs_04_48_24_extended_G, corSig == "Yes")$lncNameCol)
HiC_G <- unique(filter(CoRegPairs_04_48_24_extended_G, !loopMethod == "Neither")$lncNameCol)

CoRegPairs_04_48_24_extended_G$loopMethod_G <- "No"
CoRegPairs_04_48_24_extended_G$loopMethod_G[CoRegPairs_04_48_24_extended_G$lncNameCol %in% HiC_G] <- "Yes"
CoRegPairs_04_48_24_extended_G$corSig_G <- "No"
CoRegPairs_04_48_24_extended_G$corSig_G[CoRegPairs_04_48_24_extended_G$lncNameCol %in% Corr_G] <- "Yes"
CoRegPairs_04_48_24_extended_G$eQTL_G <- "No"
CoRegPairs_04_48_24_extended_G$eQTL_G[CoRegPairs_04_48_24_extended_G$lncNameCol %in% eQTL_G] <- "Yes"
CoRegPairs_04_48_24_extended_G$ExpectedCis_G <- "No"
CoRegPairs_04_48_24_extended_G$ExpectedCis_G[CoRegPairs_04_48_24_extended_G$EnsID %in% ControlCisLncs_exprs] <- "Yes"

CoRegPairs_04_48_24_extended_G <- unique(CoRegPairs_04_48_24_extended_G[,-c(3:6)])

compare_overlap <- data.frame("HiC" = as.numeric(CoRegPairs_04_48_24_extended_G$loopMethod_G == "Yes"), 
                              "Correlated" = as.numeric(CoRegPairs_04_48_24_extended_G$corSig_G == "Yes"), 
                              "eQTL" = as.numeric(CoRegPairs_04_48_24_extended_G$eQTL_G == "Yes"),
                              "Expected cis-acting" = as.numeric(CoRegPairs_04_48_24_extended_G$ExpectedCis_G == "Yes"))

compare_overlap$lncs <- CoRegPairs_04_48_24_extended_G$lncNameCol

compare_overlap <- filter(compare_overlap, rowSums(as.matrix(compare_overlap[,1:( length(colnames(compare_overlap)) -1) ] )) >0)

trial <- as.matrix(compare_overlap[,1:( length(colnames(compare_overlap)) -1)])
rownames(trial) <- compare_overlap$lncs
compare_overlap <- trial

compare_overlap <- compare_overlap[order(compare_overlap[,1], compare_overlap[,2], 
                                         compare_overlap[,3], compare_overlap[,4],
                                         decreasing = T),]

#note that most expected cis-acting lncs are covered by one or more layer of evidence, only one is not 
pheatmap::pheatmap(t(compare_overlap), cluster_cols = F, cluster_rows = F, angle_col = "315", legend = F, fontsize_row = 10)

#repeat for early
compare_overlap_early <- compare_overlap[rownames(compare_overlap) %in% filter(CoRegPairs_04_48_24_extended, Lnc_Timeframe == "<4hrs")$lncNameCol,]
colnames(compare_overlap_early)[4] <- "Cis-acting lncRNA lit."
pheatmap::pheatmap(t(compare_overlap_early), cluster_cols = F, cluster_rows = F, angle_col = "315", legend = F, fontsize_row = 13)

compare_overlap_earlyInduced <- compare_overlap[rownames(compare_overlap) %in% filter(CoRegPairs_04_48_24_extended, Lnc_Cluster == "Induced <4hrs")$lncNameCol,]
pheatmap::pheatmap(t(compare_overlap_earlyInduced), cluster_cols = F, cluster_rows = F, angle_col = "315", legend = F, fontsize_row = 10)

compare_overlap_earlyRepressed <- compare_overlap[rownames(compare_overlap) %in% filter(CoRegPairs_04_48_24_extended, Lnc_Cluster == "Repressed <4hrs")$lncNameCol,]
pheatmap::pheatmap(t(compare_overlap_earlyRepressed), cluster_cols = F, cluster_rows = F, angle_col = "315", legend = F, fontsize_row = 10)


#enrichment in tier 2 with corr
table(CoRegPairs_04_48_24_extended_G$ExpectedCis_G)
12/133

table(filter(CoRegPairs_04_48_24_extended_G, 
             (loopMethod_G == "Yes" |
                corSig_G == "Yes" |
                eQTL_G == "Yes"))$ExpectedCis_G)
11/94

a <- 11
b <- 94
c <- 12
d <- 133

#0.08, somewhat weak
fisher.test(data.frame("cisLnc" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")

#enrichment in tier 2 without corr
table(filter(CoRegPairs_04_48_24_extended_G, 
             (loopMethod_G == "Yes" |
                #corSig_G == "Yes" |
                eQTL_G == "Yes"))$ExpectedCis_G)

a <- 8
b <- 60
c <- 12
d <- 133

fisher.test(data.frame("cisLnc" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")


#cis lnc known: trancistor etc, Agrawal bit less discerning than others which focus on repeatable effects across Gapmer or more in-depth profiling
CisLncStrongest <- filter(ControlCisLncs, !source == "Agrawal")

filter(CoRegPairs_04_48_24_extended_G, EnsID %in% c(CisLncStrongest$EnsID.x, CisLncStrongest$EnsID.y))
7/133
filter(CoRegPairs_04_48_24_extended_G, EnsID %in% c(CisLncStrongest$EnsID.x, CisLncStrongest$EnsID.y), 
       (loopMethod_G == "Yes" |
          corSig_G == "Yes" |
          eQTL_G == "Yes"
       ))
6/94
filter(CoRegPairs_04_48_24_extended_G, lncNameCol %in% CisLncStrongest$EnsName, 
       (loopMethod_G == "Yes" |
          #corSig_G == "Yes" |
          eQTL_G == "Yes"
       ))
4/60

a <- 4
b <- 60
c <- 7
d <- 133

fisher.test(data.frame("cisLnc" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#no enrichment


#### some stats/numbers when comparing the three approaches ####

#total SCClncRNAs:
SCCPairs_04_48_24_extended <- filter(CoRegPairs_04_48_24_extended, (corSig == "Yes" | 
                                                                      !loopMethod == "Neither" | 
                                                                      eQTLvalidations >0 | #ExpectedCis == "Yes" | 
                                                                      FANTOM_eQTL == "Yes"))
length(unique(CoRegPairs_04_48_24_extended$EnsID))#95 (fits with heatmap above)
length(unique(SCCPairs_04_48_24_extended$EnsID))#95 (fits with heatmap above)

SCCPairs_04_48_24_extended_earlyInd <- filter(CoRegPairs_04_48_24_extended, 
                                              (corSig == "Yes" | !loopMethod == "Neither" | eQTLvalidations >0 | ExpectedCis == "Yes" | FANTOM_eQTL == "Yes"),
                                              EnsID %in% fpkm_allGDE_Upwithin_4$EnsID)
length(unique(SCCPairs_04_48_24_extended_earlyInd$EnsID)) #35 fits with heatmap above

SCCPairs_04_48_24_extended_earlyRep <- filter(CoRegPairs_04_48_24_extended, 
                                              (corSig == "Yes" | !loopMethod == "Neither" | eQTLvalidations >0 | ExpectedCis == "Yes" | FANTOM_eQTL == "Yes"),
                                              EnsID %in% fpkm_allGDE_Downwithin_4$EnsID)
length(unique(SCCPairs_04_48_24_extended_earlyRep$EnsID)) #25 fits with heatmap above


#eQTL supported:
unique(filter(CoRegPairs_04_48_24_extended, #eQTLvalidations > 0 | 
              FANTOM_eQTL == "Yes"
                )$lncNameCol)#37 eQTL-supported CClncRNAs, 31 from our approach, 16 from FANTOM, 6 extra from FANTOM
31/133 #23% eQTL supported
37/133 #28% with FANTOM
unique(filter(CoRegPairs_04_48_24_extended, corSig == "Yes")$lncNameCol)#71 correlated CClncRNAs
unique(filter(CoRegPairs_04_48_24_extended, corSig == "Yes", eQTLvalidations > 0 | FANTOM_eQTL == "Yes"
              )$lncNameCol)
71/133 #53% correlated
14/71 #20% correlated are eQTL-supported
#literally no enrichment of correlation amongst eQTL-supported (23% generally, 20% amongst correlated)
16/71 #22.5%, same picture with FANTOM support included

#HiC supported:
unique(filter(CoRegPairs_04_48_24_extended, !loopMethod == "Neither")$lncNameCol)
32/133
unique(filter(CoRegPairs_04_48_24_extended, corSig == "Yes")$lncNameCol)
unique(filter(CoRegPairs_04_48_24_extended, corSig == "Yes", !loopMethod == "Neither")$lncNameCol)
16/71
#literally no enrichment of correlation amongst HiC-supported (24% generally, 23% amongst correlated)

#HiC vs. eQTL:
32/133 #24% cclncs hi-c supported
unique(filter(CoRegPairs_04_48_24_extended, eQTLvalidations > 0 | FANTOM_eQTL == "Yes", !loopMethod == "Neither")$lncNameCol)
7/31 #23% eqtl cclncs hi-c supported
#literally no enrichment of eQTL-support amongst HiC-supported (24% generally, 23% amongst correlated)
8/31 #barely any movemen with FANTOM 2

#no enrichment in either of correlated, suggesting correlation not massively informative?
#or at least that the 3 find diff things

#Venn of all:
#build venn of pairs:
library(ggVennDiagram)
compare_overlapV <- list("HiC" = unique(filter(CoRegPairs_04_48_24_extended, !loopMethod == "Neither")$pairs),
                         "Correlated" = unique(filter(CoRegPairs_04_48_24_extended, corSig == "Yes")$pairs),
                         "eQTL" = unique(filter(CoRegPairs_04_48_24_extended, eQTLvalidations > 0 | FANTOM_eQTL == "Yes")$pairs)
                         #,"eQTL2" = unique(filter(CoRegPairs_04_48_24_extended, FANTOM_eQTL == "Yes")$pairs)
                         )
ggVennDiagram(compare_overlapV, label_color = "white")

#build venn of genes:
compare_overlapV <- list("HiC" = unique(filter(CoRegPairs_04_48_24_extended, !loopMethod == "Neither")$lncNameCol),
                         "Correlated" = unique(filter(CoRegPairs_04_48_24_extended, corSig == "Yes")$lncNameCol),
                         "eQTL" = unique(filter(CoRegPairs_04_48_24_extended, eQTLvalidations > 0 | FANTOM_eQTL == "Yes")$lncNameCol))
ggVennDiagram(compare_overlapV, label_color = "white")


#### GO/KEGG/REACTOME analysis same/later of 4hr lncs ####

#no strong signal could be found amongst all targets of early cclncs, tier 2 targets may refine?

#is there a collective signature amongst tier2 targets of early cclncs, that is not driven by single cluster genes like CXCL, HOX etc

#simple GO check, vs. all DEGs:

#background of vs. all DEGs may give specificity (but also reduce power)
fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE.csv", header = T)
fpkm_allGDE_filt <- filter(fpkm_allGDE, grepl("chr", chr), !grepl("artefacts", GeneClassUpdate), !is.na(GeneClassUpdate))
#remove manual fails
fpkm_allGDE <- filter(fpkm_allGDE_filt, EnsID %in% fpkm_allG_filt_manual$EnsID)
fpkm_PCGDE <- filter(fpkm_allGDE, grepl("protein_coding", EnsType))
fpkm_PCGDE$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE$EnsID)


#compare to #0-24hr targets of 4hr scclncs
#get 4hr DEGs
fpkm_allGDE_Upwithin_4 <- filter(fpkm_allGDE, (Hour0_meanFPKM >=1 | Hour4_meanFPKM >=1) &
                                   (LogFC_0_4 >= log2(1.5) & preadj_0_4 <0.05))
fpkm_allGDE_Downwithin_4 <- filter(fpkm_allGDE, (Hour0_meanFPKM >=1 | Hour4_meanFPKM >=1) &
                                     (LogFC_0_4 < -log2(1.5) & preadj_0_4 <0.05))
fpkm_allGDE_within_4 <- rbind(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Downwithin_4)

#filter the scclncs
AllTargets_T2 <- unique(filter(CoRegPairs_04_48_24_extended, (corSig == "Yes" | 
                                                                !loopMethod == "Neither" | 
                                                                eQTLvalidations >0 | ExpectedCis == "Yes" |FANTOM_eQTL == "Yes"
                                                              ),
                            EnsID %in% fpkm_allGDE_within_4$EnsID))

#all co-regulated genes vs DEGs
CoReg_DE <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", AllTargets_T2$EnsID.y)),
                     universe      = gsub("\\.[0-9]*", "", fpkm_PCGDE$EnsID),
                     keyType       = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "all",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
CoReg_DE_df <- as.data.frame(CoReg_DE)
#terms of interest, but purely driven by cxcl/hox

#select one from each (highest abundance)
AllTargets_loci1 <- filter(AllTargets_T2, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName.y))
AllTargets_loci <- filter(AllTargets_T2, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName.y))
AllTargets_loci_T2 <- rbind(AllTargets_loci, AllTargets_loci1)

fpkm_PCGDE_loci1 <- filter(fpkm_PCGDE, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName))
fpkm_PCGDE_loci <- filter(fpkm_PCGDE, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName))
fpkm_PCGDE_loci <- rbind(fpkm_PCGDE_loci, fpkm_PCGDE_loci1)

CoReg_DE <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", AllTargets_loci_T2$EnsID.y)),
                     universe      = gsub("\\.[0-9]*", "", fpkm_PCGDE_loci$EnsID),
                     keyType       = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "all",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
CoReg_DE_df_T2 <- as.data.frame(CoReg_DE)#still no terms for all, BP or MF

#KEGG
convertEnsEnt <- bitr(unique(fpkm_allG$EnsName), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

CoReg_DE <- enrichKEGG(gene = filter(convertEnsEnt, SYMBOL %in% AllTargets_loci_T2$EnsName.y)$ENTREZID,
                       universe = filter(convertEnsEnt, SYMBOL %in% fpkm_PCGDE_loci$EnsName)$ENTREZID,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
CoReg_DE_df_T2 <- data.frame(CoReg_DE)
dotplot(CoReg_DE, show = 10, font.size = 15)
#5 terms returned
#20% are cytokine-cytokine receptor interaction related, ~15% JAK/STAT and cancer transcriptional misregulation
#cancer misreg link to prolif? HMGA2, IL6, CXCL8:
filter(convertEnsEnt, ENTREZID %in% unlist(strsplit(filter(CoReg_DE_df_T2, Description == "Transcriptional misregulation in cancer")$geneID, split = "/")))

#importantly, the same test done without T2, yielded no terms:
AllTargets <- unique(filter(CoRegPairs_04_48_24_extended,
                            EnsID %in% fpkm_allGDE_within_4$EnsID))
AllTargets$EnsID_merge.y <- gsub("\\.[0-9]*", "", AllTargets$EnsID.y)
AllTargets_loci1 <- filter(AllTargets, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName.y))
AllTargets_loci <- filter(AllTargets, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName.y))
AllTargets_loci <- rbind(AllTargets_loci, AllTargets_loci1)

length(unique(AllTargets_loci$EnsID.y))
length(unique(AllTargets_loci_T2$EnsID.y))#clearly a way smaller group

CoReg_DE <- enrichKEGG(gene = filter(convertEnsEnt, SYMBOL %in% AllTargets_loci$EnsName.y)$ENTREZID,
                       universe = filter(convertEnsEnt, SYMBOL %in% fpkm_PCGDE_loci$EnsName)$ENTREZID,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
CoReg_DE_df <- data.frame(CoReg_DE)
#no terms returned, implying that T2 has narrowed to a group with more of a collective biological theme
#because unrelated genes removed? 
#does the string PPI increase?

#all early lnc targets:
#write.csv(unique(AllTargets_loci$EnsName.y), "all_earlyCCLnc_targets.csv")
#all early lnc targets with extra evidence
#write.csv(unique(AllTargets_loci_T2$EnsName.y), "T2_earlyCCLnc_targets.csv")

#p is smaller, pool size, cannae see enrichment strength tho

#REACTOME
library(ReactomePA)
CoReg_DE <- enrichPathway(gene          = unique(filter(convertEnsEnt,SYMBOL %in% AllTargets_loci_T2$EnsName.y)$ENTREZID),
                          universe      = unique(filter(convertEnsEnt,SYMBOL %in% fpkm_PCGDE_loci$EnsName)$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)
CoReg_DE_df_T2 <- as.data.frame(CoReg_DE)
#0x terms

CoReg_DE <- enrichPathway(gene          = unique(filter(convertEnsEnt,SYMBOL %in% AllTargets_loci$EnsName.y)$ENTREZID),
                          universe      = unique(filter(convertEnsEnt,SYMBOL %in% fpkm_PCGDE_loci$EnsName)$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)
CoReg_DE_df <- as.data.frame(CoReg_DE)
#again, 0x terms

#going to broader background may get more of a theme, enriched themes amongst CClnc targets vs. all ex genes
fpkm_PCG <- filter(fpkm_allG, grepl("protein_coding", EnsType))
fpkm_PCG$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCG$EnsID)
fpkm_PCG_loci1 <- filter(fpkm_PCG, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName))
fpkm_PCG_loci <- filter(fpkm_PCG, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName))
fpkm_PCG_loci <- rbind(fpkm_PCG_loci, fpkm_PCG_loci1)

#GO on broad background
CoReg_exp <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", AllTargets_loci_T2$EnsID.y)),
                     universe      = gsub("\\.[0-9]*", "", fpkm_PCG_loci$EnsID),
                     keyType       = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "all",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
CoReg_exp_df_T2 <- as.data.frame(CoReg_exp)#1x term, GO still doesn't really work

#KEGG
CoReg_exp <- enrichKEGG(gene = filter(convertEnsEnt, SYMBOL %in% AllTargets_loci_T2$EnsName.y)$ENTREZID,
                       universe = filter(convertEnsEnt, SYMBOL %in% fpkm_PCG_loci$EnsName)$ENTREZID,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
CoReg_exp_df_T2 <- data.frame(CoReg_exp)
#8x terms few more than before
dotplot(CoReg_exp, show = 6, font.size = 15)

#check the non T2 version again
CoReg_exp <- enrichKEGG(gene = filter(convertEnsEnt, SYMBOL %in% AllTargets_loci$EnsName.y)$ENTREZID,
                        universe = filter(convertEnsEnt, SYMBOL %in% fpkm_PCG_loci$EnsName)$ENTREZID,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
CoReg_exp_df <- data.frame(CoReg_exp)
#5 terms much less significance

#REACTOME
CoReg_exp <- enrichPathway(gene          = unique(filter(convertEnsEnt,SYMBOL %in% AllTargets_loci_T2$EnsName.y)$ENTREZID),
                          universe      = unique(filter(convertEnsEnt,SYMBOL %in% fpkm_PCG_loci$EnsName)$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)
CoReg_exp_df_T2 <- as.data.frame(CoReg_exp)
dotplot(CoReg_exp, show = 6)

#REACTOME + KEGG in broad agreement that cytokine response is key feature found
#KEGG gives wider flavour of possible roles (transcritpional misregulation is interesting too)

#proper plot, for KEGG:
#CoReg_exp_df_T2 <- CoReg_exp_df_T2[order(CoReg_exp_df_T2$Description),]

#CoReg_exp_df_T2$DescriptionII <- stringr::str_wrap(CoReg_exp_df_T2$Description, width = 30)

#CoReg_exp_df_T2$Description <- factor(CoReg_exp_df_T2$Description, labels = CoReg_exp_df_T2$DescriptionII)
#CoReg_exp_df_T2$Description <- factor(CoReg_exp_df_T2$Description,
#                                           levels = levels(CoReg_exp_df_T2$Description)[order(CoReg_exp_df_T2$p.adjust, decreasing = T)])

#ggplot(CoReg_exp_df_T2) + aes(y = Description, x = -log10(p.adjust)) +
#  geom_bar(stat = "identity", fill= "steelblue", color = "grey60") +
#  theme_minimal() +
#  xlab("REACTOME Term\nSignificance (-log10p)") +
#  ylab("") +
#  theme(text = element_text(size=24))


#### GO/KEGG/REACTOME analysis same/later of 4hr lncs (considering merges) ####

#no strong signal could be found amongst all targets of early cclncs, tier 2 targets may refine?

#is there a collective signature amongst tier2 targets of early cclncs, that is not driven by single cluster genes like CXCL, HOX etc

#simple GO check, vs. all DEGs:

#background of vs. all DEGs may give specificity (but also reduce power)
fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE.csv", header = T)
fpkm_allGDE_filt <- filter(fpkm_allGDE, grepl("chr", chr), !grepl("artefacts", GeneClassUpdate), !is.na(GeneClassUpdate))
#remove manual fails
fpkm_allGDE <- filter(fpkm_allGDE_filt, EnsID %in% fpkm_allG_filt_manual$EnsID)
fpkm_PCGDE <- filter(fpkm_allGDE, grepl("protein_coding", EnsType))
fpkm_PCGDE$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE$EnsID)

#for merged genes, select only the ones that are DEGs in GENv26
#expand table for merged genes:
fpkm_allG_ID <- unique(fpkm_allG[,c(1,2,3,5,27,29,31,33,35)])
fpkm_allG_ID <- filter(fpkm_allG_ID, !is.na(AllEns))

fpkm_allG_ID_long <- fpkm_allG_ID %>%
  mutate(EnsID2 = AllEns) %>%
  tidyr::separate_rows(EnsID2, sep = ", ")
fpkm_allG_ID_long$EnsID2[is.na(fpkm_allG_ID_long$EnsID2)] <- fpkm_allG_ID_long$EnsID[is.na(fpkm_allG_ID_long$EnsID2)]

fpkm_PCGDE_long <- filter(fpkm_allG_ID_long, EnsID %in% fpkm_PCGDE$EnsID)

#import GENv26 analysis:
GENv26_DEnonDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Timecourse_GENv26_DEnonDE.csv", header = T)

GENv26_DE <- filter(GENv26_DEnonDE, padj_AllTimepoints.LRT. < 0.05, fpkm_max_treatment >1, 
                    (preadj_0_4 <0.05 & (LogFC_0_4 > log2(1.5) | LogFC_0_4 < -log2(1.5))) |
                      (preadj_0_8 <0.05 & (LogFC_0_8 > log2(1.5) | LogFC_0_8 < -log2(1.5))) |
                      (preadj_0_24 <0.05 & (LogFC_0_24 > log2(1.5) | LogFC_0_24 < -log2(1.5))) |
                      (preadj_4_8 <0.05 & (LogFC_4_8 > log2(1.5) | LogFC_4_8 < -log2(1.5))) |
                      (preadj_4_24 <0.05 & (LogFC_4_24 > log2(1.5) | LogFC_4_24 < -log2(1.5))) |
                      (preadj_8_24 <0.05 & (LogFC_8_24 > log2(1.5) | LogFC_8_24 < -log2(1.5)))
                    )

#filter out PCGs that are not DE in this analysis:
fpkm_PCGDE_long_ <- filter(fpkm_PCGDE_long, !EnsID2 %in% GENv26_DE$V4)
fpkm_PCGDE_long <- filter(fpkm_PCGDE_long, EnsID2 %in% GENv26_DE$V4)

#background is a) non-merged DE PCGs:
fpkm_PCGDE_background <- unique(c(filter(fpkm_PCGDE, !EnsID %in% fpkm_allG_ID$EnsID)$EnsID,
                                  #and b) any individual DEGs from merges:
                                  fpkm_PCGDE_long$EnsID2))
#added on about 100 DE PCGs
100/4071 #2.7% so less than anticipated


#compare to #0-24hr targets of 4hr scclncs
#get 4hr DEGs
fpkm_allGDE_Upwithin_4 <- filter(fpkm_allGDE, (Hour0_meanFPKM >=1 | Hour4_meanFPKM >=1) &
                                   (LogFC_0_4 >= log2(1.5) & preadj_0_4 <0.05))
fpkm_allGDE_Downwithin_4 <- filter(fpkm_allGDE, (Hour0_meanFPKM >=1 | Hour4_meanFPKM >=1) &
                                     (LogFC_0_4 < -log2(1.5) & preadj_0_4 <0.05))
fpkm_allGDE_within_4 <- rbind(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Downwithin_4)

#filter the scclncs
AllTargets_T2 <- unique(filter(CoRegPairs_04_48_24_extended, (corSig == "Yes" | 
                                                                !loopMethod == "Neither" | 
                                                                eQTLvalidations >0 | ExpectedCis == "Yes" |FANTOM_eQTL == "Yes"
                                                              ),EnsID %in% fpkm_allGDE_within_4$EnsID))
AllTargets_T2_EnsID.y <- unique(AllTargets_T2$EnsID.y)

#16 are merged DE PCGs
mergedPCG <- AllTargets_T2_EnsID.y[AllTargets_T2_EnsID.y %in% fpkm_PCGDE_long$EnsID]
#29 DEGs hidden in the merged 16:
unique(filter(fpkm_PCGDE_long, EnsID %in% mergedPCG)$EnsID)
unique(filter(fpkm_PCGDE_long, EnsID %in% mergedPCG)$EnsID2)

#from 127 to 143
AllTargets_T2_EnsID.y <- unique(c(AllTargets_T2_EnsID.y, 
                                  filter(fpkm_PCGDE_long, EnsID %in% mergedPCG)$EnsID2)
                                )
127/143 #9% more, quite a big increase

#all co-regulated genes vs DEGs
CoReg_DE <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", AllTargets_T2_EnsID.y)),
                     universe      = gsub("\\.[0-9]*", "", fpkm_PCGDE_background),
                     keyType       = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "all",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
CoReg_DE_df <- as.data.frame(CoReg_DE)
#terms of interest, but purely driven by cxcl/hox

#select one from each of CXCL, MT1, HOXA and HOXC loci
toAdd <- filter(GENv26_DE, grepl("CXCL8|^MT1E|^HOXA10$|^HOXC10", V7))$V4
toRemove <- filter(GENv26_DE, grepl("^CXCL|^MT1|^HOXA|^HOXC", V7))$V4

AllTargets_T2_EnsID.y_loci <- AllTargets_T2_EnsID.y[!AllTargets_T2_EnsID.y %in% toRemove]
AllTargets_T2_EnsID.y_loci <- c(AllTargets_T2_EnsID.y_loci, toAdd)

fpkm_PCGDE_background_loci <- fpkm_PCGDE_background[!fpkm_PCGDE_background %in% toRemove]
fpkm_PCGDE_background_loci <- c(fpkm_PCGDE_background_loci, toAdd)


CoReg_DE <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", AllTargets_T2_EnsID.y_loci)),
                     universe      = gsub("\\.[0-9]*", "", fpkm_PCGDE_background_loci),
                     keyType       = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "all",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
CoReg_DE_df_T2 <- as.data.frame(CoReg_DE)#still no terms for all, BP or MF

#KEGG
convertEnsEnt <- bitr(gsub("\\.[0-9]*", "", fpkm_PCGDE_background_loci), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

CoReg_DE <- enrichKEGG(gene     = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", AllTargets_T2_EnsID.y_loci))$ENTREZID),
                       universe = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", fpkm_PCGDE_background_loci))$ENTREZID),
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)
CoReg_DE_df_T2 <- data.frame(CoReg_DE)
dotplot(CoReg_DE, show = 10, font.size = 15)
#2 terms returned
#much worse than previous analysis without considering the merge issue


#REACTOME
library(ReactomePA)
CoReg_DE <- enrichPathway(gene          = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", AllTargets_T2_EnsID.y_loci))$ENTREZID),
                          universe      = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", fpkm_PCGDE_background_loci))$ENTREZID),
                          organism = "human",
                          pvalueCutoff = 0.05,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)
CoReg_DE_df_T2 <- as.data.frame(CoReg_DE)
#0x terms


#going to broader background may get more of a theme, and be less complex
#enriched themes amongst CClnc targets vs. all ex genes
fpkm_PCG <- filter(fpkm_allG, grepl("protein_coding", EnsType))

fpkm_PCG_long <- filter(fpkm_allG_ID_long, EnsID %in% fpkm_PCG$EnsID)

#GENv26 expressed genes
GENv26_xp <- filter(GENv26_DEnonDE, fpkm_max_treatment >1)

#filter out PCGs that are not DE in this analysis:
fpkm_PCG_long_ <- filter(fpkm_PCG_long, !EnsID2 %in% GENv26_xp$V4)
fpkm_PCG_long <- filter(fpkm_PCG_long, EnsID2 %in% GENv26_xp$V4)

#background is non-merged DE PCGs:
fpkm_PCG_long_background <- unique(c(filter(fpkm_PCG, !EnsID %in% fpkm_allG_ID$EnsID)$EnsID,
                                  #plus any individual DEGs from merges:
                                  fpkm_PCG_long$EnsID2))

#select one from each of CXCL, MT1, HOXA and HOXC loci
toAdd_xp <- filter(GENv26_xp, grepl("CXCL8|^MT1E|^HOXA10$|^HOXC10", V7))$V4
toRemove_xp <- filter(GENv26_xp, grepl("^CXCL|^MT1|^HOXA|^HOXC", V7))$V4

fpkm_PCG_long_background_loci <- fpkm_PCG_long_background[!fpkm_PCG_long_background %in% toRemove_xp]
fpkm_PCG_long_background_loci <- c(fpkm_PCG_long_background_loci, toAdd_xp)

#GO on broad background
CoReg_exp <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", AllTargets_T2_EnsID.y_loci)),
                      universe      = gsub("\\.[0-9]*", "", fpkm_PCG_long_background_loci),
                      keyType       = "ENSEMBL",
                      OrgDb         = org.Hs.eg.db,
                      ont           = "all",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
CoReg_exp_df_T2 <- as.data.frame(CoReg_exp)#0x terms, GO still doesn't really work

#KEGG
CoReg_exp <- enrichKEGG(gene          = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", AllTargets_T2_EnsID.y_loci))$ENTREZID),
                        universe      = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", fpkm_PCG_long_background_loci))$ENTREZID),
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
CoReg_exp_df_T2 <- data.frame(CoReg_exp)
#still just 2x terms
dotplot(CoReg_exp, show = 6, font.size = 15)

#REACTOME
CoReg_exp <- enrichPathway(gene          = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", AllTargets_T2_EnsID.y_loci))$ENTREZID),
                           universe      = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", fpkm_PCG_long_background_loci))$ENTREZID),
                           organism = "human",
                           pvalueCutoff = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = TRUE)
CoReg_exp_df_T2 <- as.data.frame(CoReg_exp)
dotplot(CoReg_exp, show = 6)

#unmerging the genes has killed the analysis


#### GO/KEGG/REACTOME analysis same of 4hr lncs (considering merges) ####

#no strong signal could be found amongst all targets of early cclncs, tier 2 targets may refine?

#is there a collective signature amongst tier2 targets of early cclncs, that is not driven by single cluster genes like CXCL, HOX etc

#simple GO check, vs. all DEGs:

#look at 4hr co-regulations, so just filter the above target list by if induced in 4hrs:
GENv26_DE_4up <- filter(GENv26_DEnonDE, (Hour0_meanFPKM>1 | Hour4_meanFPKM>1), 
                    (preadj_0_4 <0.05 & LogFC_0_4 > log2(1.5)))
GENv26_DE_4down <- filter(GENv26_DEnonDE, (Hour0_meanFPKM>1 | Hour4_meanFPKM>1), 
                        (preadj_0_4 <0.05 & LogFC_0_4 < -log2(1.5)))

EarlyUpTargets_T2_EnsID.y_loci <- AllTargets_T2_EnsID.y_loci[AllTargets_T2_EnsID.y_loci %in% GENv26_DE_4up$V4]

#GO on broad background
CoReg_exp <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", EarlyUpTargets_T2_EnsID.y_loci)),
                      universe      = gsub("\\.[0-9]*", "", fpkm_PCG_long_background_loci),
                      keyType       = "ENSEMBL",
                      OrgDb         = org.Hs.eg.db,
                      ont           = "all",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
CoReg_exp_df_T2 <- as.data.frame(CoReg_exp)#3x terms, immune response

#KEGG
CoReg_exp <- enrichKEGG(gene          = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", EarlyUpTargets_T2_EnsID.y_loci))$ENTREZID),
                        universe      = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", fpkm_PCG_long_background_loci))$ENTREZID),
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
CoReg_exp_df_T2 <- data.frame(CoReg_exp)
#17x terms
dotplot(CoReg_exp, show = 6, font.size = 15)

#REACTOME
CoReg_exp <- enrichPathway(gene          = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", EarlyUpTargets_T2_EnsID.y_loci))$ENTREZID),
                           universe      = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", fpkm_PCG_long_background_loci))$ENTREZID),
                           organism = "human",
                           pvalueCutoff = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = TRUE)
CoReg_exp_df_T2 <- as.data.frame(CoReg_exp)
dotplot(CoReg_exp, show = 6)
#2x terms


#4hr repressed
GENv26_DE_4down <- filter(GENv26_DEnonDE, (Hour0_meanFPKM>1 | Hour4_meanFPKM>1), 
                          (preadj_0_4 <0.05 & LogFC_0_4 < -log2(1.5)))

EarlyDownTargets_T2_EnsID.y_loci <- AllTargets_T2_EnsID.y_loci[AllTargets_T2_EnsID.y_loci %in% GENv26_DE_4down$V4]

#GO on broad background
CoReg_exp <- enrichGO(gene         = unique(gsub("\\.[0-9]*", "", EarlyDownTargets_T2_EnsID.y_loci)),
                      universe      = gsub("\\.[0-9]*", "", fpkm_PCG_long_background_loci),
                      keyType       = "ENSEMBL",
                      OrgDb         = org.Hs.eg.db,
                      ont           = "all",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
CoReg_exp_df_T2 <- as.data.frame(CoReg_exp)#1x term

#KEGG
CoReg_exp <- enrichKEGG(gene          = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", EarlyDownTargets_T2_EnsID.y_loci))$ENTREZID),
                        universe      = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", fpkm_PCG_long_background_loci))$ENTREZID),
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
CoReg_exp_df_T2 <- data.frame(CoReg_exp)
#1 term
dotplot(CoReg_exp, show = 6, font.size = 15)

#REACTOME
CoReg_exp <- enrichPathway(gene          = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", EarlyDownTargets_T2_EnsID.y_loci))$ENTREZID),
                           universe      = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", fpkm_PCG_long_background_loci))$ENTREZID),
                           organism = "human",
                           pvalueCutoff = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = TRUE)
CoReg_exp_df_T2 <- as.data.frame(CoReg_exp)
dotplot(CoReg_exp, show = 6)
#0x terms

#4hr induced genes only, strong immune theme from KEGG
#is it still true that this beats the regular cclncRNA-targets for finding a theme?

#regular cclnc targets
AllTargets <- unique(filter(CoRegPairs_04_48_24_extended,EnsID %in% fpkm_allGDE_within_4$EnsID))

AllTargets_EnsID.y <- unique(AllTargets$EnsID.y)

#29 are merged DE PCGs
mergedPCG <- AllTargets_EnsID.y[AllTargets_EnsID.y %in% fpkm_PCGDE_long$EnsID]
#47 DEGs hidden in the merged 29:
unique(filter(fpkm_PCGDE_long, EnsID %in% mergedPCG)$EnsID)
unique(filter(fpkm_PCGDE_long, EnsID %in% mergedPCG)$EnsID2)

#from 245 to 268
AllTargets_EnsID.y <- unique(c(AllTargets_EnsID.y, 
                                  filter(fpkm_PCGDE_long, EnsID %in% mergedPCG)$EnsID2))
245/268 #9% more, quite a big increase

#sort the loci:
#select one from each of CXCL, MT1, HOXA and HOXC loci
toAdd <- filter(GENv26_DE, grepl("CXCL8|^MT1E|^HOXA10$|^HOXC10", V7))$V4
toRemove <- filter(GENv26_DE, grepl("^CXCL|^MT1|^HOXA|^HOXC", V7))$V4

AllTargets_EnsID.y <- AllTargets_EnsID.y[!AllTargets_EnsID.y %in% toRemove]
AllTargets_EnsID.y <- c(AllTargets_EnsID.y, toAdd)

#split to up/down in 4hrs
EarlyUpTargets_EnsID.y_loci <- AllTargets_EnsID.y[AllTargets_EnsID.y %in% GENv26_DE_4up$V4]

#GO on broad background
CoReg_exp <- enrichGO(gene          = unique(gsub("\\.[0-9]*", "", EarlyUpTargets_EnsID.y_loci)),
                      universe      = gsub("\\.[0-9]*", "", fpkm_PCG_long_background_loci),
                      keyType       = "ENSEMBL",
                      OrgDb         = org.Hs.eg.db,
                      ont           = "all",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
CoReg_exp_df <- as.data.frame(CoReg_exp)#1x term only (3 before)

#KEGG
CoReg_exp <- enrichKEGG(gene          = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", EarlyUpTargets_EnsID.y_loci))$ENTREZID),
                        universe      = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", fpkm_PCG_long_background_loci))$ENTREZID),
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)
CoReg_exp_df <- data.frame(CoReg_exp)
#11x terms, 17x before, somewhat weaker
dotplot(CoReg_exp, show = 6, font.size = 15)

#REACTOME
CoReg_exp <- enrichPathway(gene          = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", EarlyUpTargets_EnsID.y_loci))$ENTREZID),
                           universe      = unique(filter(convertEnsEnt, ENSEMBL %in% gsub("\\.[0-9]*", "", fpkm_PCG_long_background_loci))$ENTREZID),
                           organism = "human",
                           pvalueCutoff = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = TRUE)
CoReg_exp_df <- as.data.frame(CoReg_exp)
dotplot(CoReg_exp, show = 6)
#1x terms

#bit weaker but not massively!


#### TF enrichment ####

#LncRNAs reported to favour localisation near TFs genomically

#Can we see an enrichment of TFs amongst T2 0-4hr CClnc targets vs. other DEGs?

CoRegPairs_04_48_24_extended_T2 <- filter(CoRegPairs_04_48_24_extended, (corSig == "Yes" | !loopMethod == "Neither" | eQTLvalidations >0))

length(unique(filter(CoRegPairs_04_48_24_extended, Lnc_Timeframe == "<4hrs")$EnsName.y)) #108 T2 targets of 0-4hr cclncs
length(unique(filter(CoRegPairs_04_48_24_extended, Lnc_Timeframe == "<4hrs", 
                     grepl("TF", GeneClassUpdate.y))$EnsName.y)) #19 of which are TF CClnc targets
34/245

length(unique(filter(CoRegPairs_04_48_24_extended_T2, Lnc_Timeframe == "<4hrs")$EnsName.y)) #108 T2 targets of 0-4hr cclncs
length(unique(filter(CoRegPairs_04_48_24_extended_T2, Lnc_Timeframe == "<4hrs", 
                     grepl("TF", GeneClassUpdate.y))$EnsName.y)) #19 of which are TF CClnc targets
19/108

#amongst general DEG population:
length(unique(fpkm_PCGDE$EnsID)) #4071 0-4hr DEGs
length(unique(filter(fpkm_PCGDE, grepl("TF", GeneClassUpdate))$EnsID)) #338 TF CClnc targets
388/4071 #9.5%

#possibility for an enrichment:
a <- 19
b <- 108
c <- 388
d <- 4071

fisher.test(data.frame("cisLncTarget" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#targets of 0-4hr cclncs are enriched with TFs (probs due to Hox locus)

#remove 3 of the HOXA genes from consideration and re-test
a <- 19-3
b <- 108-3
c <- 388-3
d <- 4071-3

fisher.test(data.frame("cisLncTarget" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#still a decent enrichment
a/b
c/d
#n.b. previously much weaker, T2 highlights the TFs
(34-3)/(245-3)

a <- 34-3
b <- 245-3
c <- 388-3
d <- 4071-3

fisher.test(data.frame("cisLncTarget" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")

#### Blunter checks on TFs, IEGs, Cell cycle genes, SMC-biased genes ####

#GO/KEGG/REACTOME reveals little, check more specific lists

#TFs and CCs already in list

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

#all gene lists
trial <- split(fpkm_allG, fpkm_allG$AllEns)
triali <- lapply(trial, function(x){
  sum(gsub("\\.[0-9]", "", unlist(strsplit(x$AllNames, ", "))) %in% IEGs_hs[,4]) > 0
})
extra_IEG <- unique(filter(fpkm_allG, AllEns %in% names(triali)[triali == T])$EnsID)

trial <- split(fpkm_allG, fpkm_allG$AllEns)
triali <- lapply(trial, function(x){
  sum(gsub("\\.[0-9]", "", unlist(strsplit(x$AllEns, ", "))) %in% FANT_S10_SMC_G) > 0
})
extra_SMC <- unique(filter(fpkm_allG, AllEns %in% names(triali)[triali == T])$EnsID)

trial <- split(fpkm_allG, fpkm_allG$AllEns)
triali <- lapply(trial, function(x){
  sum(gsub("\\.[0-9]", "", unlist(strsplit(x$AllEns, ", "))) %in% FANT_S10_VSMC_G) > 0
})
extra_VSMC <- unique(filter(fpkm_allG, AllEns %in% names(triali)[triali == T])$EnsID)

trial <- split(fpkm_allG, fpkm_allG$AllEns)
triali <- lapply(trial, function(x){
  sum(gsub("\\.[0-9]", "", unlist(strsplit(x$AllNames, ", "))) %in% CRdb_dataB$CR) > 0
})
extra_CR <- unique(filter(fpkm_allG, AllEns %in% names(triali)[triali == T])$EnsID)
#add HMGA2 to this list too:
extra_CR <- c(extra_CR, "ENSG00000149948.13")

GeneLists <- list("TFs"= unique(filter(fpkm_allG, grepl("TF", GeneClassUpdate))$EnsID),
                  "CRs" = unique(filter(fpkm_allG, EnsName %in% CRdb_dataB$CR | EnsID %in% extra_CR)$EnsID),
                  "CC"= unique(filter(fpkm_allG, grepl("CC", GeneClassUpdate))$EnsID),
                  "IEGs"= unique(filter(fpkm_allG, EnsName %in% c(IEGs_hs[,4]) | EnsID %in% extra_IEG)$EnsID),
                  "SMC"= unique(filter(fpkm_allG, gsub("\\.[0-9]*", "", EnsID) %in% FANT_S10_SMC_G | EnsID %in% extra_SMC)$EnsID),
                  "VSMC"= unique(filter(fpkm_allG, gsub("\\.[0-9]*", "", EnsID) %in% FANT_S10_VSMC_G | EnsID %in% extra_VSMC)$EnsID)
)

results_list <- list()

#enrichment of TFs in lnc-targeted DE PCGs rather than PCGs
#same/later timeframe:
AllTargets_T2 <- unique(filter(CoRegPairs_04_48_24_extended, (corSig == "Yes" | !loopMethod == "Neither" | 
                                                                eQTLvalidations >0 | FANTOM_eQTL == "Yes" |
                                                              ExpectedCis == "Yes"),
                               EnsID %in% fpkm_allGDE_within_4$EnsID))
AllTargets_T2$EnsID_merge.y <- gsub("\\.[0-9]*", "", AllTargets_T2$EnsID.y)

#cut out HOXA/CXCL etc which massively bias
AllTargets_loci1 <- filter(AllTargets_T2, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName.y))
AllTargets_loci <- filter(AllTargets_T2, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName.y))
AllTargets_loci_T2 <- rbind(AllTargets_loci, AllTargets_loci1)

fpkm_PCGDE <- filter(fpkm_allGDE, grepl("protein_coding", EnsType))
fpkm_PCGDE$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE$EnsID)
fpkm_PCGDE_loci1 <- filter(fpkm_PCGDE, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName))
fpkm_PCGDE_loci <- filter(fpkm_PCGDE, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName))
fpkm_PCGDE_loci <- rbind(fpkm_PCGDE_loci, fpkm_PCGDE_loci1)

for (i in 1:length(GeneLists)){
  
  a <- length(unique(AllTargets_loci_T2$EnsID.y[ AllTargets_loci_T2$EnsID.y %in% GeneLists[[i]] ]))
  b <- length(unique(AllTargets_loci_T2$EnsID.y))
  c <- length(unique(fpkm_PCGDE_loci$EnsID[ fpkm_PCGDE_loci$EnsID %in% GeneLists[[i]]]))
  d <- length(unique(fpkm_PCGDE_loci$EnsID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate,
                                  paste(unique(AllTargets_loci_T2$EnsName.y[ AllTargets_loci_T2$EnsID.y %in% GeneLists[[i]] ]), collapse = "/"))
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:7] <- c("pval", "OR")
#cut non-vasc SMC
GeneLists_enriched <- GeneLists_enriched[-5,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")
#seems like essentially all these groups (aside from cell cycle) are enriched to quite strong degree 
#(low p but pool is low + strong OR)
GeneLists_enrichedSameLaterT_T2 <- GeneLists_enriched


#with non-T2 the result is much weaker right?
AllTargets <- unique(filter(CoRegPairs_04_48_24_extended,
                               EnsID %in% fpkm_allGDE_within_4$EnsID))
AllTargets$EnsID_merge.y <- gsub("\\.[0-9]*", "", AllTargets$EnsID.y)
#cut out HOXA/CXCL etc which massively bias
AllTargets_loci1 <- filter(AllTargets, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName.y))
AllTargets_loci <- filter(AllTargets, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName.y))
AllTargets_loci <- rbind(AllTargets_loci, AllTargets_loci1)

for (i in 1:length(GeneLists)){
  
  a <- length(unique(AllTargets_loci$EnsID.y[ AllTargets_loci$EnsID.y %in% GeneLists[[i]] ]))
  b <- length(unique(AllTargets_loci$EnsID.y))
  c <- length(unique(fpkm_PCGDE_loci$EnsID[ fpkm_PCGDE_loci$EnsID %in% GeneLists[[i]]]))
  d <- length(unique(fpkm_PCGDE_loci$EnsID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate,
                                  paste(unique(AllTargets_loci$EnsName.y[ AllTargets_loci$EnsID.y %in% GeneLists[[i]] ]), collapse = "/"))
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:7] <- c("pval", "OR")
#cut non-vasc SMC
GeneLists_enriched <- GeneLists_enriched[-5,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")
#weaker results without T2
GeneLists_enrichedSameLaterT <- GeneLists_enriched

#note FOXL1 is called an IEG (fusion to FOXC2 so not quite but kindof)


#same timeframe result:
EarlyTargets_T2 <- unique(filter(CoRegPairs_04_48_24_extended, (corSig == "Yes" | !loopMethod == "Neither" | 
                                                                  eQTLvalidations >0 | FANTOM_eQTL == "Yes" |
                                                                  ExpectedCis == "Yes"),
                               EnsID %in% fpkm_allGDE_within_4$EnsID, EnsID.y %in% fpkm_allGDE_within_4$EnsID))
EarlyTargets_T2$EnsID_merge.y <- gsub("\\.[0-9]*", "", EarlyTargets_T2$EnsID.y)

#cut out HOXA/CXCL etc which massively bias
EarlyTargets_loci1 <- filter(EarlyTargets_T2, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName.y))
EarlyTargets_loci <- filter(EarlyTargets_T2, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName.y))
EarlyTargets_loci_T2 <- rbind(EarlyTargets_loci, EarlyTargets_loci1)

fpkm_PCGDE <- filter(fpkm_allGDE_within_4, grepl("protein_coding", EnsType))
fpkm_PCGDE$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE$EnsID)
fpkm_PCGDE_loci1 <- filter(fpkm_PCGDE, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName))
fpkm_PCGDE_loci <- filter(fpkm_PCGDE, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName))
fpkm_PCGDE_loci <- rbind(fpkm_PCGDE_loci, fpkm_PCGDE_loci1)

for (i in 1:length(GeneLists)){
  
  a <- length(unique(EarlyTargets_loci_T2$EnsID.y[ EarlyTargets_loci_T2$EnsID.y %in% GeneLists[[i]] ]))
  b <- length(unique(EarlyTargets_loci_T2$EnsID.y))
  c <- length(unique(fpkm_PCGDE_loci$EnsID[ fpkm_PCGDE_loci$EnsID %in% GeneLists[[i]]]))
  d <- length(unique(fpkm_PCGDE_loci$EnsID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate,
                                  paste(unique(EarlyTargets_loci_T2$EnsName.y[ EarlyTargets_loci_T2$EnsID.y %in% GeneLists[[i]] ]), collapse = "/"))
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:7] <- c("pval", "OR")
#cut non-vasc SMC
GeneLists_enriched <- GeneLists_enriched[-5,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")
#no significance!
GeneLists_enrichedSameT_T2 <- GeneLists_enriched


#plot enrichment of interesting subtypes in same/later tier2:
GeneLists_enrichedSameLaterT_T2$selectHit <- GeneLists_enrichedSameLaterT_T2$a/GeneLists_enrichedSameLaterT_T2$b
GeneLists_enrichedSameLaterT_T2$backHit <- GeneLists_enrichedSameLaterT_T2$c/GeneLists_enrichedSameLaterT_T2$d

library(reshape2)
plot_Gclass <- melt(GeneLists_enrichedSameLaterT_T2[,c(1,10,11)])

plot_Gclass$GeneList <- c("TFs", "Chromatin Regulators", "Core S/G2M", "IEGs", "VSMC-enriched", 
                          "TFs", "Chromatin Regulators", "Core S/G2M", "IEGs", "VSMC-enriched")
plot_Gclass$GeneList <- as.factor(plot_Gclass$GeneList)
plot_Gclass$GeneList <- factor(plot_Gclass$GeneList, levels = levels(plot_Gclass$GeneList)[c(2,5,3,1,4)])

plot_Gclass$variable <- c(rep("SCClncRNA Targets", 5), rep("All DE PCGs", 5))
#plot_Gclass$variable <- as.factor(plot_Gclass$variable)
#plot_Gclass$variable <- factor(plot_Gclass$variable, levels = levels(plot_Gclass$variable)[c(2,1)])

plot_Gclass$value <- plot_Gclass$value*100

ggplot(plot_Gclass) + aes(y = GeneList, x = value, fill = variable) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c(`SCClncRNA Targets` = "mediumorchid",`All DE PCGs` = "grey60")) +
  theme_minimal() +
  ylab("") +
  xlab("%") +
  theme(text = element_text(size=24))

#is it fair? biased to early rather than all...
dim(filter(fpkm_PCGDE_loci, EnsID %in% fpkm_allGDE_within_4$EnsID))
1579/4059 #39% of PCGs early-reg
table(unique(AllTargets_loci_T2[,c(10,13)])$PCG_Cluster)
sum(table(unique(AllTargets_loci_T2[,c(10,13)])$PCG_Cluster))
59/99 #60% of lnc targets early-reg

#it is part of the explanation for sure
#the ORs for same-timeframe are still quite high, just low pool numbers
#either way, it's still showing and contextualising their potential targets

#### Blunter checks on TFs, IEGs, Cell cycle genes, SMC-biased genes (unmerged version) ####

#GO/KEGG/REACTOME reveals little, check more specific lists

#TFs and CCs already in list

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

#all gene lists
fpkm_PCG_long_background_Symbol <- unique(filter(GENv26_DEnonDE, V4 %in% fpkm_PCG_long_background))
IEGs <- fpkm_PCG_long_background[gsub("\\.[0-9]", "", fpkm_PCG_long_background) %in% IEGs_hs[,4]]

triali <- lapply(fpkm_PCG_long_background, function(x){
  sum(gsub("\\.[0-9]", "", fpkm_PCG_long_background) %in% IEGs_hs[,4]) > 0
})
extra_IEG <- unique(filter(fpkm_allG, AllEns %in% names(triali)[triali == T])$EnsID)

trial <- split(fpkm_allG, fpkm_allG$AllEns)
triali <- lapply(trial, function(x){
  sum(gsub("\\.[0-9]", "", unlist(strsplit(x$AllEns, ", "))) %in% FANT_S10_SMC_G) > 0
})
extra_SMC <- unique(filter(fpkm_allG, AllEns %in% names(triali)[triali == T])$EnsID)

trial <- split(fpkm_allG, fpkm_allG$AllEns)
triali <- lapply(trial, function(x){
  sum(gsub("\\.[0-9]", "", unlist(strsplit(x$AllEns, ", "))) %in% FANT_S10_VSMC_G) > 0
})
extra_VSMC <- unique(filter(fpkm_allG, AllEns %in% names(triali)[triali == T])$EnsID)

trial <- split(fpkm_allG, fpkm_allG$AllEns)
triali <- lapply(trial, function(x){
  sum(gsub("\\.[0-9]", "", unlist(strsplit(x$AllNames, ", "))) %in% CRdb_dataB$CR) > 0
})
extra_CR <- unique(filter(fpkm_allG, AllEns %in% names(triali)[triali == T])$EnsID)
#add HMGA2 to this list too:
extra_CR <- c(extra_CR, "ENSG00000149948.13")

GeneLists <- list("TFs"= unique(filter(fpkm_allG, grepl("TF", GeneClassUpdate))$EnsID),
                  "CRs" = unique(filter(fpkm_allG, EnsName %in% CRdb_dataB$CR | EnsID %in% extra_CR)$EnsID),
                  "CC"= unique(filter(fpkm_allG, grepl("CC", GeneClassUpdate))$EnsID),
                  "IEGs"= unique(filter(fpkm_allG, EnsName %in% c(IEGs_hs[,4]) | EnsID %in% extra_IEG)$EnsID),
                  "SMC"= unique(filter(fpkm_allG, gsub("\\.[0-9]*", "", EnsID) %in% FANT_S10_SMC_G | EnsID %in% extra_SMC)$EnsID),
                  "VSMC"= unique(filter(fpkm_allG, gsub("\\.[0-9]*", "", EnsID) %in% FANT_S10_VSMC_G | EnsID %in% extra_VSMC)$EnsID)
)

results_list <- list()

#enrichment of TFs in lnc-targeted DE PCGs rather than PCGs
#same/later timeframe:
AllTargets_T2 <- unique(filter(CoRegPairs_04_48_24_extended, (corSig == "Yes" | !loopMethod == "Neither" | 
                                                                eQTLvalidations >0 | FANTOM_eQTL == "Yes" |
                                                                ExpectedCis == "Yes"),
                               EnsID %in% fpkm_allGDE_within_4$EnsID))
AllTargets_T2$EnsID_merge.y <- gsub("\\.[0-9]*", "", AllTargets_T2$EnsID.y)

#cut out HOXA/CXCL etc which massively bias
AllTargets_loci1 <- filter(AllTargets_T2, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName.y))
AllTargets_loci <- filter(AllTargets_T2, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName.y))
AllTargets_loci_T2 <- rbind(AllTargets_loci, AllTargets_loci1)

fpkm_PCGDE <- filter(fpkm_allGDE, grepl("protein_coding", EnsType))
fpkm_PCGDE$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE$EnsID)
fpkm_PCGDE_loci1 <- filter(fpkm_PCGDE, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName))
fpkm_PCGDE_loci <- filter(fpkm_PCGDE, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName))
fpkm_PCGDE_loci <- rbind(fpkm_PCGDE_loci, fpkm_PCGDE_loci1)

for (i in 1:length(GeneLists)){
  
  a <- length(unique(AllTargets_loci_T2$EnsID.y[ AllTargets_loci_T2$EnsID.y %in% GeneLists[[i]] ]))
  b <- length(unique(AllTargets_loci_T2$EnsID.y))
  c <- length(unique(fpkm_PCGDE_loci$EnsID[ fpkm_PCGDE_loci$EnsID %in% GeneLists[[i]]]))
  d <- length(unique(fpkm_PCGDE_loci$EnsID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate,
                                  paste(unique(AllTargets_loci_T2$EnsName.y[ AllTargets_loci_T2$EnsID.y %in% GeneLists[[i]] ]), collapse = "/"))
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:7] <- c("pval", "OR")
#cut non-vasc SMC
GeneLists_enriched <- GeneLists_enriched[-5,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")
#seems like essentially all these groups (aside from cell cycle) are enriched to quite strong degree 
#(low p but pool is low + strong OR)
GeneLists_enrichedSameLaterT_T2 <- GeneLists_enriched


#with non-T2 the result is much weaker right?
AllTargets <- unique(filter(CoRegPairs_04_48_24_extended,
                            EnsID %in% fpkm_allGDE_within_4$EnsID))
AllTargets$EnsID_merge.y <- gsub("\\.[0-9]*", "", AllTargets$EnsID.y)
#cut out HOXA/CXCL etc which massively bias
AllTargets_loci1 <- filter(AllTargets, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName.y))
AllTargets_loci <- filter(AllTargets, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName.y))
AllTargets_loci <- rbind(AllTargets_loci, AllTargets_loci1)

for (i in 1:length(GeneLists)){
  
  a <- length(unique(AllTargets_loci$EnsID.y[ AllTargets_loci$EnsID.y %in% GeneLists[[i]] ]))
  b <- length(unique(AllTargets_loci$EnsID.y))
  c <- length(unique(fpkm_PCGDE_loci$EnsID[ fpkm_PCGDE_loci$EnsID %in% GeneLists[[i]]]))
  d <- length(unique(fpkm_PCGDE_loci$EnsID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate,
                                  paste(unique(AllTargets_loci$EnsName.y[ AllTargets_loci$EnsID.y %in% GeneLists[[i]] ]), collapse = "/"))
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:7] <- c("pval", "OR")
#cut non-vasc SMC
GeneLists_enriched <- GeneLists_enriched[-5,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")
#weaker results without T2
GeneLists_enrichedSameLaterT <- GeneLists_enriched

#note FOXL1 is called an IEG (fusion to FOXC2 so not quite but kindof)


#same timeframe result:
EarlyTargets_T2 <- unique(filter(CoRegPairs_04_48_24_extended, (corSig == "Yes" | !loopMethod == "Neither" | 
                                                                  eQTLvalidations >0 | FANTOM_eQTL == "Yes" |
                                                                  ExpectedCis == "Yes"),
                                 EnsID %in% fpkm_allGDE_within_4$EnsID, EnsID.y %in% fpkm_allGDE_within_4$EnsID))
EarlyTargets_T2$EnsID_merge.y <- gsub("\\.[0-9]*", "", EarlyTargets_T2$EnsID.y)

#cut out HOXA/CXCL etc which massively bias
EarlyTargets_loci1 <- filter(EarlyTargets_T2, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName.y))
EarlyTargets_loci <- filter(EarlyTargets_T2, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName.y))
EarlyTargets_loci_T2 <- rbind(EarlyTargets_loci, EarlyTargets_loci1)

fpkm_PCGDE <- filter(fpkm_allGDE_within_4, grepl("protein_coding", EnsType))
fpkm_PCGDE$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE$EnsID)
fpkm_PCGDE_loci1 <- filter(fpkm_PCGDE, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName))
fpkm_PCGDE_loci <- filter(fpkm_PCGDE, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName))
fpkm_PCGDE_loci <- rbind(fpkm_PCGDE_loci, fpkm_PCGDE_loci1)

for (i in 1:length(GeneLists)){
  
  a <- length(unique(EarlyTargets_loci_T2$EnsID.y[ EarlyTargets_loci_T2$EnsID.y %in% GeneLists[[i]] ]))
  b <- length(unique(EarlyTargets_loci_T2$EnsID.y))
  c <- length(unique(fpkm_PCGDE_loci$EnsID[ fpkm_PCGDE_loci$EnsID %in% GeneLists[[i]]]))
  d <- length(unique(fpkm_PCGDE_loci$EnsID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate,
                                  paste(unique(EarlyTargets_loci_T2$EnsName.y[ EarlyTargets_loci_T2$EnsID.y %in% GeneLists[[i]] ]), collapse = "/"))
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:7] <- c("pval", "OR")
#cut non-vasc SMC
GeneLists_enriched <- GeneLists_enriched[-5,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")
#no significance!
GeneLists_enrichedSameT_T2 <- GeneLists_enriched


#plot enrichment of interesting subtypes in same/later tier2:
GeneLists_enrichedSameLaterT_T2$selectHit <- GeneLists_enrichedSameLaterT_T2$a/GeneLists_enrichedSameLaterT_T2$b
GeneLists_enrichedSameLaterT_T2$backHit <- GeneLists_enrichedSameLaterT_T2$c/GeneLists_enrichedSameLaterT_T2$d

library(reshape2)
plot_Gclass <- melt(GeneLists_enrichedSameLaterT_T2[,c(1,10,11)])

plot_Gclass$GeneList <- c("TFs", "Chromatin Regulators", "Core S/G2M", "IEGs", "VSMC-enriched", 
                          "TFs", "Chromatin Regulators", "Core S/G2M", "IEGs", "VSMC-enriched")
plot_Gclass$GeneList <- as.factor(plot_Gclass$GeneList)
plot_Gclass$GeneList <- factor(plot_Gclass$GeneList, levels = levels(plot_Gclass$GeneList)[c(2,5,3,1,4)])

plot_Gclass$variable <- c(rep("SCClncRNA Targets", 5), rep("All DE PCGs", 5))
#plot_Gclass$variable <- as.factor(plot_Gclass$variable)
#plot_Gclass$variable <- factor(plot_Gclass$variable, levels = levels(plot_Gclass$variable)[c(2,1)])

plot_Gclass$value <- plot_Gclass$value*100

ggplot(plot_Gclass) + aes(y = GeneList, x = value, fill = variable) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c(`SCClncRNA Targets` = "mediumorchid",`All DE PCGs` = "grey60")) +
  theme_minimal() +
  ylab("") +
  xlab("%") +
  theme(text = element_text(size=24))

#is it fair? biased to early rather than all...
dim(filter(fpkm_PCGDE_loci, EnsID %in% fpkm_allGDE_within_4$EnsID))
1579/4059 #39% of PCGs early-reg
table(unique(AllTargets_loci_T2[,c(10,13)])$PCG_Cluster)
sum(table(unique(AllTargets_loci_T2[,c(10,13)])$PCG_Cluster))
59/99 #60% of lnc targets early-reg

#it is part of the explanation for sure
#the ORs for same-timeframe are still quite high, just low pool numbers
#either way, it's still showing and contextualising their potential targets


#### Fold change revisit ####
#### old code ####
#### TF enrichment ####

#LncRNAs reported to favour localisation near TFs genomically

#Can we see an enrichment of TFs amongst 0-4hr CClnc targets vs. other 0-4hr genes?

length(unique(filter(CoRegPairs_04_48_24_extended, Lnc_Timeframe == "<4hrs", PCG_Timeframe == "<4hrs")$EnsID.y)) #83 targets of 0-4hr cclncs
length(unique(filter(CoRegPairs_04_48_24_extended, Lnc_Timeframe == "<4hrs", PCG_Timeframe == "<4hrs", 
                     grepl("TF", GeneClassUpdate.y))$EnsID.y)) #32 TF CClnc targets

19/83 # 22.9% 

fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE.csv", header = T)
fpkm_allGDE_filt <- filter(fpkm_allGDE, grepl("chr", chr), !grepl("artefacts", GeneClassUpdate), !is.na(GeneClassUpdate))
#remove manual fails
fpkm_allGDE <- filter(fpkm_allGDE_filt, EnsID %in% fpkm_allG_filt_manual$EnsID)

#clusters needed for later
#regulated within 4 hours:
fpkm_allGDE_Upwithin_4 <- filter(fpkm_allGDE, (Hour0_meanFPKM >=1 | Hour4_meanFPKM >=1) &
                                   (LogFC_0_4 >= log2(1.5) & preadj_0_4 <0.05))
fpkm_allGDE_Downwithin_4 <- filter(fpkm_allGDE, (Hour0_meanFPKM >=1 | Hour4_meanFPKM >=1) &
                                     (LogFC_0_4 < -log2(1.5) & preadj_0_4 <0.05))

fpkm_allGDE_within_4 <- rbind(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Downwithin_4)

fpkm_PCGDE_within_4 <- filter(fpkm_allGDE_within_4, grepl("coding|TF|CC", GeneClassUpdate))

length(unique(fpkm_PCGDE_within_4$EnsID)) #1719 0-4hr DEGs
length(unique(filter(fpkm_PCGDE_within_4, grepl("TF", GeneClassUpdate))$EnsID)) #188 TF CClnc targets

188/1609 #11.7%


#possibility for an enrichment:
a <- 19
b <- 83
c <- 188
d <- 1609

fisher.test(data.frame("cisLncTarget" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#cclnc targets within 4hrs are enriched with TFs (probs due to Hox locus)


#amongst induced 0-4hr cclnc targets:
length(unique(filter(CoRegPairs_04_48_24_extended, Lnc_Cluster == "Induced <4hrs", PCG_Timeframe == "<4hrs")$EnsID.y)) #55 targets of 0-4hr cclncs
length(unique(filter(CoRegPairs_04_48_24_extended, Lnc_Cluster == "Induced <4hrs", PCG_Timeframe == "<4hrs", 
                     grepl("TF", GeneClassUpdate.y))$EnsID.y)) #32 TF CClnc targets

9/55 # 16.4% 

fpkm_PCGDE_Upwithin_4 <- filter(fpkm_allGDE_Upwithin_4, grepl("coding|TF|CC", GeneClassUpdate))
length(unique(fpkm_PCGDE_Upwithin_4$EnsID)) #973 0-4hr induced DEGs
length(unique(filter(fpkm_PCGDE_Upwithin_4, grepl("TF", GeneClassUpdate))$EnsID)) #188 TF CClnc targets

100/904 #10.9%

a <- 9
b <- 55
c <- 100
d <- 904

fisher.test(data.frame("cisLncTarget" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#ns


#using SCClncRNAs:
Tier2 <- filter(CoRegPairs_04_48_24_extended, !loopMethod == "Neither" | corSig == "Yes" | eQTL_supported == "Yes")

length(unique(filter(Tier2, Lnc_Timeframe == "<4hrs", PCG_Timeframe == "<4hrs")$EnsID.y)) #93 CClnc targets
length(unique(filter(Tier2, Lnc_Timeframe == "<4hrs", PCG_Timeframe == "<4hrs", grepl("TF", GeneClassUpdate.y))$EnsID.y)) #20 TF CClnc targets

14/59

a <- 14
b <- 59
c <- 188
d <- 1609

fisher.test(data.frame("cisLncTarget" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater") # p = 0.007, OR 2.38

#roughly the same T2


#### IEG enrichment ####

#involved with expected early co-ordinators?

#Can we see an enrichment of TFs amongst 0-4hr CClnc targets vs. other 0-4hr genes?

#IEGs:
IEGs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Gene lists/arner_2015_table_S5_IEGs_lit.csv")
IEGs_hs <- filter(IEGs, !Hs_symbol %in% c("", "-"))
IEGs_hs <- filter(IEGs_hs, !grepl("[a-z]", IEGs_hs$Total_count))

length(unique(filter(CoRegPairs_04_48_24_extended, Lnc_Timeframe == "<4hrs", PCG_Timeframe == "<4hrs")$EnsID.y)) #83 targets of 0-4hr cclncs
length(unique(filter(CoRegPairs_04_48_24_extended, Lnc_Timeframe == "<4hrs", PCG_Timeframe == "<4hrs", 
                     EnsName.y %in% IEGs_hs$Hs_symbol)$EnsID.y)) #32 TF CClnc targets
9/83 # 10.8%

fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE.csv", header = T)
fpkm_allGDE_filt <- filter(fpkm_allGDE, grepl("chr", chr), !grepl("artefacts", GeneClassUpdate), !is.na(GeneClassUpdate))
#remove manual fails
fpkm_allGDE <- filter(fpkm_allGDE_filt, EnsID %in% fpkm_allG_filt_manual$EnsID)

#clusters needed for later
#regulated within 4 hours:
fpkm_allGDE_Upwithin_4 <- filter(fpkm_allGDE, (Hour0_meanFPKM >=1 | Hour4_meanFPKM >=1) &
                                   (LogFC_0_4 >= log2(1.5) & preadj_0_4 <0.05))
fpkm_allGDE_Downwithin_4 <- filter(fpkm_allGDE, (Hour0_meanFPKM >=1 | Hour4_meanFPKM >=1) &
                                     (LogFC_0_4 < -log2(1.5) & preadj_0_4 <0.05))

fpkm_allGDE_within_4 <- rbind(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Downwithin_4)

fpkm_PCGDE_within_4 <- filter(fpkm_allGDE_within_4, grepl("coding|TF|CC", GeneClassUpdate))

length(unique(fpkm_PCGDE_within_4$EnsID)) #1609 0-4hr DEGs
length(unique(filter(fpkm_PCGDE_within_4, EnsName %in% IEGs_hs$Hs_symbol)$EnsID)) #188 TF CClnc targets

109/1609 #6.8%


#possibility for an enrichment:
a <- 9
b <- 83
c <- 109
d <- 1609

fisher.test(data.frame("cisLncTarget" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#cclnc targets within 4hrs are enriched with TFs (probs due to Hox locus)


#amongst induced 0-4hr cclnc targets:
length(unique(filter(CoRegPairs_04_48_24_extended, Lnc_Cluster == "Induced <4hrs", PCG_Timeframe == "<4hrs")$EnsID.y)) #55 targets of 0-4hr cclncs
length(unique(filter(CoRegPairs_04_48_24_extended, Lnc_Cluster == "Induced <4hrs", PCG_Timeframe == "<4hrs", 
                     EnsName.y %in% IEGs_hs$Hs_symbol)$EnsID.y)) #32 TF CClnc targets

8/55 # 14.5% 

length(unique(fpkm_allGDE_Upwithin_4$EnsID)) #973 0-4hr induced DEGs
length(unique(filter(fpkm_allGDE_Upwithin_4, EnsName %in% IEGs_hs$Hs_symbol)$EnsID)) #188 TF CClnc targets

92/973 #9.5%

a <- 8
b <- 55
c <- 92
d <- 973

fisher.test(data.frame("cisLncTarget" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#ns
