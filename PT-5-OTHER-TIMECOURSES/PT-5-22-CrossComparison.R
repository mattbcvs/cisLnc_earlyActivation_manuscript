library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(grid)
library(gridExtra)
library(reshape2)

#### importing gene sets - SVSMC early induced DEGs + early induced CClncRNAs ####

fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026filt.csv", header = T)

fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE_2026filt.csv", header = T)

fpkm_allGDE_Upwithin_4 <- filter(fpkm_allGDE, RegulationStart == "Induced <4hrs")
fpkm_allGDE_Downwithin_4 <- filter(fpkm_allGDE, RegulationStart == "Repressed <4hrs")

fpkm_allGDE_Upwithin_8 <- filter(fpkm_allGDE, RegulationStart == "Induced 4-8hrs")
fpkm_allGDE_Downwithin_8 <- filter(fpkm_allGDE, RegulationStart == "Repressed 4-8hrs")

fpkm_allGDE_Upwithin_24 <- filter(fpkm_allGDE, RegulationStart == "Induced 8-24hrs")
fpkm_allGDE_Downwithin_24 <- filter(fpkm_allGDE, RegulationStart == "Repressed 8-24hrs")

#above lines seperates genes into distinct buckets:
table(fpkm_allGDE$RegulationStart) 
#1105/828/821
#891/768/668

AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_info_2026.csv", header = T)
table(AllLNC_AllPCG_1$LncRNA.PCG.Relationship)#64/274/8097/3
AllLNC_AllPCG_1 <- filter(AllLNC_AllPCG_1, AbsDistLnc_PCG <250)#2564
table(AllLNC_AllPCG_1$LncRNA.PCG.Relationship)#62/274/2225/3

Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime_2026.csv", header = T)
AllLNC_AllPCG_1$GeneClassUpdate.x[AllLNC_AllPCG_1$EnsID %in% 
                                    filter(Enhancer_lociII, !is.na(EnhancerVerdict))$EnsID & 
                                    AllLNC_AllPCG_1$GeneClassUpdate.x == "Bona fide lncRNA"] <- "ELnc"

closestNeighbour <- filter(AllLNC_AllPCG_1, AbsDistLnc_PCG <250)
closestNeighbour <- split(closestNeighbour, closestNeighbour$EnsID)

#closest "surrounding" neighbours (consider up and down-stream):
#this selected for best "naive" results
closestNeighbour <- lapply(closestNeighbour, function(x){
  upstream <- filter(x, DistLnc_PCG < 0)
  downstream <- filter(x, DistLnc_PCG > 0)
  
  return(rbind(upstream[order(upstream$DistLnc_PCG, decreasing = T),][1,],
               downstream[order(downstream$DistLnc_PCG, decreasing = F),][1,]))
}
)

closestNeighbour <- bind_rows(closestNeighbour)
closestNeighbour <- filter(closestNeighbour, !is.na(pairs))

#now get the cis candidates:
AllLNC_AllPCG_1 <- closestNeighbour
CoRegPairs_04_48_24_extended_naiveSame <- filter(AllLNC_AllPCG_1,
                                                 #AllLNC_AllPCG_1,
                                                 (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                               fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                                fpkm_allGDE_Downwithin_4$EnsID)) |
                                                   (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                 fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                                  fpkm_allGDE_Downwithin_8$EnsID)) |
                                                   (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                 fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
                                                 )
#55 @250kbp - surrounding (29% bidirectional)
table(CoRegPairs_04_48_24_extended_naiveSame$LncRNA.PCG.Relationship)


#
#### importing gene sets - epiMT - early induced DEGs + early induced CClncRNAs ####

epiMT_DEGs_DE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GeneralEffect_FANTOM_ER/epiMT_DEGs_DE_2026.csv")
AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GeneralEffect_FANTOM_ER/AllLNC_AllPCG_1_epiMT_2026.csv")
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



#### importing gene sets - moLPS - early induced DEGs + early induced CClncRNAs ####

moLPS_DEGs_DE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GeneralEffect_FANTOM_ER/moLPS_DEGs_DE_2026.csv")
AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GeneralEffect_FANTOM_ER/AllLNC_AllPCG_1_moLPS_2026.csv")
GeneBiotypes <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/nature21374-s2/table11.csv")
#simplify lncRNA class:
GeneBiotypes$CAT_geneClassII <- GeneBiotypes$CAT_geneClass
GeneBiotypes$CAT_geneClassII[
  grepl("lncRNA", GeneBiotypes$CAT_geneClassII)] <- "lncRNA"

#n.b. (START FROM HERE ON RETURN)
moLPS_DEGs_DE_Upwithin_2 <- filter(moLPS_DEGs_DE, RegulationStart == "Induced <2hrs")
moLPS_DEGs_DE_Downwithin_2 <- filter(moLPS_DEGs_DE, RegulationStart == "Repressed <2hrs")

moLPS_DEGs_DE_Upwithin_4 <- filter(moLPS_DEGs_DE, RegulationStart == "Induced 2-4hrs")
moLPS_DEGs_DE_Downwithin_4 <- filter(moLPS_DEGs_DE, RegulationStart == "Repressed 2-4hrs")

moLPS_DEGs_DE_Upwithin_6 <- filter(moLPS_DEGs_DE, RegulationStart == "Induced 4-6hrs")
moLPS_DEGs_DE_Downwithin_6 <- filter(moLPS_DEGs_DE, RegulationStart == "Repressed 4-6hrs")

moLPS_DEGs_DE_Upwithin_8 <- filter(moLPS_DEGs_DE, RegulationStart == "Induced 6-8hrs")
moLPS_DEGs_DE_Downwithin_8 <- filter(moLPS_DEGs_DE, RegulationStart == "Repressed 6-8hrs")

moLPS_DEGs_DE_Upwithin_12 <- filter(moLPS_DEGs_DE, RegulationStart == "Induced 8-12hrs")
moLPS_DEGs_DE_Downwithin_12 <- filter(moLPS_DEGs_DE, RegulationStart == "Repressed 8-12hrs")

moLPS_DEGs_DE_Upwithin_16 <- filter(moLPS_DEGs_DE, RegulationStart == "Induced 12-16hrs")
moLPS_DEGs_DE_Downwithin_16 <- filter(moLPS_DEGs_DE, RegulationStart == "Repressed 12-16hrs")

moLPS_DEGs_DE_Upwithin_24 <- filter(moLPS_DEGs_DE, RegulationStart == "Induced 16-24hrs")
moLPS_DEGs_DE_Downwithin_24 <- filter(moLPS_DEGs_DE, RegulationStart == "Repressed 16-24hrs")

moLPS_DEGs_DE_Upwithin_36 <- filter(moLPS_DEGs_DE, RegulationStart == "Induced 24-36hrs")
moLPS_DEGs_DE_Downwithin_36 <- filter(moLPS_DEGs_DE, RegulationStart == "Repressed 24-36hrs")

moLPS_DEGs_DE_Upwithin_48 <- filter(moLPS_DEGs_DE, RegulationStart == "Induced 36-48hrs")
moLPS_DEGs_DE_Downwithin_48 <- filter(moLPS_DEGs_DE, RegulationStart == "Repressed 36-48hrs")

closestNeighbour_Lncs_moLPS <- split(AllLNC_AllPCG_1, AllLNC_AllPCG_1$AnchorGene)

closestNeighbour_Lncs_moLPS <- lapply(closestNeighbour_Lncs_moLPS, function(x){
  upstream <- filter(x, DistLnc_PCG < 0)
  downstream <- filter(x, DistLnc_PCG > 0)
  
  return(rbind(upstream[order(upstream$DistLnc_PCG, decreasing = T),][1,],
               downstream[order(downstream$DistLnc_PCG, decreasing = F),][1,]))
}
)

closestNeighbour_Lncs_moLPS <- bind_rows(closestNeighbour_Lncs_moLPS)
closestNeighbour_Lncs_moLPS <- filter(closestNeighbour_Lncs_moLPS, !is.na(pairs))

#now get the cis candidates:
#same timeframe only for simplicity
CoRegPairs_sameTimeframe_moLPS <- filter(closestNeighbour_Lncs_moLPS, 
                                         (AnchorGene %in% c(moLPS_DEGs_DE_Upwithin_2$CAT_geneID, moLPS_DEGs_DE_Downwithin_2$CAT_geneID) & 
                                            geneID %in% c(moLPS_DEGs_DE_Upwithin_2$CAT_geneID, moLPS_DEGs_DE_Downwithin_2$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(moLPS_DEGs_DE_Upwithin_4$CAT_geneID, moLPS_DEGs_DE_Downwithin_4$CAT_geneID) & 
                                              geneID %in% c(moLPS_DEGs_DE_Upwithin_4$CAT_geneID, moLPS_DEGs_DE_Downwithin_4$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(moLPS_DEGs_DE_Upwithin_6$CAT_geneID, moLPS_DEGs_DE_Downwithin_6$CAT_geneID) & 
                                              geneID %in% c(moLPS_DEGs_DE_Upwithin_6$CAT_geneID, moLPS_DEGs_DE_Downwithin_6$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(moLPS_DEGs_DE_Upwithin_8$CAT_geneID, moLPS_DEGs_DE_Downwithin_8$CAT_geneID) & 
                                              geneID %in% c(moLPS_DEGs_DE_Upwithin_8$CAT_geneID, moLPS_DEGs_DE_Downwithin_8$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(moLPS_DEGs_DE_Upwithin_12$CAT_geneID, moLPS_DEGs_DE_Downwithin_12$CAT_geneID) & 
                                              geneID %in% c(moLPS_DEGs_DE_Upwithin_12$CAT_geneID, moLPS_DEGs_DE_Downwithin_12$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(moLPS_DEGs_DE_Upwithin_16$CAT_geneID, moLPS_DEGs_DE_Downwithin_16$CAT_geneID) & 
                                              geneID %in% c(moLPS_DEGs_DE_Upwithin_16$CAT_geneID, moLPS_DEGs_DE_Downwithin_16$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(moLPS_DEGs_DE_Upwithin_24$CAT_geneID, moLPS_DEGs_DE_Downwithin_24$CAT_geneID) & 
                                              geneID %in% c(moLPS_DEGs_DE_Upwithin_24$CAT_geneID, moLPS_DEGs_DE_Downwithin_24$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(moLPS_DEGs_DE_Upwithin_36$CAT_geneID, moLPS_DEGs_DE_Downwithin_36$CAT_geneID) & 
                                              geneID %in% c(moLPS_DEGs_DE_Upwithin_36$CAT_geneID, moLPS_DEGs_DE_Downwithin_36$CAT_geneID)) |
                                           
                                           (AnchorGene %in% c(moLPS_DEGs_DE_Upwithin_48$CAT_geneID, moLPS_DEGs_DE_Downwithin_48$CAT_geneID) & 
                                              geneID %in% c(moLPS_DEGs_DE_Upwithin_48$CAT_geneID, moLPS_DEGs_DE_Downwithin_48$CAT_geneID))
)
#74 CCLnc-target pairs

length(unique(CoRegPairs_sameTimeframe_moLPS$AnchorGene))#71 CClncRNAs
length(unique(CoRegPairs_sameTimeframe_moLPS$geneID))#65 potential targets

#
#### for venn diagram ####

#matching on a) EnsID, b) SYMBOL and c) FANTOM ID

#SVSMC lncs matched to FANTOM:
Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime.csv", header = T)

Lnc_EnsFANT <- unique(filter(Enhancer_lociII, !is.na(CAGEvalidity))[,c(1,14)])

#
#### early phase-induced  lncRNAs ####
epiMT_earlyLncs <- filter(epiMT_DEGs_DE_Upwithin_2, CAT_geneClassII == "lncRNA")[,c(1,2)]
MoLPS_earlyLncs <- filter(moLPS_DEGs_DE_Upwithin_2, CAT_geneClassII == "lncRNA")[,c(1,2)]
SVSMC_earlyLncs <- filter(fpkm_allGDE_Upwithin_4, grepl("Lnc|fide", GeneClassUpdate))[,c(1,2)]

epiMT_earlyLncs$CAT_geneID_merge <- gsub("\\.[0-9]*", "", epiMT_earlyLncs$CAT_geneID)
MoLPS_earlyLncs$CAT_geneID_merge <- gsub("\\.[0-9]*", "", MoLPS_earlyLncs$CAT_geneID)

SVSMC_earlyLncs <- unique(merge(SVSMC_earlyLncs, Lnc_EnsFANT, by = "EnsID", all.x = T))
SVSMC_earlyLncs$EnsID_merge <- gsub("\\.[0-9]*", "", SVSMC_earlyLncs$EnsID)

#if FANTOM ID is missing, input the ENSEMBL:
SVSMC_earlyLncs$FANTOM_ID[is.na(SVSMC_earlyLncs$FANTOM_ID)] <- SVSMC_earlyLncs$EnsID_merge[is.na(SVSMC_earlyLncs$FANTOM_ID)]

#remove the suffix from all:
SVSMC_earlyLncs$FANTOM_ID_merge <- gsub("\\.[0-9]*", "", SVSMC_earlyLncs$FANTOM_ID)

#remaining "MSTRG" genes could not be mapped to a FANTOM annotation, return original name
SVSMC_earlyLncs$FANTOM_ID_merge[SVSMC_earlyLncs$FANTOM_ID_merge == "MSTRG"] <- SVSMC_earlyLncs$EnsID[SVSMC_earlyLncs$FANTOM_ID_merge == "MSTRG"]

#venn:
library(ggVennDiagram)

venn_input <- list("SVSMC" = SVSMC_earlyLncs$FANTOM_ID_merge, 
                   "EpiMT" = epiMT_earlyLncs$CAT_geneID_merge, 
                   "MoLPS" = MoLPS_earlyLncs$CAT_geneID_merge)

ggVennDiagram(venn_input, label_alpha = 0.5, label = "count", edge_size = 0) +
  ggplot2::theme_void()

common_earlyLncs <- filter(SVSMC_earlyLncs, FANTOM_ID_merge %in% Reduce(intersect, venn_input))
SVep_earlyLncs <- filter(SVSMC_earlyLncs, FANTOM_ID_merge %in% intersect(venn_input[[1]], venn_input[[2]]))
SVmo_earlyLncs <- filter(SVSMC_earlyLncs, FANTOM_ID_merge %in% intersect(venn_input[[1]], venn_input[[3]]))

#pretty clear that few if any SVSMC lncRNAs are activated in 4hrs for these other 2

#notable are RP11-221N13.3, ATP2B1-AS1, MSTRG.28277 which are naive cis lncs (latter two in molps too)


#
#### early phase PCGs ####

epiMT_earlyPCGs <- filter(epiMT_DEGs_DE_Upwithin_2, CAT_geneClassII == "coding_mRNA")[,c(1,2)]
MoLPS_earlyPCGs <- filter(moLPS_DEGs_DE_Upwithin_2, CAT_geneClassII == "coding_mRNA")[,c(1,2)]
SVSMC_earlyPCGs <- filter(fpkm_allGDE_Upwithin_4, grepl("coding|CC|TF", GeneClassUpdate))[,c(1,2)]

epiMT_earlyPCGs$CAT_geneID_merge <- gsub("\\.[0-9]*", "", epiMT_earlyPCGs$CAT_geneID)
MoLPS_earlyPCGs$CAT_geneID_merge <- gsub("\\.[0-9]*", "", MoLPS_earlyPCGs$CAT_geneID)
SVSMC_earlyPCGs$EnsID_merge <- gsub("\\.[0-9]*", "", SVSMC_earlyPCGs$EnsID)

SVSMC_earlyPCGs <- filter(SVSMC_earlyPCGs, !grepl("MSTRG|Split", EnsID))

#fine to merge on ensembl here (smc PCGs have not been matched to FANTOM annotaiton, unecessary):
length(unique(SVSMC_earlyPCGs$EnsID_merge))

#venn:
library(ggVennDiagram)

venn_input <- list("SVSMC" = SVSMC_earlyPCGs$EnsID_merge, 
                   "EpiMT" = epiMT_earlyPCGs$CAT_geneID_merge, 
                   "MoLPS" = MoLPS_earlyPCGs$CAT_geneID_merge)

ggVennDiagram(venn_input, label_alpha = 0.5, label = "count", edge_size = 0) +
  ggplot2::theme_void()

sapply(venn_input, length)

common_earlyPCGs <- filter(SVSMC_earlyPCGs, EnsID_merge %in% Reduce(intersect, venn_input))
SVep_earlyPCGs <- filter(SVSMC_earlyPCGs, EnsID_merge %in% intersect(venn_input[[1]], venn_input[[2]]))
SVmo_earlyPCGs <- filter(SVSMC_earlyPCGs, EnsID_merge %in% intersect(venn_input[[1]], venn_input[[3]]))



#
#### early phase co-reg lncs ####

epiMT_earlyLncs <- filter(epiMT_DEGs_DE_Upwithin_2, CAT_geneClassII == "lncRNA")[,c(1,2)]
MoLPS_earlyLncs <- filter(moLPS_DEGs_DE_Upwithin_2, CAT_geneClassII == "lncRNA")[,c(1,2)]
SVSMC_earlyLncs <- filter(fpkm_allGDE_Upwithin_4, grepl("Lnc|fide", GeneClassUpdate))[,c(1,2)]

epiMT_earlyLncs$CAT_geneID_merge <- gsub("\\.[0-9]*", "", epiMT_earlyLncs$CAT_geneID)
MoLPS_earlyLncs$CAT_geneID_merge <- gsub("\\.[0-9]*", "", MoLPS_earlyLncs$CAT_geneID)

SVSMC_earlyLncs <- unique(merge(SVSMC_earlyLncs, Lnc_EnsFANT, by = "EnsID", all.x = T))
SVSMC_earlyLncs$EnsID_merge <- gsub("\\.[0-9]*", "", SVSMC_earlyLncs$EnsID)

#if FANTOM ID is missing, input the ENSEMBL:
SVSMC_earlyLncs$FANTOM_ID[is.na(SVSMC_earlyLncs$FANTOM_ID)] <- SVSMC_earlyLncs$EnsID_merge[is.na(SVSMC_earlyLncs$FANTOM_ID)]

#remove the suffix from all:
SVSMC_earlyLncs$FANTOM_ID_merge <- gsub("\\.[0-9]*", "", SVSMC_earlyLncs$FANTOM_ID)

#remaining "MSTRG" genes could not be mapped, return original name
SVSMC_earlyLncs$FANTOM_ID_merge[SVSMC_earlyLncs$FANTOM_ID_merge == "MSTRG"] <- SVSMC_earlyLncs$EnsID[SVSMC_earlyLncs$FANTOM_ID_merge == "MSTRG"]

#venn:
library(ggVennDiagram)

#coReg: same up early
SVSMC_earlyLncs_sameUp <- filter(SVSMC_earlyLncs, EnsID %in% filter(CoRegPairs_04_48_24_extended_naiveSame, EnsID %in% fpkm_allGDE_Upwithin_4$EnsID, 
                                                                  EnsID.y %in% fpkm_allGDE_Upwithin_4$EnsID)$EnsID)

epiMT_earlyLncs_sameUp <- filter(epiMT_earlyLncs, CAT_geneID %in% filter(CoRegPairs_sameTimeframe_epiMT,
                                                                       AnchorGene %in% epiMT_DEGs_DE_Upwithin_2$CAT_geneID,
                                                                       geneID %in% epiMT_DEGs_DE_Upwithin_2$CAT_geneID
                                                                       )$AnchorGene)

MoLPS_earlyLncs_sameUp <- filter(MoLPS_earlyLncs, CAT_geneID %in% filter(CoRegPairs_sameTimeframe_moLPS,
                                                                       AnchorGene %in% moLPS_DEGs_DE_Upwithin_2$CAT_geneID,
                                                                       geneID %in% moLPS_DEGs_DE_Upwithin_2$CAT_geneID
                                                                       )$AnchorGene)

venn_input <- list("SVSMC" = SVSMC_earlyLncs_sameUp$FANTOM_ID_merge, 
                   "EpiMT" = epiMT_earlyLncs_sameUp$CAT_geneID_merge, 
                   "MoLPS" = MoLPS_earlyLncs_sameUp$CAT_geneID_merge)

ggVennDiagram(venn_input, label_alpha = 0.5, label = "count", edge_size = 0) +
  ggplot2::theme_void()

common_earlyLncs <- filter(SVSMC_earlyLncs, FANTOM_ID_merge %in% Reduce(intersect, venn_input))
SVep_earlyLncs <- filter(SVSMC_earlyLncs, FANTOM_ID_merge %in% intersect(venn_input[[1]], venn_input[[2]]))
SVmo_earlyLncs <- filter(SVSMC_earlyLncs, FANTOM_ID_merge %in% intersect(venn_input[[1]], venn_input[[3]]))

#largely unique with notable exceptions


#
#### fairer version, allow 4hr pairs in too ####

epiMT_earlyLncs <- rbind(filter(epiMT_DEGs_DE_Upwithin_2, CAT_geneClassII == "lncRNA")[,c(1,2)],
                         filter(epiMT_DEGs_DE_Upwithin_4, CAT_geneClassII == "lncRNA")[,c(1,2)])
MoLPS_earlyLncs <- rbind(filter(moLPS_DEGs_DE_Upwithin_2, CAT_geneClassII == "lncRNA")[,c(1,2)],
                         filter(moLPS_DEGs_DE_Upwithin_4, CAT_geneClassII == "lncRNA")[,c(1,2)])
SVSMC_earlyLncs <- filter(fpkm_allGDE_Upwithin_4, grepl("Lnc|fide", GeneClassUpdate))[,c(1,2)]

epiMT_earlyLncs$CAT_geneID_merge <- gsub("\\.[0-9]*", "", epiMT_earlyLncs$CAT_geneID)
MoLPS_earlyLncs$CAT_geneID_merge <- gsub("\\.[0-9]*", "", MoLPS_earlyLncs$CAT_geneID)

SVSMC_earlyLncs <- unique(merge(SVSMC_earlyLncs, Lnc_EnsFANT, by = "EnsID", all.x = T))
SVSMC_earlyLncs$EnsID_merge <- gsub("\\.[0-9]*", "", SVSMC_earlyLncs$EnsID)

#if FANTOM ID is missing, input the ENSEMBL:
SVSMC_earlyLncs$FANTOM_ID[is.na(SVSMC_earlyLncs$FANTOM_ID)] <- SVSMC_earlyLncs$EnsID_merge[is.na(SVSMC_earlyLncs$FANTOM_ID)]

#remove the suffix from all:
SVSMC_earlyLncs$FANTOM_ID_merge <- gsub("\\.[0-9]*", "", SVSMC_earlyLncs$FANTOM_ID)

#remaining "MSTRG" genes could not be mapped, return original name
SVSMC_earlyLncs$FANTOM_ID_merge[SVSMC_earlyLncs$FANTOM_ID_merge == "MSTRG"] <- SVSMC_earlyLncs$EnsID[SVSMC_earlyLncs$FANTOM_ID_merge == "MSTRG"]

#venn:
library(ggVennDiagram)

#coReg: same up early
SVSMC_earlyLncs_sameUp <- filter(SVSMC_earlyLncs, EnsID %in% filter(CoRegPairs_04_48_24_extended_naiveSame, EnsID %in% fpkm_allGDE_Upwithin_4$EnsID, 
                                                                    EnsID.y %in% fpkm_allGDE_Upwithin_4$EnsID)$EnsID)

#for these two can either: take same timeframe pairs in 2 or 4 hrs
epiMT_earlyLncs_sameUp <- filter(epiMT_earlyLncs, 
                                 CAT_geneID %in% filter(CoRegPairs_sameTimeframe_epiMT,
                                                        AnchorGene %in% epiMT_DEGs_DE_Upwithin_2$CAT_geneID,
                                                        geneID %in% epiMT_DEGs_DE_Upwithin_2$CAT_geneID)$AnchorGene |
                                   CAT_geneID %in% filter(CoRegPairs_sameTimeframe_epiMT,
                                                          AnchorGene %in% epiMT_DEGs_DE_Upwithin_4$CAT_geneID,
                                                          geneID %in% epiMT_DEGs_DE_Upwithin_4$CAT_geneID)$AnchorGene)

MoLPS_earlyLncs_sameUp <- filter(MoLPS_earlyLncs, 
                                 CAT_geneID %in% filter(CoRegPairs_sameTimeframe_moLPS,
                                                        AnchorGene %in% moLPS_DEGs_DE_Upwithin_2$CAT_geneID,
                                                        geneID %in% moLPS_DEGs_DE_Upwithin_2$CAT_geneID)$AnchorGene |
                                   CAT_geneID %in% filter(CoRegPairs_sameTimeframe_moLPS,
                                                          AnchorGene %in% moLPS_DEGs_DE_Upwithin_4$CAT_geneID,
                                                          geneID %in% moLPS_DEGs_DE_Upwithin_4$CAT_geneID)$AnchorGene)

venn_input <- list("SVSMC" = SVSMC_earlyLncs_sameUp$FANTOM_ID_merge, 
                   "EpiMT" = epiMT_earlyLncs_sameUp$CAT_geneID_merge, 
                   "MoLPS" = MoLPS_earlyLncs_sameUp$CAT_geneID_merge)

ggVennDiagram(venn_input, label_alpha = 0.5, label = "count", edge_size = 0) +
  ggplot2::theme_void()

common_earlyLncs <- filter(SVSMC_earlyLncs, FANTOM_ID_merge %in% Reduce(intersect, venn_input))
SVep_earlyLncs <- filter(SVSMC_earlyLncs, FANTOM_ID_merge %in% intersect(venn_input[[1]], venn_input[[2]]))
SVmo_earlyLncs <- filter(SVSMC_earlyLncs, FANTOM_ID_merge %in% intersect(venn_input[[1]], venn_input[[3]]))

#largely unique with notable exceptions


#
#### early phase co-reg lnc neighbours ####
epiMT_earlyPCGs <- filter(epiMT_DEGs_DE_Upwithin_2, CAT_geneClassII == "coding_mRNA")[,c(1,2)]
MoLPS_earlyPCGs <- filter(moLPS_DEGs_DE_Upwithin_2, CAT_geneClassII == "coding_mRNA")[,c(1,2)]
SVSMC_earlyPCGs <- filter(fpkm_allGDE_Upwithin_4, grepl("coding|CC|TF", GeneClassUpdate))[,c(1,2)]

epiMT_earlyPCGs$CAT_geneID_merge <- gsub("\\.[0-9]*", "", epiMT_earlyPCGs$CAT_geneID)
MoLPS_earlyPCGs$CAT_geneID_merge <- gsub("\\.[0-9]*", "", MoLPS_earlyPCGs$CAT_geneID)
SVSMC_earlyPCGs$EnsID_merge <- gsub("\\.[0-9]*", "", SVSMC_earlyPCGs$EnsID)

SVSMC_earlyPCGs <- filter(SVSMC_earlyPCGs, !grepl("MSTRG|Split", EnsID))

length(unique(SVSMC_earlyPCGs$EnsID_merge))

#coReg: same
SVSMC_earlyPCGs_same <- filter(SVSMC_earlyPCGs, EnsID %in% filter(CoRegPairs_04_48_24_extended_naiveSame, EnsID %in% fpkm_allGDE_Upwithin_4$EnsID, 
                                                                  EnsID.y %in% fpkm_allGDE_Upwithin_4$EnsID)$EnsID.y)

epiMT_earlyPCGs_same <- filter(epiMT_earlyPCGs, CAT_geneID %in% filter(CoRegPairs_sameTimeframe_epiMT,
                                                                       AnchorGene %in% epiMT_DEGs_DE_Upwithin_2$CAT_geneID,
                                                                       geneID %in% epiMT_DEGs_DE_Upwithin_2$CAT_geneID)$geneID)

MoLPS_earlyPCGs_same <- filter(MoLPS_earlyPCGs, CAT_geneID %in% filter(CoRegPairs_sameTimeframe_moLPS,
                                                                       AnchorGene %in% moLPS_DEGs_DE_Upwithin_2$CAT_geneID,
                                                                       geneID %in% moLPS_DEGs_DE_Upwithin_2$CAT_geneID)$geneID)

#venn:
library(ggVennDiagram)

venn_input <- list("SVSMC" = SVSMC_earlyPCGs_same$EnsID_merge, 
                   "EpiMT" = epiMT_earlyPCGs_same$CAT_geneID_merge, 
                   "MoLPS" = MoLPS_earlyPCGs_same$CAT_geneID_merge)

ggVennDiagram(venn_input, label_alpha = 0.5, label = "count", edge_size = 0) +
  ggplot2::theme_void()

sapply(venn_input, length)

#v similar picture really


#
#### fairer version - early phase co-reg lnc neighbours ####

epiMT_earlyPCGs <- rbind(filter(epiMT_DEGs_DE_Upwithin_2, CAT_geneClassII == "coding_mRNA")[,c(1,2)],
                         filter(epiMT_DEGs_DE_Upwithin_4, CAT_geneClassII == "coding_mRNA")[,c(1,2)])

MoLPS_earlyPCGs <- rbind(filter(moLPS_DEGs_DE_Upwithin_2, CAT_geneClassII == "coding_mRNA")[,c(1,2)],
                         filter(moLPS_DEGs_DE_Upwithin_4, CAT_geneClassII == "coding_mRNA")[,c(1,2)])

SVSMC_earlyPCGs <- filter(fpkm_allGDE_Upwithin_4, grepl("coding|CC|TF", GeneClassUpdate))[,c(1,2)]

epiMT_earlyPCGs$CAT_geneID_merge <- gsub("\\.[0-9]*", "", epiMT_earlyPCGs$CAT_geneID)
MoLPS_earlyPCGs$CAT_geneID_merge <- gsub("\\.[0-9]*", "", MoLPS_earlyPCGs$CAT_geneID)
SVSMC_earlyPCGs$EnsID_merge <- gsub("\\.[0-9]*", "", SVSMC_earlyPCGs$EnsID)

SVSMC_earlyPCGs <- filter(SVSMC_earlyPCGs, !grepl("MSTRG|Split", EnsID))

length(unique(SVSMC_earlyPCGs$EnsID_merge))

#coReg: same
SVSMC_earlyPCGs_same <- filter(SVSMC_earlyPCGs, EnsID %in% filter(CoRegPairs_04_48_24_extended_naiveSame, EnsID %in% fpkm_allGDE_Upwithin_4$EnsID, 
                                                                  EnsID.y %in% fpkm_allGDE_Upwithin_4$EnsID)$EnsID.y)

epiMT_earlyPCGs_same <- filter(epiMT_earlyPCGs, CAT_geneID %in% filter(CoRegPairs_sameTimeframe_epiMT,
                                                                       AnchorGene %in% epiMT_DEGs_DE_Upwithin_2$CAT_geneID,
                                                                       geneID %in% epiMT_DEGs_DE_Upwithin_2$CAT_geneID)$geneID |
                                 CAT_geneID %in% filter(CoRegPairs_sameTimeframe_epiMT,
                                                        AnchorGene %in% epiMT_DEGs_DE_Upwithin_4$CAT_geneID,
                                                        geneID %in% epiMT_DEGs_DE_Upwithin_4$CAT_geneID)$geneID)

MoLPS_earlyPCGs_same <- filter(MoLPS_earlyPCGs, CAT_geneID %in% filter(CoRegPairs_sameTimeframe_moLPS,
                                                                       AnchorGene %in% moLPS_DEGs_DE_Upwithin_2$CAT_geneID,
                                                                       geneID %in% moLPS_DEGs_DE_Upwithin_2$CAT_geneID)$geneID |
                                 CAT_geneID %in% filter(CoRegPairs_sameTimeframe_moLPS,
                                                        AnchorGene %in% moLPS_DEGs_DE_Upwithin_4$CAT_geneID,
                                                        geneID %in% moLPS_DEGs_DE_Upwithin_4$CAT_geneID)$geneID)

#venn:
library(ggVennDiagram)

venn_input <- list("SVSMC" = SVSMC_earlyPCGs_same$EnsID_merge, 
                   "EpiMT" = epiMT_earlyPCGs_same$CAT_geneID_merge, 
                   "MoLPS" = MoLPS_earlyPCGs_same$CAT_geneID_merge)

ggVennDiagram(venn_input, label_alpha = 0.5, label = "count", edge_size = 0) +
  ggplot2::theme_void()

sapply(venn_input, length)

#v similar picture really


#
#### plotting - lncRNAs ####

#View(SVSMC_earlyLncs)
venn_input <- list("SVSMC" = SVSMC_earlyLncs$FANTOM_ID_merge,
                   "EpiMT" = epiMT_earlyLncs$CAT_geneID_merge,
                   "MoLPS" = MoLPS_earlyLncs$CAT_geneID_merge)
ggVennDiagram(venn_input, label_alpha = 0.5, label = "count", edge_size = 0) +
  ggplot2::theme_void()

sapply(venn_input, length)

common_earlyLncs <- filter(SVSMC_earlyLncs, FANTOM_ID_merge %in% Reduce(intersect, venn_input))
SVep_earlyLncs <- filter(SVSMC_earlyLncs, FANTOM_ID_merge %in% intersect(venn_input[[1]], venn_input[[2]]))
#note: RP11-221N13.3
SVmo_earlyLncs <- filter(SVSMC_earlyLncs, FANTOM_ID_merge %in% intersect(venn_input[[1]], venn_input[[3]]))
#note: SNHG15

SVSMC_earlyLncs$overlap <- "SVSMC_only"
SVSMC_earlyLncs$overlap[SVSMC_earlyLncs$FANTOM_ID_merge %in% SVep_earlyLncs$FANTOM_ID_merge] <- "epiMT"
SVSMC_earlyLncs$overlap[SVSMC_earlyLncs$FANTOM_ID_merge %in% SVmo_earlyLncs$FANTOM_ID_merge] <- "MoLPS"
SVSMC_earlyLncs$overlap[SVSMC_earlyLncs$FANTOM_ID_merge %in% common_earlyLncs$FANTOM_ID_merge] <- "MoLPS/EpiMT"
table(SVSMC_earlyLncs$overlap)

#View(SVSMC_earlyLncs_same)
venn_input <- list("SVSMC" = unique(SVSMC_earlyLncs_sameUp$FANTOM_ID_merge),
                   "EpiMT" = unique(epiMT_earlyLncs_sameUp$CAT_geneID_merge),
                   "MoLPS" = unique(MoLPS_earlyLncs_sameUp$CAT_geneID_merge))
ggVennDiagram(venn_input, label_alpha = 0.5, label = "count", edge_size = 0) +
  ggplot2::theme_void()

sapply(venn_input, length)

common_earlyLncs_sameUp <- filter(SVSMC_earlyLncs_sameUp, FANTOM_ID_merge %in% Reduce(intersect, venn_input))
SVep_earlyLncs_sameUp   <- filter(SVSMC_earlyLncs_sameUp, FANTOM_ID_merge %in% intersect(venn_input[[1]], venn_input[[2]]))
SVmo_earlyLncs_sameUp   <- filter(SVSMC_earlyLncs_sameUp, FANTOM_ID_merge %in% intersect(venn_input[[1]], venn_input[[3]]))

SVSMC_earlyLncs_sameUp$overlap <- "SVSMC_only"
SVSMC_earlyLncs_sameUp$overlap[SVSMC_earlyLncs_sameUp$FANTOM_ID_merge %in% SVep_earlyLncs_sameUp$FANTOM_ID_merge] <- "epiMT"
SVSMC_earlyLncs_sameUp$overlap[SVSMC_earlyLncs_sameUp$FANTOM_ID_merge %in% SVmo_earlyLncs_sameUp$FANTOM_ID_merge] <- "MoLPS"
SVSMC_earlyLncs_sameUp$overlap[SVSMC_earlyLncs_sameUp$FANTOM_ID_merge %in% common_earlyLncs_sameUp$FANTOM_ID_merge] <- "MoLPS/EpiMT"
table(unique(SVSMC_earlyLncs_sameUp[,c(1,6)])$overlap)

SVSMC_earlyLncs$lncGroup <- "4hr-induced lncRNAs"
SVSMC_earlyLncs_sameUp$lncGroup <- "4hr-induced CClncRNAs"

trial <- rbind(SVSMC_earlyLncs[,5:7],
               SVSMC_earlyLncs_sameUp[,5:7])

ggplot(trial) + aes(x = lncGroup, fill = overlap) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) + Seurat::RotatedAxis()

trial$overlap <- factor(trial$overlap)
trial$overlap <- factor(trial$overlap, levels = levels(trial$overlap)[c(4,2,3,1)])

ggplot(trial) + aes(x = lncGroup, fill = overlap) +
  geom_bar(position = "fill") +
  coord_cartesian(ylim = c(0,0.3)) +
  scale_y_continuous(labels = scales::percent_format()) + Seurat::RotatedAxis()

triali <- as.data.frame(table(trial$overlap, trial$lncGroup))
table(trial$lncGroup)

triali$Freq_perc[1:4] <- triali$Freq[1:4]/23
triali$Freq_perc[5:8] <- triali$Freq[5:8]/69

colnames(triali)[1] <- "Found in"

ggplot(filter(triali, !`Found in` == "SVSMC_only")) + aes(x = Var2, fill = `Found in`, y = Freq_perc) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("%") +
  xlab("\nSVSMC DE lncRNA groups") +
  #coord_cartesian(ylim = c(0,0.3)) +
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0,0.3), breaks = seq(0,0.3,0.1)) + 
  theme_minimal() +
  theme(text = element_text(size=24)) +
  Seurat::RotatedAxis()

triali$Var2 <- factor(triali$Var2)
#triali$Var2 <- factor(triali$Var2, levels = levels(triali$Var2)[c(4,3,2,1)])

triali$`Found in` <- factor(triali$`Found in`, levels = levels(triali$`Found in`)[c(3,4,2,1)])

timecourse_colours <- list("timecourse" = RColorBrewer::brewer.pal(9, "Set1")[2:3])

ggplot(filter(triali, !`Found in` == "SVSMC_only")) + aes(y = Var2, fill = `Found in`, x = Freq_perc) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("%") +
  ylab("") +
  scale_fill_manual(values = c("epiMT" = timecourse_colours$timecourse[1], 
                               "MoLPS" = timecourse_colours$timecourse[2],
                               "MoLPS/EpiMT" = "black")) + 
  #coord_cartesian(ylim = c(0,0.3)) +
  scale_x_continuous(labels = scales::percent_format(),
                     limits = c(0,0.25), breaks = seq(0,0.3,0.1)) + 
  theme_minimal() +
  theme(text = element_text(size=24))


#
#### plotting - PCGs ####

#View(SVSMC_earlyPCGs)
venn_input <- list("SVSMC" = unique(SVSMC_earlyPCGs$EnsID_merge),
                   "EpiMT" = unique(epiMT_earlyPCGs$CAT_geneID_merge),
                   "MoLPS" = unique(MoLPS_earlyPCGs$CAT_geneID_merge))
ggVennDiagram(venn_input, label_alpha = 0.5, label = "count", edge_size = 0) +
  ggplot2::theme_void()

sapply(venn_input, length)

common_earlyPCGs <- filter(SVSMC_earlyPCGs, EnsID_merge %in% Reduce(intersect, venn_input))
SVep_earlyPCGs <- filter(SVSMC_earlyPCGs, EnsID_merge %in% intersect(venn_input[[1]], venn_input[[2]]))
#note: RP11-221N13.3
SVmo_earlyPCGs <- filter(SVSMC_earlyPCGs, EnsID_merge %in% intersect(venn_input[[1]], venn_input[[3]]))
#note: SNHG15

SVSMC_earlyPCGs$overlap <- "SVSMC_only"
SVSMC_earlyPCGs$overlap[SVSMC_earlyPCGs$EnsID_merge %in% SVep_earlyPCGs$EnsID_merge] <- "epiMT"
SVSMC_earlyPCGs$overlap[SVSMC_earlyPCGs$EnsID_merge %in% SVmo_earlyPCGs$EnsID_merge] <- "MoLPS"
SVSMC_earlyPCGs$overlap[SVSMC_earlyPCGs$EnsID_merge %in% common_earlyPCGs$EnsID_merge] <- "MoLPS/EpiMT"
table(SVSMC_earlyPCGs$overlap)

#View(SVSMC_earlyPCGs_same)
venn_input <- list("SVSMC" = SVSMC_earlyPCGs_same$EnsID_merge,
                   "EpiMT" = unique(epiMT_earlyPCGs_same$CAT_geneID_merge),
                   "MoLPS" = unique(MoLPS_earlyPCGs_same$CAT_geneID_merge))
ggVennDiagram(venn_input, label_alpha = 0.5, label = "count", edge_size = 0) +
  ggplot2::theme_void()

sapply(venn_input, length)

common_earlyPCGs <- filter(SVSMC_earlyPCGs_same, EnsID_merge %in% Reduce(intersect, venn_input))
SVep_earlyPCGs <- filter(SVSMC_earlyPCGs_same, EnsID_merge %in% intersect(venn_input[[1]], venn_input[[2]]))
SVmo_earlyPCGs <- filter(SVSMC_earlyPCGs_same, EnsID_merge %in% intersect(venn_input[[1]], venn_input[[3]]))

SVSMC_earlyPCGs_same$overlap <- "SVSMC_only"
SVSMC_earlyPCGs_same$overlap[SVSMC_earlyPCGs_same$EnsID_merge %in% SVep_earlyPCGs$EnsID_merge] <- "epiMT"
SVSMC_earlyPCGs_same$overlap[SVSMC_earlyPCGs_same$EnsID_merge %in% SVmo_earlyPCGs$EnsID_merge] <- "MoLPS"
SVSMC_earlyPCGs_same$overlap[SVSMC_earlyPCGs_same$EnsID_merge %in% common_earlyPCGs$EnsID_merge] <- "MoLPS/EpiMT"
table(SVSMC_earlyPCGs_same$overlap)


SVSMC_earlyPCGs$lncGroup <- "4hr-induced PCGs"
SVSMC_earlyPCGs_same$lncGroup <- "4hr-induced PCGs\nco-induced with\nnearby lncRNA"

trial <- rbind(SVSMC_earlyPCGs[,3:5],
               SVSMC_earlyPCGs_same[,3:5])

ggplot(trial) + aes(x = lncGroup, fill = overlap) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) + Seurat::RotatedAxis()

trial$overlap <- factor(trial$overlap)
trial$overlap <- factor(trial$overlap, levels = levels(trial$overlap)[c(4,2,3,1)])

ggplot(trial) + aes(x = lncGroup, fill = overlap) +
  geom_bar(position = "fill") +
  #coord_cartesian(ylim = c(0,0.3)) +
  scale_y_continuous(labels = scales::percent_format()) + Seurat::RotatedAxis()

triali <- as.data.frame(table(trial$overlap, trial$lncGroup))
table(trial$lncGroup)

triali$Freq_perc[1:4] <- triali$Freq[1:4]/1021
triali$Freq_perc[5:8] <- triali$Freq[5:8]/22

colnames(triali)[1] <- "Found in"

ggplot(filter(triali, !`Found in` == "SVSMC_only")) + aes(x = Var2, fill = `Found in`, y = Freq_perc) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("%") +
  xlab("\nSVSMC DE lncRNA groups") +
  #coord_cartesian(ylim = c(0,0.3)) +
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0,0.6), breaks = seq(0,0.6,0.2)) + 
  theme_minimal() +
  theme(text = element_text(size=24)) +
  Seurat::RotatedAxis()

triali$Var2 <- factor(triali$Var2)
triali$Var2 <- factor(triali$Var2, levels = levels(triali$Var2)[c(4,3,2,1)])

triali$`Found in` <- factor(triali$`Found in`, levels = levels(triali$`Found in`)[c(3,4,2,1)])

timecourse_colours <- list("timecourse" = RColorBrewer::brewer.pal(9, "Set1")[2:3])

ggplot(filter(triali, !`Found in` == "SVSMC_only")) + aes(y = Var2, fill = `Found in`, x = Freq_perc) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("%") +
  ylab("") +
  scale_fill_manual(values = c("epiMT" = timecourse_colours$timecourse[1], 
                               "MoLPS" = timecourse_colours$timecourse[2],
                               "MoLPS/EpiMT" = "black")) + 
  #coord_cartesian(ylim = c(0,0.3)) +
  scale_x_continuous(labels = scales::percent_format(),
                     limits = c(0,0.45), breaks = seq(0,0.6,0.2)) + 
  theme_minimal() +
  theme(text = element_text(size=24))


#### modifications ####

#could add concordant lncs/PCGs