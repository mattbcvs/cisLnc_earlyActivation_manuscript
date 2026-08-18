library(dplyr)
library(ggplot2)

#### import key data tables ####

#filter down (e.g. to just promoter regions or proximal regions of expressed genes) then re-save

fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026filt.csv", header = T)
length(unique(fpkm_allG$EnsID))

fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE_2026filt.csv", header = T)
length(unique(fpkm_allGDE$EnsID))

fpkm_allGDE_Upwithin_4 <- filter(fpkm_allGDE, RegulationStart == "Induced <4hrs")
fpkm_allGDE_Downwithin_4 <- filter(fpkm_allGDE, RegulationStart == "Repressed <4hrs")

fpkm_allGDE_Upwithin_8 <- filter(fpkm_allGDE, RegulationStart == "Induced 4-8hrs")
fpkm_allGDE_Downwithin_8 <- filter(fpkm_allGDE, RegulationStart == "Repressed 4-8hrs")

fpkm_allGDE_Upwithin_24 <- filter(fpkm_allGDE, RegulationStart == "Induced 8-24hrs")
fpkm_allGDE_Downwithin_24 <- filter(fpkm_allGDE, RegulationStart == "Repressed 8-24hrs")

#above lines seperates genes into distinct buckets:
table(fpkm_allGDE$RegulationStart) 


##### TF influence info from ISMARA: ####

#various inputs:
#0-24 hours:
TF_ISMARA <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/ISMARA_allFastQ_MotifZ.csv", header = T)

#0-4 hours:
#TF_ISMARA <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/ISMARA_firstFastQ_MotifZ.csv", header = T)

#all Timepoints/conditions:
#TF_ISMARA <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/ISMARA_allBeyond24FastQ_MotifZ.csv", header = T)

#create version with individual gene symbols
TF_ISMARA_long <- tidyr::separate_longer_delim(TF_ISMARA, cols = "Motif", delim = "_")
TF_ISMARA_long <- TF_ISMARA %>%
  mutate(Symbols = Motif) %>%
  tidyr::separate_rows(Symbols)

#annotate motifs, fpkm, DESeq2 info:
colnames(fpkm_allG)
TF_ISMARA_long <- merge(TF_ISMARA_long, unique(fpkm_allG[,c(2:3,25,27,29,31,33:36,41:44)]), by.x = "Symbols", by.y = "EnsName", all.x = T)

TF_ISMARA_long <- merge(TF_ISMARA_long, fpkm_allGDE[,c(1,45,46)], by = "EnsID", all.x = T)

#exclude miR seeds, anything with UGCA, and nothing else
TF_ISMARA_long$Symbols[grepl("^[GCUA]+$", TF_ISMARA_long$Symbols)]

TF_ISMARA_long$miR_seed <- NA
TF_ISMARA_long$miR_seed[grepl("^[GCUA]+$", TF_ISMARA_long$Symbols)] <- "miR_seed"

TF_ISMARA_long <- filter(TF_ISMARA_long, is.na(miR_seed))

#write.csv(TF_ISMARA_long, "TF_ISMARA_ranking2026.csv", row.names = F)

#
#### validation 1: increased motif influence with expression level? ####
library(ggplot2)

TF_ISMARA_long$fpkm_max_treatment2 <- TF_ISMARA_long$fpkm_max_treatment
TF_ISMARA_long$fpkm_max_treatment2[is.na(TF_ISMARA_long$fpkm_max_treatment2)] <- 0

ggplot(TF_ISMARA_long) + aes(x = fpkm_max_treatment2+1, y = Zscore) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  theme_minimal() +
  xlab("") +
  theme(axis.text.x = element_blank(),
        legend.position = "none")

#somewhat weak but v. sig correlation
cor.test(TF_ISMARA_long$fpkm_max_treatment2, TF_ISMARA_long$Zscore)

#facet by DE time
TF_ISMARA_long$RegulationStart[is.na(TF_ISMARA_long$RegulationStart)] <- "Non-DE"
TF_ISMARA_long$RegulationStart[TF_ISMARA_long$fpkm_max_treatment2 <1] <- "Non-XP"

ggplot(TF_ISMARA_long) + aes(x = fpkm_max_treatment2+1, y = Zscore) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  theme_minimal() +
  facet_wrap(~RegulationStart, scales = "free") +
  xlab("") +
  theme(legend.position = "none")

correlations <- TF_ISMARA_long %>%
  group_by(RegulationStart) %>%
  summarise(correlationP = cor.test(fpkm_max_treatment2, Zscore)$p.value,
            correlationEst =  cor.test(fpkm_max_treatment2, Zscore)$est)
print(correlations)
# no good!

# but then it's quite a messy comparison

#simpler:
library(ggbeeswarm)

#re-order:
#TF_ISMARA_long$RegulationStart2 <- as.factor(TF_ISMARA_long$RegulationStart)
#TF_ISMARA_long$RegulationStart2 <- factor(TF_ISMARA_long$RegulationStart2, levels(TF_ISMARA_long$RegulationStart2)[c(5,4,1:3,6:8)])

ggplot(TF_ISMARA_long) + aes(x = RegulationStart, y = Zscore) +
  geom_quasirandom(alpha =0.2) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  theme_minimal() +
  scale_y_log10() +
  xlab("") +
  theme(legend.position = "none",
        text = element_text(size=24)) + Seurat::RotatedAxis()

table(TF_ISMARA_long$RegulationStart)

kruskal.test(Zscore ~ RegulationStart, TF_ISMARA_long)

library(dunn.test)
dunnISMARA_Z <- dunn.test::dunn.test(TF_ISMARA_long$Zscore, TF_ISMARA_long$RegulationStart,
         method = "bonferroni", 
         kw = TRUE, 
         alt = "greater")

grep("Non-XP", dunnISMARA_Z$comparisons)
dunnISMARA_Z$comparisons
dunnISMARA_Z$altP
dunnISMARA_Z$altP.adjusted

dunnISMARA_Z$comparisons[grep("Non-XP", dunnISMARA_Z$comparisons)]
dunnISMARA_Z$altP[grep("Non-XP", dunnISMARA_Z$comparisons)]
p.adjust(dunnISMARA_Z$altP[grep("Non-XP", dunnISMARA_Z$comparisons)], method = "BH")

p.adjust(dunnISMARA_Z$altP[grep("Non-XP", dunnISMARA_Z$comparisons)], method = "BH")<0.05

#key info 1
#non XP tfs score lowly, building confidence in ISMARA

#key info 2
#tfs in induced clusters do better than those in repressed, these drive more of the overall tx'ome changes


#### top tier TFs ####

#top20% of 606 motifs tested
length(TF_ISMARA$Motif)*0.2 #~122 motifs in top 20% (many unexpressed in remainder tho)
TF_ISMARA[122,] #everything above this

#top 20:
ISMARA_influential20 <- filter(TF_ISMARA_long, fpkm_max_treatment >1 , Zscore >= TF_ISMARA$Zscore[122] )

#top 30 in use for WGCNA:
length(TF_ISMARA$Motif)*0.3 #~122 motifs in top 20% (many unexpressed in remainder tho)
TF_ISMARA[182,] #everything above this

#top 30:
ISMARA_influential30 <- filter(TF_ISMARA_long, fpkm_max_treatment >1, Zscore >= TF_ISMARA$Zscore[182] )

#how many expressed
ISMARA_exp <- filter(TF_ISMARA_long, fpkm_max_treatment >1)


#
#### compare to lists of expected TFs ####

#many unsurprising TF motifs in top including:
#RELA, FOS, MYC, NFKB2, FOXO1, KLF4, SMAD3, TEAD1
#all immune, IEG or SMC maturation/dediff or IL1a/PDGF

#check amongst top 20 etc

#expected TFs driving SMC phenos from lit knowledge:
muscle_TFs <- c("YY1", "KLF4", "SRF", "FOS", "MYOCD", "TET2", "SMAD3", "TCF21", "TEAD1", #my knowledge, TFs involved in promoting SMC mat. or dediff.
                #Miller snATAC paper, heatmap in fig2, motifs enriched in SMC marker genes in control art/athero tissue
                "MEF2A", "MEF2B", "MEF2C", "MEF2D", "TEAD", "TEAD2", "TEAD4", 
                "EBF1", "EBF", "BATF", "FRA1")
sum(muscle_TFs %in% TF_ISMARA_long$Symbols)#14 profiled
sum(muscle_TFs %in% ISMARA_influential20$Symbols)#9 in top 20%
sum(muscle_TFs %in% ISMARA_influential30$Symbols)#10 in top 30% (dim returns)

muscle_TFs[muscle_TFs %in% TF_ISMARA_long$Symbols]
muscle_TFs[muscle_TFs %in% ISMARA_influential20$Symbols]
length(unique(ISMARA_influential20$Symbols)) #182 TFs in top 20% ISMARA

a <- 10
b <- 14
c <- 168
d <- 368

fisher.test(data.frame("muscleTF" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#muscle TFs are enriched in highly influential TFs relative to other expressed TFs
#bit of a subjective list this one, but interesting


#cell cycle
SG2M_TFs <- unique(filter(fpkm_allG, grepl("TF \\+ CC", GeneClassUpdate))$EnsName)
sum(SG2M_TFs %in% TF_ISMARA_long$Symbols)#7 profiled
sum(SG2M_TFs %in% ISMARA_influential20$Symbols)#4 found
sum(SG2M_TFs %in% ISMARA_influential30$Symbols)#4 found (dim returns again)

4/7 #50% of expected for SG2M found in top 20%

SG2M_TFs[SG2M_TFs %in% ISMARA_influential$Motif]

a <- 4
b <- 7
c <- 124
d <- 368

#Sg2M TFs ns enriched in top ISMARA
fisher.test(data.frame("SG2M_TF" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#not very many profiled, but quite a few of those profiled are influential, too few for ns


#IEG
IEGs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Gene lists/arner_2015_table_S5_IEGs_lit.csv")
IEGs_hs <- filter(IEGs, !Hs_symbol %in% c("", "-"))
IEGs_hs <- filter(IEGs_hs, !grepl("[a-z]", IEGs_hs$Total_count))
IEGs_TFs <- IEGs_hs$Hs_symbol[IEGs_hs$Hs_symbol %in% filter(fpkm_allG, grepl("TF", GeneClassUpdate))$EnsName]

sum(IEGs_TFs %in% TF_ISMARA_long$Symbols)#36 profiled
sum(IEGs_TFs %in% ISMARA_influential20$Symbols)#19 in top
sum(IEGs_TFs %in% ISMARA_influential30$Symbols)#23 in top 30 (dim returns)

IEGs_TFs[IEGs_TFs %in% ISMARA_influential$Symbols] #includes NFKB1+2, SRF, MYC, FOS, FOSL1+2, KLF4, KLF6

a <- 23
b <- 36
c <- 168
d <- 368

#IEG TFs enriched in top ISMARA - yes v. much so
fisher.test(data.frame("SG2M_TF" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#n.b. this is maybe a key finding - the SMC activation genes captured are "regulated by IEGs" as well as IEGs themselves
#as "expected TFs" validation, would highlight the IEG finding, and the presence of other TFs with SMC roles in the top 20%

#key info 3
#influential TFs are enriched with key muscle and IEG TFs (top30 or top20 works)


#
#### SCClncRNA correlation to influence ####

#annotate scclnc association
SCClncRNAs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SCClncRNAs.csv")

#start simple:
TF_ISMARA_long_DE <- filter(TF_ISMARA_long, !RegulationStart %in% c("Non-DE", "Non-XP"))
TF_ISMARA_long_DE$SCClncRNA_Target <- "No"
TF_ISMARA_long_DE$SCClncRNA_Target[TF_ISMARA_long_DE$Symbols %in% SCClncRNAs$EnsName.y] <- "Yes"
table(TF_ISMARA_long_DE$SCClncRNA_Target)

ggplot(TF_ISMARA_long_DE) + aes(x = SCClncRNA_Target, y = Zscore) +
  geom_quasirandom(alpha = 0.2) +
  geom_boxplot(width = 0.5, alpha = 0.2, outlier.shape = NA) +
  facet_wrap(~RegulationStart)
#most SCClncRNA targeted TFs are from 4hrs (as established)
#those induced in 4hrs may have a slight increase vs. other TFs induced in 4hrs

TF_ISMARA_long_Up4 <- filter(TF_ISMARA_long_DE, RegulationStart == "Induced <4hrs")

ggplot(TF_ISMARA_long_Up4) + aes(x = SCClncRNA_Target, y = Zscore) +
  geom_quasirandom(alpha = 0.2) +
  geom_boxplot(width = 0.5, alpha = 0.2, outlier.shape = NA)

#not enough n
wilcox.test(Zscore ~ SCClncRNA_Target, TF_ISMARA_long_Up4, )

#but highlights 3 strong TFs in FOXL1, GMEB1 and GMEB2


#key info 3: SCClncRNA targeted TFs induced in 4hrs have decent influence score generally, not sig above median of others but 3x key examples 


#### overall ####

#ISMARA produces expected results and TFs

#highlights a group of SCClncRNA target TFs induced <4hrs, with decent predicted influence
#especially 3 at top
#1 of these FOXL1 is a CAD network TF, promotes SMC prolif
#it's motif looks a little strange
#other 2 functionally unannotated in SMC as yet (GMEB1/GMEB2)


#
#### old code ####




#trial <- filter(TF_ISMARA_long, Motif %in% filter(fpkm_allG, grepl("TF", GeneClassUpdate))$EnsName)

#606 motifs tested
length(TF_ISMARA$Motif)*0.2 #~122 motifs in top 20% (many unexpressed in remainder tho)

#20% thresh if considering only expressed TFs?
unique(TF_ISMARA$Zscore)#no duplicated scores, each score unique
unique(TF_ISMARA_long$Zscore)#284 motifs for expressed TFs
unique(TF_ISMARA_long$Motif)#284 motifs for expressed TFs

#already in order
unique(TF_ISMARA_long$Zscore)[285*0.2]
#this is a more stringent value than that used to select lab candidates (excludes GMEB2)

#Motifs above this score are in the top 20% most influential 
#(this tallies roughly with the online advice of >2 being "SUBSTANTIAL" influence level)
TF_ISMARA$Zscore[122] 

#top 20:
ISMARA_influential <- filter(TF_ISMARA_long, Zscore >= TF_ISMARA$Zscore[122] )

#CClncRNAs:
StrongCoRegPairs_04_48_24_extended <- filter(CoRegPairs_04_48_24_extended, 
                                             (corSig == "Yes" | 
                                                !loopMethod == "Neither" | 
                                                eQTLvalidations >0 | 
                                                FANTOM_eQTL == "Yes" | ExpectedCis == "Yes"
                                             ))

#ISMARA candidates
ISMARA_influential_CisLnc <- filter(ISMARA_influential, Motif %in% StrongCoRegPairs_04_48_24_extended$EnsName.y)
unique(ISMARA_influential_CisLnc$Motif)#4 genes, 5 if including expected cis
#FOXL1, GMEB1/2, NFYB... 
#HOXC8 is padj 0.052 via Spearman's (and HOTAIR expected cis-acting)
#therefore should be in the list?

