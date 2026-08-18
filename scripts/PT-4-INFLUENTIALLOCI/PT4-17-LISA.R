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


#### export test genes for LISA ####

LISA_Up <- filter(fpkm_allGDE_Upwithin_4, grepl("ENS", EnsID))[order(fpkm_allGDE_Upwithin_4$preadj_0_4, decreasing = F),]#[1:500,]
LISA_Up <- LISA_Up[order(LISA_Up$preadj_0_4, decreasing = F),][1:500,]
#write.csv(gsub("\\.[0-9]*", "", LISA_Up$EnsID), "LISA_Up.csv", row.names = F)

LISA_Down <- filter(fpkm_allGDE_Downwithin_4, grepl("ENS", EnsID))[order(fpkm_allGDE_Downwithin_4$preadj_0_4, decreasing = F),]#[1:500,]
LISA_Down <- LISA_Down[order(LISA_Down$preadj_0_4, decreasing = F),][1:500,]
#write.csv(gsub("\\.[0-9]*", "", LISA_Down$EnsID), "LISA_Down.csv", row.names = F)

LISA_Up8 <- filter(fpkm_allGDE_Upwithin_8, grepl("ENS", EnsID))[order(fpkm_allGDE_Upwithin_4$preadj_0_8, decreasing = F),]#[1:500,]
LISA_Up8 <- LISA_Up8[order(LISA_Up8$preadj_0_8, decreasing = F),][1:500,]
#write.csv(gsub("\\.[0-9]*", "", LISA_Up8$EnsID), "LISA_Up8.csv", row.names = F)

LISA_Down8 <- filter(fpkm_allGDE_Downwithin_8, grepl("ENS", EnsID))[order(fpkm_allGDE_Downwithin_8$preadj_0_8, decreasing = F),]#[1:500,]
LISA_Down8 <- LISA_Down8[order(LISA_Down8$preadj_0_8, decreasing = F),][1:500,]
#write.csv(gsub("\\.[0-9]*", "", LISA_Down8$EnsID), "LISA_Down8.csv", row.names = F)

LISA_Up24 <- filter(fpkm_allGDE_Upwithin_24, grepl("ENS", EnsID))[order(fpkm_allGDE_Upwithin_24$preadj_0_24, decreasing = F),]#[1:500,]
LISA_Up24 <- LISA_Up24[order(LISA_Up24$preadj_0_24, decreasing = F),][1:500,]
#write.csv(gsub("\\.[0-9]*", "", LISA_Up24$EnsID), "LISA_Up24.csv", row.names = F)

LISA_Down24 <- filter(fpkm_allGDE_Downwithin_24, grepl("ENS", EnsID))[order(fpkm_allGDE_Downwithin_24$preadj_0_24, decreasing = F),]#[1:500,]
LISA_Down24 <- LISA_Down24[order(LISA_Down24$preadj_0_24, decreasing = F),][1:500,]
#write.csv(gsub("\\.[0-9]*", "", LISA_Down24$EnsID), "LISA_Down24.csv", row.names = F)


#
#### output from LISA ####

Lisa_0_4up_rankTFs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/LISA_2026_4hrUp.csv")
Lisa_0_4up_rankTFs$bestP <-  as.numeric(sapply(strsplit(Lisa_0_4up_rankTFs$X1st.Sample.p.value, ";"), "[[", 2))

#0-4down
Lisa_0_4down_rankTFs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/LISA_2026_4hrDown.csv")
Lisa_0_4down_rankTFs$bestP <-  as.numeric(sapply(strsplit(Lisa_0_4down_rankTFs$X1st.Sample.p.value, ";"), "[[", 2))

#0-8up
Lisa_0_8up_rankTFs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/LISA_2026_8hrUp.csv")
Lisa_0_8up_rankTFs$bestP <-  as.numeric(sapply(strsplit(Lisa_0_8up_rankTFs$X1st.Sample.p.value, ";"), "[[", 2))

#0-8down
Lisa_0_8down_rankTFs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/LISA_2026_8hrDown.csv")
Lisa_0_8down_rankTFs$bestP <-  as.numeric(sapply(strsplit(Lisa_0_8down_rankTFs$X1st.Sample.p.value, ";"), "[[", 2))

#0-24up
Lisa_0_24up_rankTFs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/LISA_2026_24hrUp.csv")
Lisa_0_24up_rankTFs$bestP <-  as.numeric(sapply(strsplit(Lisa_0_24up_rankTFs$X1st.Sample.p.value, ";"), "[[", 2))

#0-24down
Lisa_0_24down_rankTFs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/LISA_2026_24hrDown.csv")
Lisa_0_24down_rankTFs$bestP <-  as.numeric(sapply(strsplit(Lisa_0_24down_rankTFs$X1st.Sample.p.value, ";"), "[[", 2))

#LISA annotation
Lisa_annotation <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/LISA_gene_classes.csv")

#combining all rankings, remove low significance first:
trial <- merge(Lisa_0_4up_rankTFs[,c(1,7)], Lisa_0_4down_rankTFs[,c(1,7)], by = "Transcription.Factor", all = T)
trial <- merge(trial, Lisa_0_8up_rankTFs[,c(1,7)], by = "Transcription.Factor", all = T)
trial <- merge(trial, Lisa_0_8down_rankTFs[,c(1,7)], by = "Transcription.Factor", all = T)
trial <- merge(trial, Lisa_0_24up_rankTFs[,c(1,7)], by = "Transcription.Factor", all = T)
trial <- merge(trial, Lisa_0_24down_rankTFs[,c(1,7)], by = "Transcription.Factor", all = T)
colnames(trial)[2:7] <- c("upP04", "downP04", "upP08", "downP08", "upP024", "downP024")

LISA_testInfluence <- trial
LISA_testInfluence$lowestP <- Biobase::rowMin(as.matrix(LISA_testInfluence[,2:7]))
LISA_testInfluence <- LISA_testInfluence[order(LISA_testInfluence$upP04, decreasing = F),]
rownames(LISA_testInfluence) <- NULL

#considering removing non TFs which are bit weird (hotair, il1b...)
trial <- merge(LISA_testInfluence, Lisa_annotation, by.x = "Transcription.Factor", by.y = "factor", all.x = T)
trial$factor_type[grepl("HBG1 ", trial$Transcription.Factor)] <- "other"
trial$factor_type[grepl("TEAD1 ", trial$Transcription.Factor)] <- "tf"
table(trial$factor_type)

trial <- filter(trial, factor_type %in% c("cr", "tf", "predicted chromatin regulator", "predicted transcription factor", "both predicted transcription factor and chromatin regulator"))
LISA_testInfluence <- trial

#add in expression info
LISA_testInfluence <- merge(LISA_testInfluence, unique(fpkm_allG[,c(2:3,25,27,29,31,33:36,41:44)]), by.x = "Transcription.Factor", by.y = "EnsName", all.x = T)

#add in regulation state
LISA_testInfluence <- merge(LISA_testInfluence, fpkm_allGDE[,c(2,45,46)], by.x = "Transcription.Factor", by.y = "EnsName", all.x = T)

#write.csv(LISA_testInfluence, "LISA_testInfluence2026.csv", row.names = F)


#
#### validation 1: increased motif influence with expression level? ####
library(ggplot2)

LISA_testInfluence$fpkm_max_treatment2 <- LISA_testInfluence$fpkm_max_treatment
LISA_testInfluence$fpkm_max_treatment2[is.na(LISA_testInfluence$fpkm_max_treatment2)] <- 0

ggplot(LISA_testInfluence) + aes(x = fpkm_max_treatment2+1, y = -log10(lowestP)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  theme_minimal() +
  xlab("") +
  theme(axis.text.x = element_blank(),
        legend.position = "none")

#ns, above plot misleading...
cor.test(LISA_testInfluence$fpkm_max_treatment2, -log10(LISA_testInfluence$lowestP))

#facet by DE time
LISA_testInfluence$RegulationStart[is.na(LISA_testInfluence$RegulationStart)] <- "Non-DE"
LISA_testInfluence$RegulationStart[LISA_testInfluence$fpkm_max_treatment2 <1] <- "Non-XP"

ggplot(LISA_testInfluence) + aes(x = fpkm_max_treatment2+1, y = -log10(lowestP)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10() +
  theme_minimal() +
  facet_wrap(~RegulationStart, scales = "free") +
  xlab("") +
  theme(legend.position = "none")

correlations <- LISA_testInfluence %>%
  group_by(RegulationStart) %>%
  summarise(correlationP = cor.test(fpkm_max_treatment2, -log10(lowestP))$p.value,
            correlationEst =  cor.test(fpkm_max_treatment2, -log10(lowestP))$est)
print(correlations)
# no good!

# but then it's quite a messy comparison

#simpler:
library(ggbeeswarm)

ggplot(LISA_testInfluence) + aes(x = RegulationStart, y = -log10(lowestP)) +
  geom_quasirandom(alpha =0.2) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha =0.2) +
  theme_minimal() +
  scale_y_log10() +
  xlab("") +
  theme(legend.position = "none",
        text = element_text(size=24)) + Seurat::RotatedAxis()

table(LISA_testInfluence$RegulationStart)

kruskal.test(lowestP ~ RegulationStart, LISA_testInfluence)

library(dunn.test)
dunnLISA_p <- dunn.test::dunn.test(LISA_testInfluence$lowestP, LISA_testInfluence$RegulationStart,
                                     method = "bonferroni", 
                                     kw = TRUE, 
                                     alt = "greater")

grep("Non-XP", dunnLISA_p$comparisons)
dunnLISA_p$comparisons
dunnLISA_p$altP
dunnLISA_p$altP.adjusted

dunnLISA_p$comparisons[grep("Non-XP", dunnLISA_p$comparisons)]
dunnLISA_p$altP[grep("Non-XP", dunnLISA_p$comparisons)]
p.adjust(dunnLISA_p$altP[grep("Non-XP", dunnLISA_p$comparisons)], method = "BH")

p.adjust(dunnLISA_p$altP[grep("Non-XP", dunnLISA_p$comparisons)], method = "BH")<0.05

#key info 1
#non XP tfs score lowly, building confidence in ISMARA

#key info 2
#tfs in non-early induced clusters do better than those in repressed, these drive more of the overall tx'ome changes


#### top tier TFs ####

#top20% of 606 motifs tested
LISA_testInfluence <- LISA_testInfluence[order(LISA_testInfluence$lowestP, decreasing = F),]
length(LISA_testInfluence$lowestP)*0.2 #~122 motifs in top 20% (many unexpressed in remainder tho)
LISA_testInfluence[179,] #everything above this

#top 20:
LISA_influential20 <- filter(LISA_testInfluence, fpkm_max_treatment >1 , lowestP <= LISA_testInfluence$lowestP[179] )

#top 30 in use for WGCNA:
length(LISA_testInfluence$Transcription.Factor)*0.3 #~122 motifs in top 20% (many unexpressed in remainder tho)
LISA_testInfluence[269,] #everything above this

#top 30:
LISA_influential30 <- filter(LISA_testInfluence, fpkm_max_treatment >1, lowestP <= LISA_testInfluence$lowestP[269] )

#how many expressed
LISA_exp <- filter(LISA_testInfluence, fpkm_max_treatment >1)


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
sum(muscle_TFs %in% LISA_testInfluence$Transcription.Factor)#16 profiled
sum(muscle_TFs %in% LISA_influential20$Transcription.Factor)#6 in top 20%
sum(muscle_TFs %in% LISA_influential30$Transcription.Factor)#9 in top 30%

muscle_TFs[muscle_TFs %in% LISA_testInfluence$Transcription.Factor]
muscle_TFs[muscle_TFs %in% LISA_influential30$Transcription.Factor]
length(unique(LISA_influential30$Transcription.Factor)) 

a <- 9
b <- 16
c <- 233
d <- 665

fisher.test(data.frame("muscleTF" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#muscle TFs are enriched in highly influential TFs relative to other expressed TFs
#bit of a subjective list this one, but interesting


#cell cycle
SG2M_TFs <- unique(filter(fpkm_allG, grepl("TF \\+ CC", GeneClassUpdate))$EnsName)
sum(SG2M_TFs %in% LISA_testInfluence$Transcription.Factor)#13 profiled
sum(SG2M_TFs %in% LISA_influential20$Transcription.Factor)#5 found
sum(SG2M_TFs %in% LISA_influential30$Transcription.Factor)#5 found (dim returns again)

a <- 5
b <- 13
c <- 233
d <- 665

#Sg2M TFs ns enriched in top ISMARA
fisher.test(data.frame("SG2M_TF" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#not very many profiled, but quite a few of those profiled are influential, too few for ns


#IEG
IEGs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Gene lists/arner_2015_table_S5_IEGs_lit.csv")
IEGs_hs <- filter(IEGs, !Hs_symbol %in% c("", "-"))
IEGs_hs <- filter(IEGs_hs, !grepl("[a-z]", IEGs_hs$Total_count))
IEGs_TFs <- IEGs_hs$Hs_symbol[IEGs_hs$Hs_symbol %in% filter(fpkm_allG, grepl("TF", GeneClassUpdate))$EnsName]

sum(IEGs_TFs %in% LISA_testInfluence$Transcription.Factor)#33 profiled
sum(IEGs_TFs %in% LISA_influential20$Transcription.Factor)#11 in top
sum(IEGs_TFs %in% LISA_influential30$Transcription.Factor)#19 in top 30

IEGs_TFs[IEGs_TFs %in% LISA_influential30$Transcription.Factor] #includes NFKB1+2, SRF, MYC, FOS, FOSL1+2, KLF4, KLF6

a <- 19
b <- 33
c <- 233
d <- 665

#IEG TFs enriched in top ISMARA - yes v. much so
fisher.test(data.frame("SG2M_TF" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#n.b. this is maybe a key finding - the SMC activation genes captured are "regulated by IEGs" as well as IEGs themselves
#as "expected TFs" validation, would highlight the IEG finding, and the presence of other TFs with SMC roles in the top 20%

#key info 3
#influential TFs are enriched with IEG TFs (top30 or top20 works)


#
#### SCClncRNA correlation to influence ####

#annotate scclnc association
SCClncRNAs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SCClncRNAs.csv")

#start simple:
TF_LISA_influence_DE <- filter(LISA_testInfluence, !RegulationStart %in% c("Non-DE", "Non-XP"))
TF_LISA_influence_DE$SCClncRNA_Target <- "No"
TF_LISA_influence_DE$SCClncRNA_Target[TF_LISA_influence_DE$Transcription.Factor %in% SCClncRNAs$EnsName.y] <- "Yes"
table(TF_LISA_influence_DE$SCClncRNA_Target)

ggplot(TF_LISA_influence_DE) + aes(x = SCClncRNA_Target, y = -log10(lowestP)) +
  geom_quasirandom(alpha = 0.2) +
  geom_boxplot(width = 0.5, alpha = 0.2, outlier.shape = NA) +
  facet_wrap(~RegulationStart)
#most SCClncRNA targeted TFs are from 4hrs (as established)
#those induced in 4hrs may have a slight increase vs. other TFs induced in 4hrs

TF_LISA_influence_Up4 <- filter(TF_LISA_influence_DE, RegulationStart == "Induced <4hrs")

ggplot(TF_LISA_influence_Up4) + aes(x = SCClncRNA_Target, y = -log10(lowestP)) +
  geom_quasirandom(alpha = 0.2) +
  geom_boxplot(width = 0.5, alpha = 0.2, outlier.shape = NA)

#not enough n for sig
wilcox.test(Zscore ~ SCClncRNA_Target, TF_ISMARA_long_Up4, )

#but highlights 3 strong TFs in GMEB1 and GMEB2 (again) and TET3

#or alternatively, induced lnc targets:
SCClncRNAs_4hrUp <- filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs")

TF_LISA_influence_DE$SCClncRNA_Target2 <- "No"
TF_LISA_influence_DE$SCClncRNA_Target2[TF_LISA_influence_DE$Transcription.Factor %in% SCClncRNAs_4hrUp$EnsName.y] <- "Yes"
table(TF_LISA_influence_DE$SCClncRNA_Target2)

ggplot(TF_LISA_influence_DE) + aes(x = SCClncRNA_Target2, y = -log10(lowestP)) +
  geom_quasirandom(alpha = 0.2) +
  geom_boxplot(width = 0.5, alpha = 0.2, outlier.shape = NA) +facet_wrap(~RegulationStart)


#key info 3: SCClncRNA targeted TFs induced in 4hrs have decent influence score generally, not sig above median of others but 3x key examples, 2 of which
#found in ISMARA too



