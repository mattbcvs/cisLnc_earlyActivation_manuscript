library(dplyr)
library(ggplot2)

#### import key data tables ####

fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGClassUpdate.csv",
                      header = T)
FPKM_CQV_OVERLAP_fpkm <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/FPKM_CQV_OVERLAP_fpkm.csv")
table(FPKM_CQV_OVERLAP_fpkm$IGV)#413 pass, 168 fail

#should have removed artefacts from here already
fpkm_allG_filt <- filter(fpkm_allG, grepl("chr", chr), !grepl("artefacts", GeneClassUpdate), !is.na(GeneClassUpdate))
#remove manual or Thresh4 fails
fpkm_allG_filt_manual <- filter(fpkm_allG_filt, 
                                !EnsID %in% filter(FPKM_CQV_OVERLAP_fpkm, IGV == "fail")$EnsID, #remove manual fails
)
fpkm_allG <- fpkm_allG_filt_manual


fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE_clustered.csv", header = T)
fpkm_allGDE_filt <- filter(fpkm_allGDE, grepl("chr", chr), !grepl("artefacts", GeneClassUpdate), !is.na(GeneClassUpdate))
#remove manual/Thresh4 fails
fpkm_allGDE <- filter(fpkm_allGDE_filt, EnsID %in% fpkm_allG_filt_manual$EnsID)

#CClncRNAs:
CoRegPairs_04_48_24_extended <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CoRegPairs_04_48_24_extended_SCClnc_Nov25.csv")


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

#add in gene symbols + some expression info
fpkm_allG_ID <- unique(fpkm_allG[,c(1,2,4,5,27,29,31,33,35)])

#also needs expanding!
fpkm_allG_ID_long <- fpkm_allG_ID %>%
  mutate(EnsName2 = AllNames) %>%
  tidyr::separate_rows(EnsName2, sep = ", ")

fpkm_allG_ID_long$EnsName2[is.na(fpkm_allG_ID_long$EnsName2)] <- fpkm_allG_ID_long$EnsName[is.na(fpkm_allG_ID_long$EnsName2)]

TF_ISMARA_long <- merge(TF_ISMARA_long, fpkm_allG_ID, by.x = "Symbols", by.y = "EnsName")

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


#### compare to lists of expected TFS ####

#many unsurprising TF motifs in top including:
#RELA, FOS, MYC, NFKB2, FOXO1, KLF4, SMAD3, TEAD1
#all immune or SMC maturation/dediff or IL1a/PDGF

#expected TFs driving SMC phenos from lit knowledge:
muscle_TFs <- c("YY1", "KLF4", "SRF", "FOS", "MYOCD", "TET2", "SMAD3", "TCF21", "TEAD1", #my knowledge, TFs involved in promoting SMC mat. or dediff.
                #Miller snATAC paper, heatmap in fig2, motifs enriched in SMC marker genes in control art/athero tissue
                "MEF2A", "MEF2B", "MEF2C", "MEF2D", "TEAD", "TEAD2", "TEAD4", 
                "EBF1", "EBF", "BATF", "FRA1")
sum(muscle_TFs %in% TF_ISMARA_long$Motif)#9 profiled
sum(muscle_TFs %in% ISMARA_influential$Motif)#5 in top 20%

muscle_TFs[muscle_TFs %in% TF_ISMARA_long$Motif]
muscle_TFs[muscle_TFs %in% ISMARA_influential$Motif]
length(unique(ISMARA_influential$Motif)) #120 TFs in top 20% ISMARA

a <- 6
b <- 10
c <- 120
d <- 350

#note TEAD1 and TEAD3 in there too

#"muscle" TFs ns for enriched in top ISMARA (but if went through each would probs improve)
fisher.test(data.frame("muscleTF" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#bit of a subjective list this one - wouldn't bother


#cell cycle
SG2M_TFs <- unique(filter(fpkm_allG_filt_manual, grepl("TF \\+ CC", GeneClassUpdate))$EnsName)
sum(SG2M_TFs %in% TF_ISMARA_long$Motif)#6 profiled
sum(SG2M_TFs %in% ISMARA_influential$Motif)#3 found

3/6 #50% of expected for SG2M found in top 20%

SG2M_TFs[SG2M_TFs %in% ISMARA_influential$Motif]

a <- 3
b <- 6
c <- 120
d <- 350

#Sg2M TFs ns enriched in top ISMARA
fisher.test(data.frame("SG2M_TF" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")

#not very many profiled, but quite a few of those profiled are influential, too few for ns


#IEG
IEGs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Gene lists/arner_2015_table_S5_IEGs_lit.csv")
IEGs_hs <- filter(IEGs, !Hs_symbol %in% c("", "-"))
IEGs_hs <- filter(IEGs_hs, !grepl("[a-z]", IEGs_hs$Total_count))
IEGs_TFs <- IEGs_hs$Hs_symbol[IEGs_hs$Hs_symbol %in% c(filter(fpkm_allG, grepl("TF", GeneClassUpdate))$EnsName
                                                       , unlist(strsplit(filter(fpkm_allG, grepl("TF", GeneClassUpdate), !is.na(AllNames))$AllNames, ", "))
                                                       )]


sum(IEGs_TFs %in% TF_ISMARA_long$Motif)#35 profiled
sum(IEGs_TFs %in% ISMARA_influential$Motif)#19 in top
19/35 #54% of IEG TFs in top 20
IEGs_TFs[IEGs_TFs %in% ISMARA_influential$Motif] #includes NFKB1+2, SRF, MYC, FOS, FOSL1+2, 

a <- 19
b <- 35
c <- 120
d <- 350

#IEG TFs enriched in top ISMARA - yes v. much so
fisher.test(data.frame("SG2M_TF" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")

#n.b. this is maybe a key finding - the SMC activation genes captured are "regulated by IEGs" as well as IEGs themselves

#as "expected TFs" validation, would highlight the IEG finding, and the presence of other TFs with SMC roles in the top 20%


#### more validation - expressed/non-expressed TFs ####

allEnsNames <- unique(c(fpkm_allG$EnsName
                 , unlist(strsplit(filter(fpkm_allG, !is.na(AllNames))$AllNames, ", "))
                 ))

#grep amongst motifs for expressed TFs:
trial <- lapply(allEnsNames, function(x){
  return(filter(TF_ISMARA, grepl(x, Motif)))
})

expressedTF_MARA <- unique(bind_rows(trial))

TF_ISMARA$Expressed <- "Not Expressed(<1FPKM)"
TF_ISMARA$Expressed[TF_ISMARA$Motif %in% expressedTF_MARA$Motif] <- "Expressed(>1FPKM)"

table(TF_ISMARA$Expressed)

ggplot(TF_ISMARA) + aes(x = Expressed, y = Zscore, color = Expressed) +
  geom_violin() +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  #geom_jitter(alpha = 0.2) +
  theme_minimal() +
  scale_y_log10() +
  xlab("") +
  theme(axis.text.x = element_blank(),
        legend.position = "none")

# **** association of TF expression with motif influence via 2 sample t
t.test(filter(TF_ISMARA, Expressed == "Not Expressed(<1FPKM)")$Zscore, 
       filter(TF_ISMARA, !Expressed == "Not Expressed(<1FPKM)")$Zscore, var.equal = F)

ggplot() + 
  stat_qq(aes(sample = filter(TF_ISMARA, Expressed == "Not Expressed(<1FPKM)")$Zscore), colour = "green") + 
  stat_qq(aes(sample = filter(TF_ISMARA, !Expressed == "Not Expressed(<1FPKM)")$Zscore), colour = "red") +
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2)

#more groups:
trial <- filter(fpkm_allG, fpkm_max_treatment>5)
tier2EnsNames <- unique(c(trial$EnsName
                        , unlist(strsplit(filter(trial, !is.na(AllNames))$AllNames, ", "))
))

#grep amongst motifs for expressed TFs:
trial <- lapply(tier2EnsNames, function(x){
  return(filter(TF_ISMARA, grepl(x, Motif)))
})

t2TF_MARA <- unique(bind_rows(trial))

TF_ISMARA$Expressed[TF_ISMARA$Motif %in% t2TF_MARA$Motif] <- "Expressed(>5FPKM)"

trial <- filter(fpkm_allG, fpkm_max_treatment>10)
tier3EnsNames <- unique(c(trial$EnsName
                          , unlist(strsplit(filter(trial, !is.na(AllNames))$AllNames, ", "))
))

#grep amongst motifs for expressed TFs:
trial <- lapply(tier3EnsNames, function(x){
  return(filter(TF_ISMARA, grepl(x, Motif)))
})

t3TF_MARA <- unique(bind_rows(trial))

TF_ISMARA$Expressed[TF_ISMARA$Motif %in% t3TF_MARA$Motif] <- "Expressed(>10FPKM)"

ggplot(TF_ISMARA) + aes(x = Expressed, y = Zscore, color = Expressed) +
  geom_violin() +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  #geom_jitter(alpha = 0.2) +
  theme_minimal() +
  scale_y_log10() +
  xlab("") +
  theme(axis.text.x = element_blank()#,legend.position = "none"
        )

aov_res <- aov(Zscore ~ Expressed, TF_ISMARA)
summary(aov_res)#returns corrected and non-corrected forms
TukeyHSD(aov_res)

#suggests a correlation of motif score and max FPKM would be good validation

trial <- lapply(allEnsNames, function(x){
  if(dim(filter(TF_ISMARA, grepl(x, Motif)))[1] >0 ) {
  return(data.frame(filter(TF_ISMARA, grepl(x, Motif)),
                    filter(unique(fpkm_allG[,c(4,5,35)]), grepl(x, EnsName) | grepl(x, AllNames))
                    )
         )
  }
  })
