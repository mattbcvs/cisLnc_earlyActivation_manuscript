#nuclear bias for early lncs?
library(dplyr)
library(ggplot2)

fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026filt.csv", header = T)

fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE_2026filt.csv", header = T)
table(fpkm_allGDE$RegulationStart) 

#510 expressed lncRNAs submitted to lncATLAS browser:
unique(filter(fpkm_allG, grepl("ENS", EnsID), 
              grepl("fide|Lnc", GeneClassUpdate))$EnsID)

#write.csv(unique(filter(fpkm_allG, grepl("ENS", EnsID), 
#                        grepl("fide|Lnc", GeneClassUpdate))$EnsID), "GEN_LNC_2026.csv", row.names = F)

#FPKM for all lncs:
lncATLAS_fpkm_allLncs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/2025-12-18 lncATLAS.csv")

#quick explore, CNRCI for x lncRNAs:
length(unique(filter(lncATLAS_fpkm_allLncs, Data.Type == "CNRCI")$ENSEMBL.ID))#501 genes have a nuc/cyto bias returned
length(unique(filter(lncATLAS_fpkm_allLncs, Data.Type == "cell")$ENSEMBL.ID))#501 with a cell fpkm
length(unique(filter(lncATLAS_fpkm_allLncs, Data.Type == "nucleus")$ENSEMBL.ID))#501 with a nucleus fpkm
length(unique(filter(lncATLAS_fpkm_allLncs, Data.Type == "cytoplasm")$ENSEMBL.ID))#501 with a cytoplasm fpkm

#some are <1 FPKM in both nucleus and cytoplasm in a cell, remove these as not robust?

#split by lnc and cell type
lncATLAS_fpkm_allLncs$Lnc_Cell <- paste(lncATLAS_fpkm_allLncs$ENSEMBL.ID, lncATLAS_fpkm_allLncs$Data.Source, sep = "_")

trial <- split(lncATLAS_fpkm_allLncs, lncATLAS_fpkm_allLncs$Lnc_Cell)
trial <- lapply(trial, function(x){
  filter(x, (Data.Type == "nucleus" & Value >1) | (Data.Type == "cytoplasm" & Value >1) )
})
trial <- bind_rows(trial)

#about half are probably being judged without v. strong signals
length(unique(trial$Lnc_Cell))
length(unique(lncATLAS_fpkm_allLncs$Lnc_Cell))

strong_lnc_cells <- unique(trial$Lnc_Cell)

#look at values for all, as well as a subset
lncATLAS_fpkm_allLncs_strong <-filter(lncATLAS_fpkm_allLncs, Lnc_Cell %in% strong_lnc_cells)
  


#nuc bias, CNRCI, select input here:
nucBias_all <- filter(lncATLAS_fpkm_allLncs_strong, Data.Type == "CNRCI")
#nucBias_all <- filter(lncATLAS_fpkm_allLncs, Data.Type == "CNRCI")

#calculate mean across all cells
trial <- split(nucBias_all, nucBias_all$ENSEMBL.ID)

trialiii <- sapply(trial, function(x){
  mean(x$Value, na.rm = T)
})

nucBias_CNRCI_summary <- data.frame("Ens_ID_merge" = names(trial),
                                    "MeanCNRCI" = trialiii)

#annotate with reg. state
fpkm_allGDE$Ens_ID_merge <- gsub("\\.[0-9]*", "", fpkm_allGDE$EnsID)
nucBias_CNRCI_summary_ <- unique(merge(nucBias_CNRCI_summary, fpkm_allGDE[,c(46,47)], by = "Ens_ID_merge", all.x = T))
nucBias_CNRCI_summary_$RegulationStart[is.na(nucBias_CNRCI_summary_$RegulationStart)] <- "Non-DE"

####not doing anymore: add in expressed lncs seperately as a background for the plot
#nucBias_CNRCI_summary$RegulationStart <- "Expressed "
#nucBias_CNRCI_summary <- rbind(nucBias_CNRCI_summary_, nucBias_CNRCI_summary)

#order plot
#nucBias_CNRCI_summary$RegulationStart <- as.factor(nucBias_CNRCI_summary$RegulationStart)
#nucBias_CNRCI_summary$RegulationStart <- factor(nucBias_CNRCI_summary$RegulationStart,
#                                                  levels(nucBias_CNRCI_summary$RegulationStart)[c(1,2:7)])

ggplot(nucBias_CNRCI_summary_) + aes(x = RegulationStart, y = MeanCNRCI, fill = RegulationStart) +
  #geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  xlab("") +
  ylab("Bias score (MeanCNRCI)") +
  geom_jitter(width = 0.1, alpha = 0.4, size =2) + 
  scale_fill_manual(values = c(`Expressed` = "grey60",
                               `Induced <4hrs` = "#D6604D",
                               `Induced 4-8hrs` = "#D6604D",
                               `Induced 8-24hrs` = "#D6604D",
                               `Repressed <4hrs` = "#67A9CF",
                               `Repressed 4-8hrs` = "#67A9CF",
                               `Repressed 8-24hrs` = "#67A9CF")) +
  theme_minimal() + 
  theme(text = element_text(size = 24), legend.position = "none") +
  Seurat::RotatedAxis()

#seems a visible trend
table(nucBias_CNRCI_summary_$RegulationStart)

#anova + dunnett's
library(afex)
#library(emmeans)
aov_res <- aov_ez(id = "Ens_ID_merge",
                  dv = "MeanCNRCI",
                  between = "RegulationStart",
                  data = nucBias_CNRCI_summary_)
summary(aov_res)#ANOVA is ns

#not significant, despite some visually interesting pattern, especially with 

#try with better FPKM only, still ns

#try with better FPKM and less groups:

#can still condense further:
nucBias_CNRCI_summary_$RegulationStartII <- as.character(nucBias_CNRCI_summary_$RegulationStart)
nucBias_CNRCI_summary_$RegulationStartII[grepl("Induced 4|Induced 8", nucBias_CNRCI_summary_$RegulationStartII)] <- "Induced >4hrs"
nucBias_CNRCI_summary_$RegulationStartII[grepl("Repressed 4|Repressed 8", nucBias_CNRCI_summary_$RegulationStartII)] <- "Repressed >4hrs"

ggplot(nucBias_CNRCI_summary_) + aes(x = RegulationStartII, y = MeanCNRCI, fill = RegulationStartII) +
  #geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  xlab("") +
  ylab("Bias score (MeanCNRCI)") +
  geom_jitter(width = 0.1, alpha = 0.4, size =2) + 
  scale_fill_manual(values = c(`Expressed` = "grey60",
                               `Induced <4hrs` = "#D6604D",
                               `Induced >4hrs` = "#D6604D",
                               `Repressed <4hrs` = "#67A9CF",
                               `Repressed >4hrs` = "#67A9CF")) +
  theme_minimal() + 
  theme(text = element_text(size = 24), legend.position = "none") +
  Seurat::RotatedAxis()

table(nucBias_CNRCI_summary_$RegulationStartII)
47+44+33+58#182 lncRNAs being checked if doing all
35+36+33+52#156 if including only those with fpkm >1 in nuc or cyto

#anova + dunnett's
library(afex)
#library(emmeans)
aov_res <- aov_ez(id = "Ens_ID_merge",
                  dv = "MeanCNRCI",
                  between = "RegulationStartII",
                  data = nucBias_CNRCI_summary_)
summary(aov_res)#ANOVA is * (ns if including all)

#Dunnett's no good
nucBias_CNRCI_summary_$RegulationStartII <- as.factor(nucBias_CNRCI_summary_$RegulationStartII)
nucBias_CNRCI_summary_$RegulationStartII <- relevel(nucBias_CNRCI_summary_$RegulationStartII, "Non-DE")
DescTools::DunnettTest(nucBias_CNRCI_summary_$MeanCNRCI, nucBias_CNRCI_summary_$RegulationStartII)

#simpler, early vs. later, borderline
t.test(filter(nucBias_CNRCI_summary_, RegulationStartII == "Induced <4hrs")$MeanCNRCI,
       filter(nucBias_CNRCI_summary_, RegulationStartII == "Induced >4hrs")$MeanCNRCI, var.equal = T)

#simpler, early vs. later, borderline
t.test(filter(nucBias_CNRCI_summary_, RegulationStartII == "Repressed <4hrs")$MeanCNRCI,
       filter(nucBias_CNRCI_summary_, RegulationStartII == "Repressed >4hrs")$MeanCNRCI, var.equal = T)

#even simpler, early vs. later, decent
t.test(filter(nucBias_CNRCI_summary_, grepl("<4hrs", RegulationStartII))$MeanCNRCI,
       filter(nucBias_CNRCI_summary_, grepl(">4hrs", RegulationStartII))$MeanCNRCI, var.equal = T)

#but only really interested in induced in a lot of the paper so should keep seperate
#n.b. the above is weaker if including the lower FPKM entries

ggplot(filter(nucBias_CNRCI_summary_, !RegulationStartII == "Non-DE")) + aes(x = RegulationStartII, y = MeanCNRCI, fill = RegulationStartII) +
  #geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  xlab("") +
  ylab("Bias score (MeanCNRCI)") +
  geom_jitter(width = 0.1, alpha = 0.4, size =2) + 
  scale_fill_manual(values = c(`Expressed` = "grey60",
                               `Induced <4hrs` = "#D6604D",
                               `Induced >4hrs` = "#D6604D",
                               `Repressed <4hrs` = "#67A9CF",
                               `Repressed >4hrs` = "#67A9CF")) +
  theme_minimal() + 
  theme(text = element_text(size = 24), legend.position = "none") +
  Seurat::RotatedAxis()

#what about early induced vs. everything else:
t.test(filter(nucBias_CNRCI_summary_, RegulationStartII == "Induced <4hrs")$MeanCNRCI,
       filter(nucBias_CNRCI_summary_, !RegulationStartII %in% c("Induced >4hrs", "Non-DE"))$MeanCNRCI, var.equal = T)


#### old code ####

#appropriate for t test with some exceptions

nucLocP <- c(1:6)

#focus on 0-4 induced:
triali <- filter(nucBias_CNRCI_summary, RegulationStart == "Expressed ",
                 !Ens_ID_merge %in% filter(nucBias_CNRCI_summary, RegulationStart == "Induced <4hrs")$Ens_ID_merge)
summary(triali$MeanCNRCI)
summary(filter(nucBias_CNRCI_summary, RegulationStart == "Induced <4hrs")$MeanCNRCI)

t.test(triali$MeanCNRCI, 
       filter(nucBias_CNRCI_summary, RegulationStart == "Induced <4hrs")$MeanCNRCI, var.equal = T)
wilcox.test(triali$MeanCNRCI, 
       filter(nucBias_CNRCI_summary, RegulationStart == "Induced <4hrs")$MeanCNRCI, paired = F)

nucLocP[1] <- t.test(triali$MeanCNRCI, 
                          filter(nucBias_CNRCI_summary, RegulationStart == "Induced <4hrs")$MeanCNRCI, var.equal = T)$p.value
#no sig diff
#there are 13 with nuc bias <0.05, and 10 with cyto > 0.05
#there are 9 with nuc bias <0.1 and 5 with cyto > 0.1
#a v. slight nuclear bias

#remainder:
triali <- filter(nucBias_CNRCI_summary, RegulationStart == "Expressed ",
                 !Ens_ID_merge %in% filter(nucBias_CNRCI_summary, RegulationStart == "Induced 4-8hrs")$Ens_ID_merge)
nucLocP[2] <- wilcox.test(triali$MeanCNRCI, 
                     filter(nucBias_CNRCI_summary, RegulationStart == "Induced 4-8hrs")$MeanCNRCI)$p.value

triali <- filter(nucBias_CNRCI_summary, RegulationStart == "Expressed ",
                 !Ens_ID_merge %in% filter(nucBias_CNRCI_summary, RegulationStart == "Induced 8-24hrs")$Ens_ID_merge)
nucLocP[3] <- wilcox.test(triali$MeanCNRCI, 
                          filter(nucBias_CNRCI_summary, RegulationStart == "Induced 8-24hrs")$MeanCNRCI)$p.value

#repressed
triali <- filter(nucBias_CNRCI_summary, RegulationStart == "Expressed ",
                 !Ens_ID_merge %in% filter(nucBias_CNRCI_summary, RegulationStart == "Repressed <4hrs")$Ens_ID_merge)
nucLocP[4] <- t.test(triali$MeanCNRCI, 
                          filter(nucBias_CNRCI_summary, RegulationStart == "Repressed <4hrs")$MeanCNRCI, var.equal = T)$p.value

triali <- filter(nucBias_CNRCI_summary, RegulationStart == "Expressed ",
                 !Ens_ID_merge %in% filter(nucBias_CNRCI_summary, RegulationStart == "Repressed 4-8hrs")$Ens_ID_merge)
nucLocP[5] <- t.test(triali$MeanCNRCI, 
                          filter(nucBias_CNRCI_summary, RegulationStart == "Induced 8-24hrs")$MeanCNRCI, var.equal = T)$p.value

triali <- filter(nucBias_CNRCI_summary, RegulationStart == "Expressed ",
                 !Ens_ID_merge %in% filter(nucBias_CNRCI_summary, RegulationStart == "Repressed 8-24hrs")$Ens_ID_merge)
nucLocP[6] <- wilcox.test(triali$MeanCNRCI, 
                          filter(nucBias_CNRCI_summary, RegulationStart == "Induced 8-24hrs")$MeanCNRCI)$p.value


#low n, small changes, combine induced/repressed lncs:
summary(filter(nucBias_CNRCI_summary, RegulationStart == "Induced <4hrs")$MeanCNRCI)
summary(filter(nucBias_CNRCI_summary, RegulationStart %in% c("Induced <4hrs", "Repressed <4hrs"))$MeanCNRCI)

nucBias_CNRCI_summary$RegulationStartII <- as.character(nucBias_CNRCI_summary$RegulationStart)
nucBias_CNRCI_summary$RegulationStartII[!nucBias_CNRCI_summary$RegulationStartII == "Expressed "] <- 
  sapply(strsplit(nucBias_CNRCI_summary$RegulationStartII[!nucBias_CNRCI_summary$RegulationStartII == "Expressed "], " "), "[[", 2)

nucBias_CNRCI_summary$RegulationStartII[nucBias_CNRCI_summary$RegulationStartII == "<4hrs"] <- "0-4hrs"
nucBias_CNRCI_summary$RegulationStartII <- factor(nucBias_CNRCI_summary$RegulationStartII)
nucBias_CNRCI_summary$RegulationStartII <- factor(nucBias_CNRCI_summary$RegulationStartII,
                                                  levels(nucBias_CNRCI_summary$RegulationStartII)[c(4,1:3)])

ggplot(nucBias_CNRCI_summary) + aes(x = RegulationStartII, y = MeanCNRCI) +
  #geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  xlab("") +
  ylab("Bias score (MeanCNRCI)") +
  geom_jitter(width = 0.1, alpha = 0.4, size =2) +
  theme_minimal() + 
  theme(text = element_text(size = 24), legend.position = "none") +
  Seurat::RotatedAxis()

table(nucBias_CNRCI_summary$RegulationStartII)
summary(aov(MeanCNRCI~RegulationStartII,nucBias_CNRCI_summary))#ns


#merged stats
triali <- filter(nucBias_CNRCI_summary, RegulationStart == "Expressed ",
                 !Ens_ID_merge %in% filter(nucBias_CNRCI_summary, RegulationStart %in% c("Induced <4hrs", "Repressed <4hrs"))$Ens_ID_merge)
t.test(triali$MeanCNRCI, 
       filter(nucBias_CNRCI_summary, RegulationStart %in% c("Induced <4hrs", "Repressed <4hrs"))$MeanCNRCI, var.equal = T)$p.value
wilcox.test(triali$MeanCNRCI, 
       filter(nucBias_CNRCI_summary, RegulationStart %in% c("Induced <4hrs", "Repressed <4hrs"))$MeanCNRCI)$p.value
#a small but significant shift towards nuclear bias for early regulated lncRNAs

triali <- filter(nucBias_CNRCI_summary, RegulationStart == "Expressed ",
                 !Ens_ID_merge %in% filter(nucBias_CNRCI_summary, RegulationStart %in% c("Induced 4-8hrs", "Repressed 4-8hrs"))$Ens_ID_merge)
t.test(triali$MeanCNRCI, 
       filter(nucBias_CNRCI_summary, RegulationStart %in% c("Induced 4-8hrs", "Repressed 4-8hrs"))$MeanCNRCI, var.equal = T)$p.value
wilcox.test(triali$MeanCNRCI, 
            filter(nucBias_CNRCI_summary, RegulationStart %in% c("Induced 4-8hrs", "Repressed 4-8hrs"))$MeanCNRCI)$p.value

triali <- filter(nucBias_CNRCI_summary, RegulationStart == "Expressed ",
                 !Ens_ID_merge %in% filter(nucBias_CNRCI_summary, RegulationStart %in% c("Induced 8-24hrs", "Repressed 8-24hrs"))$Ens_ID_merge)
t.test(triali$MeanCNRCI, 
       filter(nucBias_CNRCI_summary, RegulationStart %in% c("Induced 8-24hrs", "Repressed 8-24hrs"))$MeanCNRCI, var.equal = T)$p.value
wilcox.test(triali$MeanCNRCI, 
            filter(nucBias_CNRCI_summary, RegulationStart %in% c("Induced 8-24hrs", "Repressed 8-24hrs"))$MeanCNRCI)$p.value
#sig too small to withstand multiple hypo correction via T, via wilcox it is better...

#should be ANOVA/Dunnet anyway? diff format w/ no repeat values
nucBias_anova <- filter(nucBias_CNRCI_summary, RegulationStartII == "Expressed ")
nucBias_anovai <- filter(nucBias_CNRCI_summary, !RegulationStartII == "Expressed ")

nucBias_anova <- filter(nucBias_anova, !Ens_ID_merge %in% nucBias_anovai$Ens_ID_merge)

nucBias_anova$RegulationStart <- "Non-DE"

nucBias_anova <- rbind(nucBias_anova, nucBias_anovai)

#no diff amongst the 7 groups:
table(nucBias_anova$RegulationStart)
summary(aov(MeanCNRCI ~ RegulationStart, nucBias_anova))
kruskal.test(MeanCNRCI ~ RegulationStart, nucBias_anova)

#borderline if condense induce/repress:
table(nucBias_anova$RegulationStartII)
summary(aov(MeanCNRCI ~ RegulationStartII, nucBias_anova))
TukeyHSD(aov(MeanCNRCI ~ RegulationStartII, nucBias_anova))


#can still condense further:
nucBias_anova$RegulationStartIII <- as.character(nucBias_anova$RegulationStartII)
nucBias_anova$RegulationStartIII[grepl("Expressed ", nucBias_anova$RegulationStartIII)] <- "Non-DE"
nucBias_anova$RegulationStartIII[!grepl("0-4hrs", nucBias_anova$RegulationStartII) & 
                                           !grepl("Expressed ", nucBias_anova$RegulationStartII)] <- "Later (4-8hr/\n8-24hr)"
nucBias_anova$RegulationStartIII[grepl("0-4hrs", nucBias_anova$RegulationStartII) & 
                                   !grepl("Expressed ", nucBias_anova$RegulationStartII)] <- "Early (0-4hr)"

#borderline if condense induce/repress and early later:
table(nucBias_anova$RegulationStartIII)
summary(aov(MeanCNRCI ~ RegulationStartIII, nucBias_anova)) #0.0509 .

library(DescTools)
nucBias_anova$RegulationStartIII <- factor(nucBias_anova$RegulationStartIII)
nucBias_anova$RegulationStartIII <- relevel(nucBias_anova$RegulationStartIII, "Non-DE")
DunnettTest(x=nucBias_anova$MeanCNRCI, g=nucBias_anova$RegulationStartIII) # 0.0289 *
TukeyHSD(aov(MeanCNRCI ~ RegulationStartIII, nucBias_anova)) #0.04 *s

ggplot(nucBias_anova) + aes(x = RegulationStartIII, y = MeanCNRCI) +
  #geom_violin() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  xlab("") +
  ylab("Bias score (MeanCNRCI)") +
  geom_jitter(width = 0.1, alpha = 0.4, size =2) +
  theme_minimal() + 
  theme(text = element_text(size = 24), legend.position = "none") +
  Seurat::RotatedAxis()

library(rcompanion)
groupwiseMean(MeanCNRCI ~ RegulationStartIII, nucBias_anova)
groupwiseMean(MeanCNRCI ~ RegulationStartII, nucBias_anova)
groupwiseMean(MeanCNRCI ~ RegulationStart, nucBias_anova)


#induced version:
ggplot(filter(nucBias_CNRCI_summary, !grepl("Repressed", RegulationStart),
              !RegulationStartII %in% "Expressed ")) + aes(x = RegulationStartIII, y = MeanCNRCI) +
  #geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  xlab("") +
  ylab("Bias score (MeanCNRCI)") +
  geom_jitter(width = 0.1, alpha = 0.4, size =2) +
  theme_minimal() + 
  theme(text = element_text(size = 24), legend.position = "none") +
  Seurat::RotatedAxis()

#repressed version
ggplot(filter(nucBias_CNRCI_summary, !grepl("Induced", RegulationStart),
              !RegulationStartII %in% "Expressed ")) + aes(x = RegulationStartIII, y = MeanCNRCI) +
  #geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  xlab("") +
  ylab("Bias score (MeanCNRCI)") +
  geom_jitter(width = 0.1, alpha = 0.4, size =2) +
  theme_minimal() + 
  theme(text = element_text(size = 24), legend.position = "none") +
  Seurat::RotatedAxis()

t.test(filter(nucBias_CNRCI_summary, !RegulationStart %in% c("Induced <4hrs", "Repressed <4hrs"))$MeanCNRCI,
       filter(nucBias_CNRCI_summary, RegulationStart %in% c("Induced <4hrs", "Repressed <4hrs"))$MeanCNRCI, var.equal = T)


nucBias_CNRCI_summary$UpDown <- "Induced"
nucBias_CNRCI_summary$UpDown[grepl("Repressed", nucBias_CNRCI_summary$RegulationStart)] <- "Repressed"
nucBias_CNRCI_summary$UpDown[grepl("Expressed", nucBias_CNRCI_summary$RegulationStart)] <- "Expressed"

ggplot(nucBias_CNRCI_summary) + aes(x = RegulationStartII, y = MeanCNRCI) +
  #geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  xlab("") +
  ylab("Bias score (MeanCNRCI)") +
  geom_jitter(position = position_dodge(width = 0.3), 
              alpha = 0.4, size =2, aes(x = RegulationStartII, y = MeanCNRCI, color = UpDown)) +
  theme_minimal() + 
  theme(text = element_text(size = 24), legend.position = "none") +
  Seurat::RotatedAxis()
