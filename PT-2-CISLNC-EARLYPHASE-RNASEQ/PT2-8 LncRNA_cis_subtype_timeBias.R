#Key element of paper

#Which types of lncRNAs make up the bias towards early activation?

#Available sub-groups: candidate cis-acting lncRNAs (HiC/eQTL/Corr) enhancer or non-enhancer (GeneHancer/FANTOM)

#Are FCs higher near candidate cis-acting lncs than in other genomic loci? compared to those near other lncs? to those near PCG only? or near neither?

#Do co-reg pairings or CClnc pairings form more readily between lnc-PCG than PCG-PCG? (probs seperate code, going beyond would also include HiC/eQTL)


#### set-up tables etc for enrichment testing as previous ####
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(grid)
library(gridExtra)
library(reshape2)

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


#### bias of lncRNAs (done in 8 too, will go into this figure too as a comparison) ####

fpkm_allGDE$Simple <- fpkm_allGDE$GeneClassUpdate
fpkm_allGDE$Simple[grepl("fide|Lnc", fpkm_allGDE$GeneClassUpdate)] <- "LncRNA"
fpkm_allGDE$Simple[grepl("coding|TF|CC", fpkm_allGDE$GeneClassUpdate)] <- "PCG"

Cluster_biotype <- as.data.frame(table(fpkm_allGDE$Simple, fpkm_allGDE$RegulationStart))

table(fpkm_allGDE$Simple)

#bias of lncRNAs:
table(fpkm_allGDE$Simple)/length(unique(fpkm_allGDE$EnsID))*100
#221 lncs are 4.35% of all DEGs
trial <- filter(Cluster_biotype, Var1 == "LncRNA")
trial$Freq/table(fpkm_allGDE$RegulationStart)*100
#looks like biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(fpkm_allGDE$RegulationStart)

ClusterNames <- trial$Var2

LncEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- 221
  d <- length(unique(fpkm_allGDE$EnsID))
  
  LncEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative  = "greater")$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p)
}

names(LncEnrich_cluster) <- ClusterNames
triali <- bind_rows(LncEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")
#2x significant biases

#percentage plots
trial$FirstRegulation <- c("Within \n4hrs", "Within \n8hrs", "Within \n24hrs", "Within \n4hrs", "Within \n8hrs", "Within \n24hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(2,3,1)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/221*100
trial$PercBackground <- trial$selection/length(unique(fpkm_allGDE$EnsID))*100

LncWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(LncWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(LncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(LncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(LncWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory-1, label = Freq), size = 5) +
  geom_col(data = filter(LncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(LncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(LncWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory+1, label = Freq), size = 5) +
  ylab("% DE LncRNAs") +
  xlab("") +
  #scale_y_continuous(limits = c(-30,35),breaks = seq(-30,30, by = 10),
  #                   labels = (c(seq(30, 0, by = -10), seq(10,30,by=10)))) +
  theme_minimal() +
  theme(text = element_text(size=15), legend.position = "none")


#### bias of enhancer lncRNAs ####

#add in eLncs
Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime_2026.csv", header = T)
fpkm_allGDE$GeneClassUpdate[fpkm_allGDE$EnsID %in% filter(Enhancer_lociII, !is.na(EnhancerVerdict))$EnsID & fpkm_allGDE$V55 == "Bona fide lncRNA"] <- "ELnc"

fpkm_allGDE$ELnc <- fpkm_allGDE$Simple
fpkm_allGDE$ELnc[grepl("ELnc", fpkm_allGDE$GeneClassUpdate)] <- "ELnc"

Cluster_biotype <- as.data.frame(table(fpkm_allGDE$ELnc, fpkm_allGDE$RegulationStart))

table(fpkm_allGDE$ELnc)

#bias of category:
table(fpkm_allGDE$ELnc)["ELnc"]
table(fpkm_allGDE$ELnc)/length(unique(fpkm_allGDE$EnsID))*100
#54 lncs are 1.06% of all DEGs
trial <- filter(Cluster_biotype, Var1 == "ELnc")
trial$Freq/table(fpkm_allGDE$RegulationStart)*100
#looks like biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(fpkm_allGDE$RegulationStart)

ClusterNames <- trial$Var2

LncEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(fpkm_allGDE$ELnc)["ELnc"]
  d <- length(unique(fpkm_allGDE$EnsID))
  
  LncEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p)
}
names(LncEnrich_cluster) <- ClusterNames
triali <- bind_rows(LncEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")
#Induced first timepoint, rarely at second
#2x significant biases

#percentage plots
trial$FirstRegulation <- c("Within \n4hrs", "Within \n8hrs", "Within \n24hrs", "Within \n4hrs", "Within \n8hrs", "Within \n24hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(2,3,1)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(fpkm_allGDE$ELnc)["ELnc"]*100
trial$PercBackground <- trial$selection/length(unique(fpkm_allGDE$EnsID))*100

ELncWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(ELncWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(ELncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(ELncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(ELncWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory-1, label = Freq), size = 4.2) +
  geom_col(data = filter(ELncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(ELncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(ELncWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory+1, label = Freq), size = 4.2) +
  ylab("% DE ELncRNAs") +
  xlab("") +
  scale_y_continuous(limits = c(-30,60),breaks = seq(-20,60, by = 20),
                     labels = (c(seq(20, 0, by = -20), seq(20,60,by=20)))) +
  theme_minimal() +
  theme(text = element_text(size = 15), legend.position = "none")
#suggests lncRNAs at enhancer sites may have their most prominent effects early on

#all lncs background:

a <- 23
b <- 54
c <- 65
d <- 221

fisher.test(data.frame("LncRNA" = c(a,b-a),
                       "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")

#trial <- data.frame(61/198*100, 19/33*100)
#colnames(trial) <- c("All DE lncRNAs", "Enhancer-transcribed\nDE lncRNAs")

trial <- data.frame(65/221*100, 23/54*100)
colnames(trial) <- c("All DE lncRNAs", "Enhancer-transcribed\nDE lncRNAs")

#melt(trial)

ggplot(melt(trial)) + aes(x = value, y = variable) +
  geom_bar(stat = "identity", fill = "olivedrab3", color = "grey60") +
  xlab("% 0-4hr induced") +
  ylab("") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0,40,20)) +
  theme(text = element_text(size =20))


#### explore distances of co-regulations ####

AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_info_2026.csv", header = T)
#AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_info_aug25_redone.csv", header = T)

#1684 pairings of lnc-PCG expressed near to each other
#4225 at 1MB
#8438 at 1MB post-dec'25

ggplot(AllLNC_AllPCG_1) + aes(x = AbsDistLnc_PCG) +
  #geom_histogram(aes(y=after_stat(density)), bins = 20, color = "black") +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+
  geom_vline(xintercept = c(250,400, 500), linetype = "dashed") +
  xlab("LncRNA-PCG Neighbour Distance") +
  ylab("Density")

Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime_2026.csv", header = T)
AllLNC_AllPCG_1$GeneClassUpdate.x[AllLNC_AllPCG_1$EnsID %in% 
                                  filter(Enhancer_lociII, !is.na(EnhancerVerdict))$EnsID & AllLNC_AllPCG_1$GeneClassUpdate.x == "Bona fide lncRNA"] <- "ELnc"

CoRegPairs_04_48_24_extended_naive <- filter(AllLNC_AllPCG_1,
                                       #AllLNC_AllPCG_1,
                                       (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                     fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                       fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                        fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                       fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
)

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


ggplot(CoRegPairs_04_48_24_extended_naive) + aes(x = AbsDistLnc_PCG) +
  #geom_histogram(aes(y=after_stat(density)), bins = 20, color = "black") +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+
  geom_vline(xintercept = c(400,500,600), linetype = "dashed") +
  xlab("Co-regulated LncRNA-PCG\nNeighbour Distance (kbp)") +
  ylab("Density")

#pairings to well expressed PCG, shows a peak around 500kbp
ggplot(filter(CoRegPairs_04_48_24_extended_naive,
                EnsID.y %in% filter(fpkm_allG, fpkm_max_treatment >20)$EnsID)) + 
  aes(x = AbsDistLnc_PCG) +
  #geom_histogram(aes(y=after_stat(density)), bins = 20, color = "black") +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+
  geom_vline(xintercept = c(400,500,600), linetype = "dashed") +
  xlab("Co-regulated LncRNA-PCG\nNeighbour Distance (kbp)") +
  ylab("Density")

#a peak just before 500kbp develops at higher PCG

ggplot(filter(CoRegPairs_04_48_24_extended_naive, LncRNA.PCG.Relationship == "Intergenic")) + aes(x = AbsDistLnc_PCG) +
  #geom_histogram(aes(y=after_stat(density)), bins = 20, color = "black") +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+
  geom_vline(xintercept = c(400,500,600), linetype = "dashed")

#Elncs seem to have a longer range set of peaks
ggplot(filter(CoRegPairs_04_48_24_extended_naive, grepl("ELnc", GeneClassUpdate.x), 
              LncRNA.PCG.Relationship == "Intergenic")) + aes(x = AbsDistLnc_PCG) +
  #geom_histogram(aes(y=after_stat(density)), bins = 60, color = "black") +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+
  geom_vline(xintercept = c(400,500,600), linetype = "dashed")

ggplot(filter(CoRegPairs_04_48_24_extended_naive, !grepl("ELnc", GeneClassUpdate.x), 
              LncRNA.PCG.Relationship == "Intergenic")) + aes(x = AbsDistLnc_PCG) +
  #geom_histogram(aes(y=after_stat(density)), bins = 60, color = "black") +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+
  geom_vline(xintercept = c(400,500,600), linetype = "dashed")

#Elncs seem to have a longer range set of peaks - sustained if pushing the PCG exprs threshold
ggplot(filter(CoRegPairs_04_48_24_extended_naive, grepl("ELnc", GeneClassUpdate.x), 
              LncRNA.PCG.Relationship == "Intergenic",
              EnsID.y %in% filter(fpkm_allG, fpkm_max_treatment >10)$EnsID)) + aes(x = AbsDistLnc_PCG) +
  #geom_histogram(aes(y=after_stat(density)), bins = 60, color = "black") +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+
  geom_vline(xintercept = c(400,500,600), linetype = "dashed")
#250kbp? misses out on longer range
#600kbp? seems long... but hits both peaks
#longer range seem to be really biased for low expressed PCGs

#argument for 500kbp in that this seems to be a group of well-expressed PCG partners
#but conversely for enhancers, this seems to be mainly low-expressed... hmmm

#cut-off of 250kbp used initially for stats, they broadly work tho little low sig maybe
#400kbp also worked in latter iteration
#n.b. the transcriptional ripples paper also shows a bump at 250-400kbp near their IEGs, after a big peak <50kbp (mouse fibroblast activation)

#n.b. compare to all pairs:
AllLNC_AllPCG_1$CoReg <- "No"
AllLNC_AllPCG_1$CoReg[AllLNC_AllPCG_1$pairs %in% CoRegPairs_04_48_24_extended_naive$pairs] <- "Yes"

ggplot(AllLNC_AllPCG_1) + aes(x = AbsDistLnc_PCG, fill = CoReg) +
  geom_density(color="#e9ecef", alpha=0.8) +
  geom_vline(xintercept = c(400,500,600), linetype = "dashed")

ggplot(filter(AllLNC_AllPCG_1, LncRNA.PCG.Relationship == "Intergenic")) + aes(x = AbsDistLnc_PCG, fill = CoReg) +
  geom_density(color="#e9ecef", alpha=0.8)+
  geom_vline(xintercept = c(400,500,600), linetype = "dashed")
#second peak emphasis again... meaning? long-distance relationships quite specific in this range?

#really clear diffs between co-reg pair distances and other pairs - big bump in these peaking at ~325-500kbp

ggplot(filter(AllLNC_AllPCG_1, grepl("ELnc", GeneClassUpdate.x),
              LncRNA.PCG.Relationship == "Intergenic")) + aes(x = AbsDistLnc_PCG, fill = CoReg) +
  geom_density(color="#e9ecef", alpha=0.8)+
  geom_vline(xintercept = c(400,500,600), linetype = "dashed")
#stronger for enhancers... meaning? 

#600kbp will capture a lot of enhancer peaks:
View(filter(AllLNC_AllPCG_1, grepl("ELnc", GeneClassUpdate.x),
            LncRNA.PCG.Relationship == "Intergenic"))

#no indication which peaks or distances are more reliable for actual lnc-PCG relationship yet

#but going too low may bias away from enhancers...

#Genehancer paper indicates most of their connections occur within 100kbp...


#### define cclncRNAs - naive, i.e. no HiC - closest neighbour approach - enhancer association ####

#no conclusion on distance to use from above figures
#best approach: stay stringent, close pairs only (Genehancer paper)
#filter on distance, 250kbp initially done, use higher if backed up by HiC or eQTL
AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_info_2026.csv", header = T)
table(AllLNC_AllPCG_1$LncRNA.PCG.Relationship)#64/274/8097/3
AllLNC_AllPCG_1 <- filter(AllLNC_AllPCG_1, AbsDistLnc_PCG <250)#2564
table(AllLNC_AllPCG_1$LncRNA.PCG.Relationship)#62/274/2225/3

Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime_2026.csv", header = T)
AllLNC_AllPCG_1$GeneClassUpdate.x[AllLNC_AllPCG_1$EnsID %in% 
                                    filter(Enhancer_lociII, !is.na(EnhancerVerdict))$EnsID & 
                                    AllLNC_AllPCG_1$GeneClassUpdate.x == "Bona fide lncRNA"] <- "ELnc"

#previous approaches were a blanket all neighbours approach
CoRegPairs_04_48_24_extended_naive <- filter(AllLNC_AllPCG_1,
                                             #AllLNC_AllPCG_1,
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                           fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                               (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                             fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                              fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                               (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                             fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
)#282 @250kbp (8.5% bidirectional)
table(CoRegPairs_04_48_24_extended_naive$LncRNA.PCG.Relationship)

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
)#153 @250kbp (10.4% bidirectional)
table(CoRegPairs_04_48_24_extended_naiveSame$LncRNA.PCG.Relationship)


#more stringent options to try:
#closest neighbour, or closest up/down neighbours (surrounding neighbours) approach
#likely enriching for actual connections
#needs the non-abs distance (for up and downstream)

#reimporting to reset filtering:
AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_info_2026.csv", header = T)
Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime_2026.csv", header = T)
AllLNC_AllPCG_1$GeneClassUpdate.x[AllLNC_AllPCG_1$EnsID %in% 
                                    filter(Enhancer_lociII, !is.na(EnhancerVerdict))$EnsID & 
                                    AllLNC_AllPCG_1$GeneClassUpdate.x == "Bona fide lncRNA"] <- "ELnc"

closestNeighbour <- filter(AllLNC_AllPCG_1, AbsDistLnc_PCG <250)
closestNeighbour <- split(closestNeighbour, closestNeighbour$EnsID)

#some variations:
#closest neighbour
#closestNeighbour <- lapply(closestNeighbour, function(x){
#  filter(x, AbsDistLnc_PCG == min(AbsDistLnc_PCG ))
#  }
#  )

#closest x neighbours within 250kbp etc:
#closestNeighbour <- lapply(closestNeighbour, function(x){
#    x[order(x$DisLnc_PCG, decreasing = F),][1:5,]}
#  )

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
CoRegPairs_04_48_24_extended_naive <- filter(AllLNC_AllPCG_1,
                                             #AllLNC_AllPCG_1,
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                           fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                               (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                             fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                              fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                               (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                             fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
)
#58 @250kbp - closest
#96 @ 250kbp - closest surrounding (25% bidirectional)
table(CoRegPairs_04_48_24_extended_naive$LncRNA.PCG.Relationship)
length(unique(CoRegPairs_04_48_24_extended_naive$EnsID))
#81 @ 250kbp - surrounding
length(unique(filter(CoRegPairs_04_48_24_extended_naive, GeneClassUpdate.x == "ELnc")$EnsID))
#27 are elncs
27/81 #33%
54/221 #24% of all DELs

aii <- 27
bii <- 81    
cii <- 54   
dii <- 221   

fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")


#same timeframe
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
#34 @250kbp - closest
#55 @250kbp - surrounding (29% bidirectional)
table(CoRegPairs_04_48_24_extended_naiveSame$LncRNA.PCG.Relationship)
length(unique(CoRegPairs_04_48_24_extended_naiveSame$EnsID))
#51 @ 250kbp - surrounding
length(unique(filter(CoRegPairs_04_48_24_extended_naiveSame, GeneClassUpdate.x == "ELnc")$EnsID))
#21 are elncs
21/51 #41%
54/221 #24% of all DELs

aii <- 21
bii <- 51    
cii <- 54   
dii <- 221   

fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

trial <- data.frame(54/221*100, 21/51*100, 27/81*100)

colnames(trial) <- c("All DE lncRNAs", "CClncRNAs (same)",  "CClncRNAs (same/later)")

#melt(trial)

trial <- melt(trial)

trial$variable <- factor(trial$variable)
trial$variable <- factor(trial$variable, levels = levels(trial$variable)[c(1,3,2)])

ggplot(trial) + aes(x = value, y = variable) +
  geom_bar(stat = "identity", fill = "olivedrab3", color = "grey60") +
  xlab("% ElncRNAs") +
  ylab("") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0,40,20), limits = c(0,45)) +
  theme(text = element_text(size =20))

#
#### supplementary table of these cclncs ####

#add PCG and lncRNA timeframes:

table(fpkm_allGDE$RegulationStart)

trial <- merge(CoRegPairs_04_48_24_extended_naive, fpkm_allGDE[,c(1,46)], by = "EnsID")
trial <- merge(trial, fpkm_allGDE[,c(1,46)], by.x = "EnsID.y", by.y = "EnsID")

write.csv(trial, "SuppTable3_CClncRNAs.csv")

#### bias of lncs with same timeframe lnc 2D neighbours (i.e. naive approach first) ####

fpkm_allGDE$CCLncNaive <- fpkm_allGDE$Simple
fpkm_allGDE$CCLncNaive[fpkm_allGDE$EnsID %in% CoRegPairs_04_48_24_extended_naiveSame$EnsID] <- "CCLncNaive"

Cluster_biotype <- as.data.frame(table(fpkm_allGDE$CCLncNaive, fpkm_allGDE$RegulationStart))

table(fpkm_allGDE$CCLncNaive)

table(fpkm_allGDE$CCLncNaive)["CCLncNaive"]
table(fpkm_allGDE$CCLncNaive)/length(unique(fpkm_allGDE$EnsID))*100

trial <- filter(Cluster_biotype, Var1 == "CCLncNaive")
trial$Freq/table(fpkm_allGDE$RegulationStart)*100
#looks like  biases are v possible

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(fpkm_allGDE$RegulationStart)

ClusterNames <- trial$Var2

LncEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(fpkm_allGDE$CCLncNaive)["CCLncNaive"]
  d <- length(unique(fpkm_allGDE$EnsID))
  
  LncEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p)
}
names(LncEnrich_cluster) <- ClusterNames
triali <- bind_rows(LncEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")
#Induced first timepoint, rarely at second
#2x significant biases

#percentage plots
trial$FirstRegulation <- c("Within \n4hrs", "Within \n8hrs", "Within \n24hrs", "Within \n4hrs", "Within \n8hrs", "Within \n24hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(2,3,1)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(fpkm_allGDE$CCLncNaive)["CCLncNaive"]*100
trial$PercBackground <- trial$selection/length(unique(fpkm_allGDE$EnsID))*100

CCLncNaiveWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])
CCLncNaiveWaveBias$Var1 <- "CisLncSame"

ggplot(CCLncNaiveWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(CCLncNaiveWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(CCLncNaiveWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CCLncNaiveWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 4.2) +
  geom_col(data = filter(CCLncNaiveWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(CCLncNaiveWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CCLncNaiveWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory, label = Freq), size = 4.2) +
  ylab("% DE CCLncNaives") +
  xlab("") +
 # scale_y_continuous(limits = c(-30,45),breaks = seq(-20,45, by = 20),
#                     labels = (c(seq(20, 0, by = -20), seq(20,45,by=20)))) +
  theme_minimal() +
  theme(text = element_text(size = 15), legend.position = "none")

#600kbp: 2.0 OR
#500kbp: 2.0 OR
#250kbp: 2.1/0.001 OR/p
#250kbp closest: 3.2/0.0068
#250kbp surrounding: 3.2/5.3e-5

#first wave, stronger OR and p than regular

CCLncNaiveWaveBias
sum(CCLncNaiveWaveBias$Freq)

#600kbp: 2.1 
aii <- 48    
bii <- 137    
cii <- 65   
dii <- 221 
fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

#500kbp: 2.1 
aii <- 46    
bii <- 129    
cii <- 65   
dii <- 221 
fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

#250kbp version 1.8/0.038
aii <- 32    
bii <- 87    
cii <- 65   
dii <- 221
fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

#250kbp closest version 2.5/0.0141
aii <- 16    
bii <- 34    
cii <- 65   
dii <- 221
fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

#250kbp surrounding version 2.8/0.0018
aii <- 24    
bii <- 51    
cii <- 65   
dii <- 221
fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")


#
#### bias of 2d same+later CCLncRNAs ####

fpkm_allGDE$CCLnc <- fpkm_allGDE$Simple
fpkm_allGDE$CCLnc[fpkm_allGDE$EnsID %in% CoRegPairs_04_48_24_extended_naive$EnsID] <- "CCLnc"

Cluster_biotype <- as.data.frame(table(fpkm_allGDE$CCLnc, fpkm_allGDE$RegulationStart))

table(fpkm_allGDE$CCLnc)

table(fpkm_allGDE$CCLnc)["CCLnc"]
table(fpkm_allGDE$CCLnc)/length(unique(fpkm_allGDE$EnsID))*100

trial <- filter(Cluster_biotype, Var1 == "CCLnc")
trial$Freq/table(fpkm_allGDE$RegulationStart)*100
#looks like  biases are v possible

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(fpkm_allGDE$RegulationStart)

ClusterNames <- trial$Var2

LncEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(fpkm_allGDE$CCLnc)["CCLnc"]
  d <- length(unique(fpkm_allGDE$EnsID))
  
  LncEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p)
}
names(LncEnrich_cluster) <- ClusterNames
triali <- bind_rows(LncEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")
#Induced first timepoint, rarely at second
#2x significant biases

#percentage plots
trial$FirstRegulation <- c("Within \n4hrs", "Within \n8hrs", "Within \n24hrs", "Within \n4hrs", "Within \n8hrs", "Within \n24hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(2,3,1)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(fpkm_allGDE$CCLnc)["CCLnc"]*100
trial$PercBackground <- trial$selection/length(unique(fpkm_allGDE$EnsID))*100

CCLncWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])
CCLncWaveBias$Var1 <- "CCLncRNA"

ggplot(CCLncWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(CCLncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(CCLncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CCLncWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 4.2) +
  geom_col(data = filter(CCLncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(CCLncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CCLncWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory, label = Freq), size = 4.2) +
  ylab("% DE CCLncs") +
  xlab("") +
#  scale_y_continuous(limits = c(-30,45),breaks = seq(-20,45, by = 20),
 #                    labels = (c(seq(20, 0, by = -20), seq(20,45,by=20)))) +
  theme_minimal() +
  theme(text = element_text(size = 15), legend.position = "none")

#600kbp: OR 2.0
#500kbp: OR 2.2
#250kbp: OR 2.1/0.0002
#250kbp closest: 4.2/1.1e-7
#250kbp surrounding: 3.1/1.7e-6


#600kbp: 3.6 
aii <- 57    
bii <- 160    
cii <- 65   
dii <- 221 
fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

#500kbp: 4.3
aii <-  57   
bii <-  154   
cii <-  65  
dii <-  221
fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

#250kbp 2.2/0.008
aii <- 45    
bii <- 124    
cii <- 65   
dii <- 221
fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

#250kbp closest 4.3/1e-6
aii <- 31    
bii <- 58    
cii <- 65   
dii <- 221
fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

#250kbp surrounding 3.3/5.9e-5
aii <- 37    
bii <- 81    
cii <- 65   
dii <- 221
fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")


#
#### combined figure ####

AllTypesWaveBias <- rbind(LncWaveBias, 
                          ELncWaveBias,
                          CCLncNaiveWaveBias,
                          CCLncWaveBias)


AllTypesWaveBias$OR_corrected <- AllTypesWaveBias$OR
AllTypesWaveBias$OR_corrected[AllTypesWaveBias$p_adj >0.05] <- NA
AllTypesWaveBias$OR_corrected <- AllTypesWaveBias$OR_corrected + 0.001
AllTypesWaveBias$`Log2(Odds Ratio)` <- log2(AllTypesWaveBias$OR_corrected)

AllTypesWaveBias$`Log2(Odds Ratio)`[AllTypesWaveBias$`Log2(Odds Ratio)` < -5] <- -5
AllTypesWaveBias$`Log2(Odds Ratio)`[AllTypesWaveBias$`Log2(Odds Ratio)` >5] <- 5

AllTypesWaveBias$padj_corrected <- AllTypesWaveBias$p_adj
AllTypesWaveBias$padj_corrected[AllTypesWaveBias$p_adj >0.1] <- NA

AllTypesWaveBias$padj_simple <- NA
#AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.1] <- "p<0.1"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.05] <- "p<0.05"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.01] <- "p<0.01"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.001] <- "p<0.001"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.0001] <- "p<0.0001"

AllTypesWaveBias$UpDownType <- paste(AllTypesWaveBias$Var1, AllTypesWaveBias$UpDown, sep = "-")

#save a copy:
#write.csv(AllTypesWaveBias, "SVSMC_AllTypesWaveBias.csv", row.names = F)

#insert spacers for plotting:
AllTypesWaveBias <- AllTypesWaveBias[,c(15,4,12,14)]

appendSpacers <- data.frame(UpDownType = c(rep("Spacer1",3),
                                           rep("Spacer2",3),
                                           rep("Spacer3",3),
                                           rep("Spacer4",3)#,
                                           #rep("Spacer5",3),
                                           #rep("Spacer6",3),
                                           #rep("Spacer7",3),
                                           #rep("Spacer8",3)
), 
"FR" = rep(c("Within \n4hrs",
             "Within \n8hrs",
             "Within \n24hrs"), 4), 
"OR" = NA,  
`padj_simple` =NA)

colnames(appendSpacers) <- colnames(AllTypesWaveBias)

AllTypesWaveBias <- rbind(AllTypesWaveBias, appendSpacers)

AllTypesWaveBias$UpDownType <- factor(AllTypesWaveBias$UpDownType)

#original
#AllTypesWaveBias$UpDownType <- factor(AllTypesWaveBias$UpDownType, 
#                                      levels = levels(AllTypesWaveBias$UpDownType)[c(16,15,22,21,17,
#                                                                                     10,9,2,1,18,
#                                                                                     12,11,4,3,19,
#                                                                                     14,13,6,5,20,
#                                                                                     8,7)])
#for simpler:
#AllTypesWaveBias$UpDownType <- factor(AllTypesWaveBias$UpDownType, 
#                                      levels = levels(AllTypesWaveBias$UpDownType)[c(2,1,5,4,3,6,7)])

#without "remaining" lncs
#AllTypesWaveBias$UpDownType <- factor(AllTypesWaveBias$UpDownType, 
#                                      levels = levels(AllTypesWaveBias$UpDownType)[c(14,13,9,2,1,10,4,3,11,6,5,12,8,7)])

#simple3:
AllTypesWaveBias$UpDownType <- factor(AllTypesWaveBias$UpDownType, 
                                      levels = levels(AllTypesWaveBias$UpDownType)[c(2,1,9,4,3,10,6,5,11,8,7)])

#8 space version: c(16,15,17,26,25,18,10,9,19,2,1,20,12,11,21,4,3,22,14,13,23,6,5,24,8,7)

myColor <- colorRampPalette(c("steelblue", "white", "red"))(50)
myBreaks <- c(seq(min(AllTypesWaveBias$`Log2(Odds Ratio)`, na.rm = T), 0, 
                  length.out=ceiling(50/2)), 
              seq(max(AllTypesWaveBias$`Log2(Odds Ratio)`, na.rm = T)/50, 
                  max(AllTypesWaveBias$`Log2(Odds Ratio)`, na.rm = T), 
                  length.out=floor(50/2)))

AllTypesWaveBias$FirstRegulation <- as.character(AllTypesWaveBias$FirstRegulation)
AllTypesWaveBias$FirstRegulation[AllTypesWaveBias$FirstRegulation == "Within \n4hrs"] <- "0-4hrs"
AllTypesWaveBias$FirstRegulation[AllTypesWaveBias$FirstRegulation == "Within \n8hrs"] <- "4-8hrs"
AllTypesWaveBias$FirstRegulation[AllTypesWaveBias$FirstRegulation == "Within \n24hrs"] <- "8-24hrs"

AllTypesWaveBias$FirstRegulation <- factor(AllTypesWaveBias$FirstRegulation)
#AllTypesWaveBias$FirstRegulation <- factor(AllTypesWaveBias$FirstRegulation, levels = levels(AllTypesWaveBias$FirstRegulation)[c(2,3,1)])

ggplot(AllTypesWaveBias) + aes(x = FirstRegulation, y = UpDownType, size = padj_simple, fill = `Log2(Odds Ratio)`) +
  geom_point(color = "grey30", shape = 21) +
  xlab("") +
  ylab("") +
  scale_size_discrete(range = c(10,18), limits = c(#"p<0.1", 
                                                   "p<0.05", 
                                                   "p<0.01", 
                                                   "p<0.001", 
                                                   "p<0.0001")) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "red") +
  theme_minimal() +
  theme(legend.key.size = unit(1.4, "line"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        axis.text.x = element_text(size=20),
        axis.title.x = element_text(size=18),
        #axis.text.y = element_blank()
  ) + Seurat::RotatedAxis()

#write.csv(AllTypesWaveBias, "SVSMC_AllTypesWaveBias.csv", row.names = F)

#

#### combined figure 0-4hr focus ####

CCLncNaiveWaveBias
sum(CCLncNaiveWaveBias$Freq)

#recheck these if want to compare later
#bias vs de lncs -250kbp
#aii <- 33    
#bii <- 84    
#cii <- 61   
#dii <- 198   

#bias vs de lncs -150kbp
#aii <- 26    
#bii <- 67    
#cii <- 61   
#dii <- 198   

#bias vs de lncs -400kbp (improves)
#aii <- 40    
#bii <- 101    
#cii <- 61   
#dii <- 198 

#bias vs de lncs -250kbp surrounding neighbours
aii <- 24    
bii <- 51    
cii <- 65   
dii <- 221   

fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

CCLncWaveBias
sum(CCLncWaveBias$Freq)

#recheck these if using later
#bias vs de lncs -250kbp
#a <- 45
#b <- 109
#c <- 61
#d <- 198

#bias vs de lncs -150kbp
#a <- 39
#b <- 94
#c <- 61
#d <- 198

#bias vs de lncs -400kbp
#a <- 51
#b <- 129
#c <- 61
#d <- 198

#bias vs de lncs -250kbp surrounding neighbours
aii <- 37
bii <- 81    
cii <- 65   
dii <- 221   

fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

#250kbp
#trial <- data.frame(61/198*100, 33/84*100, 45/109*100)

#150kbp
#trial <- data.frame(61/198*100, 26/67*100, 39/94*100)

#400kbp
#trial <- data.frame(61/198*100, 40/101*100, 51/129*100)

#250kbp, surrounding neighbours
trial <- data.frame(65/221*100, 24/51*100, 37/81*100)

colnames(trial) <- c("All DE lncRNAs", "CClncRNAs (same)",  "CClncRNAs (same/later)")

#melt(trial)

trial <- melt(trial)

trial$variable <- factor(trial$variable)
trial$variable <- factor(trial$variable, levels = levels(trial$variable)[c(1,3,2)])

ggplot(trial) + aes(x = value, y = variable) +
  geom_bar(stat = "identity", fill = "olivedrab3", color = "grey60") +
  xlab("% 0-4hr induced") +
  ylab("") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0,50,25), limits = c(0,50)) +
  theme(text = element_text(size =20))


#
#### supplementary table for pub. naive CClncRNAs ####

#add DE clusters to the table:


#
#### enrichment of concordantly changing neighbours near lncs(same timeframe only) ####

fpkm_allGDE_within_4 <- rbind(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Downwithin_4)
fpkm_allGDE_within_8 <- rbind(fpkm_allGDE_Upwithin_8, fpkm_allGDE_Downwithin_8)
fpkm_allGDE_within_24 <- rbind(fpkm_allGDE_Upwithin_24, fpkm_allGDE_Downwithin_24)

#background of totally non-DE lncs overall whole timecourse, expressed in given timepoint:
fpkm_allG_04 <- filter(fpkm_allG, (Hour0_meanFPKM>0.8 | Hour4_meanFPKM>0.8), !EnsID %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                     fpkm_allGDE_within_24$EnsID)
                       )
fpkm_allG_48 <- filter(fpkm_allG, (Hour0_meanFPKM>0.8 | Hour4_meanFPKM>0.8 | Hour8_meanFPKM>0.8), !EnsID %in% c(fpkm_allGDE_within_4$EnsID, 
                                                                                       fpkm_allGDE_within_24$EnsID)
                       )
fpkm_allG_824 <- filter(fpkm_allG, fpkm_max_treatment >0.8, !EnsID %in% c(fpkm_allGDE_within_4$EnsID, 
                                                                                         fpkm_allGDE_within_8$EnsID)
                        )

#background is "genes expressed in that timeframe which are not DEGs at any point in the timecourse"

#simpler version: (less significant, possibly as including DEGs in other timepoints in the background)
#fpkm_allG_04 <- filter(fpkm_allG, (Hour0_meanFPKM>1 | Hour4_meanFPKM>1))
#fpkm_allG_48 <- filter(fpkm_allG, (Hour0_meanFPKM>1 | Hour4_meanFPKM>1 | Hour8_meanFPKM>1))
#fpkm_allG_824 <- filter(fpkm_allG, fpkm_max_treatment >1)

#further restrictions on defining hits, DE lncs in given timepoint, strictly those expressed in given timepoint:
#fpkm_allGDE_within_8i <- filter(fpkm_allGDE_within_8, (Hour4_meanFPKM>1 | Hour8_meanFPKM>1))
#fpkm_allGDE_within_24i <- filter(fpkm_allGDE_within_24, (Hour8_meanFPKM>1 | Hour24_meanFPKM>1))

backgrounds_list <- list(fpkm_allG_04, fpkm_allG_48, fpkm_allG_824)
results_list <- list()

#just up DELs:
selection_list <- list(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Upwithin_8, fpkm_allGDE_Upwithin_24)
#with up targets:
hits_list <- list(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Upwithin_8, fpkm_allGDE_Upwithin_24)


for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$EnsID)
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate))$EnsID))
  #all hits in background
  
  #c <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
  #                          EnsID %in% filter(closestNeighbour, EnsID.y %in% cisPotentialTargets)$EnsID)$EnsID))
  
  #or use all neighbours
  c <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
                            EnsID %in% filter(AllLNC_AllPCG_1, EnsID.y %in% cisPotentialTargets)$EnsID)$EnsID))
  
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
                            EnsID %in% selection_list[[i]]$EnsID)$EnsID))
  
  #pre-filter to closest neighbour
  #a <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
  #                          EnsID %in% filter(closestNeighbour, 
  #                                            EnsID %in% selection_list[[i]]$EnsID, 
  #                                            EnsID.y %in% cisPotentialTargets)$EnsID)$EnsID))
  
  #or take all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
                            EnsID %in% filter(AllLNC_AllPCG_1, 
                                              EnsID %in% selection_list[[i]]$EnsID, 
                                              EnsID.y %in% cisPotentialTargets)$EnsID)$EnsID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-24hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g2bi <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "#D6604D") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  xlab("") +   ylab("\n% With an induced\nneighbour") +
  theme_minimal()
g2bi

DEL_PCG_type_bi <- DEL_PCG_type
DEL_PCG_type_bi$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type_bi$OR <- unlist(results_list)[c(6,12,18)]
#write.csv(DEL_PCG_type_bi, "SVSMC_inducedDEL_EL_DEneighboursame_conc.csv", row.names = F)

#600kbp bias 2.2/0.004
#500kbp bias 2.4/0.0008
#250kbp bias 2.4/0.0015
#250kbp closest neighbour 3.16/0.0014
#250kbp closest 3 neighbours 3.1/0.0001
#               5 neighbours 2.7/0.0003
#       surrounding neighbours 3.19/0.00015

p.adjust(DEL_PCG_type_bi$p, method = "bonferroni")


#simple version for poster/GA etc
colnames(DEL_PCG_type) <- c("Cis LncRNA Candidates", "Background LncRNAs", "Timeframe")
DEL_PCG_typeM <- melt(DEL_PCG_type[,-4])

ggplot(DEL_PCG_typeM) + aes(x = Timeframe, y = value, fill = variable) +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.7), color = "grey60") +
  theme_minimal() +
  scale_fill_manual(values = c(`Cis LncRNA Candidates` = "#D6604D", `Background LncRNAs` = "grey30")) +
  xlab("") +
  ylab("% With a co-regulated\nneighbour") +
  theme(text = element_text(size =18)) + Seurat::RotatedAxis()


#just down DELs:
selection_list <- list(fpkm_allGDE_Downwithin_4, fpkm_allGDE_Downwithin_8, fpkm_allGDE_Downwithin_24)
#with down targets:
hits_list <- list(fpkm_allGDE_Downwithin_4, fpkm_allGDE_Downwithin_8, fpkm_allGDE_Downwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$EnsID)
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate))$EnsID))
  
  #option to use closest neighbour instead of all
  #c <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
  #                          EnsID %in% filter(closestNeighbour, EnsID.y %in% cisPotentialTargets)$EnsID)$EnsID))
  
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
                            EnsID %in% filter(AllLNC_AllPCG_1, EnsID.y %in% cisPotentialTargets)$EnsID
                            )$EnsID))
  
  
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
                            EnsID %in% selection_list[[i]]$EnsID
  )$EnsID))
  
  #a <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
  #                          EnsID %in% filter(closestNeighbour, 
  #                                            EnsID %in% selection_list[[i]]$EnsID, 
  #                                            EnsID.y %in% cisPotentialTargets)$EnsID)$EnsID))
                            
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
                            EnsID %in% filter(AllLNC_AllPCG_1, EnsID %in% selection_list[[i]]$EnsID, EnsID.y %in% cisPotentialTargets)$EnsID
  )$EnsID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-24hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g2bii <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "#67A9CF") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  xlab("") +   ylab("\n% With a repressed\nneighbour") +
  theme_minimal()
g2bii



DEL_PCG_type_bii <- DEL_PCG_type
DEL_PCG_type_bii$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type_bii$OR <- unlist(results_list)[c(6,12,18)]
#write.csv(DEL_PCG_type_bii, "SVSMC_repressedDEL_EL_DEneighboursame_conc.csv", row.names = F)

#600kbp bias 1.7/0.09
#500kbp bias 2.0/0.04
#250kbp bias 2.0/0.062
#250kbp closest neighbour 1.3/0.424
#250kbp closest 3 neighbours 1.35/0.33
#               5 neighbours 1.17/0.44
#               surrounding neighbours 1.25/0.410

p.adjust(DEL_PCG_type_bii$p, method = "bonferroni")

#both up and down
grid.arrange(g2bi, g2bii, ncol = 2,
             left = textGrob("", gp = gpar(fontface = 'bold', fontsize = 12), rot = 90))


#### enrichment of DISconcordantly changing neighbours near lncs(same timeframe only) ####

#keep consistent with above
backgrounds_list <- list(fpkm_allG_04, fpkm_allG_48, fpkm_allG_824)
results_list <- list()

#just up DELs:
selection_list <- list(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Upwithin_8, fpkm_allGDE_Upwithin_24)
#with DOWN targets:
hits_list <- list(fpkm_allGDE_Downwithin_4, fpkm_allGDE_Downwithin_8, fpkm_allGDE_Downwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$EnsID)
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate))$EnsID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
                            EnsID %in% filter(AllLNC_AllPCG_1, EnsID.y %in% cisPotentialTargets)$EnsID
  )$EnsID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
                            EnsID %in% selection_list[[i]]$EnsID
  )$EnsID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
                            EnsID %in% filter(AllLNC_AllPCG_1, EnsID %in% selection_list[[i]]$EnsID, EnsID.y %in% cisPotentialTargets)$EnsID
  )$EnsID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-24hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g3bi <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "#D6604D") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  xlab("") +   ylab("\n% With a repressed\nneighbour") +
  theme_minimal()
g3bi

DEL_PCG_type_bi_disc <- DEL_PCG_type
DEL_PCG_type_bi_disc$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type_bi_disc$OR <- unlist(results_list)[c(6,12,18)]
#write.csv(DEL_PCG_type, "SVSMC_inducedDEL_EL_DEneighboursame_conc.csv", row.names = F)

p.adjust(DEL_PCG_type_bi_disc$p, method = "bonferroni")


#just down DELs:
selection_list <- list(fpkm_allGDE_Downwithin_4, fpkm_allGDE_Downwithin_8, fpkm_allGDE_Downwithin_24)
#with UP targets:
hits_list <- list(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Upwithin_8, fpkm_allGDE_Upwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$EnsID)
  
  #all background lncs
  d <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate))$EnsID))
  #all hits in background
  c <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
                            EnsID %in% filter(AllLNC_AllPCG_1, EnsID.y %in% cisPotentialTargets)$EnsID
  )$EnsID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
                            EnsID %in% selection_list[[i]]$EnsID
  )$EnsID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
                            EnsID %in% filter(AllLNC_AllPCG_1, EnsID %in% selection_list[[i]]$EnsID, EnsID.y %in% cisPotentialTargets)$EnsID
  )$EnsID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-24hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g3bii <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "#67A9CF") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  xlab("") +   ylab("\n% With an induced\nneighbour") +
  theme_minimal()
g3bii

DEL_PCG_type_bii_disc <- DEL_PCG_type
DEL_PCG_type_bii_disc$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type_bii_disc$OR <- unlist(results_list)[c(6,12,18)]
#write.csv(DEL_PCG_type, "SVSMC_repressedDEL_EL_DEneighboursame_conc.csv", row.names = F)

p.adjust(DEL_PCG_type_bii$p, method = "bonferroni")

#strong for both up and down
grid.arrange(g3bi, g3bii, ncol = 2,
             left = textGrob("", gp = gpar(fontface = 'bold', fontsize = 12), rot = 90))


#### concordant/discordant comparison on one graph ####

#% induced neighbor on y axis
DEL_PCG_type_bi
DEL_PCG_type_bii_disc

g2bi 
g3bii

colnames(DEL_PCG_type_bi)[1] <- "Induced_DEL"
colnames(DEL_PCG_type_bii_disc)[1] <- "Repressed_DEL"

colnames(DEL_PCG_type_bi)[2] <- "NonDEL_pairs"
colnames(DEL_PCG_type_bii_disc)[2] <- "NonDEL_pairs"

library(reshape2)
DEL_PCG_type_plot <- unique(rbind(melt(DEL_PCG_type_bi[,1:3]),
                                  melt(DEL_PCG_type_bii_disc[,1:3])))

DEL_PCG_type_plot$variable <- as.character(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable[grepl("Induced", DEL_PCG_type_plot$variable)] <- "Induced lncRNAs"
DEL_PCG_type_plot$variable[grepl("Repressed", DEL_PCG_type_plot$variable)] <- "Repressed lncRNAs"
DEL_PCG_type_plot$variable[grepl("NonDEL", DEL_PCG_type_plot$variable)] <- "Non-DE lncRNAs"

DEL_PCG_type_plot$variable <- as.factor(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable <- factor(DEL_PCG_type_plot$variable, levels = levels(DEL_PCG_type_plot$variable)[c(1,3,2)])

ggplot(DEL_PCG_type_plot) + aes(x = RegulationBegins, y = value, fill = variable) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.7), color = "grey60") +
  theme_minimal() +
  theme(text = element_text(size = 32)) +
  scale_fill_manual(values = c(`Induced lncRNAs` = "#D6604D", `Repressed lncRNAs` = "#67A9CF", `Non-DE lncRNAs` = "grey30")) +
  xlab("") +
  ylab("% With induced\nneighbour")

#just concordant:
ggplot(filter(DEL_PCG_type_plot, !grepl("Repress", variable))) + aes(x = RegulationBegins, y = value, fill = variable) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.7), color = "grey60") +
  theme_minimal() +
  theme(text = element_text(size = 32)) + #, legend.position = "none") +
  scale_fill_manual(values = c(`Induced lncRNAs` = "#D6604D", `Repressed lncRNAs` = "#67A9CF", `Non-DE lncRNAs` = "grey30")) +
  xlab("") +
  ylab("% With induced\nneighbour")+
  theme(text = element_text(size=24)) + Seurat::RotatedAxis()

#just disconcordant:
ggplot(filter(DEL_PCG_type_plot, !grepl("Induce", variable))) + aes(x = RegulationBegins, y = value, fill = variable) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.7), color = "grey60") +
  theme_minimal() +
  theme(text = element_text(size = 32), legend.position = "none") +
  scale_fill_manual(values = c(`Induced lncRNAs` = "#D6604D", `Repressed lncRNAs` = "#67A9CF", `Non-DE lncRNAs` = "grey30")) +
  xlab("") +
  ylab("% With induced\nneighbour") +
  theme(text = element_text(size=24)) + Seurat::RotatedAxis()


#% induced neighbor on y axis
DEL_PCG_type_bii
DEL_PCG_type_bi_disc

g2bii 
g3bi

colnames(DEL_PCG_type_bii)[1] <- "Repressed_DEL"
colnames(DEL_PCG_type_bi_disc)[1] <- "Induced_DEL"

colnames(DEL_PCG_type_bii)[2] <- "NonDEL_pairs"
colnames(DEL_PCG_type_bi_disc)[2] <- "NonDEL_pairs"

library(reshape2)
DEL_PCG_type_plot <- unique(rbind(melt(DEL_PCG_type_bii[,1:3]),
                                  melt(DEL_PCG_type_bi_disc[,1:3])))

DEL_PCG_type_plot$variable <- as.character(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable[grepl("Induced", DEL_PCG_type_plot$variable)] <- "Induced lncRNAs"
DEL_PCG_type_plot$variable[grepl("Repressed", DEL_PCG_type_plot$variable)] <- "Repressed lncRNAs"
DEL_PCG_type_plot$variable[grepl("NonDEL", DEL_PCG_type_plot$variable)] <- "Non-DE lncRNAs"

DEL_PCG_type_plot$variable <- as.factor(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable <- factor(DEL_PCG_type_plot$variable, levels = levels(DEL_PCG_type_plot$variable)[c(3,1,2)])

ggplot(DEL_PCG_type_plot) + aes(x = RegulationBegins, y = value, fill = variable) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.7), color = "grey60") +
  theme_minimal() +
  scale_fill_manual(values = c(`Induced lncRNAs` = "#D6604D", `Repressed lncRNAs` = "#67A9CF", `Non-DE lncRNAs` = "grey30")) +
  xlab("") +
  ylab("% With a repressed neighbour")

#just disconcordant:
ggplot(filter(DEL_PCG_type_plot, !grepl("Repress", variable))) + aes(x = RegulationBegins, y = value, fill = variable) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.7), color = "grey60") +
  theme_minimal() +
  theme(text = element_text(size = 32), legend.position = "none") +
  scale_fill_manual(values = c(`Induced lncRNAs` = "#D6604D", `Repressed lncRNAs` = "#67A9CF", `Non-DE lncRNAs` = "grey30")) +
  xlab("") +
  ylab("% With repressed\nneighbour")+
  theme(text = element_text(size=24)) + Seurat::RotatedAxis()

#just concordant:
ggplot(filter(DEL_PCG_type_plot, !grepl("Induce", variable))) + aes(x = RegulationBegins, y = value, fill = variable) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.7), color = "grey60") +
  theme_minimal() +
  theme(text = element_text(size = 32)) + #, legend.position = "none") +
  scale_fill_manual(values = c(`Induced lncRNAs` = "#D6604D", `Repressed lncRNAs` = "#67A9CF", `Non-DE lncRNAs` = "grey30")) +
  xlab("") +
  ylab("% With repressed\nneighbour") +
  theme(text = element_text(size=24)) + Seurat::RotatedAxis()



#### (2D version - simpler) compare FC amongst groupings, 0-4hrs co-induced/repressed ####

PCGDE_Upwithin_4 <- filter(fpkm_allGDE_Upwithin_4, EnsType == "protein_coding", grepl("coding|TF|CC", GeneClassUpdate))

#now only co-induced pairs
PCGDE_Upwithin_4$DEneighbourTypeIII <- "No"

#closest neighbours only:


CoInducedLncRNA_4hr <- filter(CoRegPairs_04_48_24_extended_naive, 
                              EnsID %in% fpkm_allGDE_Upwithin_4$EnsID, 
                              EnsID.y %in% fpkm_allGDE_Upwithin_4$EnsID)

PCGDE_Upwithin_4$DEneighbourTypeIII[PCGDE_Upwithin_4$EnsID %in% CoInducedLncRNA_4hr$EnsID.y] <- "Yes"

table(PCGDE_Upwithin_4$DEneighbourTypeIII)#45 PCG coinduced with a CCLnc in 4hrs, 315 are coinduced with a PCG

PCGDE_Upwithin_4$`Concordant lncRNA\nneighbour` <- PCGDE_Upwithin_4$DEneighbourTypeIII

ggplot(PCGDE_Upwithin_4) + aes(x = `Concordant lncRNA\nneighbour`, y = abs(LogFC_0_4), color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5, color = "grey50") +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "#D6604D", "No" = "#67A9CF"))+
  xlab("Induced PCGs")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Log2FC (0-4hrs)")

#test gets to ** with 400kbp
t.test(PCGDE_Upwithin_4$LogFC_0_4 ~ PCGDE_Upwithin_4$`Concordant lncRNA\nneighbour`, PCGDE_Upwithin_4, var.equal = T)
t.test(PCGDE_Upwithin_4$LogFC_0_4 ~ PCGDE_Upwithin_4$`Concordant lncRNA\nneighbour`, PCGDE_Upwithin_4, var.equal = F)
wilcox.test(PCGDE_Upwithin_4$LogFC_0_4 ~ PCGDE_Upwithin_4$`Concordant lncRNA\nneighbour`, PCGDE_Upwithin_4)


#repressed 
PCGDE_Downwithin_4 <- filter(fpkm_allGDE_Downwithin_4, EnsType == "protein_coding", grepl("coding|TF|CC", GeneClassUpdate))

#now only co-Repressed pairs
PCGDE_Downwithin_4$DEneighbourTypeIII <- "No"

CoRepressedLncRNA_4hr <- filter(CoRegPairs_04_48_24_extended_naive, 
                                EnsID %in% fpkm_allGDE_Downwithin_4$EnsID, 
                                EnsID.y %in% fpkm_allGDE_Downwithin_4$EnsID)

PCGDE_Downwithin_4$DEneighbourTypeIII[PCGDE_Downwithin_4$EnsID %in% CoRepressedLncRNA_4hr$EnsID.y] <- "Yes"

table(PCGDE_Downwithin_4$DEneighbourTypeIII)#17 PCG coRepressed with a CCLnc in 4hrs, 164 are coRepressed with a PCG

PCGDE_Downwithin_4$`Concordant lncRNA\nneighbour` <- PCGDE_Downwithin_4$DEneighbourTypeIII

ggplot(PCGDE_Downwithin_4) + aes(x = `Concordant lncRNA\nneighbour`, y = abs(LogFC_0_4), color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5, color = "grey50") +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "#D6604D", "No" = "#67A9CF")) +
  xlab("Repressed PCGs")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Absolute Log2FC (0-4hrs)")

t.test(PCGDE_Downwithin_4$LogFC_0_4 ~ PCGDE_Downwithin_4$`Concordant lncRNA\nneighbour`, PCGDE_Downwithin_4, var.equal = T)
wilcox.test(PCGDE_Downwithin_4$LogFC_0_4 ~ PCGDE_Downwithin_4$`Concordant lncRNA\nneighbour`, PCGDE_Downwithin_4)


#### (2D version - simpler) compare FC amongst groupings, 0-8hrs co-induced/repressed ####

PCGDE_Upwithin_8 <- filter(fpkm_allGDE_Upwithin_8, grepl("coding|TF|CC", GeneClassUpdate))

#now only co-induced pairs
PCGDE_Upwithin_8$DEneighbourTypeIII <- "No"

CoInducedLncRNA_4hr <- filter(CoRegPairs_04_48_24_extended_naive, 
                              EnsID %in% fpkm_allGDE_Upwithin_8$EnsID, 
                              EnsID.y %in% fpkm_allGDE_Upwithin_8$EnsID)

PCGDE_Upwithin_8$DEneighbourTypeIII[PCGDE_Upwithin_8$EnsID %in% CoInducedLncRNA_4hr$EnsID.y] <- "Yes"

table(PCGDE_Upwithin_8$DEneighbourTypeIII)#13 PCG coinduced with a CCLnc

PCGDE_Upwithin_8$`Concordant lncRNA\nneighbour` <- PCGDE_Upwithin_8$DEneighbourTypeIII

PCGDE_Upwithin_8_assess <- filter(PCGDE_Upwithin_8, LogFC_0_8 >log2(1.5), preadj_0_8 < 0.05)

ggplot(PCGDE_Upwithin_8_assess) + aes(x = `Concordant lncRNA\nneighbour`, y = abs(LogFC_0_8), color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5, color = "grey50") +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "#D6604D", "No" = "#67A9CF"))+
  xlab("Induced PCGs")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Log2FC (0-8hrs)")

wilcox.test(PCGDE_Upwithin_8_assess$LogFC_0_8 ~ PCGDE_Upwithin_8_assess$`Concordant lncRNA\nneighbour`, PCGDE_Upwithin_8_assess)


#repressed 
PCGDE_Downwithin_8 <- filter(fpkm_allGDE_Downwithin_8, grepl("coding|TF|CC", GeneClassUpdate))

#now only co-Repressed pairs
PCGDE_Downwithin_8$DEneighbourTypeIII <- "No"

CoRepressedLncRNA_4hr <- filter(CoRegPairs_04_48_24_extended_naive, 
                                EnsID %in% fpkm_allGDE_Downwithin_8$EnsID, 
                                EnsID.y %in% fpkm_allGDE_Downwithin_8$EnsID)

PCGDE_Downwithin_8$DEneighbourTypeIII[PCGDE_Downwithin_8$EnsID %in% CoRepressedLncRNA_4hr$EnsID.y] <- "Yes"

table(PCGDE_Downwithin_8$DEneighbourTypeIII)#12 PCG coRepressed with a CCLnc 

PCGDE_Downwithin_8$`Concordant lncRNA\nneighbour` <- PCGDE_Downwithin_8$DEneighbourTypeIII

PCGDE_Downwithin_8_assess <- filter(PCGDE_Downwithin_8, LogFC_0_8 < -log2(1.5), preadj_0_8 < 0.05)

ggplot(PCGDE_Downwithin_8_assess) + aes(x = `Concordant lncRNA\nneighbour`, y = abs(LogFC_0_8), color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5, color = "grey50") +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "#D6604D", "No" = "#67A9CF")) +
  xlab("Repressed PCGs")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Absolute Log2FC (0-8hrs)")

wilcox.test(PCGDE_Downwithin_8_assess$LogFC_0_8 ~ PCGDE_Downwithin_8_assess$`Concordant lncRNA\nneighbour`, PCGDE_Downwithin_8_assess)

#### (2D version - simpler) compare FC amongst groupings, 0-24hrs co-induced/repressed ####

PCGDE_Upwithin_24 <- filter(fpkm_allGDE_Upwithin_24, grepl("coding|TF|CC", GeneClassUpdate))

#now only co-induced pairs
PCGDE_Upwithin_24$DEneighbourTypeIII <- "No"

CoInducedLncRNA_4hr <- filter(CoRegPairs_04_48_24_extended_naive, 
                              EnsID %in% fpkm_allGDE_Upwithin_24$EnsID, 
                              EnsID.y %in% fpkm_allGDE_Upwithin_24$EnsID)

PCGDE_Upwithin_24$DEneighbourTypeIII[PCGDE_Upwithin_24$EnsID %in% CoInducedLncRNA_4hr$EnsID.y] <- "Yes"

table(PCGDE_Upwithin_24$DEneighbourTypeIII)#4 PCG coinduced with a CCLnc

PCGDE_Upwithin_24$`Concordant lncRNA\nneighbour` <- PCGDE_Upwithin_24$DEneighbourTypeIII

PCGDE_Upwithin_24_assess <- filter(PCGDE_Upwithin_24, LogFC_0_24 >log2(1.5), preadj_0_24 < 0.05)

ggplot(PCGDE_Upwithin_24_assess) + aes(x = `Concordant lncRNA\nneighbour`, y = abs(LogFC_0_24), color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5, color = "grey50") +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "#D6604D", "No" = "#67A9CF"))+
  xlab("Induced PCGs")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Log2FC (0-24hrs)")

wilcox.test(PCGDE_Upwithin_24_assess$LogFC_0_24 ~ PCGDE_Upwithin_24_assess$`Concordant lncRNA\nneighbour`, PCGDE_Upwithin_24_assess, var.equal = T)


#repressed 
PCGDE_Downwithin_24 <- filter(fpkm_allGDE_Downwithin_24, grepl("coding|TF|CC", GeneClassUpdate))

#now only co-Repressed pairs
PCGDE_Downwithin_24$DEneighbourTypeIII <- "No"

CoRepressedLncRNA_4hr <- filter(CoRegPairs_04_48_24_extended_naive, 
                                EnsID %in% fpkm_allGDE_Downwithin_24$EnsID, 
                                EnsID.y %in% fpkm_allGDE_Downwithin_24$EnsID)

PCGDE_Downwithin_24$DEneighbourTypeIII[PCGDE_Downwithin_24$EnsID %in% CoRepressedLncRNA_4hr$EnsID.y] <- "Yes"

table(PCGDE_Downwithin_24$DEneighbourTypeIII)#12 PCG coRepressed with a CCLnc 

PCGDE_Downwithin_24$`Concordant lncRNA\nneighbour` <- PCGDE_Downwithin_24$DEneighbourTypeIII

PCGDE_Downwithin_24_assess <- filter(PCGDE_Downwithin_24, LogFC_0_24 < -log2(1.5), preadj_0_24 < 0.05)

ggplot(PCGDE_Downwithin_24_assess) + aes(x = `Concordant lncRNA\nneighbour`, y = abs(LogFC_0_24), color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5, color = "grey50") +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "#D6604D", "No" = "#67A9CF")) +
  xlab("Repressed PCGs")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Absolute Log2FC (0-24hrs)")

wilcox.test(PCGDE_Downwithin_24_assess$LogFC_0_24 ~ PCGDE_Downwithin_24_assess$`Concordant lncRNA\nneighbour`, PCGDE_Downwithin_24_assess)


#### how much omitting by focusing on 0hr baseline? ####

#100% for 0-4hr

dim(filter(fpkm_allGDE_within_8, grepl("coding|TF|CC", GeneClassUpdate)))
dim(filter(fpkm_allGDE_within_8, grepl("coding|TF|CC", GeneClassUpdate) , preadj_0_8 <0.05, abs(LogFC_0_8) >log2(1.5)))
1240/1335 #93% being assessed

dim(filter(fpkm_allGDE_within_24, grepl("coding|TF|CC", GeneClassUpdate)))
dim(filter(fpkm_allGDE_within_24, grepl("coding|TF|CC", GeneClassUpdate) , preadj_0_24 <0.05, abs(LogFC_0_24) >log2(1.5)))
772/1167 #66% being assessed

#quite a lot! justifies a second look at 4-24hr and 8-24hr FCs?
dim(filter(fpkm_allGDE_within_24, grepl("coding|TF|CC", GeneClassUpdate) , preadj_4_24 <0.05, abs(LogFC_4_24) >log2(1.5)))
dim(filter(fpkm_allGDE_within_24, grepl("coding|TF|CC", GeneClassUpdate) , preadj_8_24 <0.05, abs(LogFC_8_24) >log2(1.5)))
873/1167 #75%
793/1167 #68%

#not for now, 0hr is a fair baseline

#### combined plotting ####

PCGDE_Upwithin_4$time <- "0-4hr"
PCGDE_Upwithin_8_assess$time <- "8hr"
PCGDE_Upwithin_24_assess$time <- "8-24hr Induced"

PCGDE_Upwithin_4$LFC_2_compare <- PCGDE_Upwithin_4$LogFC_0_4
PCGDE_Upwithin_8_assess$LFC_2_compare <- PCGDE_Upwithin_8_assess$LogFC_0_8
PCGDE_Upwithin_24_assess$LFC_2_compare <- PCGDE_Upwithin_24_assess$LogFC_0_24

trial <- rbind(PCGDE_Upwithin_4,
               PCGDE_Upwithin_8_assess,
               PCGDE_Upwithin_24_assess)

trial$time <- factor(trial$time)
trial$time <- factor(trial$time, levels = levels(trial$time)[c(2,3,1)])

ggplot(trial) + aes(x = `Concordant lncRNA\nneighbour`, y = LFC_2_compare , color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha =0.4) +
  facet_wrap(~time) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "mediumorchid", "No" = "grey60")) +
  xlab("\nInduced LncRNA Neighbour")+
  theme_minimal() +
  theme(#axis.text.x = element_blank(),
    strip.background = element_rect(),
    legend.position = "none") +
  ylab("Log2FC vs. 0hr")

#keep distinct y axes
trial$DEneighbourTypeIII[trial$DEneighbourTypeIII == "Yes"] <- "Co-induced\n with CCLncRNA"
trial$DEneighbourTypeIII[trial$DEneighbourTypeIII == "No"] <- "Other"

ggplot(filter(trial, time == "4hr Induced")) + aes(x = DEneighbourTypeIII, 
                                                   y = LFC_2_compare ,
                                                   color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5, size = 4) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha =0.4) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "mediumorchid", "No" = "grey60")) +
  xlab("")+
  theme_minimal() +
  theme(text = element_text(size = 28),
        legend.position = "none") +
  ylab("Log2FC 0-4hr")

ggplot(filter(trial, time == "8hr Induced")) + aes(x = DEneighbourTypeIII, 
                                                   y = LFC_2_compare ,
                                                   color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5, size = 4) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha =0.4) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "mediumorchid", "No" = "grey60")) +
  xlab("")+
  theme_minimal() +
  theme(text = element_text(size = 28),
        legend.position = "none") +
  ylab("Log2FC 0-8hr")

ggplot(filter(trial, time == "24hr Induced")) + aes(x = DEneighbourTypeIII, 
                                                    y = LFC_2_compare ,
                                                    color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5, size = 4) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha =0.4) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "mediumorchid", "No" = "grey60")) +
  xlab("")+
  theme_minimal() +
  theme(text = element_text(size = 28),
        legend.position = "none") +
  ylab("Log2FC 0-24hr")


PCGDE_Downwithin_4$time <- "4hr Repressed"
PCGDE_Downwithin_8_assess$time <- "8hr Repressed"
PCGDE_Downwithin_24_assess$time <- "24hr Repressed"

PCGDE_Downwithin_4$LFC_2_compare <- PCGDE_Downwithin_4$LogFC_0_4
PCGDE_Downwithin_8_assess$LFC_2_compare <- PCGDE_Downwithin_8_assess$LogFC_0_8
PCGDE_Downwithin_24_assess$LFC_2_compare <- PCGDE_Downwithin_24_assess$LogFC_0_24

trial <- rbind(PCGDE_Downwithin_4,
               PCGDE_Downwithin_8_assess,
               PCGDE_Downwithin_24_assess)

trial$time <- factor(trial$time)
trial$time <- factor(trial$time, levels = levels(trial$time)[c(2,3,1)])

ggplot(trial) + aes(x = `Concordant lncRNA\nneighbour`, y = abs(LFC_2_compare) , color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha =0.4) +
  facet_wrap(~time) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "mediumorchid", "No" = "grey60")) +
  xlab("\nInduced LncRNA Neighbour")+
  theme_minimal() +
  theme(#axis.text.x = element_blank(),
    strip.background = element_rect(),
    legend.position = "none") +
  ylab("Absolute Log2FC vs. 0hr")


#keep distinct y axes
trial$DEneighbourTypeIII[trial$DEneighbourTypeIII == "Yes"] <- "Co-repressed\nwith CClncRNA"
trial$DEneighbourTypeIII[trial$DEneighbourTypeIII == "No"] <- "Other"

ggplot(filter(trial, time == "4hr Repressed")) + aes(x = DEneighbourTypeIII, 
                                                     y = abs(LFC_2_compare) ,
                                                     color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5, size = 4) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha =0.4) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "mediumorchid", "No" = "grey60")) +
  xlab("")+
  theme_minimal() +
  theme(text = element_text(size = 28),
        legend.position = "none") +
  ylab("Absolute Log2FC 0-4hr")

ggplot(filter(trial, time == "8hr Repressed")) + aes(x = DEneighbourTypeIII, 
                                                     y = abs(LFC_2_compare) ,
                                                     color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5, size = 4) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha =0.4) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "mediumorchid", "No" = "grey60")) +
  xlab("")+
  theme_minimal() +
  theme(text = element_text(size = 28),
        legend.position = "none") +
  ylab("Absolute Log2FC 0-8hr")

ggplot(filter(trial, time == "24hr Repressed")) + aes(x = DEneighbourTypeIII, 
                                                      y = abs(LFC_2_compare) ,
                                                      color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5, size = 4) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha =0.4) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "mediumorchid", "No" = "grey60")) +
  xlab("")+
  theme_minimal() +
  theme(text = element_text(size = 28),
        legend.position = "none") +
  ylab("Absolute Log2FC 0-24hr")


#### GO/KEGG/REACTOME for co-induced targets of 0-4hr lncs ####

library(clusterProfiler)
library(org.Hs.eg.db)

#running this with "surrounding genes" seems to provide clearest sign of unique lnc co-regulation

#first background, anything setting them apart from early up DEGs generally?
EarlyCoInduced <- filter(AllLNC_AllPCG_1,
                         (EnsID %in% fpkm_allGDE_Upwithin_4$EnsID & 
                            EnsID.y %in% fpkm_allGDE_Upwithin_4$EnsID))

EarlyCoInduced$EnsID_merge.y <- gsub("\\.[0-9]*", "", EarlyCoInduced$EnsID.y)

fpkm_PCGDE_Upwithin_4 <- filter(fpkm_allGDE_Upwithin_4, grepl("protein_coding", EnsType), grepl("coding|TF|CC", GeneClassUpdate))
fpkm_PCGDE_Upwithin_4$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE_Upwithin_4$EnsID)

#23 genes note high p thresh - 0 terms
CoUpEarly_UpEarly_GO <- enrichGO(gene          = unique(EarlyCoInduced$EnsID_merge.y),
                                 universe      = unique(fpkm_PCGDE_Upwithin_4$EnsID_merge.y),
                                 keyType       = "ENSEMBL",
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "all",    
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.1,
                                 qvalueCutoff  = 0.1,
                                 readable      = TRUE)
CoUpEarly_UpEarly_GO#0 terms

fpkm_PCGDE <- filter(fpkm_allGDE, grepl("protein_coding", EnsType), grepl("coding|TF|CC", GeneClassUpdate))
fpkm_PCGDE$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE$EnsID)

#second background, anything setting apart from DEGs generally
CoUpEarly_DEGs_GO <- enrichGO(gene          = unique(EarlyCoInduced$EnsID_merge.y),
                              universe      = unique(fpkm_PCGDE$EnsID_merge.y),
                              keyType       = "ENSEMBL",
                              OrgDb         = org.Hs.eg.db,
                              ont           = "all",    
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.1,
                              qvalueCutoff  = 0.1,
                              readable      = TRUE)

#mostly cytokine/chemokine, zinc related, inf by CXCL/MT1 only... not good to use
#because want to see collective lnc effect, not just driven by 1x loci
CoUpEarly_DEGs_GO_df <- as.data.frame(CoUpEarly_DEGs_GO)
#DNA bending: HMGA2/FOXL1
#growth factors
#cytokines
#niceeee

#third background, anything setting apart from early EGs generally
fpkm_allGDE_within_4 <- rbind(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Downwithin_4)
fpkm_allGDE_within_8 <- rbind(fpkm_allGDE_Upwithin_8, fpkm_allGDE_Downwithin_8)
fpkm_allGDE_within_24 <- rbind(fpkm_allGDE_Upwithin_24, fpkm_allGDE_Downwithin_24)

#background of totally non-DE lncs overall whole timecourse, expressed in given timepoint:
fpkm_allG_04 <- filter(fpkm_allG, (Hour0_meanFPKM>1 | Hour4_meanFPKM>1), !EnsID %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                       fpkm_allGDE_within_24$EnsID))

#adjust MT/CXCL loci
fpkm_PCG <- filter(fpkm_allG_04, grepl("protein_coding", EnsType), grepl("coding|TF|CC", GeneClassUpdate))
fpkm_PCG$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCG$EnsID)

CoUpEarly_EGs_GO <- enrichGO(gene          = unique(EarlyCoInduced$EnsID_merge.y),
                             universe      = unique(fpkm_PCG$EnsID_merge.y),
                             keyType       = "ENSEMBL",
                             OrgDb         = org.Hs.eg.db,
                             ont           = "all",    
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.1,
                             qvalueCutoff  = 0.1,
                             readable      = TRUE)
#few more interesting ones
CoUpEarly_EGs_GO_df <- as.data.frame(CoUpEarly_EGs_GO)

#try with KEGG:
convertEnsEnt <- bitr(unique(fpkm_PCG$EnsID_merge.y), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

EarlyTargets_loci_EG_KEGG <- enrichKEGG(gene          = unique(filter(convertEnsEnt,ENSEMBL %in% EarlyCoInduced$EnsID_merge.y)$ENTREZID),
                                      universe      = unique(filter(convertEnsEnt,ENSEMBL %in% fpkm_PCG$EnsID_merge.y)$ENTREZID),
                                      pAdjustMethod = "BH",
                                      pvalueCutoff  = 0.1,
                                      qvalueCutoff  = 0.1)
EarlyTargets_loci_EG_KEGG_df <- as.data.frame(EarlyTargets_loci_EG_KEGG)
#mineral absorption, calcium signalling (calcification?)
#4x transcriptional misregulaiton in cancer
#4x cytokine-cytokine receptor
#4x JAK-STAT

#try REACTOME:
library(ReactomePA)
EarlyTargets_loci_EG_PATH <- enrichPathway(gene          = unique(filter(convertEnsEnt,ENSEMBL %in% EarlyCoInduced$EnsID_merge.y)$ENTREZID),
                                           universe      = unique(filter(convertEnsEnt,ENSEMBL %in% fpkm_PCG$EnsID_merge.y)$ENTREZID),
                                           organism = "human",
                                           pvalueCutoff = 0.1,
                                           qvalueCutoff  = 0.1,
                                           readable      = TRUE)
EarlyTargets_loci_EG_PATH_df <- as.data.frame(EarlyTargets_loci_EG_PATH)
#cytokine
#senescence


#### Blunter checks on TFs, IEGs, Cell cycle genes, SMC-biased genes ####

#GO/KEGG/REACTOME reveals little

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

GeneLists <- list("TFs"= unique(filter(fpkm_allG, grepl("TF", GeneClassUpdate))$EnsID),
                  "CRs" = unique(filter(fpkm_allG, EnsName %in% CRdb_dataB$CR)$EnsID),
                  "CC"= unique(filter(fpkm_allG, grepl("CC", GeneClassUpdate))$EnsID),
                  "IEGs"= unique(filter(fpkm_allG, EnsName %in% c(IEGs_hs[,4]))$EnsID),
                  "SMC"= unique(filter(fpkm_allG, gsub("\\.[0-9]*", "", EnsID) %in% FANT_S10_SMC_G)$EnsID),
                  "VSMC"= unique(filter(fpkm_allG, gsub("\\.[0-9]*", "", EnsID) %in% FANT_S10_VSMC_G)$EnsID)
                  )

results_list <- list()

#enrichment of TFs/CRs etc in lnc-targeted DE PCGs rather than PCGs
for (i in 1:length(GeneLists)){
  
  a <- length(unique(EarlyCoInduced$EnsID.y[EarlyCoInduced$EnsID.y %in% GeneLists[[i]] ]))
  b <- length(unique(EarlyCoInduced$EnsID.y))
  c <- length(unique(fpkm_PCGDE_Upwithin_4$EnsID[ fpkm_PCGDE_Upwithin_4$EnsID %in% GeneLists[[i]]]))
  d <- length(unique(fpkm_PCGDE_Upwithin_4$EnsID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate,
                                  paste(unique(EarlyCoInduced$EnsName.y[EarlyCoInduced$EnsID.y %in% GeneLists[[i]] ]), collapse = "/"))
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:8] <- c("pval", "OR", "EnsName")
#cut non-vasc SMC
GeneLists_enriched <- GeneLists_enriched[-5,]

GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "bonferroni")
GeneLists_enrichedSameT <- GeneLists_enriched

#same/later
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
GeneLists_enriched <- GeneLists_enriched[-4,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")

GeneLists_enrichedSameLaterT <- GeneLists_enriched


#plot TF, in isolation
a <- length(unique(EarlyTargets_loci$EnsID.y[ EarlyTargets_loci$EnsID.y %in% GeneLists[[1]] ]))
b <- length(unique(EarlyTargets_loci$EnsID.y))
c <- length(unique(fpkm_PCGDE4_loci$EnsID[ fpkm_PCGDE4_loci$EnsID %in% GeneLists[[1]]]))
d <- length(unique(fpkm_PCGDE4_loci$EnsID))

ai <- length(unique(AllTargets_loci$EnsID.y[ AllTargets_loci$EnsID.y %in% GeneLists[[1]] ]))
bi <- length(unique(AllTargets_loci$EnsID.y))
ci <- length(unique(fpkm_PCGDE_loci$EnsID[ fpkm_PCGDE_loci$EnsID %in% GeneLists[[1]]]))
di <- length(unique(fpkm_PCGDE_loci$EnsID))

DEL_PCG_type <- data.frame("CClncRNA co-regulated\nneighbours" = c(a/b*100, ai/bi*100), 
                           "Background DEGs" = c(c/d*100, ci/di*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("Same", "Same/Later"))
#DEL_PCG_type$NoDEL <- c(a, ai)

DEL_PCG_type <- melt(DEL_PCG_type)

DEL_PCG_type$variable <- gsub("\\.", " ", DEL_PCG_type$variable)
DEL_PCG_type$variable <- gsub(" regulated ", "-regulated\n", DEL_PCG_type$variable)

ggplot(DEL_PCG_type) + aes(x = RegulationBegins, fill = variable, y = value) +
  geom_bar(stat= "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c(`CClncRNA co-regulated\nneighbours` = "mediumorchid",`Background DEGs` = "grey30")) +
  #geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  xlab("Timeframe") +   ylab("\n% Transcription\nFactor") +
  theme_minimal() +
  theme(text = element_text(size=20)) + Seurat::RotatedAxis()

#alt approach, by lncRNA

#enrichment of TF-targeting lncs amongst 0-4hr DE lncs vs. 0-4hr exprs lncs
fpkm_DEL_within_4 <- filter(fpkm_allGDE_within_4, grepl("fide|Lnc", GeneClassUpdate))
fpkm_EL_within_4 <- filter(fpkm_allG_04, grepl("fide|Lnc", GeneClassUpdate))
EarlyExLNC_AllPCG_1 <- filter(AllLNC_AllPCG_1, EnsID %in% fpkm_EL_within_4$EnsID)

for (i in 1:length(GeneLists)){
  
  a <- length(unique(EarlyTargets$EnsID[ EarlyTargets$EnsID.y %in% GeneLists[[i]] ]))
  b <- length(unique(fpkm_DEL_within_4$EnsID))
  c <- length(unique(EarlyExLNC_AllPCG_1$EnsID[EarlyExLNC_AllPCG_1$EnsID.y %in% fpkm_allGDE_within_4$EnsID &
                                                 EarlyExLNC_AllPCG_1$EnsID.y %in% GeneLists[[i]] ]))
  d <- length(unique(fpkm_EL_within_4$EnsID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
  
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:7] <- c("pval", "OR")
#cut non-vasc SMC
GeneLists_enriched <- GeneLists_enriched[-5,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")

GeneLists_enrichedSame_byLnc <- GeneLists_enriched

#e.g. for TFs
25/91 #27% of 0-4hr DELs have a DE TF nearby
67/312 #21% of 0-4 ELs 
#but ns

#same/later
for (i in 1:length(GeneLists)){
  
  a <- length(unique(AllTargets$EnsID[ AllTargets$EnsID.y %in% GeneLists[[i]] ]))
  b <- length(unique(fpkm_DEL_within_4$EnsID))
  c <- length(unique(EarlyExLNC_AllPCG_1$EnsID[EarlyExLNC_AllPCG_1$EnsID.y %in% fpkm_allGDE$EnsID &
                                                 EarlyExLNC_AllPCG_1$EnsID.y %in% GeneLists[[i]] ]))
  d <- length(unique(fpkm_EL_within_4$EnsID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
  
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:7] <- c("pval", "OR")
#cut non-vasc SMC
GeneLists_enriched <- GeneLists_enriched[-5,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")

GeneLists_enrichedSameLater_byLnc <- GeneLists_enriched
#ns


#### old code ####
#### GO for co-repressed targets of 0-4hr lncs ####

#first background, anything setting them apart from early Down DEGs generally?
EarlyCoRepressed <- filter(CoRegPairs_04_48_24_extended_naive,
                           (EnsID %in% fpkm_allGDE_Downwithin_4$EnsID & 
                              EnsID.y %in% fpkm_allGDE_Downwithin_4$EnsID))

#HOX cluster will dominate
EarlyCoRepressed$EnsID_merge.y <- gsub("\\.[0-9]*", "", EarlyCoRepressed$EnsID.y)
EarlyCoRepressed_loci1 <- filter(EarlyCoRepressed, grepl("^HOXA10|^HOXC6", EnsName.y))
EarlyCoRepressed_loci <- filter(EarlyCoRepressed, !grepl("^HOXA|^HOXC", EnsName.y))
EarlyCoRepressed_loci <- rbind(EarlyCoRepressed_loci, EarlyCoRepressed_loci1)

fpkm_PCGDE_Downwithin_4 <- filter(fpkm_allGDE_Downwithin_4, grepl("protein_coding", EnsType))
fpkm_PCGDE_Downwithin_4$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE_Downwithin_4$EnsID)
fpkm_PCGDE_Downwithin_4_loci1 <- filter(fpkm_PCGDE_Downwithin_4, grepl("^HOXA10|^HOXC6", EnsName))
fpkm_PCGDE_Downwithin_4_loci <- filter(fpkm_PCGDE_Downwithin_4, !grepl("^HOXA|^HOXC", EnsName))
fpkm_PCGDE_Downwithin_4_loci <- rbind(fpkm_PCGDE_Downwithin_4_loci, fpkm_PCGDE_Downwithin_4_loci1)

#note high p thresh
CoDownEarly_DownEarly_GO <- enrichGO(gene          = unique(EarlyCoRepressed_loci$EnsID_merge.y),
                                     universe      = unique(fpkm_PCGDE_Downwithin_4_loci$EnsID_merge.y),
                                     keyType       = "ENSEMBL",
                                     OrgDb         = org.Hs.eg.db,
                                     ont           = "all",    
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.1,
                                     qvalueCutoff  = 0.1,
                                     readable      = TRUE)
#0 terms
CoDownEarly_DownEarly_GO_df <- as.data.frame(CoDownEarly_DownEarly_GO)

fpkm_PCGDE <- filter(fpkm_allGDE, grepl("protein_coding", EnsType))
fpkm_PCGDE$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE$EnsID)
fpkm_PCGDE_loci1 <- filter(fpkm_PCGDE, grepl("^HOXA10|^HOXC6", EnsName))
fpkm_PCGDE_loci <- filter(fpkm_PCGDE, !grepl("^HOXA|^HOXC", EnsName))
fpkm_PCGDE_loci <- rbind(fpkm_PCGDE_loci, fpkm_PCGDE_loci1)

#second background, anything setting apart from DEGs generally
CoDownEarly_DEGs_GO <- enrichGO(gene          = unique(EarlyCoRepressed_loci$EnsID_merge.y),
                                universe      = unique(fpkm_PCGDE_loci$EnsID_merge.y),
                                keyType       = "ENSEMBL",
                                OrgDb         = org.Hs.eg.db,
                                ont           = "all",    
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.1,
                                qvalueCutoff  = 0.1,
                                readable      = TRUE)
#picking up on a developmental signature collectively being targeted, HOXA, HOXC and associated genes
CoDownEarly_DEGs_GO_df <- as.data.frame(CoDownEarly_DEGs_GO)


#third background, anything setting apart from EGs generally
fpkm_PCG <- filter(fpkm_allG_04, grepl("protein_coding", EnsType))
fpkm_PCG$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCG$EnsID)
fpkm_PCG_loci1 <- filter(fpkm_PCG, grepl("^HOXA10|^HOXC6", EnsName))
fpkm_PCG_loci <- filter(fpkm_PCG, !grepl("^HOXA|^HOXC", EnsName))
fpkm_PCG_loci <- rbind(fpkm_PCG_loci, fpkm_PCG_loci1)

CoDownEarly_EGs_GO <- enrichGO(gene          = unique(EarlyCoRepressed_loci$EnsID_merge.y),
                               universe      = unique(fpkm_PCG_loci$EnsID_merge.y),
                               keyType       = "ENSEMBL",
                               OrgDb         = org.Hs.eg.db,
                               ont           = "all",    
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.1,
                               qvalueCutoff  = 0.1,
                               readable      = TRUE)
CoDownEarly_EGs_GO_df <- as.data.frame(CoDownEarly_EGs_GO)
#as above
#no TF sig if correcting the locus issue

#no background
CoDownEarly_none_GO <- enrichGO(gene          = unique(EarlyCoRepressed_loci$EnsID_merge.y),
                               #universe      = unique(fpkm_PCG_loci$EnsID_merge.y),
                               keyType       = "ENSEMBL",
                               OrgDb         = org.Hs.eg.db,
                               ont           = "all",    
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.1,
                               qvalueCutoff  = 0.1,
                               readable      = TRUE)
CoDownEarly_none_GO_df <- as.data.frame(CoDownEarly_none_GO)
#similar effect

#### GO for concordant targets of 0-4hr lncs ####

fpkm_PCG <- filter(fpkm_allG_04, grepl("protein_coding", EnsType))
fpkm_PCG$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCG$EnsID)
fpkm_PCG_loci1 <- filter(fpkm_PCG, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName))
fpkm_PCG_loci <- filter(fpkm_PCG, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName))
fpkm_PCG_loci <- rbind(fpkm_PCG_loci, fpkm_PCG_loci1)

#joint testing of concordant early targets - increased power for e.g. a TF signal from joint up and down?
ConcEarly_EGs_GO <- enrichGO(gene          = unique(c(unique(EarlyCoRepressed_loci$EnsID_merge.y),
                                                      unique(EarlyCoInduced_loci$EnsID_merge.y))),
                             universe      = unique(fpkm_PCG_loci$EnsID_merge.y),
                             keyType       = "ENSEMBL",
                             OrgDb         = org.Hs.eg.db,
                             ont           = "all",    
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.1,
                             qvalueCutoff  = 0.1,
                             readable      = TRUE)
ConcEarly_EGs_GO_df <- as.data.frame(ConcEarly_EGs_GO)
#not worth using

ConcEarly_none_GO <- enrichGO(gene          = unique(c(unique(EarlyCoRepressed_loci$EnsID_merge.y),
                                                      unique(EarlyCoInduced_loci$EnsID_merge.y))),
                             #universe      = unique(fpkm_PCG_loci$EnsID_merge.y),
                             keyType       = "ENSEMBL",
                             OrgDb         = org.Hs.eg.db,
                             ont           = "all",    
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.1,
                             qvalueCutoff  = 0.1,
                             readable      = TRUE)
ConcEarly_none_GO_df <- as.data.frame(ConcEarly_none_GO)


#### GO/KEGG/REACTOME for same-timeframe targets of 0-4hr lncs ####

#now all early targets of early lncs - collective effect
fpkm_allGDE_within_4 <- rbind(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Downwithin_4)

#first background, anything setting them apart from early up DEGs generally?
EarlyTargets <- filter(AllLNC_AllPCG_1,
                       (EnsID %in% fpkm_allGDE_within_4$EnsID & 
                          EnsID.y %in% fpkm_allGDE_within_4$EnsID))

EarlyTargets$EnsID_merge.y <- gsub("\\.[0-9]*", "", EarlyTargets$EnsID.y)

#one per locus for MT1, CXCL, HOXA:
EarlyTargets_loci1 <- filter(EarlyTargets, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName.y))
EarlyTargets_loci <- filter(EarlyTargets, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName.y))
EarlyTargets_loci <- rbind(EarlyTargets_loci, EarlyTargets_loci1)

fpkm_PCGDE4 <- filter(fpkm_allGDE_within_4, grepl("protein_coding", EnsType))
fpkm_PCGDE4$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE4$EnsID)
fpkm_PCGDE4_loci1 <- filter(fpkm_PCGDE4, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName))
fpkm_PCGDE4_loci <- filter(fpkm_PCGDE4, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName))
fpkm_PCGDE4_loci <- rbind(fpkm_PCGDE4_loci, fpkm_PCGDE4_loci1)

#note high p thresh
EarlyTargets_GO <- enrichGO(gene          = unique(EarlyTargets_loci$EnsID_merge.y),
                            universe      = unique(fpkm_PCGDE4_loci$EnsID_merge.y),
                            keyType       = "ENSEMBL",
                            OrgDb         = org.Hs.eg.db,
                            ont           = "all",    
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.1,
                            qvalueCutoff  = 0.1,
                            readable      = TRUE)
#0 terms with locus simplification (previously TF heavy sig. but multi-target loci dominating these)
EarlyTargets_GO_df <- as.data.frame(EarlyTargets_GO)
#0 terms if all genes too


#looser background
fpkm_PCG <- filter(fpkm_allG_04, grepl("protein_coding", EnsType))
fpkm_PCG$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCG$EnsID)
fpkm_PCG_loci1 <- filter(fpkm_PCG, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName))
fpkm_PCG_loci <- filter(fpkm_PCG, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName))
fpkm_PCG_loci <- rbind(fpkm_PCG_loci, fpkm_PCG_loci1)

EarlyTargets_loci_EG_GO <- enrichGO(gene       = unique(EarlyTargets_loci$EnsID_merge.y),
                                 universe      = unique(fpkm_PCG_loci$EnsID_merge.y),
                                 keyType       = "ENSEMBL",
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "all",    
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.1,
                                 qvalueCutoff  = 0.1,
                                 readable      = TRUE)
#0 terms
EarlyTargets_loci_EG_GO_df <- as.data.frame(EarlyTargets_loci_EG_GO)

#n.b. for comparison, without locus correction here:
#EarlyTargets_loci_EG_GO <- enrichGO(gene       = unique(EarlyTargets$EnsID_merge.y),
#                                    universe      = unique(fpkm_PCG$EnsID_merge.y),
#                                    keyType       = "ENSEMBL",
#                                    OrgDb         = org.Hs.eg.db,
#                                    ont           = "all",    
#                                    pAdjustMethod = "BH",
#                                    pvalueCutoff  = 0.1,
#                                    qvalueCutoff  = 0.1,
#                                    readable      = TRUE)
#more terms, but all CXCL, MT1, HOX dominated
#EarlyTargets_loci_EG_GO_df <- as.data.frame(EarlyTargets_loci_EG_GO)


#try with KEGG:
convertEnsEnt <- bitr(unique(fpkm_PCG_loci$EnsID_merge.y), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

EarlyTargets_loci_EG_KEGG <- enrichGO(gene          = unique(filter(convertEnsEnt,ENSEMBL %in% EarlyTargets_loci$EnsID_merge.y)$ENTREZID),
                                    universe      = unique(filter(convertEnsEnt,ENSEMBL %in% fpkm_PCG_loci$EnsID_merge.y)$ENTREZID),
                                    keyType       = "ENSEMBL",
                                    OrgDb         = org.Hs.eg.db,
                                    ont           = "all",    
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.1,
                                    qvalueCutoff  = 0.1,
                                    readable      = TRUE)
EarlyTargets_loci_EG_KEGG_df <- as.data.frame(EarlyTargets_loci_EG_KEGG)

#try REACTOME:
library(ReactomePA)
EarlyTargets_loci_EG_PATH <- enrichPathway(gene          = unique(filter(convertEnsEnt,ENSEMBL %in% EarlyTargets_loci$EnsID_merge.y)$ENTREZID),
                                           universe      = unique(filter(convertEnsEnt,ENSEMBL %in% fpkm_PCG_loci$EnsID_merge.y)$ENTREZID),
                                           organism = "human",
                                           pvalueCutoff = 0.1,
                                           qvalueCutoff  = 0.1,
                                           readable      = TRUE)
EarlyTargets_loci_EG_PATH_df <- as.data.frame(EarlyTargets_loci_EG_PATH)
#2x interleukin/cytokine signalling enrichment terms


#code to get a figure left here, no longer needed but can cop/paste elsewhere
#View(CoUpEarly_EGs_GO_df[[i]])
EarlyTargets_loci_EG_PATH_df$selectHits <- as.numeric(sapply(strsplit(EarlyTargets_loci_EG_PATH_df$GeneRatio, "\\/"), "[[", 1))
EarlyTargets_loci_EG_PATH_df$select <- as.numeric(sapply(strsplit(EarlyTargets_loci_EG_PATH_df$GeneRatio, "\\/"), "[[", 2))
EarlyTargets_loci_EG_PATH_df$geneRatio <- EarlyTargets_loci_EG_PATH_df$selectHits/EarlyTargets_loci_EG_PATH_df$select*100

EarlyTargets_loci_EG_PATH_df <- EarlyTargets_loci_EG_PATH_df[order(EarlyTargets_loci_EG_PATH_df$Description),]

EarlyTargets_loci_EG_PATH_df$DescriptionII <- stringr::str_wrap(EarlyTargets_loci_EG_PATH_df$Description, width = 40)

EarlyTargets_loci_EG_PATH_df$Description <- factor(EarlyTargets_loci_EG_PATH_df$Description, labels = EarlyTargets_loci_EG_PATH_df$DescriptionII)
EarlyTargets_loci_EG_PATH_df$Description <- factor(EarlyTargets_loci_EG_PATH_df$Description,
                                          levels = levels(EarlyTargets_loci_EG_PATH_df$Description)[order(EarlyTargets_loci_EG_PATH_df$geneRatio, decreasing = F)])

EarlyTargets_loci_EG_PATH_df <- EarlyTargets_loci_EG_PATH_df[order(-EarlyTargets_loci_EG_PATH_df$geneRatio, EarlyTargets_loci_EG_PATH_df$p.adjust),]

EarlyTargets_loci_EG_PATH_df <- EarlyTargets_loci_EG_PATH_df[,c(9,2,12,15)]

colnames(EarlyTargets_loci_EG_PATH_df)[4] <- "% of co-regulated neighbours\nin REACTOME pathway"

ggplot(filter(EarlyTargets_loci_EG_PATH_df, p.adjust <0.05)) + aes(x = `% of co-regulated neighbours\nin REACTOME pathway`, 
                                                                   y = Description, 
                                                                   color = -log10(p.adjust), 
                                                                   size = Count) +
  geom_point() +
  theme_minimal()
#top terms were related to just 3x loci, MT1, MT2 and CXCL
#simplifying shows the themes bit better - including epigenetic reprogramming quite far down


#final background:
EarlyTargets_loci_none_GO <- enrichGO(gene          = unique(EarlyTargets_loci$EnsID_merge.y),
                                    #universe      = unique(fpkm_PCG_loci$EnsID_merge.y),
                                    keyType       = "ENSEMBL",
                                    OrgDb         = org.Hs.eg.db,
                                    ont           = "all",    
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.1,
                                    qvalueCutoff  = 0.1,
                                    readable      = TRUE)
#pre-ribosome again
EarlyTargets_loci_none_GO_df <- as.data.frame(EarlyTargets_loci_none_GO)


#### GO/KEGG/REACTOME for same-later timeframe targets of 0-4hr lncs ####

#wider signature, all same/later targets of early lncs?:
AllTargets <- unique(filter(CoRegPairs_04_48_24_extended_naive,
                            EnsID %in% fpkm_allGDE_within_4$EnsID))
AllTargets$EnsID_merge.y <- gsub("\\.[0-9]*", "", AllTargets$EnsID.y)
AllTargets_loci1 <- filter(AllTargets, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName.y))
AllTargets_loci <- filter(AllTargets, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName.y))
AllTargets_loci <- rbind(AllTargets_loci, AllTargets_loci1)

fpkm_PCGDE <- filter(fpkm_allGDE, grepl("protein_coding", EnsType))
fpkm_PCGDE$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE$EnsID)
fpkm_PCGDE_loci1 <- filter(fpkm_PCGDE, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName))
fpkm_PCGDE_loci <- filter(fpkm_PCGDE, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName))
fpkm_PCGDE_loci <- rbind(fpkm_PCGDE_loci, fpkm_PCGDE_loci1)

#note high p thresh
AllTargets_GO <- enrichGO(gene          = unique(AllTargets_loci$EnsID_merge.y),
                            universe      = unique(fpkm_PCGDE_loci$EnsID_merge.y),
                            keyType       = "ENSEMBL",
                            OrgDb         = org.Hs.eg.db,
                            ont           = "all",    
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.1,
                            qvalueCutoff  = 0.1,
                            readable      = TRUE)
#nothing
AllTargets_GO_df <- as.data.frame(AllTargets_GO)

#looser background
fpkm_PCG <- filter(fpkm_allG, grepl("protein_coding", EnsType))
fpkm_PCG$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCG$EnsID)
fpkm_PCG_loci1 <- filter(fpkm_PCG, grepl("CXCL8|^MT1E|^HOXA10|^HOXC6", EnsName))
fpkm_PCG_loci <- filter(fpkm_PCG, !grepl("^CXCL|^MT1|^HOXA|^HOXC", EnsName))
fpkm_PCG_loci <- rbind(fpkm_PCG_loci, fpkm_PCG_loci1)

#note high p thresh
AllTargets_GO <- enrichGO(gene          = unique(AllTargets_loci$EnsID_merge.y),
                          universe      = unique(fpkm_PCG_loci$EnsID_merge.y),
                          keyType       = "ENSEMBL",
                          OrgDb         = org.Hs.eg.db,
                          ont           = "all",    
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.15,
                          qvalueCutoff  = 0.15,
                          readable      = TRUE)
#2x dds breaks
AllTargets_GO_df <- as.data.frame(AllTargets_GO)

AllTargets_loci_EG_PATH <- enrichPathway(gene          = unique(filter(convertEnsEnt,ENSEMBL %in% AllTargets_loci$EnsID_merge.y)$ENTREZID),
                                           universe      = unique(filter(convertEnsEnt,ENSEMBL %in% fpkm_PCG_loci$EnsID_merge.y)$ENTREZID),
                                           organism = "human",
                                           pvalueCutoff = 0.1,
                                           qvalueCutoff  = 0.1,
                                           readable      = TRUE)
AllTargets_loci_EG_PATH_df <- as.data.frame(AllTargets_loci_EG_PATH)
#IL10 again


#### Blunter checks on TFs, IEGs, Cell cycle genes, SMC-biased genes ####

#GO/KEGG/REACTOME reveals little

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

GeneLists <- list("TFs"= unique(filter(fpkm_allG, grepl("TF", GeneClassUpdate))$EnsID),
                  "CRs" = unique(filter(fpkm_allG, EnsName %in% CRdb_dataB$CR | EnsID %in% extra_CR)$EnsID),
     "CC"= unique(filter(fpkm_allG, grepl("CC", GeneClassUpdate))$EnsID),
     "IEGs"= unique(filter(fpkm_allG, EnsName %in% c(IEGs_hs[,4]) | EnsID %in% extra_IEG)$EnsID),
     "SMC"= unique(filter(fpkm_allG, gsub("\\.[0-9]*", "", EnsID) %in% FANT_S10_SMC_G | EnsID %in% extra_SMC)$EnsID),
     "VSMC"= unique(filter(fpkm_allG, gsub("\\.[0-9]*", "", EnsID) %in% FANT_S10_VSMC_G | EnsID %in% extra_VSMC)$EnsID)
     )

results_list <- list()

#enrichment of TFs in lnc-targeted DE PCGs rather than PCGs
for (i in 1:length(GeneLists)){
  
  a <- length(unique(EarlyTargets_loci$EnsID.y[ EarlyTargets_loci$EnsID.y %in% GeneLists[[i]] ]))
  b <- length(unique(EarlyTargets_loci$EnsID.y))
  c <- length(unique(fpkm_PCGDE4_loci$EnsID[ fpkm_PCGDE4_loci$EnsID %in% GeneLists[[i]]]))
  d <- length(unique(fpkm_PCGDE4_loci$EnsID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate,
                         paste(unique(EarlyTargets_loci$EnsName.y[ EarlyTargets_loci$EnsID.y %in% GeneLists[[i]] ]), collapse = "/"))
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:7] <- c("pval", "OR")
#cut non-vasc SMC
GeneLists_enriched <- GeneLists_enriched[-4,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "bonferroni")

GeneLists_enrichedSameT <- GeneLists_enriched

#same/later
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
GeneLists_enriched <- GeneLists_enriched[-4,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")

GeneLists_enrichedSameLaterT <- GeneLists_enriched


#plot TF, in isolation
a <- length(unique(EarlyTargets_loci$EnsID.y[ EarlyTargets_loci$EnsID.y %in% GeneLists[[1]] ]))
b <- length(unique(EarlyTargets_loci$EnsID.y))
c <- length(unique(fpkm_PCGDE4_loci$EnsID[ fpkm_PCGDE4_loci$EnsID %in% GeneLists[[1]]]))
d <- length(unique(fpkm_PCGDE4_loci$EnsID))

ai <- length(unique(AllTargets_loci$EnsID.y[ AllTargets_loci$EnsID.y %in% GeneLists[[1]] ]))
bi <- length(unique(AllTargets_loci$EnsID.y))
ci <- length(unique(fpkm_PCGDE_loci$EnsID[ fpkm_PCGDE_loci$EnsID %in% GeneLists[[1]]]))
di <- length(unique(fpkm_PCGDE_loci$EnsID))

DEL_PCG_type <- data.frame("CClncRNA co-regulated\nneighbours" = c(a/b*100, ai/bi*100), 
                           "Background DEGs" = c(c/d*100, ci/di*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("Same", "Same/Later"))
#DEL_PCG_type$NoDEL <- c(a, ai)

DEL_PCG_type <- melt(DEL_PCG_type)

DEL_PCG_type$variable <- gsub("\\.", " ", DEL_PCG_type$variable)
DEL_PCG_type$variable <- gsub(" regulated ", "-regulated\n", DEL_PCG_type$variable)

ggplot(DEL_PCG_type) + aes(x = RegulationBegins, fill = variable, y = value) +
  geom_bar(stat= "identity", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c(`CClncRNA co-regulated\nneighbours` = "mediumorchid",`Background DEGs` = "grey30")) +
  #geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  xlab("Timeframe") +   ylab("\n% Transcription\nFactor") +
  theme_minimal() +
  theme(text = element_text(size=20)) + Seurat::RotatedAxis()

#alt approach, by lncRNA

#enrichment of TF-targeting lncs amongst 0-4hr DE lncs vs. 0-4hr exprs lncs
fpkm_DEL_within_4 <- filter(fpkm_allGDE_within_4, grepl("fide|Lnc", GeneClassUpdate))
fpkm_EL_within_4 <- filter(fpkm_allG_04, grepl("fide|Lnc", GeneClassUpdate))
EarlyExLNC_AllPCG_1 <- filter(AllLNC_AllPCG_1, EnsID %in% fpkm_EL_within_4$EnsID)

for (i in 1:length(GeneLists)){
  
  a <- length(unique(EarlyTargets$EnsID[ EarlyTargets$EnsID.y %in% GeneLists[[i]] ]))
  b <- length(unique(fpkm_DEL_within_4$EnsID))
  c <- length(unique(EarlyExLNC_AllPCG_1$EnsID[EarlyExLNC_AllPCG_1$EnsID.y %in% fpkm_allGDE_within_4$EnsID &
                                               EarlyExLNC_AllPCG_1$EnsID.y %in% GeneLists[[i]] ]))
  d <- length(unique(fpkm_EL_within_4$EnsID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)

  }

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:7] <- c("pval", "OR")
#cut non-vasc SMC
GeneLists_enriched <- GeneLists_enriched[-5,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")

GeneLists_enrichedSame_byLnc <- GeneLists_enriched

#e.g. for TFs
25/91 #27% of 0-4hr DELs have a DE TF nearby
67/312 #21% of 0-4 ELs 
#but ns

#same/later
for (i in 1:length(GeneLists)){
  
  a <- length(unique(AllTargets$EnsID[ AllTargets$EnsID.y %in% GeneLists[[i]] ]))
  b <- length(unique(fpkm_DEL_within_4$EnsID))
  c <- length(unique(EarlyExLNC_AllPCG_1$EnsID[EarlyExLNC_AllPCG_1$EnsID.y %in% fpkm_allGDE$EnsID &
                                                 EarlyExLNC_AllPCG_1$EnsID.y %in% GeneLists[[i]] ]))
  d <- length(unique(fpkm_EL_within_4$EnsID))
  
  results_list[[i]] <- data.frame(a,b,c,d, 
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                                  fisher.test(data.frame("DEL" = c(a,b-a),
                                                         "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
  
}

names(results_list) <- names(GeneLists)
GeneLists_enriched <- bind_rows(results_list, .id = "GeneList")
colnames(GeneLists_enriched)[6:7] <- c("pval", "OR")
#cut non-vasc SMC
GeneLists_enriched <- GeneLists_enriched[-5,]
GeneLists_enriched$pval_adj <- p.adjust(GeneLists_enriched$pval, method = "BH")

GeneLists_enrichedSameLater_byLnc <- GeneLists_enriched
#ns

#### STRING/enrichR export ####
#for enrichR
EarlyCoInduced_con <- bitr(unique(EarlyCoInduced$EnsID_merge.y), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
#write.csv(unique(EarlyCoInduced_con$ENTREZID), "EarlyCoInduced_entrez.csv", row.names = F)

EarlyCoRepressed_con <- bitr(unique(EarlyCoRepressed$EnsID_merge.y), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
#write.csv(unique(EarlyCoRepressed_con$ENTREZID), "EarlyCoRepressed_entrez.csv", row.names = F)
#nothing seperately or together

#STRING:
EarlyTargets <- filter(AllLNC_AllPCG_1,
                       EnsID %in% fpkm_allGDE_within_4$EnsID, 
                          EnsID.y %in% fpkm_allGDE_within_4$EnsID)

EarlyTargets$EnsID_merge.y <- gsub("\\.[0-9]*", "", EarlyTargets$EnsID.y)

write.csv(unique(EarlyTargets$EnsName.y), "earlyTargets_STRING.csv", row.names = F)

