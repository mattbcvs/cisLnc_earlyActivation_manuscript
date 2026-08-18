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


#### bias of lncRNAs ####

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


#### define cclncRNAs - closest neighbour approach - enhancer association ####

#stay stringent, close pairs only (Genehancer paper indicates most cis pairing associations in this range)
#filter on distance, 250kbp initially done, use higher if backed up by HiC or eQTL
AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_info_2026.csv", header = T)
Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime_2026.csv", header = T)
AllLNC_AllPCG_1$GeneClassUpdate.x[AllLNC_AllPCG_1$EnsID %in% 
                                    filter(Enhancer_lociII, !is.na(EnhancerVerdict))$EnsID & 
                                    AllLNC_AllPCG_1$GeneClassUpdate.x == "Bona fide lncRNA"] <- "ELnc"

#begin finding closest neighbour reguation events
closestNeighbour <- filter(AllLNC_AllPCG_1, AbsDistLnc_PCG <250)
closestNeighbour <- split(closestNeighbour, closestNeighbour$EnsID)

#closest "surrounding" neighbours (consider up and down-stream):
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

CCLncNaiveWaveBias
sum(CCLncNaiveWaveBias$Freq)
 c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

#250kbp surrounding
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

#250kbp surrounding
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
                                           rep("Spacer4",3)
), 
"FR" = rep(c("Within \n4hrs",
             "Within \n8hrs",
             "Within \n24hrs"), 4), 
"OR" = NA,  
`padj_simple` =NA)

colnames(appendSpacers) <- colnames(AllTypesWaveBias)

AllTypesWaveBias <- rbind(AllTypesWaveBias, appendSpacers)

AllTypesWaveBias$UpDownType <- factor(AllTypesWaveBias$UpDownType)
AllTypesWaveBias$UpDownType <- factor(AllTypesWaveBias$UpDownType, 
                                      levels = levels(AllTypesWaveBias$UpDownType)[c(2,1,9,4,3,10,6,5,11,8,7)])


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

#bias vs de lncs -250kbp surrounding neighbours
aii <- 24    
bii <- 51    
cii <- 65   
dii <- 221   

fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

CCLncWaveBias
sum(CCLncWaveBias$Freq)

#bias vs de lncs -250kbp surrounding neighbours
aii <- 37
bii <- 81    
cii <- 65   
dii <- 221   

fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

#250kbp, surrounding neighbours
trial <- data.frame(65/221*100, 24/51*100, 37/81*100)

colnames(trial) <- c("All DE lncRNAs", "CClncRNAs (same)",  "CClncRNAs (same/later)")

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
  
  c <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
                            EnsID %in% filter(AllLNC_AllPCG_1, EnsID.y %in% cisPotentialTargets)$EnsID)$EnsID))
  
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
                            EnsID %in% selection_list[[i]]$EnsID)$EnsID))
  
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

p.adjust(DEL_PCG_type_bi$p, method = "bonferroni")


#just down DELs:
selection_list <- list(fpkm_allGDE_Downwithin_4, fpkm_allGDE_Downwithin_8, fpkm_allGDE_Downwithin_24)
#with down targets:
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
g2bii <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "#67A9CF") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  xlab("") +   ylab("\n% With a repressed\nneighbour") +
  theme_minimal()

DEL_PCG_type_bii <- DEL_PCG_type
DEL_PCG_type_bii$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type_bii$OR <- unlist(results_list)[c(6,12,18)]
#write.csv(DEL_PCG_type_bii, "SVSMC_repressedDEL_EL_DEneighboursame_conc.csv", row.names = F)

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

DEL_PCG_type_bi_disc <- DEL_PCG_type
DEL_PCG_type_bi_disc$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type_bi_disc$OR <- unlist(results_list)[c(6,12,18)]

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


