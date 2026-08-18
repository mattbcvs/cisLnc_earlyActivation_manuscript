#### including PCG-PCG as comparison group ####
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(rcompanion)
library(ggbeeswarm)
library(rtracklayer)
library(GenomicRanges)

#seperate DE PCGs into key groups, 
#a) CCLnc target (done)
#b) co-regulated with a lncRNA (done)
#c) co-regulated with a PCG (to do here)
#d) not a neighbour of a regulated gene (any remaining)

#previously also including PCGs which are co-regulated/correlated with a co-regulated PCG as a control group
#i.e. I followed the process to predict cis-acting lncRNAs for PCGs as well, to use as a control group
#this is possible but would need to have the eQTL analysis to compare properly - big files for a large number of gene pairs probably

#look to do co-regulation only for now - revisit and discuss...

#tests will be:
#a) elevated FC (previously a strong change between groups)
#b) propensity to have a co-regulated partner
#c) propensity to have a co-regulated partner, with timing for cis-action
#d) propensity to have a co-regulated partner, with timing for cis-action and HiC/correlation back-up (similar to what worked previously)
#e) GO terms for CCLncs vs. co-regulated with a PCG
#f) optional: HiC loop contact points (requires some thought)


#### finding 2d neighbours ####

fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026filt.csv", header = T)
length(unique(fpkm_allG$EnsID))#12740

Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime_2026.csv", header = T)
Enhancer_lociII_DEsig_Enh <- Enhancer_lociII
length(unique(Enhancer_lociII_DEsig_Enh$EnsID))#597 lncs
length(unique(Enhancer_lociII_DEsig_Enh$MSTRG_Tx_ID))#1575 as expected
length(unique(filter(Enhancer_lociII_DEsig_Enh, !is.na(DiffExprs))$EnsID)) #221 DE

#get co-ords based on FANTOM TSS
Enhancer_lociII_DEsig_Enh$Enhancer_Coords <- paste(Enhancer_lociII_DEsig_Enh$chr, 
                                                   Enhancer_lociII_DEsig_Enh$Enh_Start, 
                                                   Enhancer_lociII_DEsig_Enh$Enh_Stop, sep = ",")

#Get TSS from FANTOM and TSS from GENCODE in same column
Enhancer_lociII_DEsig_Enh$TSS_FANTOM_GENCODE <- Enhancer_lociII_DEsig_Enh$BestStart

#obtain TSS co-ords using these CAGE sites or just 5' limit from GENCODE/Stringtie transcripts for others
trial <- fpkm_allG
trial$Tx_start <- as.numeric(sapply(strsplit(sapply(strsplit(trial$Tx_Locus, ":"), "[[", 2), "-"), "[[", 1))
trial$Tx_stop <- as.numeric(gsub(" [+-]", "", sapply(strsplit(sapply(strsplit(trial$Tx_Locus, ":"), "[[", 2), "-"), "[[", 2)))

trial <- unique(trial[,c(2,5,59:60,8,47:48)])
#alternate/better TSS CAGE from FANTOM for these lncRNAs:
triali <- unique(filter(Enhancer_lociII_DEsig_Enh, CAGEvalidity == "Valid CAGE")[,c(2,44)])

trial <- merge(trial, triali, by = "MSTRG_Tx_ID", all.x = T)

#for +ve strand, add CAGE to start:
trial$Tx_start[trial$str == "+" & !is.na(trial$TSS_FANTOM_GENCODE)] <- trial$TSS_FANTOM_GENCODE[trial$str == "+" & !is.na(trial$TSS_FANTOM_GENCODE)]
#for -ve, add CAGE to stop
trial$Tx_stop[trial$str == "-" & !is.na(trial$TSS_FANTOM_GENCODE)] <- trial$TSS_FANTOM_GENCODE[trial$str == "-" & !is.na(trial$TSS_FANTOM_GENCODE)]

allGB <- unique(trial[,c(1:6)])

#all genes should now have a set of tx with accurate TSS
#double check makes sense:
allGB$limitDiff <- allGB$Tx_stop - allGB$Tx_start

length(unique(fpkm_allG$EnsID))
length(unique(allGB$EnsID))#12740 apiece
length(unique(allGB$MSTRG_Tx_ID))#42511 TSS total (multiple TSS per gene now)


#isolate coding
allGB_PCGs <- filter(allGB, EnsID %in% filter(fpkm_allG, EnsType == "protein_coding", 
                                              grepl("TF|CC|coding", GeneClassUpdate))$EnsID)
length(unique(allGB_PCGs$EnsID))#11613

allGB_PCGs$TSS <- allGB_PCGs$Tx_start
allGB_PCGs$TSS[allGB_PCGs$str == "-"] <- allGB_PCGs$Tx_stop[allGB_PCGs$str == "-"]

trial <- allGB_PCGs
trial <- split(trial, trial$EnsID)

triali <- lapply(trial, function(z){
  filter(allGB_PCGs, 
         #for all PCGs, find PCGs on same chr
         (chr == unique(z$chr) & 
            #with a TSS within a given distance
            TSS > min(z$TSS)-250000 & 
            TSS < max(z$TSS)+250000)
  )
})

triali <- unique(bind_rows(triali, .id = "lnc_id"))
triali <- unique(triali)

AllPCG_AllPCG_1 <- triali

trial <- unique(merge(AllPCG_AllPCG_1[,c(1,3)], fpkm_allG[,c(2,3,58)], by.x = "lnc_id", by.y = "EnsID"))
trial <- unique(merge(trial, fpkm_allG[,c(2,3,58)], by = "EnsID"))

#remove self matches:
trial <- filter(trial, !EnsName.x == EnsName.y)

#2D pairs with TSS within set distance table:
AllPCG_AllPCG_1 <- trial

#correct colnames
colnames(AllPCG_AllPCG_1)[1:2] <- c("EnsID.y", "EnsID")
AllPCG_AllPCG_1$pairs <- paste(AllPCG_AllPCG_1$EnsID, AllPCG_AllPCG_1$EnsID.y, sep="-")

AllPCG_AllPCG_1 <- AllPCG_AllPCG_1[,c(2:4,1,5:7)]

#65938 pairs, 10764 PCGs with another nearby
length(unique(AllPCG_AllPCG_1$EnsName.x))
#write.csv(AllPCG_AllPCG_1, "AllPCG_AllPCG_1_2026_250.csv", row.names =  F)

#obsolete, makes it much easier later on to have these types of "reverse pairs"
#remove any which are same pair in reverse orientation:
#AllPCG_AllPCG_1$pairs2 <- paste(AllPCG_AllPCG_1$EnsID.y, AllPCG_AllPCG_1$EnsID, sep="-")

#to avoid duplicates in pair ID, for each row, order each pair by alphabet:
#trial <- lapply(as.list(as.data.frame(t(AllPCG_AllPCG_1[,c(1,4)]))), function(x){
#  paste(x[order(x)], collapse = "-")
#})

#AllPCG_AllPCG_1$pairs3 <- unlist(trial)
#mostly pairs, one singlet..
#table(table(AllPCG_AllPCG_1$pairs3))

#for each sorted pair label, remove any after first row:
#trial <- split(AllPCG_AllPCG_1, AllPCG_AllPCG_1$pairs3)

#triali <- lapply(trial, function(x){
#  x[1,]
#})

#triali <- bind_rows(triali)
#AllPCG_AllPCG_1 <- triali[,c(1:6,9)]

#table(table(AllPCG_AllPCG_1$EnsID))

#write.csv(AllPCG_AllPCG_1, "AllPCG_AllPCG_1_Q4.csv", row.names =  F)


#### closest neighbour filtering ####

#if pre-filtering to just closest neighbour etc:
allGB$TSS_FANTOM_GENCODE <- allGB$Tx_start
allGB$TSS_FANTOM_GENCODE[allGB$str == "-"] <- allGB$Tx_stop[allGB$str == "-"]

#n.b. multiple TSS per gene, select shortest distance to use in this section
AllPCG_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllPCG_AllPCG_1_2026_250.csv")

triali <- unique(merge(AllPCG_AllPCG_1, allGB[,c(2,8)], by.x = "EnsID", by.y = "EnsID"))
triali <- unique(merge(triali, allGB[,c(2,8)], by.x = "EnsID.y", by.y = "EnsID"))
triali$AbsDistLnc_PCG <- abs(triali$TSS_FANTOM_GENCODE.x - triali$TSS_FANTOM_GENCODE.y)/1000
triali$DistLnc_PCG <- (triali$TSS_FANTOM_GENCODE.x - triali$TSS_FANTOM_GENCODE.y)/1000

#shortest distance per pair:
trialii <- split(triali, triali$pairs)

#trialii[[25]]
trialiii <- lapply(trialii, function(x){
  unique(x[x$AbsDistLnc_PCG == min((x$AbsDistLnc_PCG)),-c(8:9)])
})

trial <- bind_rows(trialiii)

AllPCG_AllPCG_1 <- trial

closestNeighbourPCG <- AllPCG_AllPCG_1
closestNeighbourPCG <- split(closestNeighbourPCG, closestNeighbourPCG$EnsID)

#some alternatives
#closest neighbour
#closestNeighbourPCG <- lapply(closestNeighbourPCG, function(x){
#  filter(x, AbsDistLnc_PCG == min(AbsDistLnc_PCG))})

#closest neighbours
#closestNeighbourPCG <- lapply(closestNeighbourPCG, function(x){
#  x[order(x$DisLnc_PCG, decreasing = F),][1:5,]}
#)

#surrounding neighbours:
closestNeighbourPCG <- lapply(closestNeighbourPCG, function(x){
  upstream <- filter(x, DistLnc_PCG < 0)
  downstream <- filter(x, DistLnc_PCG > 0)
  
  return(rbind(upstream[order(upstream$DistLnc_PCG, decreasing = T),][1,],
               downstream[order(downstream$DistLnc_PCG, decreasing = F),][1,]))
}
)

closestNeighbourPCG <- bind_rows(closestNeighbourPCG)
closestNeighbourPCG <- filter(closestNeighbourPCG, !is.na(pairs))

#should be equal numbers in both? TSS of gene1 in range of gene2 and vice-versa
length(unique(closestNeighbourPCG$EnsID))
length(unique(closestNeighbourPCG$EnsID.y))

#24 not found in y...
missingInCol2 <- unique(closestNeighbourPCG$EnsID)[!unique(closestNeighbourPCG$EnsID) %in% unique(closestNeighbourPCG$EnsID.y)]

#to do with isoforms, keep simple and use EnsID

AllPCG_AllPCG_1 <- closestNeighbourPCG

#### (2D) CC-PCGs in same time frame ####

#AllPCG_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllPCG_AllPCG_1_2026.csv")

# like with the lncRNA figure, some of this could be driven by 0% chance for "later timeframe" CClncs to be 8-24hr 
#so run a version excluding later too

CoRegPairs_04_48_24_samePCG <- filter(AllPCG_AllPCG_1,
                                      #AllLNC_AllPCG_1,
                                      (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                    fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                     fpkm_allGDE_Downwithin_4$EnsID)) |
                                        (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                      fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID)) |
                                        (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                      fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)))
length(unique(c(CoRegPairs_04_48_24_samePCG$EnsID)))
#1100 @250kbp - surrounding neighbours

fpkm_allGDE$CCPCGSame <- fpkm_allGDE$Simple
fpkm_allGDE$CCPCGSame[fpkm_allGDE$EnsID %in% CoRegPairs_04_48_24_samePCG$EnsID] <- "CCPCGSame"

Cluster_biotype <- as.data.frame(table(fpkm_allGDE$CCPCGSame, fpkm_allGDE$RegulationStart))

table(fpkm_allGDE$CCPCGSame)

table(fpkm_allGDE$CCPCGSame)["CCPCGSame"]
table(fpkm_allGDE$CCPCGSame)/5081*100

trial <- filter(Cluster_biotype, Var1 == "CCPCGSame")
trial$Freq/table(fpkm_allGDE$RegulationStart)*100
#looks like  biases are v possible

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(fpkm_allGDE$RegulationStart)

ClusterNames <- trial$Var2

LncEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(fpkm_allGDE$CCPCGSame)["CCPCGSame"]
  d <- 5081
  
  LncEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), 
                                                              row.names = c("Cluster", "other")), alternative = "greater")$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), 
                                                              row.names = c("Cluster", "other")), alternative = "greater")$p)
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
trial$PercCategory <- trial$Freq/table(fpkm_allGDE$CCPCGSame)["CCPCGSame"]*100
trial$PercBackground <- trial$selection/5081*100

CCPCGSameWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])
CCPCGSameWaveBias$Var1 <- "CCPCGSame"

ggplot(CCPCGSameWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(CCPCGSameWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(CCPCGSameWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CCPCGSameWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 4.2) +
  geom_col(data = filter(CCPCGSameWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(CCPCGSameWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CCPCGSameWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory, label = Freq), size = 4.2) +
  ylab("% DE CCPCGs\n(same timeframe)") +
  xlab("") +
  scale_y_continuous(limits = c(-30,45),breaks = seq(-20,45, by = 20),
                     labels = (c(seq(20, 0, by = -20), seq(20,45,by=20)))) +
  theme_minimal() +
  theme(text = element_text(size = 15), legend.position = "none")


#### (2D) bias of CCPCGs (any timeframe) ####

#AllPCG_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllPCG_AllPCG_1_2026.csv")

#recreate the 2D, prior to HiC

CoRegPairs_04_48_24_extendedPCG <- filter(AllPCG_AllPCG_1,
                                          #AllLNC_AllPCG_1,
                                          (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                        fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                            (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                          fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                           fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                            (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                          fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
)
length(unique(CoRegPairs_04_48_24_extendedPCG$EnsID))
#1847 @250kbp - surrounding neighbours

fpkm_allGDE$Simple <- fpkm_allGDE$GeneClassUpdate
fpkm_allGDE$Simple[grepl("fide|Lnc", fpkm_allGDE$GeneClassUpdate)] <- "LncRNA"
fpkm_allGDE$Simple[grepl("coding|TF|CC", fpkm_allGDE$GeneClassUpdate)] <- "PCG"

fpkm_allGDE$CCPCG <- fpkm_allGDE$Simple
fpkm_allGDE$CCPCG[fpkm_allGDE$EnsID %in% c(CoRegPairs_04_48_24_extendedPCG$EnsID)] <- "CCPCG"

table(fpkm_allGDE$CCPCG)

Cluster_biotype <- as.data.frame(table(fpkm_allGDE$CCPCG, fpkm_allGDE$RegulationStart))

table(fpkm_allGDE$CCPCG)


table(fpkm_allGDE$CCPCG)["CCPCG"]
table(fpkm_allGDE$CCPCG)/5081*100

trial <- filter(Cluster_biotype, Var1 == "CCPCG")
trial$Freq/table(fpkm_allGDE$RegulationStart)*100
#looks like  biases are v possible

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(fpkm_allGDE$RegulationStart)

ClusterNames <- trial$Var2

PCGEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(fpkm_allGDE$CCPCG)["CCPCG"]
  d <- 5081
  
  PCGEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), 
                                                              row.names = c("Cluster", "other")), alternative = "greater")$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), 
                                                              row.names = c("Cluster", "other")), alternative = "greater")$p)
}
names(PCGEnrich_cluster) <- ClusterNames
triali <- bind_rows(PCGEnrich_cluster, .id = "Cluster")
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
trial$PercCategory <- trial$Freq/table(fpkm_allGDE$CCPCG)["CCPCG"]*100
trial$PercBackground <- trial$selection/5081*100

CCPCGWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])
CCPCGWaveBias$Var1 <- "CCPCGRNA"

ggplot(CCPCGWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(CCPCGWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(CCPCGWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CCPCGWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 4.2) +
  geom_col(data = filter(CCPCGWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(CCPCGWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CCPCGWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory, label = Freq), size = 4.2) +
  ylab("% DE CCPCGs") +
  xlab("") +
  scale_y_continuous(limits = c(-30,45),breaks = seq(-20,45, by = 20),
                     labels = (c(seq(20, 0, by = -20), seq(20,45,by=20)))) +
  theme_minimal() +
  theme(text = element_text(size = 15), legend.position = "none")

#first wave


#### smaller OR dotplot figure, naive approach ####

AllTypesWaveBias <- rbind(CCPCGSameWaveBias, 
                          CCPCGWaveBias)

AllTypesWaveBias$OR_corrected <- AllTypesWaveBias$OR
AllTypesWaveBias$OR_corrected[AllTypesWaveBias$p_adj >0.1] <- NA
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

appendSpacers <- data.frame(UpDownType = c(rep("Spacer1",3)#,
                                           #rep("Spacer2",3),
                                           #rep("Spacer3",3),
                                           #rep("Spacer4",3)#,
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

#for simpler:
AllTypesWaveBias$UpDownType <- factor(AllTypesWaveBias$UpDownType, 
                                      levels = levels(AllTypesWaveBias$UpDownType)[c(2,1,5,4,3)])

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
    #"p<0.05", 
    #"p<0.01", 
    #"p<0.001", 
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



#### comparing 4hr induction biases, lncs + PCGs with co-reg neighbours ####

#same/later timeframe numbers:
CCLncWaveBias
sum(CCLncWaveBias$Freq)
37/81
CCPCGWaveBias
sum(CCPCGWaveBias$Freq)
565/1847

aii <- 37
bii <- 81    
cii <- 565   
dii <- 1847   

fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

trial <- data.frame(37/81*100, 565/1847*100)

colnames(trial) <- c("CClncRNAs (same/later)",  "CCPCGs (same/later)")

#melt(trial)

trial <- melt(trial)

trial$variable <- factor(trial$variable)
trial$variable <- factor(trial$variable, levels = levels(trial$variable)[c(2,1)])

ggplot(trial) + aes(x = value, y = variable, fill = variable) +
  geom_bar(stat = "identity", color = "grey60") +
  xlab("% 0-4hr induced") +
  ylab("") +
  theme_minimal() +
  scale_fill_manual(values = c(`CClncRNAs (same/later)` = "olivedrab3", `CCPCGs (same/later)` = "mediumorchid")) +
  scale_x_continuous(breaks = seq(0,40,20), limits = c(0,50)) +
  theme(text = element_text(size =20))

#same timeframe numbers:
CCLncNaiveWaveBias
sum(CCLncNaiveWaveBias$Freq)
24/51
CCPCGSameWaveBias
sum(CCPCGSameWaveBias$Freq)
295/1100

aii <- 24
bii <- 51    
cii <- 295   
dii <- 1100   

fisher.test(data.frame("LncRNA" = c(aii,bii-aii),
                       "other" = c(cii-aii,dii-cii-(bii-aii)), row.names = c("Cluster", "other")), alternative = "greater")

trial <- data.frame(24/51*100, 295/1100*100)

colnames(trial) <- c("CClncRNAs (same/later)",  "CCPCGs (same/later)")

#melt(trial)

trial <- melt(trial)

trial$variable <- factor(trial$variable)
trial$variable <- factor(trial$variable, levels = levels(trial$variable)[c(2,1)])

ggplot(trial) + aes(x = value, y = variable, fill = variable) +
  geom_bar(stat = "identity", color = "grey60") +
  xlab("% 0-4hr induced") +
  ylab("") +
  theme_minimal() +
  scale_fill_manual(values = c(`CClncRNAs (same/later)` = "olivedrab3", `CCPCGs (same/later)` = "mediumorchid")) +
  scale_x_continuous(breaks = seq(0,40,20), limits = c(0,50)) +
  theme(text = element_text(size =20))

#
#### enrichment of concordant gene changes near lncs/PCGs (induced genes near induced, rep near rep) ####

#create lnc object in pt2-8 for this part

#up/down DELs:
backgrounds_list <- list(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Upwithin_8, fpkm_allGDE_Upwithin_24)
backgrounds_listi <- list(fpkm_allGDE_Downwithin_4, fpkm_allGDE_Downwithin_8, fpkm_allGDE_Downwithin_24)

#with concordant targets:
hits_list <- list(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Upwithin_8, fpkm_allGDE_Upwithin_24)
hits_listi <- list(fpkm_allGDE_Downwithin_4, fpkm_allGDE_Downwithin_8, fpkm_allGDE_Downwithin_24)

results_list <- list()

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$EnsID)
  cisPotentialTargetsi <- unique(hits_listi[[i]]$EnsID)
  
  #DE PCG genes
  d <- length(unique(
    c(filter(backgrounds_list[[i]], grepl("coding|TF|CC", GeneClassUpdate))$EnsID,
      filter(backgrounds_listi[[i]], grepl("coding|TF|CC", GeneClassUpdate))$EnsID)
  ))
  
  #use pre-filtered to closest neighbour, all hits in PCGs
  #c <- length(unique(
  #  c(filter(backgrounds_list[[i]], grepl("coding|TF|CC", GeneClassUpdate),
  #           EnsID %in% filter(closestNeighbourPCG, EnsID.y %in% cisPotentialTargets)$EnsID)$EnsID,
  #    filter(backgrounds_listi[[i]], grepl("coding|TF|CC", GeneClassUpdate),
  #           EnsID %in% filter(closestNeighbourPCG, EnsID.y %in% cisPotentialTargetsi)$EnsID)$EnsID)
  #))
  
  #or all hits in PCGs
  c <- length(unique(
    c(filter(backgrounds_list[[i]], grepl("coding|TF|CC", GeneClassUpdate),
             EnsID %in% filter(AllPCG_AllPCG_1, EnsID.y %in% cisPotentialTargets)$EnsID)$EnsID,
      filter(backgrounds_listi[[i]], grepl("coding|TF|CC", GeneClassUpdate),
             EnsID %in% filter(AllPCG_AllPCG_1, EnsID.y %in% cisPotentialTargetsi)$EnsID)$EnsID)
  ))
  
  #all selection lncs
  b <- length(unique(
    c(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate))$EnsID,
      filter(backgrounds_listi[[i]], grepl("fide|Lnc", GeneClassUpdate))$EnsID)
  ))
  
  #use pre-filtered to closest neighbours only - for specificity
  #a <- length(unique(
  #  c(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
  #           EnsID %in% filter(closestNeighbour, EnsID.y %in% cisPotentialTargets)$EnsID)$EnsID,
  #    filter(backgrounds_listi[[i]], grepl("fide|Lnc", GeneClassUpdate),
  #           EnsID %in% filter(closestNeighbour, EnsID.y %in% cisPotentialTargetsi)$EnsID)$EnsID)
  #))
  
  #or all hits in selection
  a <- length(unique(
    c(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
             EnsID %in% filter(AllLNC_AllPCG_1, EnsID.y %in% cisPotentialTargets)$EnsID)$EnsID,
      filter(backgrounds_listi[[i]], grepl("fide|Lnc", GeneClassUpdate),
             EnsID %in% filter(AllLNC_AllPCG_1, EnsID.y %in% cisPotentialTargetsi)$EnsID)$EnsID)
  ))
  
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "PCG" = c(c,d-c), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "PCG" = c(c,d-c), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
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
  xlab("") +   ylab("\n% With a concordant\nneighbour") +
  
  theme_minimal()
g2bi

DEL_PCG_type_conc_PCG <- DEL_PCG_type
DEL_PCG_type_conc_PCG$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type_conc_PCG$OR <- unlist(results_list)[c(6,12,18)]
#write.csv(DEL_PCG_type_bi_PCG, "SVSMC_inducedDEL_DEP_DEneighboursame_conc.csv", row.names = F)

p.adjust(DEL_PCG_type_conc_PCG$p)
#@600kbp OR 1.5/0.05
#@250kbp OR 1.3/0.105
#@250kbp closest neighbour 1.53/0.0731
#@250kbp closest 3 neighbours 1.4/0.07
#                 5 neighbours 1.3/0.13
#                 surrounding neighbours 1.67/0.020

library(reshape2)
DEL_PCG_type_plot <- unique(melt(DEL_PCG_type_conc_PCG[,1:3]))

DEL_PCG_type_plot$variable <- as.character(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable[grepl("^DEL", DEL_PCG_type_plot$variable)] <- "DE lncRNAs"
DEL_PCG_type_plot$variable[grepl("DEP", DEL_PCG_type_plot$variable)] <- "DE PCGs"

DEL_PCG_type_plot$variable <- as.factor(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable <- factor(DEL_PCG_type_plot$variable, levels = levels(DEL_PCG_type_plot$variable)[c(1,3,2)])

ggplot(DEL_PCG_type_plot) + aes(x = RegulationBegins, y = value, fill = variable) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.7), color = "grey60") +
  theme_minimal() +
  theme(text = element_text(size = 32)) +
  scale_fill_manual(values = c(`DE lncRNAs` = "olivedrab3", `DE PCGs` = "mediumorchid")) +
  xlab("") +
  ylab("% With a concordant\nneighbour")


#just up DELs?
backgrounds_list <- list(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Upwithin_8, fpkm_allGDE_Upwithin_24)
#with up targets:
hits_list <- list(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Upwithin_8, fpkm_allGDE_Upwithin_24)

results_list <- list()

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$EnsID)
  
  #PCG genes
  d <- length(unique(filter(backgrounds_list[[i]], grepl("coding|TF|CC", GeneClassUpdate))$EnsID))
  
  #prefiltered to closest neighbour etc
  #c <- length(unique(filter(backgrounds_list[[i]], grepl("coding|TF|CC", GeneClassUpdate),
  #                          EnsID %in% filter(closestNeighbourPCG, EnsID.y %in% cisPotentialTargets)$EnsID
  #                          )$EnsID))
  
  #all hits in PCGs
  c <- length(unique(filter(backgrounds_list[[i]], grepl("coding|TF|CC", GeneClassUpdate),
                            EnsID %in% filter(AllPCG_AllPCG_1, EnsID.y %in% cisPotentialTargets)$EnsID
                            )$EnsID))
  
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate))$EnsID))
 
  #prefiltered to closest neighbour etc
  #a <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
  #                          EnsID %in% filter(closestNeighbour, EnsID.y %in% cisPotentialTargets)$EnsID
  #                          )$EnsID))
   
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
                            EnsID %in% filter(AllLNC_AllPCG_1, EnsID.y %in% cisPotentialTargets)$EnsID
                            )$EnsID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "PCG" = c(c,d-c), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "PCG" = c(c,d-c), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, 
                                           results_list[[2]][1]/results_list[[2]][2]*100, 
                                           results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, 
                                           results_list[[2]][3]/results_list[[2]][4]*100, 
                                           results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-24hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g2bi <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "#D6604D") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", linetype = "dashed", size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  xlab("") +   ylab("\n% With an induced \nneighbour") +
  
  theme_minimal()
g2bi

DEL_PCG_type_bi_PCG <- DEL_PCG_type
DEL_PCG_type_bi_PCG$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type_bi_PCG$OR <- unlist(results_list)[c(6,12,18)]
#write.csv(DEL_PCG_type_bi_PCG, "SVSMC_inducedDEL_DEP_DEneighboursame_conc.csv", row.names = F)

p.adjust(DEL_PCG_type_bi_PCG$p)
#@250kbp - all neighbours 1.3/0.188
#@250kbp - closest neighbour 1.72/0.0598
#@250kbp - closest 3 neighbours 1.6/0.051
#          closest 5 neighbours 1.6/0.051
#          surrounding neighbours 1.92/0.0127

library(reshape2)
DEL_PCG_type_plot <- unique(melt(DEL_PCG_type_bi_PCG[,1:3]))

DEL_PCG_type_plot$variable <- as.character(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable[grepl("^DEL", DEL_PCG_type_plot$variable)] <- "Induced lncRNAs"
DEL_PCG_type_plot$variable[grepl("DEP", DEL_PCG_type_plot$variable)] <- "Induced PCGs"

DEL_PCG_type_plot$variable <- as.factor(DEL_PCG_type_plot$variable)
DEL_PCG_type_plot$variable <- factor(DEL_PCG_type_plot$variable, levels = levels(DEL_PCG_type_plot$variable)[c(1,3,2)])

ggplot(DEL_PCG_type_plot) + aes(x = RegulationBegins, y = value, fill = variable) + 
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.7), color = "grey60") +
  theme_minimal() +
  theme(text = element_text(size = 24)) + #, legend.position = "none") +
  scale_fill_manual(values = c(`Induced lncRNAs` = "#D6604D", `Induced PCGs` = "grey80")) +
  xlab("") +
  ylab("% With induced\nneighbour") + Seurat::RotatedAxis()


#just down DELs:
backgrounds_list <- list(fpkm_allGDE_Downwithin_4, fpkm_allGDE_Downwithin_8, fpkm_allGDE_Downwithin_24)
#with down targets:
hits_list <- list(fpkm_allGDE_Downwithin_4, fpkm_allGDE_Downwithin_8, fpkm_allGDE_Downwithin_24)

for (i in (1:length(backgrounds_list))){
  
  #for finding hits
  cisPotentialTargets <- unique(hits_list[[i]]$EnsID)
  
  #PCG genes
  d <- length(unique(filter(backgrounds_list[[i]], grepl("coding|TF|CC", GeneClassUpdate))$EnsID))
  #all hits in PCGs
  c <- length(unique(filter(backgrounds_list[[i]], grepl("coding|TF|CC", GeneClassUpdate),
                            EnsID %in% filter(closestNeighbourPCG, EnsID.y %in% cisPotentialTargets)$EnsID
  )$EnsID))
  #all selection lncs
  b <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate))$EnsID))
  #all hits in selection
  a <- length(unique(filter(backgrounds_list[[i]], grepl("fide|Lnc", GeneClassUpdate),
                            EnsID %in% filter(closestNeighbour, EnsID.y %in% cisPotentialTargets)$EnsID
  )$EnsID))
  
  #background genes contains selection genes:
  results_list[[i]] <- c(a,b,c,d,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "PCG" = c(c,d-c), row.names = c("CC-target", "Other")), alternative = "greater")$p,
                         fisher.test(data.frame("DEL" = c(a,b-a),
                                                "PCG" = c(c,d-c), row.names = c("CC-target", "Other")), alternative = "greater")$estimate)
}

DEL_PCG_type <- data.frame("DEL_pairs" = c(results_list[[1]][1]/results_list[[1]][2]*100, results_list[[2]][1]/results_list[[2]][2]*100, results_list[[3]][1]/results_list[[3]][2]*100), 
                           "DEP_pairs" = c(results_list[[1]][3]/results_list[[1]][4]*100, results_list[[2]][3]/results_list[[2]][4]*100, results_list[[3]][3]/results_list[[3]][4]*100))
DEL_PCG_type$RegulationBegins <- as.factor(c("0-4hrs", "4-8hrs", "8-24hrs"))
DEL_PCG_type$NoDEL <- c(results_list[[1]][1], results_list[[2]][1], results_list[[3]][1])
g2bii <- ggplot(DEL_PCG_type) +
  #scale_y_continuous(limits = c(0,86)) +
  geom_col(aes(x = RegulationBegins, y = DEL_pairs), fill = "#67A9CF") +
  geom_col(aes(x = RegulationBegins, y = DEP_pairs), fill = NA, color = "grey30", #linetype = "dashed", 
           size = 1.2) +
  geom_label(aes(x = RegulationBegins, y = DEL_pairs, label = NoDEL), size = 3) +  
  xlab("") +   ylab("\n% With a repressed\nneighbour") +
  
  theme_minimal()
g2bii

DEL_PCG_type_bii_PCG <- DEL_PCG_type
DEL_PCG_type_bii_PCG$p <- unlist(results_list)[c(5,11,17)]
DEL_PCG_type_bii_PCG$OR <- unlist(results_list)[c(6,12,18)]
#write.csv(DEL_PCG_type_bii_PCG, "SVSMC_repressedDEL_DEP_DEneighboursame_conc.csv", row.names = F)

p.adjust(DEL_PCG_type_bii_PCG$p)#ns

#for both up and down
grid.arrange(g2bi, g2bii, ncol = 2,
             left = textGrob("", gp = gpar(fontface = 'bold', fontsize = 12), rot = 90))

#lack of significance induced, show the concordant graph instead


#### (2D version) compare FC amongst groupings, 0-4hrs co-induced ####

AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_info_2026.csv", header = T)
AllLNC_AllPCG_1 <- filter(AllLNC_AllPCG_1, DisLnc_PCG <250)

#reminder: up and down FCs are different, don't combine effect on induced/repressed in same graph:
summary(fpkm_allGDE_Upwithin_4$LogFC_0_4)
summary(abs(fpkm_allGDE_Downwithin_4$LogFC_0_4))

PCGDE_Upwithin_4 <- filter(fpkm_allGDE_Upwithin_4, EnsType == "protein_coding", grepl("coding|TF|CC", GeneClassUpdate))

#now only co-induced pairs
PCGDE_Upwithin_4$DEneighbourTypeIII <- "None"

#identifies all PCG neighbour pairs with a change in 4hrs
#option to insert a "closest neighbours" object
CoInducedPCG_4hr <- filter(closestNeighbourPCG, EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID),
                           EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID))

CoInducedLncRNA_4hr <- filter(closestNeighbour, EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID),
                              EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID))

PCGDE_Upwithin_4$DEneighbourTypeIII <- "No PCG or CClncRNA neighbour"
PCGDE_Upwithin_4$DEneighbourTypeIII[PCGDE_Upwithin_4$EnsID %in% c(CoInducedPCG_4hr$EnsID)] <- "PCG neighbour only"
PCGDE_Upwithin_4$DEneighbourTypeIII[PCGDE_Upwithin_4$EnsID %in% CoInducedLncRNA_4hr$EnsID.y] <- "CCLncRNA neighbour"

table(PCGDE_Upwithin_4$DEneighbourTypeIII)
#@250kbp: 50 PCG coinduced with a CCLnc in 4hrs, 386 are coinduced with a PCG

PCGDE_Upwithin_4$`Concordant with:` <- PCGDE_Upwithin_4$DEneighbourTypeIII

ggplot(PCGDE_Upwithin_4) + aes(x = `Concordant with:`, y = abs(LogFC_0_4), color = `Concordant with:`) +
  geom_quasirandom(alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c(`CCLncRNA neighbour` = "olivedrab3", `PCG neighbour only` = "mediumorchid", 
                               `No PCG or CClncRNA neighbour` = "grey")) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect(),
        axis.title.y = element_text(size=20),
        legend.position = "bottom",
        axis.text.y = element_text(size=18)) +
  ylab("Log2FC (0-4hrs)")

#anova first
summary(aov(abs(LogFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Upwithin_4))
TukeyHSD(aov(abs(LogFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Upwithin_4))
#Dunnet's for a check too
DescTools::DunnettTest(x = PCGDE_Upwithin_4$LogFC_0_4, PCGDE_Upwithin_4$DEneighbourTypeIII)#same

kruskal.test(abs(LogFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Upwithin_4)#**
#Dunn's post-hoc recommended for KW
FSA::dunnTest(PCGDE_Upwithin_4$LogFC_0_4, PCGDE_Upwithin_4$`Concordant with:`, method = "bonferroni")

summary(filter(PCGDE_Upwithin_4, DEneighbourTypeIII == "CCLncRNA")$LogFC_0_4)
summary(filter(PCGDE_Upwithin_4, DEneighbourTypeIII == "None")$LogFC_0_4)
summary(filter(PCGDE_Upwithin_4, DEneighbourTypeIII == "PCG_only")$LogFC_0_4)
#establishes that no strong diffs in FC between these groups, not worth including


#repressed 
PCGDE_Downwithin_4 <- filter(fpkm_allGDE_Downwithin_4, grepl("coding|TF|CC", GeneClassUpdate))

#now only co-Repressed pairs
PCGDE_Downwithin_4$DEneighbourTypeIII <- "None"

#identifies all PCG neighbour pairs with a change in 4hrs
CoRepressedPCG_4hr <- filter(AllPCG_AllPCG_1, EnsID %in% c(fpkm_allGDE_Downwithin_4$EnsID),
                             EnsID.y %in% c(fpkm_allGDE_Downwithin_4$EnsID))

CoRepressedLncRNA_4hr <- filter(CoRegPairs_04_48_24_extended_naiveSame, EnsID %in% c(fpkm_allGDE_Downwithin_4$EnsID),
                                EnsID.y %in% c(fpkm_allGDE_Downwithin_4$EnsID))

PCGDE_Downwithin_4$DEneighbourTypeIII <- "No PCG or CClncRNA neighbour"
PCGDE_Downwithin_4$DEneighbourTypeIII[PCGDE_Downwithin_4$EnsID %in% c(CoRepressedPCG_4hr$EnsID)] <- "PCG neighbour only"
PCGDE_Downwithin_4$DEneighbourTypeIII[PCGDE_Downwithin_4$EnsID %in% CoRepressedLncRNA_4hr$EnsID.y] <- "CCLncRNA neighbour"

table(PCGDE_Downwithin_4$DEneighbourTypeIII)#17 PCG coRepressed with a CCLnc in 4hrs, 164 are coRepressed with a PCG

PCGDE_Downwithin_4$`Concordant with:` <- PCGDE_Downwithin_4$DEneighbourTypeIII

PCGDE_Downwithin_4_assess <- filter(PCGDE_Downwithin_4, preadj_0_4 <0.05, LogFC_0_4 < -log2(1.5))

ggplot(PCGDE_Downwithin_4_assess) + aes(x = `Concordant with:`, y = abs(LogFC_0_4), color = `Concordant with:`) +
  geom_quasirandom(alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c(`CCLncRNA neighbour` = "olivedrab3", `PCG neighbour only` = "mediumorchid", 
                                `No PCG or CClncRNA neighbour` = "grey")) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect(),
        axis.title.y = element_text(size=20),
        legend.position = "bottom",
        axis.text.y = element_text(size=18)) +
  ylab("Absolute Log2FC (0-4hrs)")

#anova first - pool size bit on the edge for a parametric
#summary(aov(abs(LogFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Downwithin_4_assess))# **
kruskal.test(abs(LogFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Downwithin_4_assess)#**

#parametric looks great, pool is just about ok in size, relatively normal dist too
TukeyHSD(aov(abs(LogFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Downwithin_4_assess))

#non-para for reference
#Dunn's post-hoc recommended for KW, just about
FSA::dunnTest(PCGDE_Downwithin_4_assess$LogFC_0_4, PCGDE_Downwithin_4_assess$`Concordant with:`, method = "bonferroni")

summary(filter(PCGDE_Downwithin_4, DEneighbourTypeIII == "CCLncRNA")$LogFC_0_4)
summary(filter(PCGDE_Downwithin_4, DEneighbourTypeIII == "None")$LogFC_0_4)
summary(filter(PCGDE_Downwithin_4, DEneighbourTypeIII == "PCG_only")$LogFC_0_4)

#same problem as above, driven by 4x HOXA?
#take average value across CXCL and re-check stats (is it driven by just this locus?)
PCGDE_Downwithin_4_LociAv <- PCGDE_Downwithin_4

PCGDE_Downwithin_4_LociAv$EnsName[grepl("HOXA", PCGDE_Downwithin_4_LociAv$EnsName)] <- "HOXA_locus"
PCGDE_Downwithin_4_LociAv$LogFC_0_4[grepl("HOXA", PCGDE_Downwithin_4_LociAv$EnsName)] <- min(PCGDE_Downwithin_4_LociAv$LogFC_0_4[grepl("HOXA", PCGDE_Downwithin_4_LociAv$EnsName)])

PCGDE_Downwithin_4_LociAv <- unique(PCGDE_Downwithin_4_LociAv[,c(50,51,34,2)])

ggplot(PCGDE_Downwithin_4_LociAv) + aes(x = `Concordant with:`, y = abs(LogFC_0_4), color = `Concordant with:`) +
  geom_quasirandom(alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c(`CCLncRNA neighbour` = "olivedrab3", `PCG neighbour only` = "mediumorchid", 
                                `No PCG or CClncRNA neighbour` = "grey")) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect(),
        axis.title.y = element_text(size=20),
        #legend.position = "bottom",
        axis.text.y = element_text(size=18)) +
  ylab("Absolute Log2FC (0-4hrs)")

kruskal.test(abs(LogFC_0_4) ~ DEneighbourTypeIII, data = PCGDE_Downwithin_4_LociAv)#n.s.


#### 0-8hr ####

PCGDE_Upwithin_8 <- filter(fpkm_allGDE_Upwithin_8, grepl("coding|TF|CC", GeneClassUpdate))

#now only co-induced pairs
PCGDE_Upwithin_8$DEneighbourTypeIII <- "None"

#identifies all PCG neighbour pairs with a change in 4hrs
CoInducedPCG_8hr <- filter(AllPCG_AllPCG_1, EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID),
                           EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID))

CoInducedLncRNA_8hr <- filter(AllLNC_AllPCG_1, EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID),
                              EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID))

PCGDE_Upwithin_8$DEneighbourTypeIII <- "No PCG or CClncRNA neighbour"
PCGDE_Upwithin_8$DEneighbourTypeIII[PCGDE_Upwithin_8$EnsID %in% c(CoInducedPCG_8hr$EnsID)] <- "PCG neighbour only"
PCGDE_Upwithin_8$DEneighbourTypeIII[PCGDE_Upwithin_8$EnsID %in% CoInducedLncRNA_8hr$EnsID.y] <- "CCLncRNA neighbour"

table(PCGDE_Upwithin_8$DEneighbourTypeIII)
#@400kbp: 15 PCG coinduced with a lnc, 328 PCG only

PCGDE_Upwithin_8$`Concordant with:` <- PCGDE_Upwithin_8$DEneighbourTypeIII

PCGDE_Upwithin_8_assess <- filter(PCGDE_Upwithin_8, LogFC_0_8 >log2(1.5), preadj_0_8 < 0.05)

ggplot(PCGDE_Upwithin_8_assess) + aes(x = `Concordant with:`, y = abs(LogFC_0_8), color = `Concordant with:`) +
  geom_quasirandom(alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c(`CCLncRNA neighbour` = "olivedrab3", `PCG neighbour only` = "mediumorchid", 
                                `No PCG or CClncRNA neighbour` = "grey")) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect(),
        axis.title.y = element_text(size=20),
        legend.position = "bottom",
        axis.text.y = element_text(size=18)) +
  ylab("Log2FC (0-8hrs)")

kruskal.test(abs(LogFC_0_8) ~ DEneighbourTypeIII, data = PCGDE_Upwithin_8_assess)#*

#non-para for reference
#Dunn's post-hoc recommended for KW, ns
FSA::dunnTest(PCGDE_Upwithin_8_assess$LogFC_0_8, PCGDE_Upwithin_8_assess$`Concordant with:`, method = "bonferroni")

summary(filter(PCGDE_Upwithin_8, DEneighbourTypeIII == "CCLncRNA")$LogFC_0_8)
summary(filter(PCGDE_Upwithin_8, DEneighbourTypeIII == "None")$LogFC_0_8)
summary(filter(PCGDE_Upwithin_8, DEneighbourTypeIII == "PCG_only")$LogFC_0_8)



#repressed 
PCGDE_Downwithin_8 <- filter(fpkm_allGDE_Downwithin_8, grepl("coding|TF|CC", GeneClassUpdate))

#now only co-Repressed pairs
PCGDE_Downwithin_8$DEneighbourTypeIII <- "None"

#identifies all PCG neighbour pairs with a change in 4hrs
CoRepressedPCG_8hr <- filter(AllPCG_AllPCG_1, EnsID %in% c(fpkm_allGDE_Downwithin_8$EnsID),
                             EnsID.y %in% c(fpkm_allGDE_Downwithin_8$EnsID))

CoRepressedLncRNA_8hr <- filter(CoRegPairs_04_48_24_extended_naiveSame, EnsID %in% c(fpkm_allGDE_Downwithin_8$EnsID),
                                EnsID.y %in% c(fpkm_allGDE_Downwithin_8$EnsID))

PCGDE_Downwithin_8$DEneighbourTypeIII <- "No PCG or CClncRNA neighbour"
PCGDE_Downwithin_8$DEneighbourTypeIII[PCGDE_Downwithin_8$EnsID %in% c(CoRepressedPCG_8hr$EnsID)] <- "PCG neighbour only"
PCGDE_Downwithin_8$DEneighbourTypeIII[PCGDE_Downwithin_8$EnsID %in% CoRepressedLncRNA_8hr$EnsID.y] <- "CCLncRNA neighbour"

table(PCGDE_Downwithin_8$DEneighbourTypeIII)#16 PCG coRepressed with a CCLnc in 8hrs

PCGDE_Downwithin_8$`Concordant with:` <- PCGDE_Downwithin_8$DEneighbourTypeIII

PCGDE_Downwithin_8_assess <- filter(PCGDE_Downwithin_8, preadj_0_8 <0.05, LogFC_0_8 < -log2(1.5))

ggplot(PCGDE_Downwithin_8_assess) + aes(x = `Concordant with:`, y = abs(LogFC_0_8), color = `Concordant with:`) +
  geom_quasirandom(alpha = 1) +
  #geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c(`CCLncRNA neighbour` = "olivedrab3", `PCG neighbour only` = "mediumorchid", 
                                `No PCG or CClncRNA neighbour` = "grey")) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect(),
        axis.title.y = element_text(size=20),
        legend.position = "bottom",
        legend.text = element_text(size=15),
        axis.text.y = element_text(size=18)) +
  ylab("Absolute Log2FC (0-8hrs)")

#anova first - pool size bit on the edge for a parametric
kruskal.test(abs(LogFC_0_8) ~ DEneighbourTypeIII, data = PCGDE_Downwithin_8_assess)#ns

summary(filter(PCGDE_Downwithin_8, DEneighbourTypeIII == "CCLncRNA")$LogFC_0_8)
summary(filter(PCGDE_Downwithin_8, DEneighbourTypeIII == "None")$LogFC_0_8)
summary(filter(PCGDE_Downwithin_8, DEneighbourTypeIII == "PCG_only")$LogFC_0_8)


#### 0-24hr ####

PCGDE_Upwithin_24 <- filter(fpkm_allGDE_Upwithin_24, grepl("coding|TF|CC", GeneClassUpdate))

#now only co-induced pairs
PCGDE_Upwithin_24$DEneighbourTypeIII <- "None"

#identifies all PCG neighbour pairs with a change in 4hrs
CoInducedPCG_24hr <- filter(AllPCG_AllPCG_1, EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID),
                            EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID))

CoInducedLncRNA_24hr <- filter(AllLNC_AllPCG_1, EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID),
                               EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID))

PCGDE_Upwithin_24$DEneighbourTypeIII <- "No PCG or CClncRNA neighbour"
PCGDE_Upwithin_24$DEneighbourTypeIII[PCGDE_Upwithin_24$EnsID %in% c(CoInducedPCG_24hr$EnsID)] <- "PCG neighbour only"
PCGDE_Upwithin_24$DEneighbourTypeIII[PCGDE_Upwithin_24$EnsID %in% CoInducedLncRNA_24hr$EnsID.y] <- "CCLncRNA neighbour"

table(PCGDE_Upwithin_24$DEneighbourTypeIII)
#@400kbp: 5 PCG coinduced with a lnc, 178 PCG only

PCGDE_Upwithin_24$`Concordant with:` <- PCGDE_Upwithin_24$DEneighbourTypeIII

PCGDE_Upwithin_24_assess <- filter(PCGDE_Upwithin_24, LogFC_0_24 >log2(1.5), preadj_0_24 < 0.05)

ggplot(PCGDE_Upwithin_24_assess) + aes(x = `Concordant with:`, y = abs(LogFC_0_24), color = `Concordant with:`) +
  geom_quasirandom(alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c(`CCLncRNA neighbour` = "olivedrab3", `PCG neighbour only` = "mediumorchid", 
                                `No PCG or CClncRNA neighbour` = "grey")) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect(),
        axis.title.y = element_text(size=20),
        legend.position = "bottom",
        axis.text.y = element_text(size=18)) +
  ylab("Log2FC (0-24hrs)")

kruskal.test(abs(LogFC_0_24) ~ DEneighbourTypeIII, data = PCGDE_Upwithin_24_assess)#ns

summary(filter(PCGDE_Upwithin_24, DEneighbourTypeIII == "CCLncRNA")$LogFC_0_24)
summary(filter(PCGDE_Upwithin_24, DEneighbourTypeIII == "None")$LogFC_0_24)
summary(filter(PCGDE_Upwithin_24, DEneighbourTypeIII == "PCG_only")$LogFC_0_24)



#repressed 
PCGDE_Downwithin_24 <- filter(fpkm_allGDE_Downwithin_24, grepl("coding|TF|CC", GeneClassUpdate))

#now only co-Repressed pairs
PCGDE_Downwithin_24$DEneighbourTypeIII <- "None"

#identifies all PCG neighbour pairs with a change in 4hrs
CoRepressedPCG_24hr <- filter(AllPCG_AllPCG_1, EnsID %in% c(fpkm_allGDE_Downwithin_24$EnsID),
                              EnsID.y %in% c(fpkm_allGDE_Downwithin_24$EnsID))

CoRepressedLncRNA_24hr <- filter(CoRegPairs_04_48_24_extended_naiveSame, EnsID %in% c(fpkm_allGDE_Downwithin_24$EnsID),
                                 EnsID.y %in% c(fpkm_allGDE_Downwithin_24$EnsID))

PCGDE_Downwithin_24$DEneighbourTypeIII <- "No PCG or CClncRNA neighbour"
PCGDE_Downwithin_24$DEneighbourTypeIII[PCGDE_Downwithin_24$EnsID %in% c(CoRepressedPCG_24hr$EnsID)] <- "PCG neighbour only"
PCGDE_Downwithin_24$DEneighbourTypeIII[PCGDE_Downwithin_24$EnsID %in% CoRepressedLncRNA_24hr$EnsID.y] <- "CCLncRNA neighbour"

table(PCGDE_Downwithin_24$DEneighbourTypeIII)#14 PCG coRepressed with a CCLnc in 24hrs

PCGDE_Downwithin_24$`Concordant with:` <- PCGDE_Downwithin_24$DEneighbourTypeIII

PCGDE_Downwithin_24_assess <- filter(PCGDE_Downwithin_24, preadj_0_24 <0.05, LogFC_0_24 < -log2(1.5))

ggplot(PCGDE_Downwithin_24_assess) + aes(x = `Concordant with:`, y = abs(LogFC_0_24), color = `Concordant with:`) +
  geom_quasirandom(alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c(`CCLncRNA neighbour` = "olivedrab3", `PCG neighbour only` = "mediumorchid", 
                                `No PCG or CClncRNA neighbour` = "grey")) +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  xlab("")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect(),
        axis.title.y = element_text(size=20),
        legend.position = "bottom",
        axis.text.y = element_text(size=18)) +
  ylab("Absolute Log2FC (0-24hrs)")

kruskal.test(abs(LogFC_0_24) ~ DEneighbourTypeIII, data = PCGDE_Downwithin_24_assess)#**

summary(filter(PCGDE_Downwithin_24, DEneighbourTypeIII == "CCLncRNA")$LogFC_0_24)
summary(filter(PCGDE_Downwithin_24, DEneighbourTypeIII == "None")$LogFC_0_24)
summary(filter(PCGDE_Downwithin_24, DEneighbourTypeIII == "PCG_only")$LogFC_0_24)

