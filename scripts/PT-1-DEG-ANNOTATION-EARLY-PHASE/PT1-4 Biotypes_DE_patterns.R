library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(ggbeeswarm)

#### import expressed G and DE G ####

fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026filt.csv")

fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE_2026filt.csv")


#### expected/positive controls ####

#expectations, early immune response and IEGs, cell cycle induction later, loss of SMC markers

#heatmap

#SMC markers
SMC_phenoSwitch <- c("ACTA2", "CNN1", "VIM", "DES", "TAGLN", "CARMN", 
                     "YY1", "KLF4", "MYOCD", "TET2", "SMAD3", "TCF21")

#IEGs (inc. lots of immune)
IEGs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Gene lists/arner_2015_table_S5_IEGs_lit.csv")
IEGs_hs <- filter(IEGs, !Hs_symbol %in% c("", "-"))
IEGs_hs <- filter(IEGs_hs, !grepl("[a-z]", IEGs_hs$Total_count))


plot_expectedGenes <- filter(fpkm_allGDE, (grepl("CC", GeneClassUpdate)) | EnsName %in% c(IEGs_hs[,4]))

plot_expectedGenes <- filter(fpkm_allGDE, (grepl("CC", GeneClassUpdate)))

plot_expectedGenes <- filter(fpkm_allGDE, EnsName %in% c(IEGs_hs[,4]))

plot_expectedGenes <- filter(fpkm_allGDE, EnsName %in% SMC_phenoSwitch)

mat <- plot_expectedGenes[,4:19]
rownames(mat) <- plot_expectedGenes[,2]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
mat <- t(apply(mat, 1, cal_z_score))

myColor <- colorRampPalette(c("steelblue", "white", "red"))(50)
myBreaks <- c(seq(min(mat), 0, 
                  length.out=ceiling(50/2)), 
              seq(max(mat)/50, 
                  max(mat), 
                  length.out=floor(50/2)))


condition<-c(rep("0hr",4),rep("4hr",4),rep("8hr",4),rep("24hr",4))
patient<-rep(c("1","2","3","4"),4)
data_cols<-data.frame(condition=condition,patient=patient)

data_colsHeat <- data.frame("Hours" = data_cols[,1], stringsAsFactors = T)
rownames(data_colsHeat) <- colnames(mat)
data_colsHeat$Hours <- factor(data_colsHeat$Hours, levels(data_colsHeat$Hours)[c(1,3,4,2)])


library(pheatmap)
p <- pheatmap(mat,clustering_method = "complete",annotation_legend = F,
              annotation_col = data_colsHeat,
              show_colnames = F, 
              show_rownames = F, 
              cluster_cols = F,
              cluster_rows = T,
              #cutree_rows = 2,
              treeheight_col = 0, 
              treeheight_row = 45,
              #legend = F,
              color = myColor, 
              breaks = myBreaks,
              border_color = NA)

rownames(mat)[p$tree_row$order]#first 23 are down

#SMC pheno
p <- pheatmap(mat,clustering_method = "complete",annotation_legend = F,
              annotation_col = data_colsHeat,
              show_colnames = F, 
              show_rownames = T, 
              cluster_cols = F,
              cluster_rows = T,
              #cutree_rows = 2,
              treeheight_col = 0, 
              treeheight_row = 45,
              #legend = F,
              color = myColor, 
              breaks = myBreaks,
              border_color = NA)


#### temporal clustering ####

#seperate genes by the timing of their first regulation

length(unique(fpkm_allGDE$EnsID))

#regulated within 4 hours:
fpkm_allGDE_Upwithin_4 <- filter(fpkm_allGDE, (Hour0_meanFPKM >=0.8 | Hour4_meanFPKM >=0.8) &
                                   (LogFC_0_4 >= log2(1.5) & preadj_0_4 <0.05))
fpkm_allGDE_Downwithin_4 <- filter(fpkm_allGDE, (Hour0_meanFPKM >=0.8 | Hour4_meanFPKM >=0.8) &
                                     (LogFC_0_4 < -log2(1.5) & preadj_0_4 <0.05))
#of remaining, within 8 hours:
fpkm_allGDE_Upwithin_8 <- filter(fpkm_allGDE, 
                                 !EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                                 (Hour0_meanFPKM >=0.8 | Hour4_meanFPKM >=0.8 | Hour8_meanFPKM >=0.8),
                                 ((LogFC_0_8 >= log2(1.5) & preadj_0_8 <0.05) |
                                   (LogFC_4_8 >= log2(1.5) & preadj_4_8 <0.05))
                                 )
fpkm_allGDE_Downwithin_8 <- filter(fpkm_allGDE, 
                                 !EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                                 (Hour0_meanFPKM >=0.8 | Hour4_meanFPKM >=0.8 | Hour8_meanFPKM >=0.8),
                                 ((LogFC_0_8 <= -log2(1.5) & preadj_0_8 <0.05) |
                                    (LogFC_4_8 <= -log2(1.5) & preadj_4_8 <0.05))
                                 )

#of remaining, within 24 hours:
fpkm_allGDE_Upwithin_24 <- filter(fpkm_allGDE, 
                                 !EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID,
                                               fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID),
                                 (fpkm_max_treatment >=0.8),
                                 ((LogFC_0_24 >= log2(1.5) & preadj_0_24 <0.05) |
                                    (LogFC_4_24 >= log2(1.5) & preadj_4_24 <0.05) |
                                    (LogFC_8_24 >= log2(1.5) & preadj_8_24 <0.05)
                                  )
                                 )
fpkm_allGDE_Downwithin_24 <- filter(fpkm_allGDE, 
                                  !EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID,
                                                fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID),
                                  (fpkm_max_treatment >=0.8),
                                  ((LogFC_0_24 <= -log2(1.5) & preadj_0_24 <0.05) |
                                     (LogFC_4_24 <= -log2(1.5) & preadj_4_24 <0.05) |
                                     (LogFC_8_24 <= -log2(1.5) & preadj_8_24 <0.05)
                                  )
                                  )

#for figures, assign each gene to wave, should have all the genes seperated by the above:
973+746+749+652+666+559#4345 as expected
1105+828+821+891+768+668#

fpkm_allGDE$RegulationStart[fpkm_allGDE$EnsID %in% fpkm_allGDE_Upwithin_4$EnsID] <- "Induced <4hrs"
fpkm_allGDE$RegulationStart[fpkm_allGDE$EnsID %in% fpkm_allGDE_Downwithin_4$EnsID] <- "Repressed <4hrs"

fpkm_allGDE$RegulationStart[fpkm_allGDE$EnsID %in% fpkm_allGDE_Upwithin_8$EnsID] <- "Induced 4-8hrs"
fpkm_allGDE$RegulationStart[fpkm_allGDE$EnsID %in% fpkm_allGDE_Downwithin_8$EnsID] <- "Repressed 4-8hrs"

fpkm_allGDE$RegulationStart[fpkm_allGDE$EnsID %in% fpkm_allGDE_Upwithin_24$EnsID] <- "Induced 8-24hrs"
fpkm_allGDE$RegulationStart[fpkm_allGDE$EnsID %in% fpkm_allGDE_Downwithin_24$EnsID] <- "Repressed 8-24hrs"

#remove others (only 1):
fpkm_allGDE <- filter(fpkm_allGDE, !is.na(RegulationStart))

#write.csv(fpkm_allGDE, "fpkm_allGDE_2026filt.csv", row.names = F)


#### bar plots showing size of waves ####

#size of waves
trial <- as.data.frame(table(fpkm_allGDE$RegulationStart))
trial$Freq/length(fpkm_allGDE$EnsID)
trial$FirstRegulation <- c("0-4hrs", "4-8hrs", "8-24hrs", "0-4hrs", "4-8hrs", "8-24hrs")
trial$UpDown <- sapply(sapply(as.character(trial$Var1), strsplit, " "),"[[" , 1)
trial$PercAllLncs <- trial$Freq/length(fpkm_allGDE$EnsID)

trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation))#[c(2,3,1)])

AllGWaveBias <- trial

ggplot(trial, aes(x = FirstRegulation)) +
  geom_bar(stat = "identity", data = filter(trial, grepl("Induced", Var1)), 
           aes(y = PercAllLncs, fill = UpDown), color = "black") +
  geom_label(data = filter(trial, grepl("Induced", Var1)), 
             aes(y = PercAllLncs-0.05, label = Freq), size = 4.2) +
  geom_bar(stat = "identity", data = filter(trial, grepl("Repressed", Var1)), 
           aes(y = -PercAllLncs-0.05, fill = UpDown), color = "black") +
  geom_label(data = filter(trial, grepl("Repressed", Var1)), 
             aes(y = -PercAllLncs, label = Freq), size = 4.2) +
  ylab("% all DEGs") +
  xlab("") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(-0.25,0.25), breaks = seq(-0.20,0.20, by = 0.20)) +
  scale_fill_manual(values = c("Induced" = "#D6604D", "Repressed" = "#67A9CF", "Notable" = "purple")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 18, angle = 45, hjust = 1.2, vjust = 1.4),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 18))


#### clusters enriched GO terms ####

#build support for clustering - biological theming
fpkm_allGDE$Ens_ID_merge <- sapply(sapply(fpkm_allGDE$EnsID, strsplit, "\\."), "[[", 1)
fpkm_allG$Ens_ID_merge <- sapply(sapply(fpkm_allG$EnsID, strsplit, "\\."), "[[", 1)

trial <- split(fpkm_allGDE, fpkm_allGDE$RegulationStart)
sapply(trial, dim)

library(clusterProfiler)
library(org.Hs.eg.db)
trial <- lapply(trial, function(x){
  enrichGO(gene          = unique(x$Ens_ID_merge),
           universe      = unique(fpkm_allG$Ens_ID_merge),
           keyType       = "ENSEMBL",
           OrgDb         = org.Hs.eg.db,
           ont           = "ALL",
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05,
           readable      = TRUE)
})

saveRDS(trial, "GO_list_SVSMC_clusters_Dec25.rds")

#write.csv(as.data.frame(trial[[1]]), "Up_4.csv", row.names = F)
#write.csv(as.data.frame(trial[[2]]), "Up_8.csv", row.names = F)
#write.csv(as.data.frame(trial[[3]]), "Up_24.csv", row.names = F)
#write.csv(as.data.frame(trial[[4]]), "Down_4.csv", row.names = F)
#write.csv(as.data.frame(trial[[5]]), "Down_8.csv", row.names = F)
#write.csv(as.data.frame(trial[[6]]), "Down_24.csv", row.names = F)


GO_list_SVSMC_clusters25.rds <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GO_list_SVSMC_clusters25.rds")

sapply(trial, dim)#only lost within 8hrs is weakly themed

edox2 <- enrichplot::pairwise_termsim(trial[[1]])
enrichplot::treeplot(edox2, showCategory = 30, #fontSize =2, #extend = 0.1, #hilight =F #nWords = 3,
                     cluster.params = list(n = 5), 
                     #label_format_tiplab = function(x) stringr::str_wrap(x, width=40),
                     label_format = function(x) stringr::str_wrap(x, width=25))
View(data.frame(trial[[1]]))

triali <- lapply(GO_list_SVSMC_clusters25.rds, simplify)

dotplot(triali[[1]], showCategory = 20, label_format = function(x) stringr::str_wrap(x, width=37))
dotplot(triali[[2]], showCategory = 20, label_format = function(x) stringr::str_wrap(x, width=37))
dotplot(triali[[3]], showCategory = 20, label_format = function(x) stringr::str_wrap(x, width=37))

dotplot(triali[[4]], showCategory = 20, label_format = function(x) stringr::str_wrap(x, width=37))
dotplot(triali[[5]], showCategory = 20, label_format = function(x) stringr::str_wrap(x, width=37))
dotplot(triali[[6]], showCategory = 20, label_format = function(x) stringr::str_wrap(x, width=37))


#### bias of gene sets to times of first initiation/repression ####

fpkm_allGDE$Simple <- fpkm_allGDE$GeneClassUpdate
fpkm_allGDE$Simple[grepl("fide|Lnc", fpkm_allGDE$GeneClassUpdate)] <- "LncRNA"

Cluster_biotype <- as.data.frame(table(fpkm_allGDE$Simple, fpkm_allGDE$RegulationStart))

table(fpkm_allGDE$Simple)

#bias of lncRNAs:
table(fpkm_allGDE$Simple)/length(fpkm_allGDE$EnsID)*100
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
  d <- length(fpkm_allGDE$EnsID)
  
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
#2x significant biases

#percentage plots
trial$FirstRegulation <- c("Within \n4hrs", "Within \n8hrs", "Within \n24hrs", "Within \n4hrs", "Within \n8hrs", "Within \n24hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(2,3,1)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/221*100
trial$PercBackground <- trial$selection/5081*100

LncWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(LncWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(LncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(LncWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(LncWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 5) +
  geom_col(data = filter(LncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(LncWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(LncWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory, label = Freq), size = 5) +
  ylab("% DE LncRNAs") +
  xlab("") +
  scale_y_continuous(limits = c(-30,35),breaks = seq(-30,30, by = 10),
                     labels = (c(seq(30, 0, by = -10), seq(10,30,by=10)))) +
  theme_minimal() +
  theme(text = element_text(size=15))


#Cell cycle - Giotti S/G2M
fpkm_allGDE$CC <- fpkm_allGDE$GeneClassUpdate
fpkm_allGDE$CC[grepl("CC", fpkm_allGDE$GeneClassUpdate)] <- "CC"

Selected <- "CC"
table(fpkm_allGDE$CC)[Selected]
table(fpkm_allGDE$CC)/length(fpkm_allGDE$EnsID)*100

Cluster_biotype <- as.data.frame(table(fpkm_allGDE$CC, fpkm_allGDE$RegulationStart))
trial <- filter(Cluster_biotype, grepl(Selected, Var1))
trial$Freq/table(fpkm_allGDE$RegulationStart)*100
#looks like definite biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(fpkm_allGDE$RegulationStart)

ClusterNames <- trial$Var2

CCEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(fpkm_allGDE$CC)[Selected]
  d <- length(fpkm_allGDE$EnsID)
  
  CCEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                             "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est,
                                      fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                             "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p)
}
names(CCEnrich_cluster) <- ClusterNames
triali <- bind_rows(CCEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")
#all significant biases

#percentage plots
trial$FirstRegulation <- c("Within \n4hrs", "Within \n8hrs", "Within \n24hrs", "Within \n4hrs", "Within \n8hrs", "Within \n24hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(2,3,1)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(fpkm_allGDE$CC)[Selected]*100
trial$PercBackground <- trial$selection/length(fpkm_allGDE$EnsID)*100

CCWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(CCWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(CCWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(CCWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CCWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 4.2) +
  geom_col(data = filter(CCWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(CCWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(CCWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory, label = Freq), size = 4.2) +
  ylab("% DE S/G2M") +
  xlab("") +
  #scale_y_continuous(limits = c(-30,35),breaks = seq(-30,30, by = 10),
  #                   labels = (c(seq(30, 0, by = -10), seq(10,30,by=10)))) +
  theme_minimal() +
  theme(text = element_text(size = 15))


#IEGs
#inserting new gene categories:
IEGs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Gene lists/arner_2015_table_S5_IEGs_lit.csv")
IEGs_hs <- filter(IEGs, !Hs_symbol %in% c("", "-"))
IEGs_hs <- filter(IEGs_hs, !grepl("[a-z]", IEGs_hs$Total_count))

fpkm_allGDE$IEG <- fpkm_allGDE$GeneClassUpdate
fpkm_allGDE$IEG[fpkm_allGDE$EnsName %in% IEGs_hs$Hs_symbol] <- "IEG"

Selected <- "IEG"
table(fpkm_allGDE$IEG)[Selected]
table(fpkm_allGDE$IEG)/length(fpkm_allGDE$EnsID)*100
#141 IEGs are 2.8% of all DEGs

Cluster_biotype <- as.data.frame(table(fpkm_allGDE$IEG, fpkm_allGDE$RegulationStart))
trial <- filter(Cluster_biotype, grepl(Selected, Var1))
trial$Freq/table(fpkm_allGDE$RegulationStart)*100
#looks like definite biases

#fisher test: background DE PCGs, selection cluster, hit lncRNA
trial$selection <- table(fpkm_allGDE$RegulationStart)

ClusterNames <- trial$Var2

IEGEnrich_cluster <- list()

for (i in 1:length(ClusterNames)){
  a <- trial[i,3]
  b <- trial[i,4]
  c <- table(fpkm_allGDE$IEG)[Selected]
  d <- length(fpkm_allGDE$EnsID)
  
  IEGEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p)
}
names(IEGEnrich_cluster) <- ClusterNames
triali <- bind_rows(IEGEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")
#all significant biases

#percentage plots
trial$FirstRegulation <- c("Within \n4hrs", "Within \n8hrs", "Within \n24hrs", "Within \n4hrs", "Within \n8hrs", "Within \n24hrs")
trial$FirstRegulation <- as.factor(trial$FirstRegulation)
trial$FirstRegulation <- factor(trial$FirstRegulation, levels = levels(trial$FirstRegulation)[c(2,3,1)])
trial$UpDown <- sapply(sapply(as.character(trial$Var2), strsplit, " "),"[[" , 1)
trial$PercCategory <- trial$Freq/table(fpkm_allGDE$IEG)[Selected]*100
trial$PercBackground <- trial$selection/length(fpkm_allGDE$EnsID)*100

IEGWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(IEGWaveBias, aes(x = FirstRegulation)) +
  geom_col(data = filter(IEGWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown)) +
  geom_col(data = filter(IEGWaveBias, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(IEGWaveBias, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 4.2) +
  geom_col(data = filter(IEGWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown)) +
  geom_col(data = filter(IEGWaveBias, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed") +
  geom_label(data = filter(IEGWaveBias, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory, label = Freq), size = 4.2) +
  ylab("% DE IEGs") +
  xlab("") +
  #scale_y_continuous(limits = c(-30,35),breaks = seq(-30,30, by = 10),
  #                   labels = (c(seq(30, 0, by = -10), seq(10,30,by=10)))) +
  theme_minimal() +
  theme(text = element_text(size = 15))


#combined figure:
AllTypesWaveBias <- rbind(PCGWaveBias, 
  CCWaveBias, IEGWaveBias, LncWaveBias
  )

#correct timeframe labels:
AllTypesWaveBias$FirstRegulation <- as.character(AllTypesWaveBias$FirstRegulation)
AllTypesWaveBias$FirstRegulation[AllTypesWaveBias$FirstRegulation == "Within \n4hrs"] <- "Within \n0-4hrs"
AllTypesWaveBias$FirstRegulation[AllTypesWaveBias$FirstRegulation == "Within \n8hrs"] <- "Within \n4-8hrs"
AllTypesWaveBias$FirstRegulation[AllTypesWaveBias$FirstRegulation == "Within \n24hrs"] <- "Within \n8-24hrs"

AllTypesWaveBias$OR_corrected <- AllTypesWaveBias$OR
AllTypesWaveBias$OR_corrected[AllTypesWaveBias$p_adj >0.1] <- NA
AllTypesWaveBias$OR_corrected <- AllTypesWaveBias$OR_corrected + 0.001
AllTypesWaveBias$`Log2(Odds Ratio)` <- log2(AllTypesWaveBias$OR_corrected)

AllTypesWaveBias$`Log2(Odds Ratio)`[AllTypesWaveBias$`Log2(Odds Ratio)` < -5] <- -5
AllTypesWaveBias$`Log2(Odds Ratio)`[AllTypesWaveBias$`Log2(Odds Ratio)` >5] <- 5

AllTypesWaveBias$padj_corrected <- AllTypesWaveBias$p_adj
AllTypesWaveBias$padj_corrected[AllTypesWaveBias$p_adj >0.1] <- NA

AllTypesWaveBias$padj_simple <- NA
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.1] <- "p<0.1"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.05] <- "p<0.05"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.01] <- "p<0.01"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.001] <- "p<0.001"
AllTypesWaveBias$padj_simple[AllTypesWaveBias$padj_corrected < 0.0001] <- "p<0.0001"

AllTypesWaveBias$UpDownType <- paste(AllTypesWaveBias$Var1, AllTypesWaveBias$UpDown, sep = "-")

AllTypesWaveBias <- AllTypesWaveBias[,c(4,15,14,12)]

#insert spacer between each biotype:
spacers <- data.frame(
           "FirstRegulation" = c(rep("Within \n0-4hrs", 4), rep("Within \n4-8hrs", 4), rep("Within \n8-24hrs", 4)),
           "UpDownType" = rep(c("Space1", "Space2", "Space3", "Space4"), 3),
           "padj_simple" = NA,
           "Log2(Odds Ratio)" = NA
           )
colnames(spacers) <- colnames(AllTypesWaveBias)

AllTypesWaveBias <- rbind(AllTypesWaveBias, spacers)

AllTypesWaveBias$UpDownType <- factor(AllTypesWaveBias$UpDownType)
AllTypesWaveBias$UpDownType <- factor(AllTypesWaveBias$UpDownType, 
                                      levels = levels(AllTypesWaveBias$UpDownType)[c(6,5,9,2,1,10,4,3,11,8,7)])

myColor <- colorRampPalette(c("steelblue", "white", "red"))(50)
myBreaks <- c(seq(min(AllTypesWaveBias$`Log2(Odds Ratio)`, na.rm = T), 0, 
                  length.out=ceiling(50/2)), 
              seq(max(AllTypesWaveBias$`Log2(Odds Ratio)`, na.rm = T)/50, 
                  max(AllTypesWaveBias$`Log2(Odds Ratio)`, na.rm = T), 
                  length.out=floor(50/2)))

ggplot(AllTypesWaveBias[-c(1:6),]) + aes(x = FirstRegulation, y = UpDownType, size = padj_simple, fill = `Log2(Odds Ratio)`) +
  geom_point(color = "grey30", shape = 21) +
  xlab("") +
  #ylab("") +
  scale_size_discrete(range = c(6,12), limits = c("p<0.05", #"p<0.01","p<0.001", 
                                                  "p<0.0001"#, "p<0.000001", "p<0.00000001"
                                                  )) +
  scale_fill_gradient2(low = "steelblue", mid = "lightsalmon", high = "red") +
  theme_minimal() +
  theme(legend.key.size = unit(1.4, "line"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        axis.text.x = element_text(size=22),
        #axis.text.y = element_blank()
  )

    
#### % of DEGs per biotype ####

#miRNA expression from small RNAseq, 10RPMM cut off for expression:
rpm_allmiRs_annotated <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/rpm_allmiRs_annotated_2025.csv", header = T)
rpm_allmiRs_annotated <-  filter(rpm_allmiRs_annotated, rpm_max_treatment >10)
rpm_allmiRs_annotated_DE <-  filter(rpm_allmiRs_annotated, rpm_max_treatment >10, DE_consensus == "DE")
length(unique(rpm_allmiRs_annotated_DE$nameStars))#22 unique DE miRNAs
length(unique(rpm_allmiRs_annotated$nameStars))#358unique miRNAs expressed
22/338*100 #6.5% of miRNAs change over time


#lnc/PCG equivalent
table(unique(fpkm_allG[,c(2,58)])$GeneClassUpdate)
table(fpkm_allGDE$Simple)

221/597#37% lncRNAs

length(unique(filter(fpkm_allG, EnsType == "protein_coding", grepl("coding|TF|CC", GeneClassUpdate))$EnsID))
#11613
length(unique(filter(fpkm_allGDE, EnsType == "protein_coding", grepl("coding|TF|CC", GeneClassUpdate))$EnsID))
#4701

4701/11613

#graph of this
plot_DEType <- data.frame("DE" = c(22,221,4701), "Stable" = c(338-22,597-221,11613-4701), "SpeciesType" = c("miRNA", "LncRNA", "PCG"))
plot_DEType <- reshape2::melt(plot_DEType)

plot_DEType$SpeciesType <- as.factor(plot_DEType$SpeciesType)
plot_DEType$SpeciesType <- factor(plot_DEType$SpeciesType, levels(plot_DEType$SpeciesType)[c(3,1,2)])

ggplot(plot_DEType) + aes(x = SpeciesType, y = value, fill = variable) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(labels = c("DE" = "DE", "Stable" = "No Change"), values = c("DE" = "grey30", "Stable" = "grey90")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 18, angle = 45, vjust =1.05, hjust = 0.95),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(10, 10, 10, 100),
        legend.title = element_blank())

#miRNAs are far more IP-invariant than other gene classes on the whole
#lncRNAs on a par with PCG, TF
#curated list of IEGs and core S/G2Ms v. likely to be DE in the timecourse

#simplified:
ggplot(filter(plot_typeDE, SpeciesType %in% c("PCG", "LncRNA", "miRNA"))) + aes(x = SpeciesType, fill = DE_consensus) +
  geom_bar(position = "fill", color = "black") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(labels = c("DE" = "DE", "Stable" = "No Change"),
                     values = c("DE" = "grey30", "Stable" = "grey90", "Notable" = "purple")) +
  theme_minimal() +
  theme(text = element_text(size = 25),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        #plot.margin = margin(10, 10, 10, 100),
        legend.title = element_blank()) + Seurat::RotatedAxis()


#### furthers early bias depiction ####

#line plot of FCs across all genes (including ns...)

fpkm_lncGDE <- filter(fpkm_allGDE, grepl("Lnc|fide", GeneClassUpdate))

#baseline 0hr
colnames(fpkm_lncGDE)
fpkm_lncGDE <- fpkm_lncGDE[,c(1,30,36,38)]

colnames(fpkm_lncGDE) <- c("EnsID", "4", "8", "24")

fpkm_lncGDE$`0` <- 0

library(reshape2)

trial <- melt(fpkm_lncGDE)

#trial$variable <- factor(trial$variable, levels(trial$variable)[c(4,1,2,3)])

trial$variable <- as.numeric(as.character(trial$variable))

trial$geneSet <- NA
trial$geneSet[trial$EnsID %in% fpkm_allGDE_Downwithin_24$EnsID] <- "Repressed 0-24hr"
trial$geneSet[trial$EnsID %in% fpkm_allGDE_Upwithin_24$EnsID] <- "Induced 0-24hr"
trial$geneSet[trial$EnsID %in% fpkm_allGDE_Downwithin_8$EnsID] <- "Repressed 4-8hr"
trial$geneSet[trial$EnsID %in% fpkm_allGDE_Upwithin_8$EnsID] <- "Induced 4-8hr"
trial$geneSet[trial$EnsID %in% fpkm_allGDE_Downwithin_4$EnsID] <- "Repressed 0-4hr"
trial$geneSet[trial$EnsID %in% fpkm_allGDE_Upwithin_4$EnsID] <- "Induced 0-4hr"

table(unique(trial[,c(1,4)])$geneSet)

RColorBrewer::brewer.pal(n = 6, name = "RdBu")

trial$geneSet <- as.factor(trial$geneSet)
trial$geneSet <- factor(trial$geneSet, levels(trial$geneSet)[c(2,3,1,5,6,4)])

ggplot(trial) + aes(x = variable, y = value, color = geneSetII) +
  geom_line(aes(x = variable, y = value, group = EnsID, color = geneSet)) +
  #geom_line(data = filter(trial, !is.na(geneSetIII)), size = 1.2,
   #         aes(x = variable, y = value, group = EnsID, color = geneSetIII, alpha = 0.5)) +
  #geom_point() + 
  scale_color_manual(values = c("Induced 0-24hr" = "#FDDBC7", "Repressed 0-24hr" = "#D1E5F0", 
                                "Induced 4-8hr" = "#EF8A62", "Repressed 4-8hr" = "#67A9CF", 
                                "Induced 0-4hr" = "#B2182B", "Repressed 0-4hr" = "#2166AC", "Notable" = "purple"))+
  scale_x_continuous(breaks = c(0,4,8,24)) +
  xlab("\nHours IL-1a/PDGF-BB Exposure") +
  ylab("LogFC vs. 0hrs") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size =15))

trial$geneSetII <- "Repressed"
trial$geneSetII[grepl("Induced", trial$geneSet)] <- "Induced"

#mark notable lncs:
NotableLncs <- filter(fpkm_allGDE, EnsName %in% c("SMILR", "MIR3142HG") | EnsID %in% c("MSTRG.24277"))$EnsID

trial$geneSetIII[trial$EnsID %in% NotableLncs] <- "Notable"

RColorBrewer::brewer.pal(n = 8, name = "RdBu")

ggplot(trial) + aes(x = variable, y = value, color = geneSetII) + 
  geom_line(aes(x = variable, y = value, group = EnsID, color = geneSetII, alpha = 0.5)) +
  geom_line(data = filter(trial, !is.na(geneSetIII)), size = 0.8,
            aes(x = variable, y = value, group = EnsID, color = geneSetIII, alpha = 0.5)) +
  scale_color_manual(values = c("Induced" = "#D6604D", "Repressed" = "#67A9CF", "Notable" = "purple")) +
  scale_x_continuous(breaks = c(0,4,8,24)) +
  xlab("\nHours IL-1a/PDGF-BB\nExposure") +
  ylab("LogFC vs. 0hrs") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size =15))


#IEGs/CCs etc:
fpkm_IEGDE <- filter(fpkm_allGDE, IEG == "IEG")

#baseline 0hr
colnames(fpkm_IEGDE)
fpkm_IEGDE <- fpkm_IEGDE[,c(1,30,36,38)]

colnames(fpkm_IEGDE) <- c("EnsID", "4", "8", "24")

fpkm_IEGDE$`0` <- 0

library(reshape2)

trial <- melt(fpkm_IEGDE)

trial$variable <- as.numeric(as.character(trial$variable))

trial$geneSetII <- "Repressed"
trial$geneSetII[trial$EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Upwithin_24$EnsID)] <- "Induced"

ggplot(trial) + aes(x = variable, y = value, color = geneSetII) + 
  geom_line(aes(x = variable, y = value, group = EnsID, color = geneSetII, alpha = 0.5)) +
  scale_color_manual(values = c("Induced" = "#D6604D", "Repressed" = "#67A9CF", "Notable" = "purple")) +
  scale_x_continuous(breaks = c(0,4,8,24)) +
  xlab("\nHours IL-1a/PDGF-BB\nExposure") +
  ylab("LogFC vs. 0hrs") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size =15))


#IEGs/CCs etc:
fpkm_IEGDE <- filter(fpkm_allGDE, grepl("CC", GeneClassUpdate))

#baseline 0hr
colnames(fpkm_IEGDE)
fpkm_IEGDE <- fpkm_IEGDE[,c(1,30,36,38)]

colnames(fpkm_IEGDE) <- c("EnsID", "4", "8", "24")

fpkm_IEGDE$`0` <- 0

library(reshape2)

trial <- melt(fpkm_IEGDE)

trial$variable <- as.numeric(as.character(trial$variable))

trial$geneSetII <- "Repressed"
trial$geneSetII[trial$EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Upwithin_24$EnsID)] <- "Induced"

ggplot(trial) + aes(x = variable, y = value, color = geneSetII) + 
  geom_line(aes(x = variable, y = value, group = EnsID, color = geneSetII, alpha = 0.5)) +
  scale_color_manual(values = c("Induced" = "#D6604D", "Repressed" = "#67A9CF", "Notable" = "purple")) +
  scale_x_continuous(breaks = c(0,4,8,24)) +
  xlab("\nHours IL-1a/PDGF-BB\nExposure") +
  ylab("LogFC vs. 0hrs") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size =15))
