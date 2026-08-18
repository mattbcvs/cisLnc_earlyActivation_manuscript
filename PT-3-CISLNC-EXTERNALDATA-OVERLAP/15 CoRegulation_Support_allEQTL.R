#### overview of all p vals from all tests ####
library(dplyr)
library(ggplot2)

#import p GTEX:
GTEX_csv <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/400kbp (outOfDate)/", pattern = "GTEX_eQTL*", full.names = T)
GTEX_csv_names <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/400kbp (outOfDate)/", pattern = "GTEX_eQTL*", full.names = F)

GTEX_csv_names <- paste(sapply(strsplit(GTEX_csv_names, "_"), "[[", 3), sapply(strsplit(GTEX_csv_names, "_"), "[[", 4), sep = "_")
GTEX_csv_names <- gsub("_df.csv", "", GTEX_csv_names)

GTEX_csv <- lapply(GTEX_csv, read.csv)
names(GTEX_csv) <- GTEX_csv_names

GTEX_runs <- bind_rows(GTEX_csv, .id = "Tissue_pair")

#import p shu:
GTEX_csv <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Shu_SMC_Biobank/", pattern = "Shu_*", full.names = T)
GTEX_csv_names <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Shu_SMC_Biobank/", pattern = "Shu_*", full.names = F)

GTEX_csv_names <- paste(sapply(strsplit(GTEX_csv_names, "_"), "[[", 5), sapply(strsplit(GTEX_csv_names, "_"), "[[", 4), sep = "_")
GTEX_csv_names <- gsub(".csv_df", "", GTEX_csv_names)
GTEX_csv_names <- gsub("df", "", GTEX_csv_names)
GTEX_csv_names <- gsub(".csv_", "", GTEX_csv_names)

GTEX_csv <- lapply(GTEX_csv, read.csv)
names(GTEX_csv) <- GTEX_csv_names

Shu_runs <- bind_rows(GTEX_csv, .id = "Tissue_pair")


#remove <5 eQTL supported pairs in background from considerations
#Shu_runs <- filter(Shu_runs, V5 > 5)
#GTEX_runs <- filter(GTEX_runs, c > 5)

#combine to tissue, pairs, eQTL p thresh, fisher's p
Shu_runs$Tissue <- "SMC_Biobank"
colnames(Shu_runs)

Shu_runs_ <- Shu_runs[,c(13,1,10,3,5:8)]
Shu_runs_$CClncRNA_type <- "Same/later"
Shu_runs_$CClncRNA_type[grepl("same", Shu_runs_$Tissue_pair)] <- "Same"

Shu_runs_$CClncRNA_eQTL_overlap <- "Promoter"
Shu_runs_$CClncRNA_eQTL_overlap[grepl("exon", Shu_runs_$Tissue_pair)] <- "Exon"
Shu_runs_$CClncRNA_eQTL_overlap[grepl("promExon", Shu_runs_$Tissue_pair)] <- "Promoter/Exon"
Shu_runs_$CClncRNA_eQTL_overlap[grepl("locus", Shu_runs_$Tissue_pair)] <- "Whole locus"
Shu_runs_$CClncRNA_eQTL_overlap[grepl("intron", Shu_runs_$Tissue_pair)] <- "Intron"
Shu_runs_$CClncRNA_eQTL_overlap[grepl("splice", Shu_runs_$Tissue_pair)] <- "Splice junction"
Shu_runs_$CClncRNA_eQTL_overlap[grepl("TTS", Shu_runs_$Tissue_pair)] <- "TTS"


colnames(GTEX_runs)
GTEX_runs_ <- GTEX_runs[,c(9,1,2,3,5:8)]
GTEX_runs_$CClncRNA_type <- "Same/later"
GTEX_runs_$CClncRNA_type[grepl("same", GTEX_runs_$Tissue_pair)] <- "Same"

GTEX_runs_$CClncRNA_eQTL_overlap <- "Promoter"
GTEX_runs_$CClncRNA_eQTL_overlap[grepl("exon", GTEX_runs_$Tissue_pair)] <- "Exon"
GTEX_runs_$CClncRNA_eQTL_overlap[grepl("promoterExon", GTEX_runs_$Tissue_pair)] <- "Promoter/Exon"
GTEX_runs_$CClncRNA_eQTL_overlap[grepl("locus", GTEX_runs_$Tissue_pair)] <- "Whole locus"
GTEX_runs_$CClncRNA_eQTL_overlap[grepl("Intron", GTEX_runs_$Tissue_pair)] <- "Intron"
GTEX_runs_$CClncRNA_eQTL_overlap[grepl("Splice", GTEX_runs_$Tissue_pair)] <- "Splice junction"
GTEX_runs_$CClncRNA_eQTL_overlap[grepl("TTS", GTEX_runs_$Tissue_pair)] <- "TTS"

GTEX_runs_$ThreshP <- as.numeric(sapply(sapply(GTEX_runs_$Run, strsplit, "_"), "[[", 2))
Shu_runs_$ThreshP <- as.numeric(Shu_runs_$ThreshP)

colnames(GTEX_runs_)
GTEX_runs_ <- GTEX_runs_[,c(1,9,10,11,4,5:8)]

colnames(Shu_runs_)
Shu_runs_ <- Shu_runs_[,c(1,9,10,3,4,5:8)]

colnames(GTEX_runs_) <- colnames(Shu_runs_)

all_runs <- rbind(Shu_runs_, GTEX_runs_)

#linegraph, all tests, eQTL p vs. Fisher's p:
#is overall trend for increasing Fisher's with more stringency
all_runs$Tissue_pair_overlap <- paste(all_runs$Tissue, all_runs$CClncRNA_type, all_runs$CClncRNA_eQTL_overlap, sep = "_")

ggplot(all_runs) + aes(x = -log10(ThreshP), y = -log10(V1), color = Tissue_pair_overlap) +
  geom_line(alpha = 0.2) +
  theme(legend.position = "none")

#via boxplots is maybe easier to communicate:
all_runs$ThreshP_char <- as.factor(as.character(all_runs$ThreshP))
all_runs$ThreshP_char <- factor(all_runs$ThreshP_char, levels = levels(all_runs$ThreshP_char)[c(4,3,2,1,5,10,6,7,8,9)])

ggplot(filter(all_runs, !Tissue %in% "SMC_Biobank")) + aes(x = ThreshP_char, y = -log10(V1)) +
  geom_beeswarm(alpha = 0.5, size = 4) +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.5) +
  xlab("eQTL-eGene Confidence Threshold (p)") +
  ylab("eQTL-eGene Enrichment in\nCClncRNA-Targets\n(-log10p)") +
  theme_minimal() + 
  theme(axis.title.y = element_text(angle =0, vjust = 0.5),
        legend.position = "none",
        text = element_text(size =24))

#just tests where greater than 10 found in all tests:
trial <- split(all_runs, all_runs$Tissue_pair_overlap)
triali <- sapply(trial, function(x){
  sum(x$V5 >10) == 6
})

table(triali)

triali <- data.frame("test" = names(triali), triali)
triali <- filter(triali, triali == TRUE)

trialii <- filter(all_runs, Tissue_pair_overlap %in% triali$test)

#just GTEX or SMC biobank will stick out
ggplot(filter(trialii, !Tissue %in% "SMC_Biobank")) + aes(x = -log10(ThreshP), y = -log10(V1), color = Tissue_pair_overlap) +
  geom_line(alpha = 0.4, linewidth = 1) +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(filter(trialii, !Tissue %in% "SMC_Biobank")) + aes(x = ThreshP_char, y = -log10(V1)) +
  geom_beeswarm(alpha = 0.5, size = 4) +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.5) +
  xlab("eQTL-eGene Confidence Threshold (p)") +
  ylab("eQTL-eGene Enrichment in\nCClncRNA-Targets\n(-log10p)") +
  theme_minimal() + 
  theme(axis.title.y = element_text(angle =0, vjust = 0.5),
        legend.position = "none",
        text = element_text(size =24))

#split to same and same/later for write-up
ggplot(filter(trialii, CClncRNA_type == "Same/later", !Tissue %in% "SMC_Biobank")) + aes(x = ThreshP_char, y = -log10(V1)) +
  geom_beeswarm(alpha = 0.5, size = 4) +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.5) +
  xlab("eQTL-eGene Confidence Threshold (p)") +
  ylab("eQTL-eGene Enrichment in\nCClncRNA-Targets\n(-log10p)") +
  theme_minimal() + 
  theme(axis.title.y = element_text(angle =0, vjust = 0.5),
        legend.position = "none",
        text = element_text(size =24))
ggplot(filter(trialii, CClncRNA_type == "Same", !Tissue %in% "SMC_Biobank")) + aes(x = ThreshP_char, y = -log10(V1)) +
  geom_beeswarm(alpha = 0.5, size = 4) +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.5) +
  xlab("eQTL-eGene Confidence Threshold (p)") +
  ylab("eQTL-eGene Enrichment in\nCClncRNA-Targets\n(-log10p)") +
  theme_minimal() + 
  theme(axis.title.y = element_text(angle =0, vjust = 0.5),
        legend.position = "none",
        text = element_text(size =24))

#just tests where something is significant:
trial <- split(trialii, trialii$Tissue_pair_overlap)
triali <- sapply(trial, function(x){
  min(x$V1) < 0.05
})

table(triali)

triali <- data.frame("test" = names(triali), triali)
triali <- filter(triali, triali == TRUE)

trialiii <- filter(trialii, Tissue_pair_overlap %in% triali$test)

ggplot(trialiii) + aes(x = -log10(ThreshP), y = -log10(V1), color = Tissue_pair_overlap) +
  geom_line(alpha = 0.4, linewidth = 1) +
  theme_minimal() + theme(legend.position = "none")

#just GTEX or SMC biobank will stick out
ggplot(filter(trialiii, !Tissue %in% "SMC_Biobank")) + aes(x = -log10(ThreshP), y = -log10(V1), color = Tissue_pair_overlap) +
  geom_line(alpha = 0.4, linewidth = 1) +
  theme_minimal() + theme(legend.position = "none")

#grey only
ggplot(filter(trialiii, !Tissue %in% "SMC_Biobank")) + aes(x = -log10(ThreshP), y = -log10(V1), group = Tissue_pair_overlap) +
  geom_line(alpha = 0.6, linewidth = 1, color = "grey50") +
  theme_minimal() + theme(legend.position = "none")

#via boxplots is maybe easier to communicate:
trialiii$ThreshP_char <- as.factor(as.character(trialiii$ThreshP))
trialiii$ThreshP_char <- factor(trialiii$ThreshP_char, levels = levels(trialiii$ThreshP_char)[c(4,3,2,1,5,10,6,7,8,9)])

ggplot(filter(trialiii, !Tissue %in% "SMC_Biobank")) + aes(x = ThreshP_char, y = -log10(V1)) +
  geom_beeswarm(alpha = 0.5, size = 4) +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.5) +
  xlab("eQTL-eGene Confidence Threshold (p)") +
  ylab("eQTL-eGene Enrichment in\nCClncRNA-Targets\n(-log10p)") +
  theme_minimal() + 
  theme(axis.title.y = element_text(angle =0, vjust = 0.5),
        legend.position = "none",
        text = element_text(size =24))

#select best p per Tissue_pair_overlap:
all_runs$Tissue_pair_overlap <- paste(all_runs$Tissue, all_runs$CClncRNA_type, all_runs$CClncRNA_eQTL_overlap, sep = "_")
trial <- split(all_runs, all_runs$Tissue_pair_overlap)

triali <- lapply(trial, function(x){
  x[order(x$V1, decreasing = F),][1,]
})

triali <- bind_rows(triali)

best_runs <- triali

#consider those with low eQTL-support in background as n.s.
best_runs$V1[best_runs$V5 <6] <- 1
best_runs$V1[is.na(best_runs$V1)] <- 1

#transform p:
best_runs$pTransform <- -log10(best_runs$V1)

#plot all
best_runs <- unique(best_runs[,-c(4:10)])

trial <- tidyr::pivot_wider(best_runs, names_from = "Tissue", values_from = "pTransform", values_fill = 0)

triali <- as.matrix(trial[,3:52])
rownames(triali) <- paste(trial$CClncRNA_eQTL_overlap, trial$CClncRNA_type, sep = "_")

#order by overlap type:
rownames(triali)[c(3,10,
                   1,8,
                   4,11,
                   2,9,
                   5,12,
                   6,13,
                   7,14)]
triali <- triali[c(3,10,
                   1,8,
                   4,11,
                   2,9,
                   5,12,
                   6,13,
                   7,14),]

#order cols by strongest fisher's value:
triali <- triali[,match(unique(best_runs$Tissue[order(best_runs$pTransform, decreasing = T)]), colnames(triali))]

pheatmap::pheatmap(triali, gaps_row = c(2,4,6,8,10,12),
                   cluster_cols = F, cluster_rows = F, scale = "none", angle_col = 315, 
                   colorRampPalette(c("white", "red"))(500))

pheatmap::pheatmap(triali[c(1,3,5,7,9,11,13),],
                   cluster_cols = F, cluster_rows = F, scale = "none", angle_col = 315, 
                   colorRampPalette(c("white", "red"))(500))

pheatmap::pheatmap(triali[c(1,3,7,9,11,13),],
                   cluster_cols = F, cluster_rows = F, scale = "none", angle_col = 315, 
                   colorRampPalette(c("white", "red"))(500))


#plot simpler p vals - the above implies that lots of strong associations may be found
#whereas most do not reach above 0.05 significance
#also hard to tell that exon and promoter do better than intron or TTS
best_runs_ <- best_runs
best_runs_$pTransformII <- "n.s"
best_runs_$pTransformII[best_runs_$pTransform > -log10(0.1)] <- "0.1"
best_runs_$pTransformII[best_runs_$pTransform > -log10(0.05)] <- "0.05"
best_runs_$pTransformII[best_runs_$pTransform > -log10(0.01)] <- "0.01"
best_runs_$pTransformII[best_runs_$pTransform > -log10(0.005)] <- "0.005"
best_runs_$pTransformII[best_runs_$pTransform > -log10(0.001)] <- "0.001"

best_runs_ <- best_runs_[,-4]

best_runs_$rowNames <- paste(best_runs_$CClncRNA_type, best_runs_$CClncRNA_eQTL_overlap, sep = "_")

#order by overlap type:
best_runs_$rowNames <- as.factor(best_runs_$rowNames)
best_runs_$rowNames <- factor(best_runs_$rowNames, 
                              levels = levels(best_runs_$rowNames)[rev(c(3,10,
                                                                         1,8,
                                                                         4,11,
                                                                         2,9,
                                                                         5,12,
                                                                         6,13,
                                                                         7,14))])

best_runs_$Tissue <- as.factor(best_runs_$Tissue)
best_runs_$Tissue <- factor(best_runs_$Tissue, 
                            levels = levels(best_runs_$Tissue)[match(unique(best_runs$Tissue[order(best_runs$pTransform, decreasing = T)]), 
                                                                     levels(best_runs_$Tissue))])
#best_runs_ <- filter(best_runs_, !is.na(pTransformII))

ggplot(best_runs_) + aes(y = rowNames, x = Tissue, fill = pTransformII) +
  geom_tile(colour = "grey50"
  ) +
  scale_fill_manual(values = c("n.s" = "white", 
                               "0.1" = "lightsteelblue1",
                               "0.05" = "steelblue1",
                               "0.01" = "royalblue1",
                               "0.005" = "royalblue3",
                               "0.001" = "royalblue4")) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  theme(#panel.grid.major =  element_blank(),
    axis.text.x = element_text(size=9, angle =315, hjust = 0))


#just same/later:
ggplot(filter(best_runs_, CClncRNA_type == "Same/later")) + aes(y = rowNames, x = Tissue, fill = pTransformII) +
  geom_tile(colour = "grey80"
  ) +
  scale_fill_manual(values = c("n.s" = "white", 
                               "0.1" = "lightsteelblue1",
                               "0.05" = "steelblue1",
                               "0.01" = "royalblue1",
                               "0.005" = "royalblue3",
                               "0.001" = "royalblue4")) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  theme(panel.grid.major =  element_blank(),
        axis.text.x = element_text(size=9, angle =315, hjust = 0))

#just same/later, no promoter/exon:
ggplot(filter(best_runs_, !CClncRNA_eQTL_overlap == "Promoter/Exon", CClncRNA_type == "Same/later")) + aes(y = rowNames, x = Tissue, fill = pTransformII) +
  geom_tile(colour = "grey80"
  ) +
  scale_fill_manual(values = c("n.s" = "white", 
                               "0.1" = "#C6DBEF",
                               "0.05" = "#9ECAE1",
                               "0.01" = "#6BAED6",
                               "0.005" = "#3182BD",
                               "0.001" = "#08519C")) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  theme(panel.grid.major =  element_blank(),
        axis.text.x = element_text(size=9, angle =315, hjust = 0))

#just same, no promoter/exon:
ggplot(filter(best_runs_, !CClncRNA_eQTL_overlap == "Promoter/Exon", CClncRNA_type == "Same")) + aes(y = rowNames, x = Tissue, fill = pTransformII) +
  geom_tile(colour = "grey80"
  ) +
  scale_fill_manual(values = c("n.s" = "white", 
                               "0.1" = "#C6DBEF",
                               "0.05" = "#9ECAE1",
                               "0.01" = "#6BAED6",
                               "0.005" = "#3182BD",
                               "0.001" = "#08519C")) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  theme(panel.grid.major =  element_blank(),
        axis.text.x = element_text(size=9, angle =315, hjust = 0))


#### p adjust format exploring: ####

#mass p adjust first:
summary(p.adjust(all_runs$V1, method = "BH"))
#on only viable runs:
summary(p.adjust(filter(all_runs, V5 >5)$V1, method = "BH"))
#on only viable runs, per each set of cclnc-target pairs:
summary(p.adjust(filter(all_runs, CClncRNA_type == "Same", V5 >5)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, CClncRNA_type == "Same/later", V5 >5)$V1, method = "BH"))
#on only more viable runs, per each set of cclnc-target pairs:
summary(p.adjust(filter(all_runs, CClncRNA_type == "Same", V5 >10)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, CClncRNA_type == "Same/later", V5 >10)$V1, method = "BH"))
#on only more viable runs, per each set of cclnc-target pairs:
summary(p.adjust(filter(all_runs, CClncRNA_type == "Same", V5 >20)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, CClncRNA_type == "Same/later", V5 >20)$V1, method = "BH"))

#focusing on obvious overlap sites of interest:
summary(p.adjust(filter(all_runs, CClncRNA_eQTL_overlap %in% c("Promoter", "Exon"), 
                        CClncRNA_type == "Same", V5 >5)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, CClncRNA_eQTL_overlap %in% c("Promoter", "Exon"), 
                        CClncRNA_type == "Same/later", V5 >5)$V1, method = "BH"))
#more viable
summary(p.adjust(filter(all_runs, CClncRNA_eQTL_overlap %in% c("Promoter", "Exon"), 
                        CClncRNA_type == "Same", V5 >10)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, CClncRNA_eQTL_overlap %in% c("Promoter", "Exon"), 
                        CClncRNA_type == "Same/later", V5 >10)$V1, method = "BH"))
#even more viable
summary(p.adjust(filter(all_runs, CClncRNA_eQTL_overlap %in% c("Promoter", "Exon"), 
                        CClncRNA_type == "Same", V5 >20)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, CClncRNA_eQTL_overlap %in% c("Promoter", "Exon"), 
                        CClncRNA_type == "Same/later", V5 >20)$V1, method = "BH"))

#less pThresh
summary(p.adjust(filter(all_runs, ThreshP %in% c(0.05, 1e-06, 1e-08, 0.001, 0.00001), CClncRNA_eQTL_overlap %in% c("Promoter", "Exon"), 
                        CClncRNA_type == "Same", V5 >5)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, ThreshP %in% c(0.05, 1e-06, 1e-08, 0.001, 0.00001), CClncRNA_eQTL_overlap %in% c("Promoter", "Exon"), 
                        CClncRNA_type == "Same/later", V5 >5)$V1, method = "BH"))
#more viable
summary(p.adjust(filter(all_runs, ThreshP %in% c(0.05, 1e-06, 1e-08, 0.001, 0.00001), CClncRNA_eQTL_overlap %in% c("Promoter", "Exon"), 
                        CClncRNA_type == "Same", V5 >10)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, ThreshP %in% c(0.05, 1e-06, 1e-08, 0.001, 0.00001), CClncRNA_eQTL_overlap %in% c("Promoter", "Exon"), 
                        CClncRNA_type == "Same/later", V5 >10)$V1, method = "BH"))

#no pThresh
summary(p.adjust(filter(all_runs, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Promoter", "Exon"), 
                        CClncRNA_type == "Same", V5 >5)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Promoter", "Exon"), 
                        CClncRNA_type == "Same/later", V5 >5)$V1, method = "BH"))
#more viable
summary(p.adjust(filter(all_runs, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Promoter", "Exon"), 
                        CClncRNA_type == "Same", V5 >10)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Promoter", "Exon"), 
                        CClncRNA_type == "Same/later", V5 >10)$V1, method = "BH"))

#per just promoter:
summary(p.adjust(filter(all_runs, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Promoter"), 
                        CClncRNA_type == "Same", V5 >5)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Promoter"), 
                        CClncRNA_type == "Same/later", V5 >5)$V1, method = "BH"))
#more viable
summary(p.adjust(filter(all_runs, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Promoter"), 
                        CClncRNA_type == "Same", V5 >10)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Promoter"), 
                        CClncRNA_type == "Same/later", V5 >10)$V1, method = "BH"))
#per just exon
summary(p.adjust(filter(all_runs, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Exon"), 
                        CClncRNA_type == "Same", V5 >5)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Exon"), 
                        CClncRNA_type == "Same/later", V5 >5)$V1, method = "BH"))
#more viable
summary(p.adjust(filter(all_runs, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Exon"), 
                        CClncRNA_type == "Same", V5 >10)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Exon"), 
                        CClncRNA_type == "Same/later", V5 >10)$V1, method = "BH"))

#all of the above are done across all tissues, doing per tissue is justifiable as each can be considered a seperate dataset
#with diff power/rep (as will be shown in the supplement)
# this is last resort, ideally will adjust across all viable tissues, some may be obviously non viable:

#could filter datasets to just more powerful? so those with high number of eGenes/genes or expressed SVSMC PCGs near lncs:

#number of expressed SVSCM PCGs near lncs here (any of GTEX tables)
GTEX_eQTL_TTSFish_df <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/GTEX_eQTL_TTSFish_df.csv")

pairedPCG_eGenes_PerTissue <- unique(filter(GTEX_eQTL_TTSFish_df, Run == "Run_0.05")[,c(8,7)])
pairedPCG_eGenes_PerTissue$d[is.na(pairedPCG_eGenes_PerTissue$d)] <- 0
#add SMC biobank:
filter(Shu_runs, Tissue_pair == "exon", ThreshP == "0.05")
pairedPCG_eGenes_PerTissue <- rbind(pairedPCG_eGenes_PerTissue, data.frame("tissue" = "SMC_Biobank", "d" = 836))

pairedPCG_eGenes_PerTissue <- pairedPCG_eGenes_PerTissue[order(pairedPCG_eGenes_PerTissue$tissue),]
pairedPCG_eGenes_PerTissue$tissue <- factor(pairedPCG_eGenes_PerTissue$tissue)
pairedPCG_eGenes_PerTissue$tissue <- factor(pairedPCG_eGenes_PerTissue$tissue, 
                                            levels = levels(pairedPCG_eGenes_PerTissue$tissue)[order(pairedPCG_eGenes_PerTissue$d)])

ggplot(pairedPCG_eGenes_PerTissue) + aes(x = tissue, y = d) +
  geom_bar(stat = "identity") + Seurat::RotatedAxis() +
  ylab("PCGs near SVSMC lncRNAs\ndetected >1 TPM")
#3x tissues where no CClncRNA targets are eGenes:
#SMC biobank obvious outlier

#data from GTEX:
GTEX_info <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/eQTL_dataset_info.csv")
GTEX_info$Tissue <- gsub(" - ", "_", GTEX_info$Tissue)
GTEX_info$Tissue <- gsub(" ", "_", GTEX_info$Tissue)
GTEX_info$Tissue <- gsub("\\(", "", GTEX_info$Tissue)
GTEX_info$Tissue <- gsub("\\)", "", GTEX_info$Tissue)
GTEX_info$Tissue <- gsub("_q", "_", GTEX_info$Tissue)

GTEX_info <- merge(GTEX_info, pairedPCG_eGenes_PerTissue, by.x = "Tissue", by.y = "tissue", all.y = T)

colnames(GTEX_info) <- gsub("X\\.\\.", "", colnames(GTEX_info))

#fill in blanks for Shu dataset:
GTEX_info$RNASeq.and.Genotyped.samples[GTEX_info$Tissue == "SMC_Biobank"] <- 1486
#total no. eGenes? from Table S3: 42,257 total
#GTEX is likely just reporting PCG and lncRNA (GENCODE)
16672+12154
GTEX_info$eGenes[GTEX_info$Tissue == "SMC_Biobank"] <- 28826
#52,271 or for PCG/lncs:
19549+15207
GTEX_info$eGenes.Total.Expressed.Genes[GTEX_info$Tissue == "SMC_Biobank"] <- 28826/34756

#eGenes/total genes indicates how powerful the dataset is at finding eQTLs:
#presumably gets better with sample no.:
ggplot(GTEX_info) + aes(x = RNASeq.and.Genotyped.samples, y = eGenes.Total.Expressed.Genes) +
  geom_point()

#yes v. much so

#datasets with higher power also get more of the CClnc targets?
ggplot(GTEX_info) + aes(x = d, y = eGenes.Total.Expressed.Genes) +
  geom_point() +
  ylab("eQTL dataset power\n(eGenes/Expressed genes)") +
  xlab("No. PCGs near lncRNAs in SVSMC timecourse\ndetected >1 TPM (per each eQTL dataset)")

#so leading on from this, best results from best powered tissues?
unique(all_runs$Tissue[order(all_runs$V1, decreasing = F)])
bestTissues <- unique(filter(all_runs, V1 <0.01)$Tissue)

#simpler to just label tissues which work well in final heatmap after padjust:
bestTissues <- unique(all_runs$Tissue)[c(1,2,3,5,41,43,36,40,21)]
GTEX_info$label <- NA
GTEX_info$label[GTEX_info$Tissue %in% bestTissues] <- "Predictive value for\nSVSMC CClncRNA-targets"

ggplot(GTEX_info) + aes(x = d, y = eGenes.Total.Expressed.Genes, color = label) +
  geom_point(size = 3, alpha =0.8) +
  ggrepel::geom_label_repel(data = filter(GTEX_info, !is.na(label)),
                            aes(x = d, y = eGenes.Total.Expressed.Genes, color = label, label = Tissue), 
                            nudge_x = 120, nudge_y = -0.1, force = 50) +
  theme_minimal() +
  scale_x_continuous(limits = c(0,1010)) +
  theme(text = element_text(size =21),
        legend.position = "none") +
  ylab("eQTL detection power\n(eGenes/Expressed genes)") +
  xlab("No. CClncRNA neighbours\ndetected as eGene")

ggplot(GTEX_info) + aes(x = RNASeq.and.Genotyped.samples, y = eGenes.Total.Expressed.Genes, color = label) +
  geom_point(size = 3, alpha = 0.8) +
  ggrepel::geom_label_repel(data = filter(GTEX_info, !is.na(label)),
                            aes(x = RNASeq.and.Genotyped.samples, y = eGenes.Total.Expressed.Genes, color = label, label = Tissue), 
                            nudge_x = 500, nudge_y = -0.1, force = 10) +
  theme_minimal() +
  theme(text = element_text(size =21),
        legend.position = "none") +
  ylab("eQTL detection power\n(eGenes/Expressed genes)") +
  xlab("No. samples\n")

powerTissues <- filter(GTEX_info, eGenes.Total.Expressed.Genes >0.5)$Tissue

#on only viable runs, per each set of cclnc-target pairs:
summary(p.adjust(filter(all_runs, Tissue %in% powerTissues, CClncRNA_type == "Same", V5 >5)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, Tissue %in% powerTissues, CClncRNA_type == "Same/later", V5 >5)$V1, method = "BH"))
#on only more viable runs, per each set of cclnc-target pairs:
summary(p.adjust(filter(all_runs, Tissue %in% powerTissues, CClncRNA_type == "Same", V5 >10)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, Tissue %in% powerTissues, CClncRNA_type == "Same/later", V5 >10)$V1, method = "BH"))
#on only more viable runs, per each set of cclnc-target pairs:
summary(p.adjust(filter(all_runs, Tissue %in% powerTissues, CClncRNA_type == "Same", V5 >20)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, Tissue %in% powerTissues, CClncRNA_type == "Same/later", V5 >20)$V1, method = "BH"))


#per just promoter:
summary(p.adjust(filter(all_runs, Tissue %in% powerTissues, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Promoter"), 
                        CClncRNA_type == "Same", V5 >5)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, Tissue %in% powerTissues, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Promoter"), 
                        CClncRNA_type == "Same/later", V5 >5)$V1, method = "BH"))
#more viable
summary(p.adjust(filter(all_runs, Tissue %in% powerTissues, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Promoter"), 
                        CClncRNA_type == "Same", V5 >10)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, Tissue %in% powerTissues, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Promoter"), 
                        CClncRNA_type == "Same/later", V5 >10)$V1, method = "BH"))
#per just exon
summary(p.adjust(filter(all_runs, Tissue %in% powerTissues, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Exon"), 
                        CClncRNA_type == "Same", V5 >5)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, Tissue %in% powerTissues, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Exon"), 
                        CClncRNA_type == "Same/later", V5 >5)$V1, method = "BH"))
#more viable
summary(p.adjust(filter(all_runs, Tissue %in% powerTissues, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Exon"), 
                        CClncRNA_type == "Same", V5 >10)$V1, method = "BH"))
summary(p.adjust(filter(all_runs, Tissue %in% powerTissues, ThreshP == 0.05, CClncRNA_eQTL_overlap %in% c("Exon"), 
                        CClncRNA_type == "Same/later", V5 >10)$V1, method = "BH"))

#so basically, the adjustment HAS to be per tissue/overlap/pair - only option left
#cannot be a mass p adjust across several tissues
#follow up by comparing to FANTOM eQTL required
#edit, the good part is that this lines up v. nicely with FANTOM eQTL-PCG pairs


#### p adjust per overlap/tissue/pair ####

#select best p per Tissue_pair_overlap:
all_runs$Tissue_pair_overlap <- paste(all_runs$Tissue, all_runs$CClncRNA_type, all_runs$CClncRNA_eQTL_overlap, sep = "_")
trial <- split(all_runs, all_runs$Tissue_pair_overlap)

triali <- lapply(trial, function(x){
  y <- filter(x, V5 >10)
  y$V1_padj <- p.adjust(y$V1, method = "BH")
  return(y[order(y$V1, decreasing = F),][1,])
})

triali <- bind_rows(triali)

best_runs <- filter(triali, !is.na(V1_padj))

#transform p:
best_runs$pTransform <- -log10(best_runs$V1_padj)

#plot all, save a version
best_runs_full <- best_runs
best_runs <- unique(best_runs[,-c(4:11)])

trial <- tidyr::pivot_wider(best_runs, names_from = "Tissue", values_from = "pTransform", values_fill = 0)

triali <- as.matrix(trial[,3:47])
rownames(triali) <- paste(trial$CClncRNA_eQTL_overlap, trial$CClncRNA_type, sep = "_")

#order by overlap type:
rownames(triali)[c(3,10,
                   1,8,
                   4,11,
                   2,9,
                   5,14,
                   6,12,
                   7,13)]
triali <- triali[c(3,10,
                   1,8,
                   4,11,
                   2,9,
                   5,14,
                   6,12,
                   7,13),]

#order cols by strongest fisher's value:
triali <- triali[,match(unique(best_runs$Tissue[order(best_runs$pTransform, decreasing = T)]), colnames(triali))]

pheatmap::pheatmap(triali, gaps_row = c(2,4,6,8,10,12),
                   cluster_cols = F, cluster_rows = F, scale = "none", angle_col = 315, 
                   colorRampPalette(c("white", "red"))(500))

pheatmap::pheatmap(triali[c(1,3,5,7,9,11,13),],
                   cluster_cols = F, cluster_rows = F, scale = "none", angle_col = 315, 
                   colorRampPalette(c("white", "red"))(500))

pheatmap::pheatmap(triali[c(1,3,7,9,11,13),],
                   cluster_cols = F, cluster_rows = F, scale = "none", angle_col = 315, 
                   colorRampPalette(c("white", "red"))(500))


#plot simpler p vals - the above implies that lots of strong associations may be found
#whereas most do not reach above 0.05 significance
#also hard to tell that exon and promoter do better than intron or TTS
best_runs_ <- best_runs
best_runs_$pTransformII <- NA
#best_runs_$pTransformII[best_runs_$pTransform > -log10(0.1)] <- "0.1"
best_runs_$pTransformII[best_runs_$pTransform > -log10(0.05)] <- "<0.05"
best_runs_$pTransformII[best_runs_$pTransform > -log10(0.01)] <- "<0.01"
best_runs_$pTransformII[best_runs_$pTransform > -log10(0.005)] <- "<0.005"
best_runs_$pTransformII[best_runs_$pTransform > -log10(0.001)] <- "<0.001"

best_runs_ <- best_runs_[,-4]

best_runs_$rowNames <- paste(best_runs_$CClncRNA_type, best_runs_$CClncRNA_eQTL_overlap, sep = "_")

#remove tissues with no predictive value:
goodTissues <- unique(filter(best_runs_, !is.na(pTransformII))$Tissue)

best_runs_ <- filter(best_runs_, Tissue %in% goodTissues)

#make wider to put NA values back in where needed
trial <- tidyr::pivot_wider(best_runs_, names_from = "Tissue", values_from = "pTransformII")
#melt to get back in right format
trial <- reshape2::melt(trial, id.vars = colnames(trial)[1:3])
trial$value[is.na(trial$value)] <- "n.s"

best_runs_ <- trial
colnames(best_runs_)[4:5] <- c("Tissue", "pTransformII")

#order by overlap type:
best_runs_$rowNames <- as.factor(best_runs_$rowNames)
best_runs_$rowNames <- factor(best_runs_$rowNames, 
                              levels = levels(best_runs_$rowNames)[rev(c(3,10,
                                                                         1,8,
                                                                         4,11,
                                                                         2,9,
                                                                         5,12,
                                                                         6,13,
                                                                         7,14))])

best_runs_$Tissue <- as.factor(best_runs_$Tissue)
best_runs_orderPlot <- filter(best_runs, Tissue %in% goodTissues)
best_runs_$Tissue <- factor(best_runs_$Tissue, 
                            levels = levels(best_runs_$Tissue)[match(unique(best_runs_orderPlot$Tissue[order(best_runs_orderPlot$pTransform, decreasing = T)]), 
                                                                     levels(best_runs_$Tissue))])
#best_runs_ <- filter(best_runs_, !is.na(pTransformII))

ggplot(best_runs_) + aes(y = rowNames, x = Tissue, fill = pTransformII) +
  geom_tile(colour = "grey50"
  ) +
  scale_fill_manual(values = c("n.s" = "white", 
                               "0.1" = "#C6DBEF",
                               "0.05" = "#9ECAE1",
                               "0.01" = "#6BAED6",
                               "0.005" = "#3182BD",
                               "0.001" = "#08519C")) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  theme(panel.grid.major =  element_blank(),
        axis.text.x = element_text(size=9, angle =315, hjust = 0))


#just same/later:
ggplot(filter(best_runs_, CClncRNA_type == "Same/later")) + aes(y = rowNames, x = Tissue, fill = pTransformII) +
  geom_tile(colour = "grey80"
  ) +
  scale_fill_manual(values = c("n.s" = "white", 
                               "0.1" = "#C6DBEF",
                               "0.05" = "#9ECAE1",
                               "0.01" = "#6BAED6",
                               "0.005" = "#3182BD",
                               "0.001" = "#08519C")) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  theme(panel.grid.major =  element_blank(),
        axis.text.x = element_text(size=9, angle =315, hjust = 0))

#just same/later, no promoter/exon:
best_runs_$CClncRNA_eQTL_overlap <- as.factor(best_runs_$CClncRNA_eQTL_overlap)
best_runs_$CClncRNA_eQTL_overlap <- factor(best_runs_$CClncRNA_eQTL_overlap, 
                                           levels = levels(best_runs_$CClncRNA_eQTL_overlap)[rev(c(3,1,4,2,5,6,7))])

ggplot(filter(best_runs_, !CClncRNA_eQTL_overlap == "Promoter/Exon", CClncRNA_type == "Same/later")) + 
  aes(y = CClncRNA_eQTL_overlap, x = Tissue, fill = pTransformII) +
  geom_tile(colour = "grey80"
  ) +
  scale_fill_manual(values = c("n.s" = "white", 
                               #"0.1" = "#C6DBEF",
                               "<0.05" = "#9ECAE1",
                               "<0.01" = "#6BAED6",
                               "<0.005" = "#3182BD",
                               "<0.001" = "#08519C")) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  theme(panel.grid.major =  element_blank(),
        axis.text.x = element_text(size=18, angle =315, hjust = 0),
        axis.text.y = element_text(size=18),
        legend.text = element_text(size=18),
        axis.ticks.x.bottom = element_line())

#just same, no promoter/exon:
ggplot(filter(best_runs_, !CClncRNA_eQTL_overlap == "Promoter/Exon", CClncRNA_type == "Same")) + 
  aes(y = CClncRNA_eQTL_overlap, x = Tissue, fill = pTransformII) +
  geom_tile(colour = "grey80"
  ) +
  scale_fill_manual(values = c("n.s" = "white", 
                               #"0.1" = "#C6DBEF",
                               "<0.05" = "#9ECAE1",
                               "<0.01" = "#6BAED6",
                               "<0.005" = "#3182BD",
                               "<0.001" = "#08519C")) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  theme(panel.grid.major =  element_blank(),
        axis.text.x = element_text(size=18, angle =315, hjust = 0),
        axis.text.y = element_text(size=18),
        legend.text = element_text(size=18),
        axis.ticks.x.bottom = element_line())


#### total eQTL-supported cclncRNA-target pairs from Shu/GTEX ####

#15x runs to use from GTEX or Shu
best_runs_toSelect_GTEX <- filter(best_runs_full, !Tissue == "SMC_Biobank", !CClncRNA_eQTL_overlap %in% c("Promoter/Exon"), V1_padj<0.05)

eQTLconfirmed_pairs <- list()

Pairs2Test_list <- list("Same/later" = CoRegPairs_04_48_24_extended,
                        "Same" = CoRegPairs_04_48_24_extendedSame)

for (i in 1:length(best_runs_toSelect_GTEX$Tissue)){
  
  triali <- AllLNC_AllPCG_2d3d
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(best_runs_toSelect$Tissue[i], tissueType), 
                                                        OverlapType == best_runs_toSelect$CClncRNA_eQTL_overlap[i], 
                                                        pval_nominal < best_runs_toSelect$ThreshP[i],
                                                        TotalOverlaps == 1)$pairs] <- "Yes"
  
  #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
  triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[best_runs_toSelect$Tissue[i]]], 
                                               pval_nominal < best_runs_toSelect$ThreshP[i])$gene_id)
  
  #remove anything un-expressed >1 TPM in this tissue:
  findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", best_runs_toSelect$Tissue[i])), 
                                                      gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
  table(findExprs)
  triali <- filter(triali, EnsID.y %in% GTEX_exprs[findExprs, 1])
  
  triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% Pairs2Test_list[[best_runs_toSelect$CClncRNA_type[i]]]$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  DELNC_DEPCG_1_eQTL$tissueName <- best_runs_toSelect$Tissue[i]
  DELNC_DEPCG_1_eQTL$overlapType <- best_runs_toSelect$CClncRNA_eQTL_overlap[i]
  DELNC_DEPCG_1_eQTL$pairType <- best_runs_toSelect$CClncRNA_type[i]
  eQTLconfirmed_pairs[[i]] <- unique(DELNC_DEPCG_1_eQTL[,c(3,20:22)])
}


#for Shu Biobank:
Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]

triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                 EnsName.y %in% filter(Shu_allVar_p, pvalue < best_runs_toSelect$ThreshP[15])$GeneName)

triali$eQTL_validated_tissue <- "No"
triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < best_runs_toSelect$ThreshP[15],
                                                            OverlapType == best_runs_toSelect$CClncRNA_eQTL_overlap[15], 
                                                            TotalOverlaps == 1)$pairs] <- "Yes"

#all potential targets for potential cisLncs
DELNC_DEPCG_1 <- filter(triali, pairs %in% Pairs2Test_list[[best_runs_toSelect$CClncRNA_type[15]]]$pairs)
DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
DELNC_DEPCG_1_eQTL$tissueName <- best_runs_toSelect$Tissue[15]
DELNC_DEPCG_1_eQTL$overlapType <- best_runs_toSelect$CClncRNA_eQTL_overlap[15]
DELNC_DEPCG_1_eQTL$pairType <- best_runs_toSelect$CClncRNA_type[15]
eQTLconfirmed_pairs[[15]] <- unique(DELNC_DEPCG_1_eQTL[,c(3,20:22)])

#add columns to table:
trial <- CoRegPairs_04_48_24_extended

for (i in 1:length(eQTLconfirmed_pairs)){
  trial$newCol <- FALSE
  trial$newCol[trial$pairs %in% eQTLconfirmed_pairs[[i]]$pairs] <- TRUE
  colnames(trial)[length(colnames(trial))] <- paste(eQTLconfirmed_pairs[[i]][1,2],
                                                    eQTLconfirmed_pairs[[i]][1,3],
                                                    eQTLconfirmed_pairs[[i]][1,4],sep = "_")
}

CoRegPairs_04_48_24_extended_eQTL <- trial

CoRegPairs_04_48_24_extended_eQTL$eQTLvalidations <- rowSums(CoRegPairs_04_48_24_extended_eQTL[,c(18:32)])

#write.csv(CoRegPairs_04_48_24_extended_eQTL, "CoRegPairs_04_48_24_extended_eQTL25.csv",row.names = F)


#### compare to previous version ####

CoRegPairs_04_48_24_extended_eQTL <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CoRegPairs_04_48_24_extended_eQTL_Oct25.csv")
eQTL_validated <- filter(CoRegPairs_04_48_24_extended_eQTL, eQTLvalidations>0)

#compare to 250kbp (i.e. did the TPM levels add much?):
CoRegPairs_04_48_24_extended250 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CoRegPairs_04_48_24_extendedIII.csv")
colnames(CoRegPairs_04_48_24_extended250)

CoRegPairs_04_48_24_extended250$eQTLvalidationsGTEX <- rowSums(CoRegPairs_04_48_24_extended250[,c(23,24,25,26:30)])

eQTL_validated250 <- filter(CoRegPairs_04_48_24_extended250, eQTLvalidationsGTEX>0)

#similar numbers, but surprisingly little overlap
sum(eQTL_validated$pairs %in% eQTL_validated250$pairs)
23/45
23/45

#TPM levels add different pairs...

#which does best on expected cis lncs:
ControlCisLncs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/ControlCisLncs.csv")[,-1]

#no MSTRG.IDs available for the CATG IDs:
ControlCisLncs$Ens_ID[grepl("CATG", ControlCisLncs$Ens_ID)] %in% gsub("\\.[0-9]*", "", Enhancer_lociII_DEsig_Enh$FANTOM_ID)

#how many lncs with expected cis function expressed in timecourse:
fpkm_allG$EnsID_merge <- gsub("\\.[0-9]*", "", fpkm_allG$EnsID)
#23x found expressed and still called lncs
ControlCisLncs_timecourse <- filter(ControlCisLncs, Ens_ID %in% filter(fpkm_allG, grepl("Lnc|fide", GeneClassUpdate))$EnsID_merge)

eQTL_validated250$EnsID_merge <- gsub("\\.[0-9]*", "", eQTL_validated250$EnsID)
eQTL_validated$EnsID_merge <- gsub("\\.[0-9]*", "", eQTL_validated$EnsID)

#same for both:
sum(unique(ControlCisLncs_timecourse$Ens_ID) %in% eQTL_validated$EnsID_merge)
sum(unique(ControlCisLncs_timecourse$Ens_ID) %in% eQTL_validated250$EnsID_merge)
ControlCisLncseQTLval <- filter(GTEX_validated, EnsID_merge %in% ControlCisLncs_timecourse$Ens_ID)

#quick check of HiC whilst here:
HiC_validated <- filter(CoRegPairs_04_48_24_extended, !loopMethod == "Neither")
HiC_validated250 <- filter(CoRegPairs_04_48_24_extended250, !loopMethod == "Neither")
#same pairs found, makese sense, expanding to 400kbp wouldn't change this
HiC_validated$pairs %in% HiC_validated250$pairs

HiC_validated$EnsID_merge <- gsub("\\.[0-9]*", "", HiC_validated$EnsID)

#3x found:
sum(ControlCisLncs_timecourse$Ens_ID %in% HiC_validated$EnsID_merge)
length(unique(filter(fpkm_allG, grepl("Lnc|fide", GeneClassUpdate))$EnsID))
length(unique(HiC_validated$EnsID))
23/428
3/32
#somewhat of an enrichment but low pool...


#### overlap identified eQTL-supported pairs with those from FANTOM analysis ####

CoRegPairs_04_48_24_extended_eQTL <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CoRegPairs_04_48_24_extended_eQTL_Oct25.csv")
eQTL_validated <- filter(CoRegPairs_04_48_24_extended_eQTL, eQTLvalidations>0)

#compare these to Hon 2017
#slight diff method - they looked for any eQTL-linked pairs which were significantly more co-expressed 
#than background shuffled pairs
#within similar distance, with same orientation
#correlation across all 1.8k samples - non-specific

#we compared to a set of eQTL linked lncRNA-mRNA from previous analysis
#which was a) more generic, not focused on confirming co-expression in one cell type
#b) looking only for linear correlation (not timing based synchronising like done here)

FANTOM_eQTL_pairs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/nature21374-s2/table16.csv")
FANTOM_eQTL_pairs_cis <- filter(FANTOM_eQTL_pairs, cis_correlated_candidate == "yes")
#5264 pairs as reported in Hon et al 2017

#pairs found as neighbours:
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

#includes SNHG15-TBRG4, no other lab candidates validated by FANTOM approach

#confirmed by eQTLs from any tissue at any lnc site:
triali <- filter(trial, pairs %in% eQTL_validated$pairs, 
                 FANTpairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff |
                   pairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff)
length(unique(triali$pairs))#11 found both expressed in SVSMC and in a co-reg pairing and supported by this study's eQTL approach

11/45 # 24% 11 of the 45 eQTL validated pairs are found in FANTOM approach too
18/364 # 5% of co-regs generally

#nice overlap

a <- 11
b <- 45
c <- 18
d <- 364

fisher.test(data.frame("cisLnc" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#p = 5e-7, OR = 14.2

#with obvious advantage that I have found (45-18 =) 27 more LOL

DEL_PCG_type <- rbind("eQTL-supported CClncRNA-\n-target pairs (this study)" = c(a,b, a/b),
                      "All CClncRNA-target pairs" = c(c,d, c/d))
colnames(DEL_PCG_type) <- c("eQTL-supported (FANTOM)", "Other", "%")

library(reshape2)
DEL_PCG_typei <- melt(DEL_PCG_type)[5:6,]

DEL_PCG_typei$value <- DEL_PCG_typei$value*100

DEL_PCG_typei$Var1 <- as.factor(DEL_PCG_typei$Var1)
DEL_PCG_typei$Var1 <- factor(DEL_PCG_typei$Var1, levels = levels(DEL_PCG_typei$Var1)[c(2,1)])

ggplot(DEL_PCG_typei) + aes(y = Var1, x = value) + 
  geom_bar(stat = "identity", color = "grey60") +
  xlab("% eQTL-supported\n(FANTOM, Hon et al. 2017)") +
  ylab("") +
  theme_minimal() +
  theme(text = element_text(size=24))# + Seurat::RotatedAxis()

DEL_PCG_type$Type <- as.factor(DEL_PCG_type$Type)
#DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,3,1)])
DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,1,3)])

colnames(DEL_PCG_type)[c(2,5)] <- c("LncRNA-PCG Pairings", "noPairs")

ggplot(DEL_PCG_type) + aes(y = Type, x = value, fill = `LncRNA-PCG Pairings`) + #, label = noPairs) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "grey60") +
  scale_fill_manual(values = c(`CCLncRNA targets\n(same timeframe)` = "mediumorchid", `All neighbours` = "grey60")) +
  #geom_label(position = position_dodge(width = 0.7), color = "white", size =5) +
  #ylab("\nEarly CClncRNAs\n(with target in same timeframe)") +
  xlab("% HiC Connected") +
  scale_x_continuous(breaks = seq(0,15,5), limits = c(0,15)) +
  theme_minimal() +
  #ggtitle("LncRNAs with a DE PCG neighbour", subtitle = "(same timeframe)") +
  theme(text = element_text(size=24))# + Seurat::RotatedAxis()


#SMC biobank specific:
#confirmed by Shu SMC biobank eQTLs at the lnc promoter:
triali <- filter(trial, pairs %in% ShuBioBankSupported$pairs, 
                 FANTpairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff |
                   pairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff)
length(unique(triali$pairs))#6 found both expressed in SVSMC and in a co-reg pairing and supported by the biobank

6/21 #28% of the shu are found
18/364 #5% of the co-regs generally

#nice overlap

a <- 6
b <- 21
c <- 18
d <- 364

fisher.test(data.frame("cisLnc" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")
#p = 0.0002, OR = 10.9


#### no. early induced/repressed ####

CoRegPairs_04_48_24_extended_eQTL <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CoRegPairs_04_48_24_extended_eQTL_Oct25.csv")

#any additional?
CoRegPairs_04_48_24_extended_eQTL$FANTOM_eQTL <- "No"
CoRegPairs_04_48_24_extended_eQTL$FANTOM_eQTL[CoRegPairs_04_48_24_extended$pairs %in% triali$pairs] <- "Yes"

eQTL_validated <- filter(CoRegPairs_04_48_24_extended_eQTL, eQTLvalidations>0 | FANTOM_eQTL == "Yes"
                         )

#37x lncRNAs
length(unique(eQTL_validated$EnsID))

eQTL_val_time <- as.data.frame(table(unique(eQTL_validated[,c(2,9)])$Lnc_Cluster))
eQTL_val_time$selection <- table(fpkm_allGDE$RegulationStart)

sum(eQTL_val_time$Freq)
11/37
973/4345

LncEnrich_cluster <- list()
for (i in 1:length(eQTL_val_time$Var1)){
  a <- eQTL_val_time[i,2]
  b <- eQTL_val_time[i,3]
  c <- sum(eQTL_val_time$Freq)
  d <- 4345
  
  LncEnrich_cluster[[i]] <- data.frame(fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")))$est,
                                       fisher.test(data.frame("LncRNA" = c(a,b-a),
                                                              "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")))$p)
}
names(LncEnrich_cluster) <- eQTL_val_time$Var1
triali <- bind_rows(LncEnrich_cluster, .id = "Cluster")
rownames(triali) <- NULL
colnames(triali) <- c("Cluster", "OR", "p")
triali$p_adj <- p.adjust(triali$p, method = "BH")

#ns

22/37
(973+746)/4345

a <- 22
b <- 973+746
c <- 37
d <- 4345

#* somewhat towards early regulation
fisher.test(data.frame("LncRNA" = c(a,b-a),
                       "other" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")))

#percentage plots
eQTL_val_time$FirstRegulation <- c("Within \n4hrs", "Within \n8hrs", "Within \n24hrs", "Within \n4hrs", "Within \n8hrs", "Within \n24hrs")
eQTL_val_time$FirstRegulation <- as.factor(eQTL_val_time$FirstRegulation)
eQTL_val_time$FirstRegulation <- factor(eQTL_val_time$FirstRegulation, levels = levels(eQTL_val_time$FirstRegulation)[c(2,3,1)])

eQTL_val_time$UpDown <- sapply(sapply(as.character(eQTL_val_time$Var1), strsplit, " "),"[[" , 1)

eQTL_val_time$PercCategory <- eQTL_val_time$Freq/sum(eQTL_val_time$Freq)*100
eQTL_val_time$PercBackground <- table(fpkm_allGDE$RegulationStart)/4345*100

#ELncWaveBias <- cbind(trial[,-c(2)], triali[,-c(1)])

ggplot(eQTL_val_time, aes(x = FirstRegulation)) +
  geom_col(data = filter(eQTL_val_time, grepl("Induced", UpDown)), 
           aes(y = PercCategory, fill = UpDown), color = "black") +
  geom_col(data = filter(eQTL_val_time, grepl("Induced", UpDown)), 
           aes(y = PercBackground), fill = NA, color = "grey30", linetype = "dashed", linewidth = 1) +
  geom_label(data = filter(eQTL_val_time, grepl("Induced", UpDown)), 
             aes(y = PercCategory, label = Freq), size = 4.2) +
  geom_col(data = filter(eQTL_val_time, grepl("Repressed", UpDown)), 
           aes(y = -PercCategory, fill = UpDown), color = "black") +
  geom_col(data = filter(eQTL_val_time, grepl("Repressed", UpDown)), 
           aes(y = -PercBackground), fill = NA, color = "grey30", linetype = "dashed", linewidth = 1) +
  geom_label(data = filter(eQTL_val_time, grepl("Repressed", UpDown)), 
             aes(y = -PercCategory+1, label = Freq), size = 4.2) +
  ylab("% eQTL-supported\nCClncRNAs") +
  xlab("") +
  scale_fill_manual(values = c("Induced" = "#D6604D", "Repressed" = "#67A9CF")) +
  #scale_y_continuous(limits = c(-30,60),breaks = seq(-20,60, by = 20),
  #                   labels = (c(seq(20, 0, by = -20), seq(20,60,by=20)))) +
  theme_minimal() +
  theme(text = element_text(size = 15), legend.position = "none")


#vs. all lncRNAs too?
trial <- data.frame(61/198*100, 22/37*100)
colnames(trial) <- c("All DE lncRNAs", "Enhancer-transcribed\nDE lncRNAs")

#melt(trial)

ggplot(melt(trial)) + aes(x = value, y = variable) +
  geom_bar(stat = "identity", fill = "olivedrab3", color = "grey60") +
  xlab("% 0-4hr induced") +
  ylab("") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0,50,25)) +
  theme(text = element_text(size =24))
