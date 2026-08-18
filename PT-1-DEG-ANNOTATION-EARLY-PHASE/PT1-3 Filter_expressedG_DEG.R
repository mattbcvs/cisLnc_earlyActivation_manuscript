#Filtering to expressed genes and transcripts
library(dplyr)
PLAR_Timecourse_4Timepoints_DEnonDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/PLAR_Timecourse_4Timepoints_DEnonDE_2026.csv", header = T)


#######remove complexity by taking out low exprs tx#############
trial <- data.frame(PLAR_Timecourse_4Timepoints_DEnonDE, 
                    "Ens_ID_Factor" = as.factor(PLAR_Timecourse_4Timepoints_DEnonDE$EnsID), stringsAsFactors = F)

#could also do via the iso_pct column of RSEM output, but actually more stringent to include more tx for lnc annotation (chance that low exprs may be PCG)
#probs useful for later tests to get rid of v. low isos across the board
trialii <- split(trial, trial$Ens_ID_Factor)
trialiii <- lapply(trialii, function(x){
  x[(x$Tx_Max_Average/sum(x$Tx_Max_Average) >= 0.05), 47]
}) #return name of tx which provide over 5% of fpkm in each gene
trialiv <- sapply(trialiii, length)
table(trialiv) #no genes returning >7 transcripts, no need to trim further

#can now subset by the returned transcripts, use max tx fpkm to detect whether tx is expressed, average probs not particularly accurate so may as well make inclusive
fpkm_lowCompTx <- unlist(unlist(trialiii))
length(unique(fpkm_lowCompTx)) #48962 (one is NA) tx in total (~33k are <5%)

fpkm_allTx <- filter(PLAR_Timecourse_4Timepoints_DEnonDE, MSTRG_Tx_ID %in% fpkm_lowCompTx) 
length(unique(fpkm_allTx$MSTRG_Tx_ID))#48962 Tx
length(unique(fpkm_allTx$EnsID))#18672 Genes


##### assign tx/gene classes ####
#simplify PLAR pre-filter classes (any coding ORF in databases put into one group)
fpkm_allTx$linc_pred_level[grep("Fragment", fpkm_allTx$PLAR_Prefilter)] <- "Protein coding"
fpkm_allTx$linc_pred_level[grep("Full", fpkm_allTx$PLAR_Prefilter)] <- "Protein coding"

table(fpkm_allTx$linc_pred_level)#40879(29626) discarded through known ORFs

fpkm_allTx$linc_pred_level[fpkm_allTx$linc_pred_level == "potential_linc"] <- "Putative lncRNA"
fpkm_allTx$linc_pred_level[fpkm_allTx$linc_pred_level == "clean_linc"] <- "Putative lncRNA"
fpkm_allTx$linc_pred_level[fpkm_allTx$linc_pred_level == "pure_linc"] <- "Bona fide lncRNA"

fpkm_allTx$linc_pred_level[is.na(fpkm_allTx$linc_pred_level)] <- "Repetitive/rare artefacts"
fpkm_allTx$linc_pred_level[fpkm_allTx$PLAR_Prefilter == "Single exon stringent"] <- "Single exon artefacts"
fpkm_allTx$linc_pred_level[grepl("Too low", fpkm_allTx$PLAR_Prefilter)] <- "Too low"

table(fpkm_allTx$linc_pred_level)
#2932(869) bona fide lncRNA transcripts expressed

#extrapolate to gene level:
fpkm_allG <- filter(fpkm_allTx, fpkm_max_treatment >0.8)
length(unique(fpkm_allG$MSTRG_Tx_ID))#44225 Tx
length(unique(fpkm_allG$EnsID))#14235 Genes
table(fpkm_allG$linc_pred_level)

#assign classes to each gene
fpkm_remnantsG <- fpkm_allG[is.na(fpkm_allG$linc_pred_level),]
fpkm_allG[fpkm_allG$EnsID %in% fpkm_remnantsG$EnsID, 55] <- "Remnants"
fpkm_repG <- fpkm_allG$EnsID[fpkm_allG$PLAR_Prefilter == "Repetitive"]
fpkm_allG[fpkm_allG$EnsID %in% fpkm_repG, 55] <- "Repetitive/rare artefacts"
fpkm_lowG <- fpkm_allG$EnsID[fpkm_allG$linc_pred_level == "Too low"]
fpkm_allG[fpkm_allG$EnsID %in% fpkm_lowG, 55] <- "Too low"
fpkm_1exG <- fpkm_allG$EnsID[fpkm_allG$PLAR_Prefilter == "Single exon stringent"]
fpkm_allG[fpkm_allG$EnsID %in% fpkm_1exG, 55] <- "Single exon artefacts"

#simple method of lncRNA classification
fpkm_bonaG <- (filter(fpkm_allG, (grepl("Bona", linc_pred_level))))$EnsID
fpkm_allG[fpkm_allG$EnsID %in% fpkm_bonaG, 55] <- "Bona fide lncRNA"
fpkm_putG <- (filter(fpkm_allG, (grepl("Putative", linc_pred_level))))$EnsID
fpkm_allG[fpkm_allG$EnsID %in% fpkm_putG, 55] <- "Putative lncRNA"
fpkm_PCG <- (filter(fpkm_allG, (grepl("Full", PLAR_Prefilter) | grepl("Fragment", PLAR_Prefilter))))$EnsID
fpkm_allG[fpkm_allG$EnsID %in% fpkm_PCG, 55] <- "Protein coding"

table(unique(fpkm_allG[,c(2,55)])$V55)
length(unique(fpkm_allG$EnsID))#14235 expressed genes:
#847 bona fide lncs
#11887 PCGs (increases relative to previous)

#sanity check - search for duplicates
ref_ids <- unique(fpkm_allG[,1:2])
ref_ids <- ref_ids[duplicated(ref_ids$EnsID),1:2]
trial <- fpkm_allG[fpkm_allG$EnsID %in% ref_ids$EnsID, 1:2]
#none found

#deal with rounding error duplicates
fpkm_allG[,c(25:33,51:52)] <- round(fpkm_allG[,c(25:33,51:52)], 2)
fpkm_allG[,34:46] <- as.numeric(sapply(fpkm_allG[,34:46], formatC, format = "e", digits = 2))
fpkm_allG <- unique(fpkm_allG)

table(unique(fpkm_allG[,c(2,55)])$V55)
length(unique(fpkm_allG$EnsID))#14235
length(unique(fpkm_allG$MSTRG_Tx_ID))#44225 txs

#write.csv(fpkm_allG, "fpkm_allG_2026.csv", row.names = F)

#
##### DEGs ####

#start matrix to find DEGs
fpkm_allGDE_lrt <- unique(filter(fpkm_allG, padj < 0.05)[,-c(1,47:54)])
fpkm_allGDE_pairs <-  unique(filter(fpkm_allG, (preadj_0_4 <0.05 & (LogFC_0_4 > log2(1.5) | LogFC_0_4 < -log2(1.5))) |
                                      (preadj_0_8 <0.05 & (LogFC_0_8 > log2(1.5) | LogFC_0_8 < -log2(1.5))) |
                                      (preadj_0_24 <0.05 & (LogFC_0_24 > log2(1.5) | LogFC_0_24 < -log2(1.5))) |
                                      (preadj_4_8 <0.05 & (LogFC_4_8 > log2(1.5) | LogFC_4_8 < -log2(1.5))) |
                                      (preadj_4_24 <0.05 & (LogFC_4_24 > log2(1.5) | LogFC_4_24 < -log2(1.5))) |
                                      (preadj_8_24 <0.05 & (LogFC_8_24 > log2(1.5) | LogFC_8_24 < -log2(1.5))))[,-c(1,47:54)])

length(unique(filter(fpkm_allGDE_lrt, EnsID %in% fpkm_allGDE_pairs$EnsID)$EnsID)) #5177 overlap
length(unique(fpkm_allGDE_lrt$EnsID)) #5412 LRT DEGs 4636/4865 
length(unique(fpkm_allGDE_pairs$EnsID)) #5880 pairwise DEGs 4636/5226

#Overlap probs best, kicks out low fold changing genes (LRT only) or inconsistent changing genes (Pairwise only)
fpkm_allGDE <- filter(fpkm_allGDE_lrt, EnsID %in% fpkm_allGDE_pairs$EnsID)

#no duplicates
dupEnsID <- fpkm_allGDE$EnsID[duplicated(fpkm_allGDE[,1])]

table(fpkm_allGDE$V55)#285 lncs

#now using this option to filter out artefacts:
fpkm_allGDE <- filter(fpkm_allGDE, !grepl("artefact", V55))
length(unique(fpkm_allGDE$EnsID))
table(unique(fpkm_allGDE[,c(1,46)])$V55)
#4508 DEGs, 262 lncs, 4579 PCGs 

#gene checking
filter(fpkm_allG, EnsID == "MSTRG.24277")

filter(fpkm_allGDE, EnsID == "MSTRG.12915")
filter(fpkm_allGDE, EnsName == "AC002480.4")
