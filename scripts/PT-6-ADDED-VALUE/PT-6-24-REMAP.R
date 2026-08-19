library(dplyr)
library(GenomicRanges)
library(ggplot2)

#
#### import table ####

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

#table of all lncRNAs + CAGE sites if available + TSS based on CAGE if available:
Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime_2026.csv", header = T)
Enhancer_lociII_DEsig_Enh <- Enhancer_lociII
length(unique(Enhancer_lociII_DEsig_Enh$EnsID))#s597
length(unique(Enhancer_lociII_DEsig_Enh$MSTRG_Tx_ID))#1575
length(unique(filter(Enhancer_lociII_DEsig_Enh, !is.na(DiffExprs))$EnsID))#221

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

#start with just promoters...
allGB_GR <- makeGRangesFromDataFrame(allGB[,-7], 
                                     start.field = "Tx_start", 
                                     end.field = "Tx_stop", 
                                     seqnames.field = "chr", 
                                     strand.field = "str", keep.extra.columns = T)

#consider shortening later 2000bp is a lot
allGB_GR_promoters <- promoters(allGB_GR)

#flatten to aid next steps
allGB_GR_promoters_reduce <- reduce(allGB_GR_promoters)
length(allGB_GR_promoters_reduce)#18838 regions

#bed format:
allGB_GR_promoters_reduce_bed <- as.data.frame(allGB_GR_promoters_reduce)

#write.table(allGB_GR_promoters_reduce_bed[,c(1:3)], "allGB_GR_promoters_reduce_bed.bed", row.names = F, col.names = F, quote = F)


#
#### remap overlay ####
remapCatalog2022hg38 <- read.delim("heavyFiles4R/remap2022_nr_macs2_hg38_v1_0.bed", nrows = 5000000,
                                   header = F)
dim(remapCatalog2022hg38) #5M as expected

trial <- makeGRangesFromDataFrame(remapCatalog2022hg38, 
                                  start.field = "V2", 
                                  end.field = "V3", seqnames.field = "V1", keep.extra.columns = T)
rm(remapCatalog2022hg38)

#remove non standards
seqlevels(trial)
table(as.character(seqnames(trial)))
trial <- keepStandardChromosomes(trial, pruning.mode = "coarse")
seqlevels(trial)
table(as.character(seqnames(trial)))

#query all genes granges with remap:
remapindex <- findOverlaps(query = trial, subject = allGB_GR_promoters_reduce)

#batch1, 825,243 NR promoter peaks:
batch1_SVSMC_match <- unique(trial[queryHits(remapindex)])
write.table(batch1_SVSMC_match, "batch1_SVSMC_matchi_extraPCG.bed", row.names = F, col.names = F, quote = F)
dim(batch1_SVSMC_match)


#next 5mil-20mil
remapCatalog2022hg38_2 <- read.delim("heavyFiles4R/remap2022_nr_macs2_hg38_v1_0.bed", skip = 5000000, nrows = 15000000,
                                     header = F)
dim(remapCatalog2022hg38_2)

trial <- makeGRangesFromDataFrame(remapCatalog2022hg38_2, 
                                  start.field = "V2", 
                                  end.field = "V3", seqnames.field = "V1", keep.extra.columns = T)
rm(remapCatalog2022hg38_2)

#remove non standards
seqlevels(trial)
table(as.character(seqnames(trial)))
trial <- keepStandardChromosomes(trial, pruning.mode = "coarse")
seqlevels(trial)
table(as.character(seqnames(trial)))

#query all genes granges with remap:
remapindex <- findOverlaps(query = trial, subject = allGB_GR_promoters_reduce)

#batch2, 2,163,643
batch2_SVSMC_match <- unique(trial[queryHits(remapindex)])
write.table(batch2_SVSMC_match, "batch2_SVSMC_matchi_extraPCG.bed", row.names = F, col.names = F, quote = F)


#next 20mil-55mil
remapCatalog2022hg38_3 <- read.delim("heavyFiles4R/remap2022_nr_macs2_hg38_v1_0.bed", 
                                     skip  = 20000000, 
                                     nrows = 35000000,
                                     header = F)
dim(remapCatalog2022hg38_3)

trial <- makeGRangesFromDataFrame(remapCatalog2022hg38_3, 
                                  start.field = "V2", 
                                  end.field = "V3", seqnames.field = "V1", keep.extra.columns = T)
rm(remapCatalog2022hg38_3)

#remove non standards
seqlevels(trial)
table(as.character(seqnames(trial)))
trial <- keepStandardChromosomes(trial, pruning.mode = "coarse")
seqlevels(trial)
table(as.character(seqnames(trial)))

#query all genes granges with remap:
remapindex <- findOverlaps(query = trial, subject = allGB_GR_promoters_reduce)

#batch3 5,312,220
batch3_SVSMC_match <- unique(trial[queryHits(remapindex)])
write.table(batch3_SVSMC_match, "batch3_SVSMC_matchi_extraPCG.bed", row.names = F, col.names = F, quote = F)


remapCatalog2022hg38_4 <- read.delim("heavyFiles4R/remap2022_nr_macs2_hg38_v1_0.bed", 
                                     skip  = 55000000,
                                     header = F)
dim(remapCatalog2022hg38_4)

trial <- makeGRangesFromDataFrame(remapCatalog2022hg38_4, 
                                  start.field = "V2", 
                                  end.field = "V3", seqnames.field = "V1", keep.extra.columns = T)
rm(remapCatalog2022hg38_4)

#remove non standards
seqlevels(trial)
table(as.character(seqnames(trial)))
trial <- keepStandardChromosomes(trial, pruning.mode = "coarse")
seqlevels(trial)
table(as.character(seqnames(trial)))

#query all genes granges with remap:
remapindex <- findOverlaps(query = trial, subject = allGB_GR_promoters_reduce)

#batch4 1,733,473
batch4_SVSMC_match <- unique(trial[queryHits(remapindex)])
write.table(batch4_SVSMC_match, "batch4_SVSMC_matchi_extraPCG.bed", row.names = F, col.names = F, quote = F)


#
#### import NR ReMap peaks overlapping SVSMC promoter regions (10M peaks overlapping 2kbp regions) ####

batch1_SVSMC_match <- read.table("batch1_SVSMC_matchi_extraPCG.bed")
batch2_SVSMC_match <- read.table("batch2_SVSMC_matchi_extraPCG.bed")
batch3_SVSMC_match <- read.table("batch3_SVSMC_matchi_extraPCG.bed")
batch4_SVSMC_match <- read.table("batch4_SVSMC_matchi_extraPCG.bed")

#colnames(batch1_SVSMC_match) <- c("seqnames", "start", "end", "width", "V5", "V6", "V7", "V8", "peak_start", "peak_end", "V11")
#colnames(batch2_SVSMC_match) <- c("seqnames", "start", "end", "width", "V5", "V6", "V7", "V8", "peak_start", "peak_end", "V11")
#colnames(batch3_SVSMC_match) <- c("seqnames", "start", "end", "width", "V5", "V6", "V7", "V8", "peak_start", "peak_end", "V11")
#colnames(batch4_SVSMC_match) <- c("seqnames", "start", "end", "width", "V5", "V6", "V7", "V8", "peak_start", "peak_end", "V11")

remap_SVSMC_prom <- unique(rbind(batch1_SVSMC_match, batch2_SVSMC_match, batch3_SVSMC_match, batch4_SVSMC_match))

remap_SVSMC_prom_GR <- makeGRangesFromDataFrame(remap_SVSMC_prom, start.field = "V2", end.field = "V3", seqnames = "V1", keep.extra.columns = T)

rm(remap_SVSMC_prom)
rm(batch1_SVSMC_match)
rm(batch2_SVSMC_match)
rm(batch3_SVSMC_match)
rm(batch4_SVSMC_match)
gc()

#saveRDS(remap_SVSMC_prom_GR, "remap_SVSMC_prom_GR_extraPCG_2026.rds")

remap_SVSMC_prom_GR <- readRDS("remap_SVSMC_prom_GR_extraPCG_2026.rds")

#
#### 200bp promoters, find all TFs bound to each gene ####

#per gene yields a handful of enriched TFs

#switching to smaller promoters, to a) fit with AME b) see if increases TP signal
allGB_GR_promoters200 <- promoters(allGB_GR, upstream = 200)

#flatten to aid next steps
allGB_GR_promoters_reduce200 <- reduce(allGB_GR_promoters200)
length(allGB_GR_promoters_reduce200)#21805 regions (about 3k more than 2000)

#repeat assignment of remap peaks
remap_SVSMC_prom_GR$peakID <- paste(as.character(seqnames(remap_SVSMC_prom_GR)),
                                    start(remap_SVSMC_prom_GR), 
                                    end(remap_SVSMC_prom_GR), 
                                    remap_SVSMC_prom_GR$V4, sep = "_")

remap_SVSMC_prom_GR$peakTF <- sapply(strsplit(remap_SVSMC_prom_GR$V6, ":"), "[[", 1)


#prom ID for promoter co-ordinates:
allGB_GR_promoters_reduce200$promID <- paste(as.character(seqnames(allGB_GR_promoters_reduce200)),
                                          start(allGB_GR_promoters_reduce200), end(allGB_GR_promoters_reduce200), 
                                          #allGB_GR_promoters_reduce200$MSTRG_Tx_ID, 
                                          sep = "_")

#add in gene ID per reduced prom:
remapindex <- findOverlaps(query = allGB_GR_promoters, subject = allGB_GR_promoters_reduce200)

#return matched peaks from ReMap2022:
allPromRed_EnsID_200 <- unique(data.frame("EnsID" = allGB_GR_promoters$EnsID[queryHits(remapindex)],
                                      "promID" = allGB_GR_promoters_reduce200$promID[subjectHits(remapindex)]))

#now make match of reduced proms to remap peaks
remapindex <- findOverlaps(query = remap_SVSMC_prom_GR, subject = allGB_GR_promoters_reduce200)

#return matched peaks from ReMap2022:
allProm_remap2022_200 <- data.frame("peakID" = remap_SVSMC_prom_GR$peakID[queryHits(remapindex)],
                                "peakTF" = remap_SVSMC_prom_GR$peakTF[queryHits(remapindex)],
                                "peakTF_cell" = remap_SVSMC_prom_GR$V6[queryHits(remapindex)],
                                "promID" = allGB_GR_promoters_reduce200$promID[subjectHits(remapindex)])

#now can combine all above to get for each gene ID, add any bound TFs in ReMap 
trial <- unique(merge(unique(allProm_remap2022_200[,c(2:4)]), allPromRed_EnsID_200, by = "promID", all.y = T))

length(unique(trial$promID))#21805 promoters
length(unique(allProm_remap2022_200$promID))#20981
length(unique(allGB_GR_promoters_reduce200$promID))#21805 (probs inc some without a remap overlap)

length(unique(trial$EnsID))#12740 genes -i.e. even if a gene is not bound by a TF, should be in here
length(unique(allGB_GR_promoters$EnsID))#12740 genes
length(unique(fpkm_allG$EnsID))#12740 genes
length(unique(allPromRed_EnsID_200$EnsID))#12740 genes

allG_remap2022_TFs_200 <- trial
rm(trial)
gc()

#write.csv(allG_remap2022_TFs_200, "allG_remap2022_200_TFs.csv", row.names = F)

allG_remap2022_TFs_200 <- read.csv("allG_remap2022_200_TFs.csv")
#\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/

#
#### 0-4hr run - 200bp ####

#use this to find ReMap regulators which are expressed / with high potential to be active in each window:
fpkm_thresh <- 5
ReMaps_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >fpkm_thresh | Hour4_meanFPKM >fpkm_thresh))

#potential background, promoters expressed per window:
expr_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)$EnsID)

#test any TFs expressed in window of interest, with a peak in target gene set (here up4):
allG_remap2022_TFs_window <- filter(allG_remap2022_TFs_200,
                                    EnsID %in% fpkm_allGDE_Upwithin_4$EnsID,
                                    peakTF %in% ReMaps_04$EnsName)

TFs_inTargetProms <- unique(allG_remap2022_TFs_window$peakTF)
TFs_inTargetProms <- TFs_inTargetProms[!is.na(TFs_inTargetProms)]

#promoter regions for lncRNAs in given cluster of interest 
TargetGenes <- allG_remap2022_TFs_window[allG_remap2022_TFs_window$EnsID %in% fpkm_allGDE_Upwithin_4$EnsID,]

TFs_in_4hrDEGProms_all_Fish <- list()

#restrict (less sig TFs, more confidence) or expand background (more sig TFs, requires harsh filtering):
BackgroundGenes <- filter(allG_remap2022_TFs_200, 
                          #hash in or out
                          EnsID %in% expr_04
                          )

for(i in 1:length(TFs_inTargetProms)){
  
  a <- length(unique(filter(TargetGenes, peakTF == TFs_inTargetProms[i])$EnsID))
  b <- length(unique(TargetGenes$EnsID))
  c <- length(unique(filter(BackgroundGenes, peakTF == TFs_inTargetProms[i])$EnsID))
  d <- length(unique(BackgroundGenes$EnsID))
  
  TFs_in_4hrDEGProms_all_Fish[[i]] <- c(
    a,b,c,d,
    fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                           "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p,
    fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                           "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est
  )
}

names(TFs_in_4hrDEGProms_all_Fish) <- TFs_inTargetProms

TFs_in_4hrDEGProms_all_Fish_df <- as.data.frame(t(bind_rows(TFs_in_4hrDEGProms_all_Fish, .id = "TFs")))

TFs_in_4hrDEGProms_all_Fish_df$peakTF <- TFs_inTargetProms

TFs_in_4hrDEGProms_all_Fish_df$pval_BH <- p.adjust(TFs_in_4hrDEGProms_all_Fish_df$V5, method = "BH")

#store some variants here
Prom0.2k_TF5fpkm_Prom0.8 <- TFs_in_4hrDEGProms_all_Fish_df

ggplot(Prom0.2k_TF5fpkm_Prom0.8) + aes(x = -log10(pval_BH), y = V6) +
  geom_point() +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")

ggplot(Prom0.2k_TF5fpkm_PromAll) + aes(x = -log10(pval_BH), y = V6) +
  geom_point() +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")



#### 0-4hr run - 200bp - split experiments ####

#use this to find ReMap regulators which are expressed / with high potential to be active in each window:
#better to be stricter than previous (?)
fpkm_thresh_TFs <- 5
ReMaps_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >fpkm_thresh_TFs | Hour4_meanFPKM >fpkm_thresh_TFs))

#set background genes, promoters expressed per window:
expr_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)$EnsID)

#select test set, PCGs vs lncs vs cclncs vs scclncs...
fpkm_PCGDE_Upwithin_4 <- filter(fpkm_allGDE_Upwithin_4, EnsType == "protein_coding", grepl("coding", GeneClassUpdate))
SCClncRNAs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SCClncRNAs.csv")
SCClncRNAs_04 <- filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs")
#from PT2-8:
CClncRNAs_04 <- filter(CoRegPairs_04_48_24_extended_naiveSame, EnsID %in% fpkm_allGDE_Upwithin_4$EnsID)
#probs no point for these ones
#non scclnc early induced:
fpkm_nonSCClncDE_Upwithin_4 <- filter(fpkm_allGDE_Upwithin_4, grepl("fide|Lnc", fpkm_allGDE_Upwithin_4$GeneClassUpdate),
                                !EnsID %in% SCClncRNAs_04$EnsID)

testSet <- fpkm_nonSCClncDE_Upwithin_4

#test any TFs expressed in window of interest, with a peak in target gene set (here up4):
allG_remap2022_TFs_window <- filter(allG_remap2022_TFs_200,
                                    #test genes
                                    EnsID %in% testSet$EnsID,
                                    #with a TF of interest peak in promoter 
                                    peakTF %in% ReMaps_04$EnsName)

#expand so can split and test TFs per cell type too
allG_remap2022_TFs_window <- allG_remap2022_TFs_window %>%
  tidyr::separate_rows(peakTF_cell, sep = ",")

allG_remap2022_TFs_window$cell <- sub("^[^:]*:\\s*", "", allG_remap2022_TFs_window$peakTF_cell)
allG_remap2022_TFs_window$peakTF_cell <- paste(allG_remap2022_TFs_window$peakTF,
                                               allG_remap2022_TFs_window$cell, sep = "_")
  
TFs_inTargetProms <- unique(allG_remap2022_TFs_window$peakTF_cell)
TFs_inTargetProms <- TFs_inTargetProms[!is.na(TFs_inTargetProms)]
#note: looks like ~0.5-1k more than ENCODE4

#promoter regions and bound TFs in a given experiment for genes of interest 
TargetGenes <- allG_remap2022_TFs_window[allG_remap2022_TFs_window$EnsID %in% testSet$EnsID,]
#cut unecessary cols
TargetGenes <- unique(TargetGenes[,c(3,4)])

#restrict (less sig TFs, more confidence) or expand background (more sig TFs, requires harsh filtering):
BackgroundGenes <- filter(allG_remap2022_TFs_200,
                          EnsID %in% expr_04
                          )

#also needs to be expanded:
BackgroundGenes <- BackgroundGenes %>%
  tidyr::separate_rows(peakTF_cell, sep = "\\s*,\\s*")
BackgroundGenes$cell <- sub("^[^:]*:\\s*", "", BackgroundGenes$peakTF_cell)
BackgroundGenes$peakTF_cell <- paste(BackgroundGenes$peakTF,
                                     BackgroundGenes$cell, sep = "_")

#cut unecessary cols
BackgroundGenes <- unique(BackgroundGenes[,c(3,4)])
gc()

#store results in here
TFs_in_4hrDEGProms_all_Fish <- list()

for(i in 1:length(TFs_inTargetProms)){
  
  a <- length(unique(filter(TargetGenes, peakTF_cell == TFs_inTargetProms[i])$EnsID))
  b <- length(unique(TargetGenes$EnsID))
  c <- length(unique(filter(BackgroundGenes, peakTF_cell == TFs_inTargetProms[i])$EnsID))
  d <- length(unique(BackgroundGenes$EnsID))
  
  TFs_in_4hrDEGProms_all_Fish[[i]] <- c(
    a,b,c,d,
    fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                           "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p,
    fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                           "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est
  )
}

names(TFs_in_4hrDEGProms_all_Fish) <- TFs_inTargetProms#[1:400]

TFs_in_4hrDEGProms_all_Fish_df <- as.data.frame(t(bind_rows(TFs_in_4hrDEGProms_all_Fish, .id = "TFs")))

TFs_in_4hrDEGProms_all_Fish_df$peakTF <- TFs_inTargetProms#[1:400]

TFs_in_4hrDEGProms_all_Fish_dfi <- filter(TFs_in_4hrDEGProms_all_Fish_df, V3 >20)

TFs_in_4hrDEGProms_all_Fish_dfi$pval_BH <- p.adjust(TFs_in_4hrDEGProms_all_Fish_dfi$V5, method = "BH")

#store PCG version (keep all <20 for now in saved object)
#Prom0.2k_remapSplit_TF5fpkm_PCGProm0.8 <- TFs_in_4hrDEGProms_all_Fish_df
#write.csv(Prom0.2k_remapSplit_TF5fpkm_PCGProm0.8, "ReMap_enrichedTFs_Up4hr_PCG.csv", row.names = F)

#plotting
Prom0.2k_remapSplit_TF5fpkm_PCGProm0.8$TF <- sapply(strsplit(Prom0.2k_remapSplit_TF5fpkm_PCGProm0.8$peakTF, "_"), "[[", 1)
Prom0.2k_remapSplit_TF5fpkm_PCGProm0.8$TF[Prom0.2k_remapSplit_TF5fpkm_PCGProm0.8$pval_BH >0.001 | 
                                            Prom0.2k_remapSplit_TF5fpkm_PCGProm0.8$V6 < 1.5] <- "Other"
Prom0.2k_remapSplit_TF5fpkm_PCGProm0.8$HitNo <- Prom0.2k_remapSplit_TF5fpkm_PCGProm0.8$V1
#Prom0.2k_remapSplit_TF5fpkm_PCGProm0.8$HitNo[Prom0.2k_remapSplit_TF5fpkm_PCGProm0.8$pval_BH >0.001 | 
#                                            Prom0.2k_remapSplit_TF5fpkm_PCGProm0.8$V6 < 1.5] <- 1

ggplot(filter(Prom0.2k_remapSplit_TF5fpkm_PCGProm0.8, pval_BH <0.05, V6 > 1.4)) + aes(x = -log10(pval_BH), y = V6, 
                                                                                      color = TF, 
                                                                                      #size = HitNo
                                                                                      ) +
  geom_point(alpha = 0.65, shape = 19, size = 5) +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
  scale_size(range = c(1,7)) +
  theme_minimal() +
  theme(text = element_text(size=15)) +
  xlab("Significance (-log10 adjusted p)") +
  ylab("Strength (odds ratio)")

#Prom0.2k_remapSplit_TF5fpkm_SCClncProm0.8 <- TFs_in_4hrDEGProms_all_Fish_df


#store scclncRNA version here:
Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8 <- TFs_in_4hrDEGProms_all_Fish_df
#write.csv(Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8, "ReMap_enrichedTFs_Up4hr_SCCL.csv")

#plotting
Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8$TF <- sapply(strsplit(Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8$peakTF, "_"), "[[", 1)
#Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8$TF[Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8$V5 >0.05 | 
#                                            Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8$V6 < 1.5] <- "Other"
Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8$HitNo <- Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8$V1
#Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8$HitNo[Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8$pval_BH >0.001 | 
#                                            Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8$V6 < 1.5] <- 1

ggplot(filter(Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8, V5 <0.025, V6 > 1.5)) + aes(x = -log10(V5), y = V6, 
                                                                                      color = TF, 
                                                                                      size = HitNo
                                                                                      ) +
  geom_point(alpha = 0.65, shape = 19, size = 5
             ) +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
  scale_size(range = c(3,7)) +
  theme_minimal() +
  theme(text = element_text(size=15)) +
  xlab("Significance (-log10 p)") +
  ylab("Strength (odds ratio)")


#store non scclncRNA version here
Prom0.2k_remapSplit_TF5fpkm_otherLProm0.8 <- TFs_in_4hrDEGProms_all_Fish_df
#write.csv(Prom0.2k_remapSplit_TF5fpkm_otherLProm0.8, "ReMap_enrichedTFs_Up4hr_otherL.csv", row.names = F)

#plotting
Prom0.2k_remapSplit_TF5fpkm_otherLProm0.8$TF <- sapply(strsplit(Prom0.2k_remapSplit_TF5fpkm_otherLProm0.8$peakTF, "_"), "[[", 1)
#Prom0.2k_remapSplit_TF5fpkm_otherLProm0.8$TF[Prom0.2k_remapSplit_TF5fpkm_otherLProm0.8$V3 >0.001 | 
#                                               Prom0.2k_remapSplit_TF5fpkm_otherLProm0.8$V6 < 1.5] <- "Other"
Prom0.2k_remapSplit_TF5fpkm_otherLProm0.8$HitNo <- Prom0.2k_remapSplit_TF5fpkm_otherLProm0.8$V1
#Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8$HitNo[Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8$pval_BH >0.001 | 
#                                            Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8$V6 < 1.5] <- 1

ggplot(filter(Prom0.2k_remapSplit_TF5fpkm_otherLProm0.8, V5 <0.005, V6 > 1.5)) + aes(x = -log10(V5), y = V6, 
                                                                                   color = TF, 
                                                                                   size = HitNo
                                                                                   ) +
  geom_point(alpha = 0.65, shape = 19, size = 5) +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
  scale_size(range = c(3,7)) +
  theme_minimal() +
  theme(text = element_text(size=15)) +
  xlab("Significance (-log10 p)") +
  ylab("Strength (odds ratio)")


#
#### prepare plotting ####

Prom0.2k_remapSplit_TF5fpkm_PCGProm0.8 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/ReMap_enrichedTFs_Up4hr_PCG.csv")
Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/ReMap_enrichedTFs_Up4hr_SCCL.csv")
Prom0.2k_remapSplit_TF5fpkm_otherLProm0.8 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/ReMap_enrichedTFs_Up4hr_OtherL.csv")

#number of datasets tested:
length(unique(c(Prom0.2k_remapSplit_TF5fpkm_PCGProm0.8$peakTF, 
         Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8$peakTF, 
         Prom0.2k_remapSplit_TF5fpkm_otherLProm0.8$peakTF)))

#correct p vals:

#plot across all 3:
#enrichments in all
ReMap_Up4_PCG <- filter(Prom0.2k_remapSplit_TF5fpkm_PCGProm0.8)#, V5 < 0.05, V6 > 1.2)
ReMap_Up4_SCCL <- filter(Prom0.2k_remapSplit_TF5fpkm_SCCLProm0.8)#, V5 < 0.05, V6 > 1.2)
ReMap_Up4_OtherL <- filter(Prom0.2k_remapSplit_TF5fpkm_otherLProm0.8)#, V5 < 0.05, V6 > 1.2)

#hit rates
ReMap_Up4_PCG$HitRate <- ReMap_Up4_PCG$V1/ReMap_Up4_PCG$V2
ReMap_Up4_SCCL$HitRate <- ReMap_Up4_SCCL$V1/ReMap_Up4_SCCL$V2
ReMap_Up4_OtherL$HitRate <- ReMap_Up4_OtherL$V1/ReMap_Up4_OtherL$V2
#back rates:
ReMap_Up4_PCG$BGRate <- ReMap_Up4_PCG$V3/ReMap_Up4_PCG$V4
ReMap_Up4_SCCL$BGRate <- ReMap_Up4_SCCL$V3/ReMap_Up4_SCCL$V4
ReMap_Up4_OtherL$BGRate <- ReMap_Up4_OtherL$V3/ReMap_Up4_OtherL$V4
#fold increase from background
ReMap_Up4_PCG$FoldEnr <- ReMap_Up4_PCG$HitRate/ReMap_Up4_PCG$BGRate
ReMap_Up4_SCCL$FoldEnr <- ReMap_Up4_SCCL$HitRate/ReMap_Up4_SCCL$BGRate
ReMap_Up4_OtherL$FoldEnr <- ReMap_Up4_OtherL$HitRate/ReMap_Up4_OtherL$BGRate

#TF column
ReMap_Up4_PCG$TF <- sapply(strsplit(ReMap_Up4_PCG$peakTF, "_"), "[[", 1)
ReMap_Up4_SCCL$TF <- sapply(strsplit(ReMap_Up4_SCCL$peakTF, "_"), "[[", 1)
ReMap_Up4_OtherL$TF <- sapply(strsplit(ReMap_Up4_OtherL$peakTF, "_"), "[[", 1)

#option to remove small numbers in this could be here

#option to plot diffs in hit rates could be here

#
#### compare SCClnc and PCG ####

#rbind plots to compare PCG and SCClnc:
ReMap_Up4_PCG$GeneType <- "PCG"
ReMap_Up4_SCCL$GeneType <- "SCCLnc"

#enriched experiments in one or other:
enrichedTFs <- unique(c(filter(ReMap_Up4_PCG,  V5 < 0.01, V6 > 1.5, V1>=4)$peakTF,
                 filter(ReMap_Up4_SCCL,  V5 < 0.01, V6 > 1.5, V1>=4)$peakTF))

#create a TF ordering by strongest in SCCL vs strongest in DE PCG:
trial <- filter(ReMap_Up4_PCG,  peakTF %in% enrichedTFs)
triali <- filter(ReMap_Up4_SCCL,  peakTF %in% enrichedTFs)

ReMap_Up4_PCG_v_SCCL <- merge(trial[,c(7,11, 5, 6,10)], triali[,c(7,11, 5, 6,10)], by = "peakTF", all = T)

ReMap_Up4_PCG_v_SCCL$peakTF <- as.factor(ReMap_Up4_PCG_v_SCCL$peakTF)
orderTFs <- ReMap_Up4_PCG_v_SCCL$peakTF[order(-ReMap_Up4_PCG_v_SCCL$V6.y, ReMap_Up4_PCG_v_SCCL$V6.x)]

colnames(ReMap_Up4_PCG_v_SCCL)[c(4,8)] <- c("PCGs", "SCClncs")
colnames(ReMap_Up4_PCG_v_SCCL)[c(3,7)] <- c("PCGs.p", "SCClncs.p")

trialiii <- reshape2::melt(ReMap_Up4_PCG_v_SCCL[,c(1,4,8)])

trialiii$peakTF <- factor(trialiii$peakTF)
trialiii$peakTF <- factor(trialiii$peakTF,
                                    levels = levels(trialiii$peakTF)[
                                      match(orderTFs, levels(trialiii$peakTF))
                                    ])

trialiii$p_val <- reshape2::melt(ReMap_Up4_PCG_v_SCCL[,c(1,3,7)])[,3]
trialiii$p_val_sig <- "NotSig"
trialiii$p_val_sig[trialiii$p_val < 0.01] <- "Sig"

trialiii$variableii <- paste(trialiii$variable, trialiii$p_val_sig, sep = "_")

ggplot(filter(trialiii)) + aes(x = peakTF, y = value, fill = variableii) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "grey40") +
  ylab("Odds ratio hits in selected\n vs. background promoters") +
  scale_fill_manual(values = c("PCGs_Sig" = "mediumorchid1", "PCGs_NotSig" = "plum1", 
                               "SCClncs_Sig" = "olivedrab3", "SCClncs_NotSig" = "lemonchiffon")) +
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  Seurat::RotatedAxis()

#simplfying, display per TF:
trialiii$TF <- sapply(strsplit(as.character(trialiii$peakTF), "_"), "[[", 1)

ggplot(filter(trialiii, TF == "RELA")) + aes(x = TF, y = value, fill = variableii) +
  geom_boxplot(position = position_dodge(width = 0.9), color = "grey40") +
  ylab("Odds ratio hits in selected\n vs. background promoters") +
  scale_fill_manual(values = c("PCGs_Sig" = "mediumorchid1", "PCGs_NotSig" = "plum1", 
                               "SCClncs_Sig" = "olivedrab3", "SCClncs_NotSig" = "lemonchiffon")) +
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  Seurat::RotatedAxis()

#or take the max per TF
trial$TF <- sapply(strsplit(trial$peakTF, "_"), "[[", 1)
triali$TF <- sapply(strsplit(triali$peakTF, "_"), "[[", 1)

#take best result out of all:
trialsplit <- split(trial, trial$TF)
trialspliti <- lapply(trialsplit, function(x){
  x[which(x$V6 == max(x$V6)),]
})
trial_perTF <- bind_rows(trialspliti)

trialsplit <- split(triali, triali$TF)
trialspliti <- lapply(trialsplit, function(x){
  x[which(x$V6 == max(x$V6)),]
})
triali_perTF <- bind_rows(trialspliti)

ReMap_Up4_PCG_v_SCCL_max <- merge(trial_perTF[,c(8,12, 5,6, 10,11)], triali_perTF[,c(8,12, 5,6, 10,11)], by = "TF", all = T)

ReMap_Up4_PCG_v_SCCL$peakTF <- as.factor(ReMap_Up4_PCG_v_SCCL$TF)
orderTFs <- ReMap_Up4_PCG_v_SCCL$peakTF[order(-ReMap_Up4_PCG_v_SCCL$V6.y, ReMap_Up4_PCG_v_SCCL$V6.x)]

colnames(ReMap_Up4_PCG_v_SCCL)[c(5,10)] <- c("PCGs", "SCClncs")
colnames(ReMap_Up4_PCG_v_SCCL)[c(4,9)] <- c("PCGs.p", "SCClncs.p")

trialiii <- reshape2::melt(ReMap_Up4_PCG_v_SCCL[,c(1,5,10)])

trialiii$peakTF <- factor(trialiii$TF)
trialiii$peakTF <- factor(trialiii$peakTF,
                          levels = levels(trialiii$peakTF)[
                            match(orderTFs, levels(trialiii$peakTF))
                          ])

trialiii$p_val <- reshape2::melt(ReMap_Up4_PCG_v_SCCL[,c(1,4,9)])[,3]
trialiii$p_val_sig <- "NotSig"
trialiii$p_val_sig[trialiii$p_val < 0.01] <- "Sig"

trialiii$variableii <- paste(trialiii$variable, trialiii$p_val_sig, sep = "_")

ggplot(filter(trialiii)) + aes(x = peakTF, y = value, fill = variableii) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "grey40") +
  ylab("Odds ratio - hits in selected\n vs. background promoters") +
  scale_fill_manual(values = c("PCGs_Sig" = "mediumorchid1", "PCGs_NotSig" = "plum1", 
                               "SCClncs_Sig" = "olivedrab3", "SCClncs_NotSig" = "lemonchiffon")) +
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  Seurat::RotatedAxis()

#TFs without enrichment or with much lesser enrichment in PCG:
SCClnc_TFs <- filter(trialiii, variableii == "SCClncs_Sig")
PCG_TFs <- filter(trialiii, variableii == "PCGs_Sig")

SCClnc_TFs_exPCG <- filter(SCClnc_TFs, !TF %in% PCG_TFs$TF)

#or also allow TFs with smaller enrichment in PCG
ReMap_Up4_PCG_v_SCCL$OR_diff <- ReMap_Up4_PCG_v_SCCL$SCClncs/ReMap_Up4_PCG_v_SCCL$PCGs

SCClnc_TFs_exPCG2 <- filter(SCClnc_TFs, TF %in% filter(ReMap_Up4_PCG_v_SCCL, OR_diff >1.5)$TF)


#
#### compare SCClnc and other lncs ####

#rbind plots to compare PCG and SCClnc:
ReMap_Up4_OtherL$GeneType <- "OtherLnc"
ReMap_Up4_SCCL$GeneType <- "SCCLnc"

#enriched experiments in one or other:
enrichedTFs <- unique(c(filter(ReMap_Up4_OtherL,  V5 < 0.01, V6 > 1.5, V1>=4)$peakTF,
                        filter(ReMap_Up4_SCCL,  V5 < 0.01, V6 > 1.5, V1>=4)$peakTF))

#create a TF ordering by strongest in SCCL vs strongest in DE PCG:
trial <- filter(ReMap_Up4_OtherL,  peakTF %in% enrichedTFs)
triali <- filter(ReMap_Up4_SCCL,  peakTF %in% enrichedTFs)

ReMap_Up4_OtherL_v_SCCL <- merge(trial[,c(7,11, 5, 6,10)], triali[,c(7,11, 5, 6,10)], by = "peakTF", all = T)

ReMap_Up4_OtherL_v_SCCL$peakTF <- as.factor(ReMap_Up4_OtherL_v_SCCL$peakTF)
orderTFs <- ReMap_Up4_OtherL_v_SCCL$peakTF[order(-ReMap_Up4_OtherL_v_SCCL$V6.y, ReMap_Up4_OtherL_v_SCCL$V6.x)]

colnames(ReMap_Up4_OtherL_v_SCCL)[c(4,8)] <- c("OtherLncs", "SCClncs")
colnames(ReMap_Up4_OtherL_v_SCCL)[c(3,7)] <- c("OtherLncs.p", "SCClncs.p")

trialiii <- reshape2::melt(ReMap_Up4_OtherL_v_SCCL[,c(1,4,8)])

trialiii$peakTF <- factor(trialiii$peakTF)
trialiii$peakTF <- factor(trialiii$peakTF,
                          levels = levels(trialiii$peakTF)[
                            match(orderTFs, levels(trialiii$peakTF))
                          ])

trialiii$p_val <- reshape2::melt(ReMap_Up4_OtherL_v_SCCL[,c(1,3,7)])[,3]
trialiii$p_val_sig <- "NotSig"
trialiii$p_val_sig[trialiii$p_val < 0.01] <- "Sig"

trialiii$variableii <- paste(trialiii$variable, trialiii$p_val_sig, sep = "_")

ggplot(filter(trialiii)) + aes(x = peakTF, y = value, fill = variableii) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "grey40") +
  ylab("Odds ratio - hits in selected\n vs. background promoters") +
  scale_fill_manual(values = c("OtherLncs_Sig" = "mediumorchid1", "OtherLncs_NotSig" = "plum1", 
                               "SCClncs_Sig" = "olivedrab3", "SCClncs_NotSig" = "lemonchiffon")) +
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  Seurat::RotatedAxis()

#simplfying, taking the max per TF

#take best result out of all:
trialsplit <- split(trial, trial$TF)
trialspliti <- lapply(trialsplit, function(x){
  x[which(x$V6 == max(x$V6)),]
})
trial_perTF <- bind_rows(trialspliti)

trialsplit <- split(triali, triali$TF)
trialspliti <- lapply(trialsplit, function(x){
  x[which(x$V6 == max(x$V6)),]
})
triali_perTF <- bind_rows(trialspliti)

ReMap_Up4_OtherL_v_SCCL <- merge(trial_perTF[,c(8,12, 5,6, 10,11)], triali_perTF[,c(8,12, 5,6, 10,11)], by = "TF", all = T)

ReMap_Up4_OtherL_v_SCCL$peakTF <- as.factor(ReMap_Up4_OtherL_v_SCCL$TF)
orderTFs <- ReMap_Up4_OtherL_v_SCCL$peakTF[order(-ReMap_Up4_OtherL_v_SCCL$FoldEnr.y, ReMap_Up4_OtherL_v_SCCL$FoldEnr.x)]

colnames(ReMap_Up4_OtherL_v_SCCL)[c(5,10)] <- c("OtherLncs", "SCClncs")
colnames(ReMap_Up4_OtherL_v_SCCL)[c(4,9)] <- c("OtherLncs.p", "SCClncs.p")

trialiii <- reshape2::melt(ReMap_Up4_OtherL_v_SCCL[,c(1,5,10)])

trialiii$peakTF <- factor(trialiii$TF)
trialiii$peakTF <- factor(trialiii$peakTF,
                          levels = levels(trialiii$peakTF)[
                            match(orderTFs, levels(trialiii$peakTF))
                          ])

trialiii$p_val <- reshape2::melt(ReMap_Up4_OtherL_v_SCCL[,c(1,4,9)])[,3]
trialiii$p_val_sig <- "(Fisher's p >0.05)"
trialiii$p_val_sig[trialiii$p_val < 0.01] <- "(Fisher's p <0.05)"

trialiii$variableii <- paste(trialiii$variable, trialiii$p_val_sig, sep = " ")

ggplot(filter(trialiii)) + aes(x = peakTF, y = value, fill = variableii) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "grey40") +
  ylab("Odds ratio - hits in selected\n vs. background promoters") +
  scale_fill_manual(values = c("OtherLncs (Fisher's p <0.05)" = "mediumorchid1", "OtherLncs (Fisher's p >0.05)" = "plum1", 
                               "SCClncs (Fisher's p <0.05)" = "olivedrab3", "SCClncs (Fisher's p >0.05)" = "lemonchiffon")) +
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  Seurat::RotatedAxis()

#TFs without enrichment or with much lesser enrichment in PCG:
SCClnc_TFs <- filter(trialiii, variableii == "SCClncs_Sig")
Otherlnc_TFs <- filter(trialiii, variableii == "OtherLncs_Sig")

SCClnc_TFs_exOlnc <- filter(SCClnc_TFs, !TF %in% Otherlnc_TFs$TF)

#or significant with a much stronger OR than in other lncs
ReMap_Up4_OtherL_v_SCCL$OR_diff <- ReMap_Up4_OtherL_v_SCCL$SCClncs/ReMap_Up4_OtherL_v_SCCL$OtherLncs

SCClnc_TFs_exOlnc2 <- filter(SCClnc_TFs, TF %in% filter(ReMap_Up4_OtherL_v_SCCL, OR_diff >1.5)$TF)


#