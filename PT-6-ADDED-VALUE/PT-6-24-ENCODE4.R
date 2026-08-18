library(dplyr)
library(GenomicRanges)


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
#### import ENCODE .csvs ####

ENC4_overlap_exprs_files <- list.files("C://Users/mbennet5/Downloads/", pattern = ".csv", full.names = T)

ENC4_overlap_exprs_files <- ENC4_overlap_exprs_files[grepl("ENC4", ENC4_overlap_exprs_files)]

trial <- lapply(ENC4_overlap_exprs_files, read.csv)

head(trial[[2]])

keyCOls <- c(1,2,3,4,5,13,14,15,16)

trial <- lapply(trial, function(x){
  x[,keyCOls]
})

sapply(trial, dim)

#just expressed TFs well in 0-4hr:
fpkm_thresh <- 5
ReMaps_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >fpkm_thresh | Hour4_meanFPKM >fpkm_thresh))

trial <- lapply(trial, function(x){
  filter(x, factor %in% ReMaps_04$EnsName)
})

#full table:
triali <- bind_rows(trial)

length(unique(triali$factor))

#HMGB2 missing
sum(unique(triali$factor) %in% "HMGB2")
#RELA present obviously
sum(unique(triali$factor) %in% "RELA")
#others
sum(unique(triali$factor) %in% "KLF4")
sum(unique(triali$factor) %in% "HMGA2")
sum(unique(triali$factor) %in% "FOXL1")
sum(unique(triali$factor) %in% "FOXC2")
sum(unique(triali$factor) %in% "FOXC1")
sum(unique(triali$factor) %in% "AR")
sum(unique(triali$factor) %in% "GMEB2")
sum(unique(triali$factor) %in% "GMEB1")
sum(unique(triali$factor) %in% "KLF9")

ENC4_overlap_exprs_df <- triali
rm(triali)
ENC4_overlap_exprs_df <- unique(ENC4_overlap_exprs_df)

#check scores per TF, scaled to 1000 now (for bed)
trial <- split(ENC4_overlap_exprs_df, ENC4_overlap_exprs_df$factor)

#443 TFs expressed in early timeframe with some ChIPseq found, expecting big mix of high/low confidence variation...

#take 5 randoms:
samplerTFs <- bind_rows(trial[c(1:10)])
table(samplerTFs$factor)
samplerTFs <- filter(samplerTFs, score >0)
table(samplerTFs$factor)

ggplot(samplerTFs) + aes(x = factor, y = score) +
  ggbeeswarm::geom_quasirandom(alpha = 0.2) +
  geom_boxplot(width = 0.25) +
  scale_y_log10() + Seurat::RotatedAxis()

#some of the TFs found quite ubiquitously:
topTFs_allBg <- c("JUN", "RELA", "FOS", "SMAD3", "SMC1A", "TEAD4", "BRD2", "BRD9")

bigTFs <- bind_rows(trial[topTFs_allBg])
bigTFs <- filter(bigTFs, score >0)
table(bigTFs$factor)

ggplot(bigTFs) + aes(x = factor, y = score) +
  ggbeeswarm::geom_quasirandom(alpha = 0.2) +
  geom_boxplot(width = 0.25) +
  scale_y_log10() + Seurat::RotatedAxis()

#takeaway, missing HMGB2, but a useful comparison
#conf scoring seems a little more documented, could do a straight cut (<10 or <50 etc)
#but just try a normal analysis as done for ReMap2022 first

#
#### (skip)find all TFs bound to each gene ####

enc4_SVSMC_prom_GR <- makeGRangesFromDataFrame(ENC4_overlap_exprs_df, 
                                               start.field = "chromStart", 
                                               end.field = "chromEnd", 
                                               seqnames = "X.chrom", keep.extra.columns = T)

#(get single yes or no per tf per gene)
enc4_SVSMC_prom_GR$peakID <- paste(as.character(seqnames(enc4_SVSMC_prom_GR)),
                                    start(enc4_SVSMC_prom_GR), 
                                    end(enc4_SVSMC_prom_GR), 
                                    enc4_SVSMC_prom_GR$name, sep = "_")

#prom ID for promoter co-ordinates:
allGB_GR_promoters_reduce$promID <- paste(as.character(seqnames(allGB_GR_promoters_reduce)),
                                          start(allGB_GR_promoters_reduce), end(allGB_GR_promoters_reduce), 
                                          #allGB_GR_promoters_reduce$MSTRG_Tx_ID, 
                                          sep = "_")

#add in gene ID per reduced prom:
enc4index <- findOverlaps(query = allGB_GR_promoters, subject = allGB_GR_promoters_reduce)

#return matched peaks from enc4:
allPromRed_EnsID <- unique(data.frame("EnsID" = allGB_GR_promoters$EnsID[queryHits(enc4index)],
                                      "promID" = allGB_GR_promoters_reduce$promID[subjectHits(enc4index)]))

#now make match of reduced proms to enc4 peaks
enc4index <- findOverlaps(query = enc4_SVSMC_prom_GR, subject = allGB_GR_promoters_reduce)

#return matched peaks from enc4:
allProm_enc4 <- data.frame("peakID" = enc4_SVSMC_prom_GR$peakID[queryHits(enc4index)],
                           "peakTF" = enc4_SVSMC_prom_GR$factor[queryHits(enc4index)],
                           "peakScore" = enc4_SVSMC_prom_GR$score[queryHits(enc4index)],
                                #"peakTF_cell" = enc4_SVSMC_prom_GR$[queryHits(enc4index)],
                                "promID" = allGB_GR_promoters_reduce$promID[subjectHits(enc4index)]
                                )

#now can combine all above to get for each gene ID, add any bound TFs in enc4 
trial <- unique(merge(unique(allProm_enc4[,c(2:4)]), allPromRed_EnsID, by = "promID", all.y = T))

length(unique(trial$promID))#18838 promoters
length(unique(allProm_enc4$promID))#17241
length(unique(allGB_GR_promoters_reduce$promID))#18838 (probs inc some without a remap overlap)

length(unique(trial$EnsID))#12740 genes
length(unique(allGB_GR_promoters$EnsID))#12740 genes
length(unique(fpkm_allG$EnsID))#12740 genes
#length(unique(filter(gencode_v26_gtf_table_extra, !gene_id %in% allGB$EnsID)$gene_id))#extra 5939 genes in background (if using)
length(unique(allPromRed_EnsID$EnsID))#12740 genes

allG_enc4_TFs <- trial
rm(trial)
gc()

#write.csv(allG_enc4_TFs, "allG_enc4_TFs.csv", row.names = F)


#
#### (skip)0-4hr run ####

#use this to find ReMap regulators which are expressed / with high potential to be active in each window:
fpkm_thresh <- 5
ReMaps_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >fpkm_thresh | Hour4_meanFPKM >fpkm_thresh))

#potential background, promoters expressed per window:
expr_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)$EnsID)

#test any TFs expressed in window of interest, with a peak in target gene set (here up4):
allG_enc4_TFs_window <- filter(allG_enc4_TFs, peakScore > 80,
                                    EnsID %in% fpkm_allGDE_Upwithin_4$EnsID,
                                    peakTF %in% ReMaps_04$EnsName)

table(allG_enc4_TFs_window$peakTF)

TFs_inTargetProms <- unique(allG_enc4_TFs_window$peakTF)
TFs_inTargetProms <- TFs_inTargetProms[!is.na(TFs_inTargetProms)]

#promoter regions for lncRNAs in given cluster of interest 
TargetGenes <- allG_enc4_TFs_window[allG_enc4_TFs_window$EnsID %in% fpkm_allGDE_Upwithin_4$EnsID,]

TFs_in_4hrDEGProms_all_Fish <- list()

#restrict (less sig TFs, more confidence) or expand background (more sig TFs, requires harsh filtering):
BackgroundGenes <- filter(allG_enc4_TFs, peakScore > 80,
                          #hash in or out
                          EnsID %in% expr_04
                          )

for(i in 1:length(TFs_inTargetProms)){
  
  #optional: remove bottom 20% lowest scoring peaks per TF (rather than simple cut off per TF done above)
  
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
#Prom2k_enc4_TF5fpkm_Prom0.8 <- TFs_in_4hrDEGProms_all_Fish_df
#Prom2k_enc4_TF5fpkm_Prom0.8_conf50 <- TFs_in_4hrDEGProms_all_Fish_df
Prom2k_enc4_TF5fpkm_Prom0.8_conf80 <- TFs_in_4hrDEGProms_all_Fish_df

#Prom2k_enc4_TF5fpkm_Prom2 <- TFs_in_4hrDEGProms_all_Fish_df
library(ggplot2)

ggplot(Prom2k_enc4_TF5fpkm_Prom0.8) + aes(x = -log10(pval_BH), y = V6) +
  geom_point() +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")

ggplot(Prom2k_enc4_TF5fpkm_Prom0.8_conf50) + aes(x = -log10(pval_BH), y = V6) +
  geom_point() +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")

#less power than ReMap2022


#### 200bp promoters, find all TFs bound to each gene ####

#per gene yields a handful of enriched TFs
#switching to smaller promoters, to a) fit with AME b) see if increases TP signal
allGB_GR_promoters200 <- promoters(allGB_GR, upstream = 200)

#flatten to aid next steps
allGB_GR_promoters_reduce200 <- reduce(allGB_GR_promoters200)
length(allGB_GR_promoters_reduce200)#21805 regions (about 3k more than 2000)


enc4_SVSMC_prom_GR <- makeGRangesFromDataFrame(ENC4_overlap_exprs_df, 
                                               start.field = "chromStart", 
                                               end.field = "chromEnd", 
                                               seqnames = "X.chrom", keep.extra.columns = T)

#repeat assignment of remap peaks
enc4_SVSMC_prom_GR$peakID <- paste(as.character(seqnames(enc4_SVSMC_prom_GR)),
                                    start(enc4_SVSMC_prom_GR), 
                                    end(enc4_SVSMC_prom_GR), 
                                    enc4_SVSMC_prom_GR$V4, sep = "_")

#prom ID for promoter co-ordinates:
allGB_GR_promoters_reduce200$promID <- paste(as.character(seqnames(allGB_GR_promoters_reduce200)),
                                             start(allGB_GR_promoters_reduce200), end(allGB_GR_promoters_reduce200), 
                                             #allGB_GR_promoters_reduce200$MSTRG_Tx_ID, 
                                             sep = "_")

#add in gene ID per reduced prom:
enc4index <- findOverlaps(query = allGB_GR_promoters, subject = allGB_GR_promoters_reduce200)

#return matched peaks from ReMap2022:
allPromRed_EnsID_200 <- unique(data.frame("EnsID" = allGB_GR_promoters$EnsID[queryHits(enc4index)],
                                          "promID" = allGB_GR_promoters_reduce200$promID[subjectHits(enc4index)]))

#now make match of reduced proms to remap peaks
enc4index <- findOverlaps(query = enc4_SVSMC_prom_GR, subject = allGB_GR_promoters_reduce200)

#return matched peaks from ReMap2022:
allProm_enc4_200 <- data.frame("peakID" = enc4_SVSMC_prom_GR$peakID[queryHits(enc4index)],
                                    "peakTF" = enc4_SVSMC_prom_GR$factor[queryHits(enc4index)],
                               "peakScore" = enc4_SVSMC_prom_GR$score[queryHits(enc4index)],
                               "Experiments" = enc4_SVSMC_prom_GR$exp[queryHits(enc4index)],
                                    #"peakTF_cell" = enc4_SVSMC_prom_GR$V6[queryHits(enc4index)],
                                    "promID" = allGB_GR_promoters_reduce200$promID[subjectHits(enc4index)])

#now can combine all above to get for each gene ID, add any bound TFs in ReMap 
trial <- unique(merge(unique(allProm_enc4_200[,c(2:5)]), allPromRed_EnsID_200, by = "promID", all.y = T))

length(unique(trial$promID))#21805 promoters
length(unique(allProm_remap2022_200$promID))#31751
length(unique(allGB_GR_promoters_reduce200$promID))#21805 (probs inc some without a remap overlap)

length(unique(trial$EnsID))#12740 genes -i.e. even if a gene is not bound by a TF, should be in here
length(unique(allGB_GR_promoters$EnsID))#12740 genes
length(unique(fpkm_allG$EnsID))#12740 genes
length(unique(allPromRed_EnsID_200$EnsID))#12740 genes

allG_enc4_TFs_200 <- trial
rm(trial)
gc()

#write.csv(allG_enc4_TFs_200, "allG_enc4_200_TFs.csv", row.names = F)


#
#### 0-4hr run - 200bp ####

#use this to find ReMap regulators which are expressed / with high potential to be active in each window:
fpkm_thresh <- 5
ReMaps_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >fpkm_thresh | Hour4_meanFPKM >fpkm_thresh))

#potential background, promoters expressed per window:
expr_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)$EnsID)

#alter DE strength:
fpkm_allGDE_Upwithin_4hi <- filter(fpkm_allGDE_Upwithin_4, preadj_0_4 < 0.01, LogFC_0_4> log2(2))

#test any TFs expressed in window of interest, with a peak in target gene set (here up4):
allG_enc4_TFs_window <- filter(allG_enc4_TFs_200, #peakScore >50,
                                    EnsID %in% fpkm_allGDE_Upwithin_4$EnsID,
                                    peakTF %in% ReMaps_04$EnsName)

TFs_inTargetProms <- unique(allG_enc4_TFs_window$peakTF)
TFs_inTargetProms <- TFs_inTargetProms[!is.na(TFs_inTargetProms)]

#promoter regions for lncRNAs in given cluster of interest 
TargetGenes <- allG_enc4_TFs_window[allG_enc4_TFs_window$EnsID %in% fpkm_allGDE_Upwithin_4$EnsID,]

TFs_in_4hrDEGProms_all_Fish <- list()

#restrict (less sig TFs, more confidence) or expand background (more sig TFs, requires harsh filtering):
BackgroundGenes <- filter(allG_enc4_TFs_200, #peakScore >50,
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
#Prom0.2k_enc4_TF5fpkm_Prom0.8 <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_enc4_TF5fpkm_Prom0.8_bigDEGs <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_enc4_TF5fpkm_Prom0.8_conf50 <- TFs_in_4hrDEGProms_all_Fish_df
Prom0.2k_enc4_TF5fpkm_Prom0.8_conf50_bigDEGs <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_enc4_TF5fpkm_Prom0.8_conf80 <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_enc4_TF5fpkm_Prom0.8_conf80_bigDEGs <- TFs_in_4hrDEGProms_all_Fish_df

#still HDAC6, but more sig tho

ggplot(Prom0.2k_enc4_TF5fpkm_Prom0.8) + aes(x = -log10(pval_BH), y = V6) +
  geom_point() +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")

#loss of HDAC6
ggplot(Prom0.2k_enc4_TF5fpkm_Prom0.8_conf50) + aes(x = -log10(pval_BH), y = V6) +
  geom_point() +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")

#
ggplot(Prom0.2k_enc4_TF5fpkm_Prom0.8_conf80) + aes(x = -log10(pval_BH), y = V6) +
  geom_point() +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")

#single result is not v convincing

#
#### 

#### 0-4hr run - 200bp - by experiment ####

#breaking down into individual experiments could work better
#as cell type diffs would not be merged, inappropriate experiments could be lost

#use this to find ReMap regulators which are expressed / with high potential to be active in each window:
fpkm_thresh <- 5
ReMaps_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >fpkm_thresh | Hour4_meanFPKM >fpkm_thresh))

#potential background, promoters expressed per window:
expr_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)$EnsID)

#alter DE strength:
fpkm_allGDE_Upwithin_4hi <- filter(fpkm_allGDE_Upwithin_4, preadj_0_4 < 0.01, LogFC_0_4> log2(2))

#test any TFs expressed in window of interest, with a peak in target gene set (here up4):
allG_enc4_TFs_window <- filter(allG_enc4_TFs_200, #peakScore >50,
                               EnsID %in% fpkm_allGDE_Upwithin_4$EnsID,
                               peakTF %in% ReMaps_04$EnsName)

#per experiment:
allG_enc4_TFs_window_exp <- allG_enc4_TFs_window %>%
  tidyr::separate_rows(Experiments, sep = "\\s*,\\s*")

table(allG_enc4_TFs_window_exp$Experiments)
length(unique(allG_enc4_TFs_window_exp$Experiments))

#each exp should be for one TF?
allG_enc4_TFs_window_exp$peakTF_experiment <- paste(allG_enc4_TFs_window_exp$peakTF, 
                                                    allG_enc4_TFs_window_exp$Experiments, sep = "_")
length(unique(allG_enc4_TFs_window_exp$peakTF_experiment))
#yes

TFs_inTargetProms <- unique(allG_enc4_TFs_window_exp$peakTF_experiment)
TFs_inTargetProms <- TFs_inTargetProms[!is.na(TFs_inTargetProms)]

#promoter regions for lncRNAs in given cluster of interest 
TargetGenes <- allG_enc4_TFs_window_exp[allG_enc4_TFs_window_exp$EnsID %in% fpkm_allGDE_Upwithin_4$EnsID,]

TFs_in_4hrDEGProms_all_Fish <- list()

#restrict (less sig TFs, more confidence) or expand background (more sig TFs, requires harsh filtering):
BackgroundGenes <- filter(allG_enc4_TFs_200, #peakScore >50,
                          #hash in or out
                          EnsID %in% expr_04,
                          )

#also needs to be expanded:
BackgroundGenes_exp <- BackgroundGenes %>%
  tidyr::separate_rows(Experiments, sep = "\\s*,\\s*")
BackgroundGenes_exp$peakTF_experiment <- paste(BackgroundGenes_exp$peakTF, 
                                               BackgroundGenes_exp$Experiments, sep = "_")

#now running it for 1642 experiments...
#may mean there is far more multiple hypothesis issues
#may mean there is far more sig tests tho too...
for(i in 1201:1642){
  
  a <- length(unique(filter(TargetGenes, peakTF_experiment == TFs_inTargetProms[i])$EnsID))
  b <- length(unique(TargetGenes$EnsID))
  c <- length(unique(filter(BackgroundGenes_exp, peakTF_experiment == TFs_inTargetProms[i])$EnsID))
  d <- length(unique(BackgroundGenes_exp$EnsID))
  
  TFs_in_4hrDEGProms_all_Fish[[i]] <- c(
    a,b,c,d,
    fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                           "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p,
    fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                           "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est
  )
}

names(TFs_in_4hrDEGProms_all_Fish) <- TFs_inTargetProms#[1:1200]

TFs_in_4hrDEGProms_all_Fish_df <- as.data.frame(t(bind_rows(TFs_in_4hrDEGProms_all_Fish, .id = "TFs")))

TFs_in_4hrDEGProms_all_Fish_df$peakTF <- TFs_inTargetProms#[1:1200]

TFs_in_4hrDEGProms_all_Fish_df$pval_BH <- p.adjust(TFs_in_4hrDEGProms_all_Fish_df$V5, method = "BH")

TFs_in_4hrDEGProms_all_Fish_df$Target_perc <- TFs_in_4hrDEGProms_all_Fish_df$V1/TFs_in_4hrDEGProms_all_Fish_df$V2 *100
TFs_in_4hrDEGProms_all_Fish_df$Bg_perc <- TFs_in_4hrDEGProms_all_Fish_df$V3/TFs_in_4hrDEGProms_all_Fish_df$V4 *100
TFs_in_4hrDEGProms_all_Fish_df$FoldChange <- TFs_in_4hrDEGProms_all_Fish_df$Target_perc/TFs_in_4hrDEGProms_all_Fish_df$Bg_perc

#store some variants here
Prom0.2k_enc4split_TF5fpkm_Prom0.8 <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_enc4split_TF5fpkm_Prom0.8_bigDEGs <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_enc4split_TF5fpkm_Prom0.8_conf50 <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_enc4split_TF5fpkm_Prom0.8_conf50_bigDEGs <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_enc4split_TF5fpkm_Prom0.8_conf80 <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_enc4split_TF5fpkm_Prom0.8_conf80_bigDEGs <- TFs_in_4hrDEGProms_all_Fish_df

#still HDAC6, but more sig tho

ggplot(Prom0.2k_enc4_TF5fpkm_Prom0.8) + aes(x = -log10(pval_BH), y = V6) +
  geom_point() +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")

#loss of HDAC6
ggplot(Prom0.2k_enc4split_TF5fpkm_Prom0.8) + aes(x = -log10(pval_BH), y = V6) +
  geom_point() +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")

#
ggplot(Prom0.2k_enc4split_TF5fpkm_Prom0.8_conf80) + aes(x = -log10(pval_BH), y = V6) +
  geom_point() +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")

#single result is not v convincing

#
#### 

#### 0-4hr scclncRNA run - 200bp - by experiment ####

#breaking down into individual experiments could work better
#as cell type diffs would not be merged, inappropriate experiments could be lost

#use this to find ReMap regulators which are expressed / with high potential to be active in each window:
fpkm_thresh <- 20
ReMaps_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >fpkm_thresh | Hour4_meanFPKM >fpkm_thresh))

SCClncRNAs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SCClncRNAs.csv")
SCClncRNAs_04 <- filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs")

#potential background, promoters expressed per window:
expr_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)$EnsID)

#alter DE strength:
#fpkm_allGDE_Upwithin_4hi <- filter(fpkm_allGDE_Upwithin_4, preadj_0_4 < 0.01, LogFC_0_4> log2(2))

#test any TFs expressed in window of interest, with a peak in target gene set (here up4):
allG_enc4_TFs_window <- filter(allG_enc4_TFs_200, #peakScore >50,
                               EnsID %in% SCClncRNAs_04$EnsID,
                               peakTF %in% ReMaps_04$EnsName)

#per experiment:
allG_enc4_TFs_window_exp <- allG_enc4_TFs_window %>%
  tidyr::separate_rows(Experiments, sep = "\\s*,\\s*")

table(allG_enc4_TFs_window_exp$Experiments)
length(unique(allG_enc4_TFs_window_exp$Experiments))

#each exp should be for one TF?
allG_enc4_TFs_window_exp$peakTF_experiment <- paste(allG_enc4_TFs_window_exp$peakTF, 
                                                    allG_enc4_TFs_window_exp$Experiments, sep = "_")
length(unique(allG_enc4_TFs_window_exp$peakTF_experiment))
#yes

TFs_inTargetProms <- unique(allG_enc4_TFs_window_exp$peakTF_experiment)
TFs_inTargetProms <- TFs_inTargetProms[!is.na(TFs_inTargetProms)]

#promoter regions for lncRNAs in given cluster of interest 
TargetGenes <- allG_enc4_TFs_window_exp[allG_enc4_TFs_window_exp$EnsID %in% fpkm_allGDE_Upwithin_4$EnsID,]

TFs_in_4hrDEGProms_all_Fish <- list()

#restrict (less sig TFs, more confidence) or expand background (more sig TFs, requires harsh filtering):
BackgroundGenes <- filter(allG_enc4_TFs_200, #peakScore >50,
                          #hash in or out
                          EnsID %in% expr_04,
)

#also needs to be expanded:
BackgroundGenes_exp <- BackgroundGenes %>%
  tidyr::separate_rows(Experiments, sep = "\\s*,\\s*")
BackgroundGenes_exp$peakTF_experiment <- paste(BackgroundGenes_exp$peakTF, 
                                               BackgroundGenes_exp$Experiments, sep = "_")

#now running it for 1642 experiments...
#may mean there is far more multiple hypothesis issues
#may mean there is far more sig tests tho too...
for(i in 1:343){
  
  a <- length(unique(filter(TargetGenes, peakTF_experiment == TFs_inTargetProms[i])$EnsID))
  b <- length(unique(TargetGenes$EnsID))
  c <- length(unique(filter(BackgroundGenes_exp, peakTF_experiment == TFs_inTargetProms[i])$EnsID))
  d <- length(unique(BackgroundGenes_exp$EnsID))
  
  TFs_in_4hrDEGProms_all_Fish[[i]] <- c(
    a,b,c,d,
    fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                           "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p,
    fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                           "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est
  )
}

names(TFs_in_4hrDEGProms_all_Fish) <- TFs_inTargetProms#[1:800]

TFs_in_4hrDEGProms_all_Fish_df <- as.data.frame(t(bind_rows(TFs_in_4hrDEGProms_all_Fish, .id = "TFs")))

TFs_in_4hrDEGProms_all_Fish_df$peakTF <- TFs_inTargetProms#[1:800]

TFs_in_4hrDEGProms_all_Fish_df$pval_BH <- p.adjust(TFs_in_4hrDEGProms_all_Fish_df$V5, method = "BH")

TFs_in_4hrDEGProms_all_Fish_df$Target_perc <- TFs_in_4hrDEGProms_all_Fish_df$V1/TFs_in_4hrDEGProms_all_Fish_df$V2 *100
TFs_in_4hrDEGProms_all_Fish_df$Bg_perc <- TFs_in_4hrDEGProms_all_Fish_df$V3/TFs_in_4hrDEGProms_all_Fish_df$V4 *100
TFs_in_4hrDEGProms_all_Fish_df$FoldChange <- TFs_in_4hrDEGProms_all_Fish_df$Target_perc/TFs_in_4hrDEGProms_all_Fish_df$Bg_perc

#store some variants here
Prom0.2k_enc4split_TF5fpkm_SCClncProm0.8 <- TFs_in_4hrDEGProms_all_Fish_df
Prom0.2k_enc4split_TF20fpkm_SCClncProm0.8 <- TFs_in_4hrDEGProms_all_Fish_df

#nothing convincing
