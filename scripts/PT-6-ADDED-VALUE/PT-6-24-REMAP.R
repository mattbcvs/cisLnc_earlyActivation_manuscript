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
#### (not recommended) alternative input ####

#not just promoters expressed in the SVSMC, but all promoters in GENCODE

#more background, to elevate significance of hits

#note: this is analogous to probs most GO analyses, which do not set the background to be just expressed genes

#additional GENv26 genes:
gencode_v26_gtf_table <- read.delim("heavyFiles4R/gencode.v26.primary_assembly.annotation.gtf", 
                                    header=FALSE, stringsAsFactors=FALSE, sep = "\t", skip = 5)
#prom regions needs tx
# (stick to well supported PCGs for more reliable)
gencode_v26_gtf_table <- filter(gencode_v26_gtf_table, V3 == "transcript", grepl("protein_coding", V9), grepl("support_level 1", V9))

gencode_v26_gtf_table$tx_id <- gsub("transcript_id ", "", sapply(strsplit(gencode_v26_gtf_table$V9, "; "), "[[", 2))
gencode_v26_gtf_table$gene_id <- gsub("gene_id ", "", sapply(strsplit(gencode_v26_gtf_table$V9, "; "), "[[", 1))

#additional 23k strong tx for allGB_GR:
gencode_v26_gtf_table_extra <- filter(gencode_v26_gtf_table, !tx_id %in% allGB$MSTRG_Tx_ID)

#extra 23k transcripts to use
head(allGB)
head(gencode_v26_gtf_table_extra)

appendtoAll <- gencode_v26_gtf_table_extra[,c(10,11,1,4,5,7)]

head(allGB)
head(appendtoAll)

allGB_append <- allGB[,1:6]

colnames(appendtoAll) <- colnames(allGB_append)

allGB_append <- rbind(allGB_append, appendtoAll)

allGB_GR <- makeGRangesFromDataFrame(allGB_append, 
                                     start.field = "Tx_start", 
                                     end.field = "Tx_stop", 
                                     seqnames.field = "chr", 
                                     strand.field = "str", keep.extra.columns = T)

#consider shortening later 2000bp is a lot
allGB_GR_promoters <- promoters(allGB_GR)

#flatten to aid next steps
allGB_GR_promoters_reduce <- reduce(allGB_GR_promoters)
length(allGB_GR_promoters_reduce)#28716 regions

#bed format:
allGB_GR_promoters_reduce_bed <- as.data.frame(allGB_GR_promoters_reduce)

write.table(allGB_GR_promoters_reduce_bed[,c(1:3,6)], "allGB_GR_promoters_reduce_bed.bed", row.names = F, col.names = F, quote = F)


#
#### remap explore ####

#sampler, what input is available for non-redundant
remapCatalog2022hg38 <- read.delim("heavyFiles4R/remap2022_nr_macs2_hg38_v1_0.bed", nrows = 50000,
                                   header = F)
head(remapCatalog2022hg38, 30)
table(remapCatalog2022hg38$V5)
#looks likely to be number of peaks at a given nr locus
#no obvious confidence scoring? which is a bit worrying

#sampler of peaks
remapCatalog2022hg38_peaks <- read.delim("heavyFiles4R/remap2022_all_macs2_hg38_v1_0.bed", nrows = 50000,
                                   header = F)
#number of peaks in nr region chr1:9902-10328 for ZBTB40 expected to be 4 based on number in nr peaks table:
filter(remapCatalog2022hg38_peaks, V2 > 9801, V3 < 10429, V1 == "chr1", grepl("ZBTB40", V4))
#matches

#low confidence peaks rate? include a few more:
remapCatalog2022hg38_peaks <- read.delim("heavyFiles4R/remap2022_all_macs2_hg38_v1_0.bed", skip = 4000000, nrows = 500000,
                                         header = F)
remapCatalog2022hg38_peaks$TF <- sapply(strsplit(remapCatalog2022hg38_peaks$V4, "\\."), "[[", 2)

trial <- split(remapCatalog2022hg38_peaks, remapCatalog2022hg38_peaks$TF)

#1163 TFs found, expecting big mix of high/low confidence variation...

#take 5 randoms:
samplerTFs <- bind_rows(trial[c(1:10)])
table(samplerTFs$TF)
samplerTFs <- filter(samplerTFs, V5 >0)
table(samplerTFs$TF)

ggplot(samplerTFs) + aes(x = TF, y = V5) +
  ggbeeswarm::geom_quasirandom(alpha = 0.2) +
  geom_boxplot(width = 0.25) +
  scale_y_log10() + Seurat::RotatedAxis()

#some of the TFs found quite ubiquitously:
topTFs_allBg <- c("JUN", "RELA", "FOS", "SMAD3", "SMC1A", "TEAD4", "BRD2", "BRD9")

bigTFs <- bind_rows(trial[topTFs_allBg])
bigTFs <- filter(bigTFs, V5 >0)
table(bigTFs$TF)

ggplot(bigTFs) + aes(x = TF, y = V5) +
  ggbeeswarm::geom_quasirandom(alpha = 0.2) +
  geom_boxplot(width = 0.25) +
  scale_y_log10() + Seurat::RotatedAxis()

notableTFs_allBg <- c("NFKB2", "HMGB2", "GPS2")

notableTFs <- bind_rows(trial[notableTFs_allBg])
notableTFs <- filter(notableTFs, V5 >0)
table(notableTFs$TF)

ggplot(notableTFs) + aes(x = TF, y = V5) +
  ggbeeswarm::geom_quasirandom(alpha = 0.2) +
  geom_boxplot(width = 0.25) +
  scale_y_log10() + Seurat::RotatedAxis()


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
#### (skip) find all peaks within early up promoters ####

#subset the allGB object for prom coords for early induced lncs
ER_GR_promoters <- allGB_GR_promoters[allGB_GR_promoters$EnsID %in% fpkm_allGDE_Upwithin_4$EnsID]

#match up these genes to regions in the remap object
remapindex <- findOverlaps(query = remap_SVSMC_prom_GR, subject = ER_GR_promoters)

#make a col to match on
remap_SVSMC_prom_GR$peakID <- paste(as.character(seqnames(remap_SVSMC_prom_GR)),
                                    start(remap_SVSMC_prom_GR), end(remap_SVSMC_prom_GR), remap_SVSMC_prom_GR$V4, sep = "_")

#return nr remap peaks found in early induced promoters:
ER_remap2022 <- data.frame("peakID" = remap_SVSMC_prom_GR$peakID[queryHits(remapindex)],
                              "peakTF_cell" = remap_SVSMC_prom_GR$V6[queryHits(remapindex)],
                              "LncRNA" = ER_GR_promoters$MSTRG_Tx_ID[subjectHits(remapindex)])

#if needed matched peaks as a GR:
ER_remap2022_GR <- remap_SVSMC_prom_GR[remap_SVSMC_prom_GR$peakID %in% ER_remap2022$peakID]


#column format, just the TF
ER_remap2022_GR$peakTF <- sapply(strsplit(ER_remap2022_GR$V6, ":"), "[[", 1)
remap_SVSMC_prom_GR$peakTF <- sapply(strsplit(remap_SVSMC_prom_GR$V6, ":"), "[[", 1)


#### alternatively, preferable, find all TFs bound to each gene ####

#option 1 could be driven by a few loci with lots of motifs
#per gene is fairer at this stage

#second option (single yes or no per tf per gene) looks like this
remap_SVSMC_prom_GR$peakID <- paste(as.character(seqnames(remap_SVSMC_prom_GR)),
                                    start(remap_SVSMC_prom_GR), 
                                    end(remap_SVSMC_prom_GR), 
                                    remap_SVSMC_prom_GR$V4, sep = "_")

remap_SVSMC_prom_GR$peakTF <- sapply(strsplit(remap_SVSMC_prom_GR$V6, ":"), "[[", 1)

#prom ID for promoter co-ordinates:
allGB_GR_promoters_reduce$promID <- paste(as.character(seqnames(allGB_GR_promoters_reduce)),
                                          start(allGB_GR_promoters_reduce), end(allGB_GR_promoters_reduce), 
                                          #allGB_GR_promoters_reduce$MSTRG_Tx_ID, 
                                          sep = "_")

#add in gene ID per reduced prom:
remapindex <- findOverlaps(query = allGB_GR_promoters, subject = allGB_GR_promoters_reduce)

#return matched peaks from ReMap2022:
allPromRed_EnsID <- unique(data.frame("EnsID" = allGB_GR_promoters$EnsID[queryHits(remapindex)],
                                      "promID" = allGB_GR_promoters_reduce$promID[subjectHits(remapindex)]))

#now make match of reduced proms to remap peaks
remapindex <- findOverlaps(query = remap_SVSMC_prom_GR, subject = allGB_GR_promoters_reduce)

#return matched peaks from ReMap2022:
allProm_remap2022 <- data.frame("peakID" = remap_SVSMC_prom_GR$peakID[queryHits(remapindex)],
                                "peakTF" = remap_SVSMC_prom_GR$peakTF[queryHits(remapindex)],
                                "peakTF_cell" = remap_SVSMC_prom_GR$V6[queryHits(remapindex)],
                                "promID" = allGB_GR_promoters_reduce$promID[subjectHits(remapindex)])

#now can combine all above to get for each gene ID, add any bound TFs in ReMap 
trial <- unique(merge(unique(allProm_remap2022[,c(2,4)]), allPromRed_EnsID, by = "promID", all.y = T))

length(unique(trial$promID))#18838 promoters
length(unique(allProm_remap2022$promID))#18376
length(unique(allGB_GR_promoters_reduce$promID))#18838 (probs inc some without a remap overlap)

length(unique(trial$EnsID))#12740 genes
length(unique(allGB_GR_promoters$EnsID))#12740 genes
length(unique(fpkm_allG$EnsID))#12740 genes
length(unique(filter(gencode_v26_gtf_table_extra, !gene_id %in% allGB$EnsID)$gene_id))#extra 5939 genes in background (if using)
length(unique(allPromRed_EnsID$EnsID))#12740 genes

allG_remap2022_TFs <- trial
rm(trial)
gc()

#write.csv(allG_remap2022_TFs, "allG_remap2022_TFs_extraPCGs.csv", row.names = F)


#
#### 0-4hr run ####

#use this to find ReMap regulators which are expressed / with high potential to be active in each window:
fpkm_thresh <- 5
ReMaps_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >fpkm_thresh | Hour4_meanFPKM >fpkm_thresh))

#potential background, promoters expressed per window:
expr_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)$EnsID)

#test any TFs expressed in window of interest, with a peak in target gene set (here up4):
allG_remap2022_TFs_window <- filter(allG_remap2022_TFs,
                                    EnsID %in% fpkm_allGDE_Upwithin_4$EnsID,
                                    peakTF %in% ReMaps_04$EnsName)

TFs_inTargetProms <- unique(allG_remap2022_TFs_window$peakTF)
TFs_inTargetProms <- TFs_inTargetProms[!is.na(TFs_inTargetProms)]

#promoter regions for lncRNAs in given cluster of interest 
TargetGenes <- allG_remap2022_TFs_window[allG_remap2022_TFs_window$EnsID %in% fpkm_allGDE_Upwithin_4$EnsID,]

TFs_in_4hrDEGProms_all_Fish <- list()

#restrict (less sig TFs, more confidence) or expand background (more sig TFs, requires harsh filtering):
BackgroundGenes <- filter(allG_remap2022_TFs, 
                          #hash in or out
                          #EnsID %in% expr_04
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
Prom2k_TF5fpkm_Prom0.8 <- TFs_in_4hrDEGProms_all_Fish_df
#Prom2k_TF5fpkm_PromAll <- TFs_in_4hrDEGProms_all_Fish_df

#Prom2k_TF5fpkm_Prom2 <- TFs_in_4hrDEGProms_all_Fish_df
library(ggplot2)

ggplot(Prom2k_TF5fpkm_Prom0.8) + aes(x = -log10(pval_BH), y = V6) +
  geom_point() +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")

ggplot(Prom2k_TF5fpkm_PromAll) + aes(x = -log10(pval_BH), y = V6) +
  geom_point() +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")

#same genes well covered in both?


#bit weakening if including less promoters

#NFKB2 - up4 - non-canonical NFkB pathway (p100 precursor processed into p52, slower + more persistent than canonical)
#JARID2 - up4 - PRC2 recruitment, gene silencer
#SUPT16H - not DE - FACT complex, H2A/B interactor
#TSCCD4 - not DE (see TSCCD1)
#BNC2 - not DE (see BNC1 tho)

#RELA, JUN, quite ubiquitous, suggest a 200bp only promoter version:


#### (best iteration) 200bp promoters, find all TFs bound to each gene ####

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
#Prom0.2k_TF5fpkm_PromAll <- TFs_in_4hrDEGProms_all_Fish_df

#no more TFs, but change in emphasis:

#JARID2 - up4 - fine-tuning, PRC2, need for proper development of heart/liver/spleen, multi-lineage differentiation
#GPS2 - down8 - co-repressor antagonises ERa, maybe pioneer, maybe anti-inflamm
#TSC22D4 - notDE (see 1 and 3) - interacts Akt-1 during nutrient deprivation, repress senescence
#BNC2 - notDE (see 1) - myofibroblast activation
#NFKB2 - up4 - non-canonical NFkB pathway

ggplot(Prom0.2k_TF5fpkm_Prom0.8) + aes(x = -log10(pval_BH), y = V6) +
  geom_point() +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")

ggplot(Prom0.2k_TF5fpkm_PromAll) + aes(x = -log10(pval_BH), y = V6) +
  geom_point() +
  geom_hline(yintercept = 1.1, linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")



#### (best iteration) 0-4hr run - 200bp - split experiments ####

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
#### identify TFs with affinity to SCClncs over PCGs or other lncs ####

ReMap_Up4_PCG_enr <- filter(ReMap_Up4_PCG)
trialsplit <- split(ReMap_Up4_PCG_enr, ReMap_Up4_PCG_enr$TF)
trialspliti <- lapply(trialsplit, function(x){
  x[which(x$V6 == max(x$V6)),]
})
ReMap_Up4_PCG_enr <- bind_rows(trialspliti)

ReMap_Up4_SCCL_enr <- filter(ReMap_Up4_SCCL)
trialsplit <- split(ReMap_Up4_SCCL_enr, ReMap_Up4_SCCL_enr$TF)
trialspliti <- lapply(trialsplit, function(x){
  x[which(x$V6 == max(x$V6)),]
})
ReMap_Up4_SCCL_enr <- bind_rows(trialspliti)

ReMap_Up4_OtherL_enr <- filter(ReMap_Up4_OtherL)
trialsplit <- split(ReMap_Up4_OtherL_enr, ReMap_Up4_OtherL_enr$TF)
trialspliti <- lapply(trialsplit, function(x){
  x[which(x$V6 == max(x$V6)),]
})
ReMap_Up4_OtherL_enr <- bind_rows(trialspliti)

#tile map of these ones:
trial <- merge(ReMap_Up4_PCG_enr[,c(11,5,6)], ReMap_Up4_SCCL_enr[,c(11,5,6)], by = "TF", all = T)
trial <- merge(trial, ReMap_Up4_OtherL_enr[,c(11,5,6)], by = "TF", all = T)
colnames(trial) <- c("TF", "PCG.p", "PCG.or", "SCCL.p", "SCCL.or", "lnc.p", "lnc.or")

triali <- reshape2::melt(trial[,c(1,3,5,7)])

enrichedTFs <- unique(c(filter(ReMap_Up4_OtherL_enr,  V5 < 0.01, V6 > 1.5, V1>=4)$TF,
                        filter(ReMap_Up4_SCCL_enr,  V5 < 0.01, V6 > 1.5, V1>=4)$TF,
                        filter(ReMap_Up4_PCG_enr,  V5 < 0.01, V6 > 1.5, V1>=4)$TF))

triali <- filter(triali, TF %in% enrichedTFs)

triali$TF <- as.factor(triali$TF)

#order by SCClncs:
trialii <- filter(trial, TF %in% enrichedTFs)
orderTFs <- trialii$TF[order(-trialii$SCCL.or, trialii$PCG.or)]
triali$TF <- factor(triali$TF,
                          levels = levels(triali$TF)[
                            match(orderTFs, levels(triali$TF))
                          ])

ggplot(triali) + aes(x = TF, y = variable, fill = log2(value)) +
  geom_tile() + Seurat::RotatedAxis()


#
#### compare to ISMARA/LISA ####

#quite a blunt method compared to others, do we find same TFs? 
#highly ranked TFs in ISMARA:

#write.csv(TF_ISMARA_long, "TF_ISMARA_ranking2026.csv", row.names = F)
TF_ISMARA_ranking2026 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/TF_ISMARA_ranking2026.csv")

#significantly enriched in 0-4hr up PCG promoters:
#ReMap_Up4_PCG$TF <- sapply(strsplit(ReMap_Up4_PCG$peakTF, "_"), "[[", 1)
#ReMap_Up4_PCG_enr <- filter(ReMap_Up4_PCG, V5 < 0.05, V6 > 1.25)

#take max fold enr (odds ratio doesn't work as well for this)
trialsplit <- split(ReMap_Up4_PCG, ReMap_Up4_PCG$TF)
trialspliti <- lapply(trialsplit, function(x){
  x[which(x$FoldEnr == max(x$FoldEnr)),]
})
trial_perTF <- bind_rows(trialspliti)

#higher fold enrichment should mean higher ISMARA ranking:
#n.b., only TFs >5 were profiled in ReMap2022:
TF_ISMARA_ranking2026_5FPKM <- filter(TF_ISMARA_ranking2026, Symbols %in% ReMaps_04$EnsName)

ISMARA_v_ReMap <- merge(TF_ISMARA_ranking2026_5FPKM[,c(2,4)], trial_perTF[,c(11,10,5)], by.x = "Symbols", by.y = "TF", all = T)

ISMARA_v_ReMap$Sig <- "Non_sig"
ISMARA_v_ReMap$Sig[ISMARA_v_ReMap$V5 < 0.05] <- "p<0.05"
ISMARA_v_ReMap$Sig[ISMARA_v_ReMap$V5 < 0.01] <- "p<0.01"
ISMARA_v_ReMap$Sig[ISMARA_v_ReMap$V5 < 0.001] <- "p<0.001"

ggplot(ISMARA_v_ReMap) + aes(x = Zscore, y = FoldEnr) +
  geom_point(aes(x = Zscore, y = FoldEnr, color = Sig), alpha = 0.6, size = 2) +
  theme_minimal() +
  scale_color_manual(values = c("p<0.05" = "grey80", "p<0.01" = "grey60", "p<0.001" = "grey10")) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "lm", color = "grey60") +
  xlab("ISMARA (Zscore)") +
  ylab("ReMap (Fold Enrichment)") +
  theme(text = element_text(size = 15))

cor.test(ISMARA_v_ReMap$Zscore, ISMARA_v_ReMap$FoldEnr) #nice correlation between the two


#
#### look for SCClncRNA targets amongst enriched TFs ####

SCClncRNAs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SCClncRNAs.csv")
SCClncRNAs_04 <- filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs")

#5x TFs in ReMap
filter(SCClncRNAs_04, EnsName.y %in% ReMap_Up4_PCG$TF)
#non sig
filter(SCClncRNAs_04, EnsName.y %in% filter(ReMap_Up4_PCG, V5 < 0.05)$TF)
#non with signs of enrich
filter(SCClncRNAs_04, EnsName.y %in% filter(ReMap_Up4_PCG, FoldEnr > 1.25)$TF)


#
#### look for SCClncRNAs bound by same TF as target ####

#this may signal a prospective mechanism linking the two loci
#e.g. lncRNA recruits a TF that binds the PCG
#particularly notable if the same TF comes up a lot
SCClncRNAs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SCClncRNAs.csv")
SCClncRNAs_04 <- filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs")

#can use any TF (even low exprs fine) for this one, but optionally:
allG_remap2022_TFs_200_exprs <- filter(allG_remap2022_TFs_200, peakTF %in% ReMaps_04$EnsName)

#per pair, identify co-binding TFs:
trial <- split(SCClncRNAs_04[,c(1,5,3,7)], SCClncRNAs_04$pairs)

TFoverlap_per_pair_numbers <- list()
TFoverlap_per_pair <- list()

for (i in 1:length(unique(SCClncRNAs_04$pairs))){
  lncTFs <- filter(allG_remap2022_TFs_200_exprs, EnsID %in% trial[[i]]$EnsID)
  PCGTFs <- filter(allG_remap2022_TFs_200_exprs, EnsID %in% trial[[i]]$EnsID.y)
  JointTFs <- merge(lncTFs, PCGTFs, by = "peakTF")
  #TFoverlap_per_pair[[i]] <- JointTFs
  
  TFoverlap_per_pair_numbers[[i]] <- c(length(unique(lncTFs$peakTF)),
                                       length(unique(PCGTFs$peakTF)),
                                       length(unique(JointTFs$peakTF)))
}

TFoverlap_per_pair_numbersi <- t(data.frame(TFoverlap_per_pair_numbers))
TFoverlap_per_pair_numbersi <- as.data.frame(TFoverlap_per_pair_numbersi)
TFoverlap_per_pair_numbersi$pairs <- names(trial)
TFoverlap_per_pair_numbersi <- merge(SCClncRNAs_04[,c(2,5,7)], TFoverlap_per_pair_numbersi)

TFoverlap_per_pair_numbersi$perc_LncTFs_inTarget <- TFoverlap_per_pair_numbersi$V3/TFoverlap_per_pair_numbersi$V1*100

#very common to get good TF overlap
summary(TFoverlap_per_pair_numbersi$perc_LncTFs_inTarget)

#what about shuffling pairs?
y0 <- SCClncRNAs_04$EnsID.y
y_new <- sample(y0)

while(any(y_new == y0)) {
  y_new <- sample(y0)
}

SCClncRNAs_04_bg <- transform(SCClncRNAs_04, EnsID.y = y_new)
SCClncRNAs_04_bg$pairs <- paste(SCClncRNAs_04_bg$EnsID, SCClncRNAs_04_bg$EnsID.y, sep = "-")

trial <- split(SCClncRNAs_04_bg[,c(1,5,3,7)], SCClncRNAs_04_bg$pairs)

TFoverlap_per_pair_numbers_bg <- list()

for (i in 1:length(unique(SCClncRNAs_04$pairs))){
  lncTFs <- filter(allG_remap2022_TFs_200_exprs, EnsID %in% trial[[i]]$EnsID)
  PCGTFs <- filter(allG_remap2022_TFs_200_exprs, EnsID %in% trial[[i]]$EnsID.y)
  JointTFs <- merge(lncTFs, PCGTFs, by = "peakTF")
  TFoverlap_per_pair[[i]] <- JointTFs
  
  TFoverlap_per_pair_numbers_bg[[i]] <- c(length(unique(lncTFs$peakTF)),
                                          length(unique(PCGTFs$peakTF)),
                                          length(unique(JointTFs$peakTF)))
}

TFoverlap_per_pair_numbersii <- t(data.frame(TFoverlap_per_pair_numbers_bg))
TFoverlap_per_pair_numbersii <- as.data.frame(TFoverlap_per_pair_numbersii)
TFoverlap_per_pair_numbersii$pairs <- names(trial)
TFoverlap_per_pair_numbersii$perc_LncTFs_inTarget <- TFoverlap_per_pair_numbersii$V3/TFoverlap_per_pair_numbersii$V1*100

#TFs >5FPKM
#very common to get good TF overlap
summary(TFoverlap_per_pair_numbersi$perc_LncTFs_inTarget)#68-90
#shuffled background show that this is probably elevated above random
summary(TFoverlap_per_pair_numbersii$perc_LncTFs_inTarget)#50-86, 50-87, 50-88


#all TFs:
#very common to get good TF overlap
summary(TFoverlap_per_pair_numbersi$perc_LncTFs_inTarget)#63-87
#shuffled background show that this is probably elevated above random
summary(TFoverlap_per_pair_numbersii$perc_LncTFs_inTarget)#53-81, 52-81, 50-81

boxplot(TFoverlap_per_pair_numbersi$perc_LncTFs_inTarget,
        TFoverlap_per_pair_numbersii$perc_LncTFs_inTarget)

t.test(TFoverlap_per_pair_numbersi$perc_LncTFs_inTarget,
       TFoverlap_per_pair_numbersii$perc_LncTFs_inTarget, var.equal = T)


#
####################
#### older code ####
####################
#### SCClnc specific view ####

#unlikely to get significance, however if they have a differing set of TFs then possible

#but at least helps confirm if patterns from all up4 apply also to these lncs

SCClncRNAs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SCClncRNAs.csv")
SCClncRNAs_04 <- filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs")

#use this to find ReMap regulators which are expressed / with high potential to be active in each window:
fpkm_thresh <- 5
ReMaps_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >fpkm_thresh | Hour4_meanFPKM >fpkm_thresh))

#potential background, promoters expressed per window:
expr_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)$EnsID)

#test any TFs expressed in window of interest, with a peak in target gene set (here up4):
allG_remap2022_TFs_window <- filter(allG_remap2022_TFs_200,
                                    EnsID %in% SCClncRNAs_04$EnsID,
                                    peakTF %in% c(fpkm_allGDE_Upwithin_4$EnsName, fpkm_allGDE_Downwithin_4$EnsName),
                                    peakTF %in% ReMaps_04$EnsName)

TFs_inTargetProms <- unique(allG_remap2022_TFs_window$peakTF)
TFs_inTargetProms <- TFs_inTargetProms[!is.na(TFs_inTargetProms)]

#promoter regions for lncRNAs in given cluster of interest 
TargetGenes <- allG_remap2022_TFs_window[allG_remap2022_TFs_window$EnsID %in% SCClncRNAs_04$EnsID,]

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
Prom0.2k_TF5fpkm_SCCLncProm0.8 <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_TF5fpkm_SCCLncPromAll <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_TF15fpkm_SCCLncPromAll <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_TF15fpkmDE_SCCLncPromAll <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_TF10fpkmDE_SCCLncPromAll <- TFs_in_4hrDEGProms_all_Fish_df

#Prom0.2k_TF10fpkm_SCCLncProm0.8 <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_TF10fpkm4DE_SCCLncProm0.8 <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_TF10fpkm4DE0.01_SCCLncProm0.8 <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_TF15fpkm4DE0.01_SCCLncProm0.8 <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_TF15fpkm4DE0.0001_SCCLncProm0.8 <- TFs_in_4hrDEGProms_all_Fish_df
#Prom0.2k_TF15fpkm4DE0.0001pFC1.5_SCCLncProm0.8 <- TFs_in_4hrDEGProms_all_Fish_df

#no sig TFs, can get closer with v. strongly DE TFs

#regardless HMGB2 is strongest signal, highly downregulated within 4hrs, far weaker amongst all 4hr up genes

#also DAXX - up4 - apoptosis, H3 chaperone, oncogene + suppressor.... 

#DAXX in particular has a quite good % (42% of SCClncs, 25% of general)
#HMGB2 bit weaker but good 19% vs 6%

#from previous set: NFKB2 is still around (14% vs. 9%) and in the highly DE set
#whilst SUPT16H (28% vs 24%), GPS2 (8% vs 4%)
#but JARID2 is lost (possible depletion? 33% vs 52%)


#### top DE ReMap TFs in 0-4hrs ####

#which are in top 10 most dynamic:

topEarlyTFs <- filter(fpkm_allGDE, EnsName %in% allG_remap2022_TFs_200$peakTF)

#strongest p val
topEarlyTFs <- topEarlyTFs[order(topEarlyTFs$preadj_0_4, decreasing = F),]
topEarlyTFs$EarlyP_rank <- 1:length(topEarlyTFs$EnsID)

#strongest max fpkm
topEarlyTFs$earlyFPKM <- Biobase::rowMax(as.matrix(topEarlyTFs[,c(20,22)]))
topEarlyTFs <- topEarlyTFs[order(topEarlyTFs$earlyFPKM, decreasing = T),]
topEarlyTFs$EarlyFPKM_rank <- 1:length(topEarlyTFs$EnsID)

#make a score:
topEarlyTFs$EarlyP_FPKM_score <- topEarlyTFs$EarlyP_rank * topEarlyTFs$EarlyFPKM_rank
topEarlyTFs <- topEarlyTFs[order(topEarlyTFs$EarlyP_FPKM_score, decreasing = F),]

#top 30 includes HMGB2 + DAXX + NFKB2 


#
#### visualisations - bar plot? jaccard similarity? ####

#ReMap TFs in background, Up PCGs, Up SCClncs:

#all TFs >5? just DE? just top DE?

#all TFs >5 to start:

#use this to find ReMap regulators which are expressed / with high potential to be active in each window:
fpkm_thresh <- 5
ReMaps_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >fpkm_thresh | Hour4_meanFPKM >fpkm_thresh))

#background, promoters per window:
expr_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)$EnsID)

#TFs expressed above the TF thresh, and with peak in an expressed gene in target set:
allG_remap2022_TFs_window <- filter(allG_remap2022_TFs_200, #EnsID %in% SCClncRNAs_04$EnsID, 
                                    #peakTF %in% c(filter(fpkm_allGDE_Upwithin_4, abs(LogFC_0_4) >1.25, preadj_0_4 < 0.0001)$EnsName, 
                                    #              filter(fpkm_allGDE_Downwithin_4, abs(LogFC_0_4) >1.25, preadj_0_4 < 0.0001)$EnsName),
                                    peakTF %in% ReMaps_04$EnsName)

TFs_inTargetProms <- unique(allG_remap2022_TFs_window$peakTF)
TFs_inTargetProms <- TFs_inTargetProms[!is.na(TFs_inTargetProms)]

#promoter regions for lncRNAs in given cluster of interest 
TargetGenes1 <- allG_remap2022_TFs_200[allG_remap2022_TFs_200$EnsID %in% SCClncRNAs_04$EnsID,]
TargetGenes2 <- allG_remap2022_TFs_200[allG_remap2022_TFs_200$EnsID %in% filter(fpkm_allGDE_Upwithin_4, 
                                                                                grepl("coding|TF|CC", GeneClassUpdate),
                                                                                EnsType == "protein_coding")$EnsID,]
TargetGenes3 <- allG_remap2022_TFs_200[allG_remap2022_TFs_200$EnsID %in% filter(fpkm_allGDE_Upwithin_4, 
                                                                                grepl("fide", GeneClassUpdate))$EnsID,]


TFs_in_4hrDEGProms_all_Fish <- list()

for(i in 1:length(TFs_inTargetProms)){
  
  a <- length(unique(filter(TargetGenes1, peakTF == TFs_inTargetProms[i])$EnsID))
  b <- length(unique(TargetGenes1$EnsID))
  ai <- length(unique(filter(TargetGenes2, peakTF == TFs_inTargetProms[i])$EnsID))
  bi <- length(unique(TargetGenes2$EnsID))
  aii <- length(unique(filter(TargetGenes3, peakTF == TFs_inTargetProms[i])$EnsID))
  bii <- length(unique(TargetGenes3$EnsID))
  c <- length(unique(filter(allG_remap2022_TFs_200, peakTF == TFs_inTargetProms[i])$EnsID))
  d <- length(unique(allG_remap2022_TFs_200$EnsID))
  
  TFs_in_4hrDEGProms_all_Fish[[i]] <- c(
    a,b,ai,bi,aii,bii,c,d)
}

names(TFs_in_4hrDEGProms_all_Fish) <- TFs_inTargetProms

TFs_in_4hrDEGProms_all_Fish_df <- as.data.frame(t(bind_rows(TFs_in_4hrDEGProms_all_Fish, .id = "TFs")))

TFs_in_4hrDEGProms_all_Fish_df$peakTF <- TFs_inTargetProms

TFs_in_4hrDEGProms_all_Fish_df$SCCLHitRate <- TFs_in_4hrDEGProms_all_Fish_df$V1/TFs_in_4hrDEGProms_all_Fish_df$V2
TFs_in_4hrDEGProms_all_Fish_df$Up4HitRate <- TFs_in_4hrDEGProms_all_Fish_df$V3/TFs_in_4hrDEGProms_all_Fish_df$V4
TFs_in_4hrDEGProms_all_Fish_df$LncHitRate <- TFs_in_4hrDEGProms_all_Fish_df$V5/TFs_in_4hrDEGProms_all_Fish_df$V6
TFs_in_4hrDEGProms_all_Fish_df$Expr4HitRate <- TFs_in_4hrDEGProms_all_Fish_df$V7/TFs_in_4hrDEGProms_all_Fish_df$V8

#don't include small numbers in this:
TFs_in_4hrDEGProms_all_Fish_df$SCCLHitRate[TFs_in_4hrDEGProms_all_Fish_df$V1 < 6] <- 0
TFs_in_4hrDEGProms_all_Fish_df$Up4HitRate[TFs_in_4hrDEGProms_all_Fish_df$V3 < 6] <- 0
TFs_in_4hrDEGProms_all_Fish_df$LncHitRate[TFs_in_4hrDEGProms_all_Fish_df$V5 < 6] <- 0

#fold increase from background
TFs_in_4hrDEGProms_all_Fish_df$SCCL_fold <- TFs_in_4hrDEGProms_all_Fish_df$SCCLHitRate/TFs_in_4hrDEGProms_all_Fish_df$Expr4HitRate
TFs_in_4hrDEGProms_all_Fish_df$Up4_fold <- TFs_in_4hrDEGProms_all_Fish_df$Up4HitRate/TFs_in_4hrDEGProms_all_Fish_df$Expr4HitRate
TFs_in_4hrDEGProms_all_Fish_df$Lnc_fold <- TFs_in_4hrDEGProms_all_Fish_df$LncHitRate/TFs_in_4hrDEGProms_all_Fish_df$Expr4HitRate

#fold increase from lnc to PCG
TFs_in_4hrDEGProms_all_Fish_df$SCCL_v_Up4 <- TFs_in_4hrDEGProms_all_Fish_df$SCCLHitRate/TFs_in_4hrDEGProms_all_Fish_df$Up4HitRate
TFs_in_4hrDEGProms_all_Fish_df$Lnc_v_Up4 <- TFs_in_4hrDEGProms_all_Fish_df$LncHitRate/TFs_in_4hrDEGProms_all_Fish_df$Up4HitRate

#plot TFs which have enrichment potential in one of UP DE PCGs, DE lncs or SCCLNCs
plot_remap_earlyUp <- filter(TFs_in_4hrDEGProms_all_Fish_df, 
                             #at least 10 hits in background
                             V7 >10, 
                             #found in at least slightly higher proportions in either sccl or de pcg
                             (SCCL_fold > 1.1) | (Up4_fold >1.1) | (Lnc_fold > 1.1))
plot_remap_earlyUp <- melt(plot_remap_earlyUp[,c(9,14:16)])

#order TFs by strongest in SCCL vs strongest in DE PCG:
orderTFs <- filter(TFs_in_4hrDEGProms_all_Fish_df, peakTF %in% plot_remap_earlyUp$peakTF)
orderTFs <- orderTFs[order(-orderTFs$Lnc_fold, orderTFs$Up4_fold),]

plot_remap_earlyUp$peakTF <- factor(plot_remap_earlyUp$peakTF)
plot_remap_earlyUp$peakTF <- factor(plot_remap_earlyUp$peakTF,
                                    levels = levels(plot_remap_earlyUp$peakTF)[
                                      match(orderTFs$peakTF, levels(plot_remap_earlyUp$peakTF))
                                    ])

ggplot(filter(plot_remap_earlyUp)) + aes(x = peakTF, y = value, fill = variable) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  ylab("Fold change promoter hits\nin early DEGs vs. expressed") +
  Seurat::RotatedAxis()

#so basically, the TFs with enrichment potential in the SCClncs and the TFs with enrichment potential in DE PCGs are v. diff
#there is basically only one SUPT16H, and it is v weak
#note that HMGB2 is only 7 TFs, DAXX only 5

#adding in lncRNAs maybe just adds confusion...

#going back a layer to just % hits in promoters (and adding in the % in background as dashed line) may help make clearer

#TBC, the HMGB2 thing is v interesting for cis lncs


#
#### takeaways ####

#seems like different TF binding pattern at scclncs

#not the best dataset, caveats
# - some of this may be driven by lack of peaks in lnc sites, more specific so not as widely captured as PCG
# - this shouldn't apply to e.g. HMGB2 though as in this case the inverse is true, found despite lower chance of profiling lncs in database

#but downsampling may be fairer way to do this (downsample the PCG peaks so similar to the lncs)

#but only 2 TFs close to significance
#  - v strongly regulated relative to other TFs, and more enriched than others in scclncs
#  - v decent odds ratios too, the p is lacking after adjustment
#  - even if being v stringent on TF inclusion

#HMGB2 lit. search done previously was of interest
#weaker now than it was before (14 before, probalby from doing all 0-4hr lncs)

