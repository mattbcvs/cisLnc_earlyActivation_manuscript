#### Promoter analysis ####
library(dplyr)
library(GenomicRanges)

#HOMER

#build .bed files for promoters:
fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026filt.csv", header = T)
length(unique(fpkm_allG$EnsID))

fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE_2026filt.csv", header = T)
length(unique(fpkm_allGDE$EnsID))

fpkm_allGDE_Upwithin_4 <- filter(fpkm_allGDE, RegulationStart == "Induced <4hrs")
fpkm_allGDE_Downwithin_4 <- filter(fpkm_allGDE, RegulationStart == "Repressed <4hrs")

fpkm_allGDE_Upwithin_8 <- filter(fpkm_allGDE, RegulationStart == "Induced 4-8hrs")
fpkm_allGDE_Downwithin_8 <- filter(fpkm_allGDE, RegulationStart == "Repressed 4-8hrs")

fpkm_allGDE_Upwithin_24 <- filter(fpkm_allGDE, RegulationStart == "Induced 8-24hrs")
fpkm_allGDE_Downwithin_24 <- filter(fpkm_allGDE, RegulationStart == "Repressed 8-24hrs")

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

#build genomic ranges:
#start with just promoters...
allGB_GR <- makeGRangesFromDataFrame(allGB[,-7], 
                                     start.field = "Tx_start", 
                                     end.field = "Tx_stop", 
                                     seqnames.field = "chr", 
                                     strand.field = "str", keep.extra.columns = T)

#200bp for HOMER:
allGB_GR_promoters <- promoters(allGB_GR, upstream = 200, downstream = 20)


#bed for all 4hr up DEGs
EarlyInduced_GR_promoters <- allGB_GR_promoters[allGB_GR_promoters$EnsID %in% fpkm_allGDE_Upwithin_4$EnsID]
#flatten overlaps
EarlyInduced_GR_promoters_reduce <- GenomicRanges::reduce(EarlyInduced_GR_promoters)
#put gene IDs back:
reduceGeneIDindex <- findOverlaps(query = EarlyInduced_GR_promoters_reduce, subject = allGB_GR_promoters)

EarlyInduced_GR_promoters_reduce$peakID <- paste(as.character(seqnames(EarlyInduced_GR_promoters_reduce)),
                                                 start(EarlyInduced_GR_promoters_reduce), end(EarlyInduced_GR_promoters_reduce), 
                                                 strand(EarlyInduced_GR_promoters_reduce), sep = "_")

#return matched peaks:
reduceGeneIDindex_df <- unique(data.frame("Reduced" = EarlyInduced_GR_promoters_reduce$peakID[queryHits(reduceGeneIDindex)],
                                          "GeneID" = allGB_GR_promoters$EnsID[subjectHits(reduceGeneIDindex)]))

#build .bed:
reduceGeneIDindex_df$prom_chr <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 1)
reduceGeneIDindex_df$prom_start <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 2)
reduceGeneIDindex_df$prom_end <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 3)
reduceGeneIDindex_df$prom_str <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 4)

#concatenate multiple geneIDs per peak:
length((reduceGeneIDindex_df$Reduced))
length(unique(reduceGeneIDindex_df$Reduced))

trial <- reduceGeneIDindex_df  %>%
  group_by(Reduced) %>%
  summarise(across(c(GeneID, prom_chr, prom_start, prom_end, prom_str), ~ paste(unique(.), collapse = ", ")))

EarlyInduced_GR_promoters_reduce_bed <- trial[,c(3,4,5,1,2,6)]

#bed format
#write.table(EarlyInduced_GR_promoters_reduce_bed, "EarlyInduced_GR_promoters_reduce_bed.bed", row.names = F, col.names = F, quote = F, sep = "\t")


#bed for 4hr up PCGs (can skip for AME)
EarlyInduced_GR_promoters <- allGB_GR_promoters[allGB_GR_promoters$EnsID %in% 
                                                  filter(fpkm_allGDE_Upwithin_4, 
                                                         EnsType == "protein_coding", 
                                                         grepl("coding", GeneClassUpdate))$EnsID]

#flatten overlaps
EarlyInduced_GR_promoters_reduce <- GenomicRanges::reduce(EarlyInduced_GR_promoters)
#put gene IDs back:
reduceGeneIDindex <- findOverlaps(query = EarlyInduced_GR_promoters_reduce, subject = allGB_GR_promoters)

EarlyInduced_GR_promoters_reduce$peakID <- paste(as.character(seqnames(EarlyInduced_GR_promoters_reduce)),
                                                 start(EarlyInduced_GR_promoters_reduce), 
                                                 end(EarlyInduced_GR_promoters_reduce), 
                                                 strand(EarlyInduced_GR_promoters_reduce), sep = "_")

#return matched peaks:
reduceGeneIDindex_df <- unique(data.frame("Reduced" = EarlyInduced_GR_promoters_reduce$peakID[queryHits(reduceGeneIDindex)],
                                   "GeneID" = allGB_GR_promoters$EnsID[subjectHits(reduceGeneIDindex)]))

#build .bed:
reduceGeneIDindex_df$prom_chr <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 1)
reduceGeneIDindex_df$prom_start <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 2)
reduceGeneIDindex_df$prom_end <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 3)
reduceGeneIDindex_df$prom_str <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 4)

#concatenate multiple geneIDs per peak:
length((reduceGeneIDindex_df$Reduced))
length(unique(reduceGeneIDindex_df$Reduced))

trial <- reduceGeneIDindex_df  %>%
  group_by(Reduced) %>%
  summarise(across(c(GeneID, prom_chr, prom_start, prom_end, prom_str), ~ paste(unique(.), collapse = ", ")))

EarlyInduced_GR_promoters_reduce_bed <- trial[,c(3,4,5,1,2,6)]

#bed format
#write.table(EarlyInduced_GR_promoters_reduce_bed, "EarlyInduced_PCG_GR_promoters_reduce_bed.bed", row.names = F, col.names = F, quote = F, sep = "\t")


#bed for 4hr up CClncRNAs
EarlyInduced_GR_promoters <- allGB_GR_promoters[allGB_GR_promoters$EnsID %in% fpkm_allGDE_Upwithin_4$EnsID &
                                                  allGB_GR_promoters$EnsID %in% CoRegPairs_04_48_24_extended_naiveSame$EnsID]

#flatten overlaps
EarlyInduced_GR_promoters_reduce <- GenomicRanges::reduce(EarlyInduced_GR_promoters)
#put gene IDs back:
reduceGeneIDindex <- findOverlaps(query = EarlyInduced_GR_promoters_reduce, subject = allGB_GR_promoters)

EarlyInduced_GR_promoters_reduce$peakID <- paste(as.character(seqnames(EarlyInduced_GR_promoters_reduce)),
                                                 start(EarlyInduced_GR_promoters_reduce), end(EarlyInduced_GR_promoters_reduce), 
                                                 strand(EarlyInduced_GR_promoters_reduce), sep = "_")

#return matched peaks:
reduceGeneIDindex_df <- unique(data.frame("Reduced" = EarlyInduced_GR_promoters_reduce$peakID[queryHits(reduceGeneIDindex)],
                                          "GeneID" = allGB_GR_promoters$EnsID[subjectHits(reduceGeneIDindex)]))

#build .bed:
reduceGeneIDindex_df$prom_chr <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 1)
reduceGeneIDindex_df$prom_start <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 2)
reduceGeneIDindex_df$prom_end <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 3)
reduceGeneIDindex_df$prom_str <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 4)

#concatenate multiple geneIDs per peak:
length((reduceGeneIDindex_df$Reduced))
length(unique(reduceGeneIDindex_df$Reduced))

trial <- reduceGeneIDindex_df  %>%
  group_by(Reduced) %>%
  summarise(across(c(GeneID, prom_chr, prom_start, prom_end, prom_str), ~ paste(unique(.), collapse = ", ")))

EarlyInduced_GR_promoters_reduce_bed <- trial[,c(3,4,5,1,2,6)]

#bed format
#write.table(EarlyInduced_GR_promoters_reduce_bed, "EarlyInduced_CClnc_GR_promoters_reduce_bed.bed", row.names = F, col.names = F, quote = F, sep = "\t")


#bed for 4hr up strong CClncRNAs
SCClncRNAs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SCClncRNAs.csv")
SCClncRNAs_04 <- filter(SCClncRNAs, Lnc_Cluster == "Induced <4hrs")
EarlyInduced_GR_promoters <- allGB_GR_promoters[allGB_GR_promoters$EnsID %in% 
                                                  SCClncRNAs_04$EnsID]
#flatten overlaps
EarlyInduced_GR_promoters_reduce <- GenomicRanges::reduce(EarlyInduced_GR_promoters)
#put gene IDs back:
reduceGeneIDindex <- findOverlaps(query = EarlyInduced_GR_promoters_reduce, subject = allGB_GR_promoters)

EarlyInduced_GR_promoters_reduce$peakID <- paste(as.character(seqnames(EarlyInduced_GR_promoters_reduce)),
                                                 start(EarlyInduced_GR_promoters_reduce), end(EarlyInduced_GR_promoters_reduce), 
                                                 strand(EarlyInduced_GR_promoters_reduce), sep = "_")

#return matched peaks:
reduceGeneIDindex_df <- unique(data.frame("Reduced" = EarlyInduced_GR_promoters_reduce$peakID[queryHits(reduceGeneIDindex)],
                                          "GeneID" = allGB_GR_promoters$EnsID[subjectHits(reduceGeneIDindex)]))

#build .bed:
reduceGeneIDindex_df$prom_chr <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 1)
reduceGeneIDindex_df$prom_start <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 2)
reduceGeneIDindex_df$prom_end <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 3)
reduceGeneIDindex_df$prom_str <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 4)

#concatenate multiple geneIDs per peak:
length((reduceGeneIDindex_df$Reduced))
length(unique(reduceGeneIDindex_df$Reduced))

trial <- reduceGeneIDindex_df  %>%
  group_by(Reduced) %>%
  summarise(across(c(GeneID, prom_chr, prom_start, prom_end, prom_str), ~ paste(unique(.), collapse = ", ")))

EarlyInduced_GR_promoters_reduce_bed <- trial[,c(3,4,5,1,2,6)]

#bed format
#write.table(EarlyInduced_GR_promoters_reduce_bed, "EarlyInduced_SCClnc_GR_promoters_reduce_bed.bed", row.names = F, col.names = F, quote = F, sep = "\t")


#bed for ALL CClncRNAs
EarlyInduced_GR_promoters <- allGB_GR_promoters[allGB_GR_promoters$EnsID %in% CoRegPairs_04_48_24_extended_naiveSame$EnsID]

#flatten overlaps
EarlyInduced_GR_promoters_reduce <- GenomicRanges::reduce(EarlyInduced_GR_promoters)
#put gene IDs back:
reduceGeneIDindex <- findOverlaps(query = EarlyInduced_GR_promoters_reduce, subject = allGB_GR_promoters)

EarlyInduced_GR_promoters_reduce$peakID <- paste(as.character(seqnames(EarlyInduced_GR_promoters_reduce)),
                                                 start(EarlyInduced_GR_promoters_reduce), end(EarlyInduced_GR_promoters_reduce), 
                                                 strand(EarlyInduced_GR_promoters_reduce), sep = "_")

#return matched peaks:
reduceGeneIDindex_df <- unique(data.frame("Reduced" = EarlyInduced_GR_promoters_reduce$peakID[queryHits(reduceGeneIDindex)],
                                          "GeneID" = allGB_GR_promoters$EnsID[subjectHits(reduceGeneIDindex)]))

#build .bed:
reduceGeneIDindex_df$prom_chr <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 1)
reduceGeneIDindex_df$prom_start <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 2)
reduceGeneIDindex_df$prom_end <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 3)
reduceGeneIDindex_df$prom_str <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 4)

#concatenate multiple geneIDs per peak:
length((reduceGeneIDindex_df$Reduced))
length(unique(reduceGeneIDindex_df$Reduced))

trial <- reduceGeneIDindex_df  %>%
  group_by(Reduced) %>%
  summarise(across(c(GeneID, prom_chr, prom_start, prom_end, prom_str), ~ paste(unique(.), collapse = ", ")))

EarlyInduced_GR_promoters_reduce_bed <- trial[,c(3,4,5,1,2,6)]

#bed format
#write.table(EarlyInduced_GR_promoters_reduce_bed, "All_CClnc_GR_promoters_reduce_bed.bed", row.names = F, col.names = F, quote = F, sep = "\t")


#bed for ALLstrong CClncRNAs
SCClncRNAs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SCClncRNAs.csv")
EarlyInduced_GR_promoters <- allGB_GR_promoters[allGB_GR_promoters$EnsID %in% 
                                                  SCClncRNAs$EnsID]
#flatten overlaps
EarlyInduced_GR_promoters_reduce <- GenomicRanges::reduce(EarlyInduced_GR_promoters)
#put gene IDs back:
reduceGeneIDindex <- findOverlaps(query = EarlyInduced_GR_promoters_reduce, subject = allGB_GR_promoters)

EarlyInduced_GR_promoters_reduce$peakID <- paste(as.character(seqnames(EarlyInduced_GR_promoters_reduce)),
                                                 start(EarlyInduced_GR_promoters_reduce), end(EarlyInduced_GR_promoters_reduce), 
                                                 strand(EarlyInduced_GR_promoters_reduce), sep = "_")

#return matched peaks:
reduceGeneIDindex_df <- unique(data.frame("Reduced" = EarlyInduced_GR_promoters_reduce$peakID[queryHits(reduceGeneIDindex)],
                                          "GeneID" = allGB_GR_promoters$EnsID[subjectHits(reduceGeneIDindex)]))

#build .bed:
reduceGeneIDindex_df$prom_chr <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 1)
reduceGeneIDindex_df$prom_start <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 2)
reduceGeneIDindex_df$prom_end <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 3)
reduceGeneIDindex_df$prom_str <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 4)

#concatenate multiple geneIDs per peak:
length((reduceGeneIDindex_df$Reduced))
length(unique(reduceGeneIDindex_df$Reduced))

trial <- reduceGeneIDindex_df  %>%
  group_by(Reduced) %>%
  summarise(across(c(GeneID, prom_chr, prom_start, prom_end, prom_str), ~ paste(unique(.), collapse = ", ")))

EarlyInduced_GR_promoters_reduce_bed <- trial[,c(3,4,5,1,2,6)]

#bed format
#write.table(EarlyInduced_GR_promoters_reduce_bed, "All_SCClnc_GR_promoters_reduce_bed.bed", row.names = F, col.names = F, quote = F, sep = "\t")


#background promoters for 4hr up genes
EarlyBg_GR_promoters <- allGB_GR_promoters[allGB_GR_promoters$EnsID %in% 
                                             filter(fpkm_allG, Hour0_meanFPKM > 1 | 
                                                      Hour4_meanFPKM >1)$EnsID]

#flatten overlaps
EarlyBg_GR_promoters_reduce <- GenomicRanges::reduce(EarlyBg_GR_promoters)
EarlyBg_GR_promoters_reduce$peakID <- paste(as.character(seqnames(EarlyBg_GR_promoters_reduce)),
                                            start(EarlyBg_GR_promoters_reduce), 
                                            end(EarlyBg_GR_promoters_reduce), 
                                            strand(EarlyBg_GR_promoters_reduce), sep = "_")

#put gene IDs back:
reduceGeneIDindex <- findOverlaps(query = EarlyBg_GR_promoters_reduce, subject = allGB_GR_promoters)

#return matched peaks:
reduceGeneIDindex_df <- unique(data.frame("Reduced" = EarlyBg_GR_promoters_reduce$peakID[queryHits(reduceGeneIDindex)],
                                          "GeneID" = allGB_GR_promoters$EnsID[subjectHits(reduceGeneIDindex)]))

#build .bed:
reduceGeneIDindex_df$prom_chr <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 1)
reduceGeneIDindex_df$prom_start <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 2)
reduceGeneIDindex_df$prom_end <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 3)
reduceGeneIDindex_df$prom_str <- sapply(strsplit(reduceGeneIDindex_df$Reduced, "_"), "[[", 4)

#concatenate multiple geneIDs per peak:
length((reduceGeneIDindex_df$Reduced))
length(unique(reduceGeneIDindex_df$Reduced))

trial <- reduceGeneIDindex_df  %>%
  group_by(Reduced) %>%
  summarise(across(c(GeneID, prom_chr, prom_start, prom_end, prom_str), ~ paste(unique(.), collapse = ", ")))

EarlyExprs_GR_promoters_reduce_bed <- trial[,c(3,4,5,1,2,6)]

#bed format
#write.table(EarlyExprs_GR_promoters_reduce_bed, "EarlyExprs_GR_promoters_reduce_bed.bed", row.names = F, col.names = F, quote = F, sep = "\t")


#
#### check the output, known motifs, 0-4hr induced CClncRNAs vs PCG ####

#whole genome background to help significance along:
CCLncs_04Up_known_noBg <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/HOMER/EarlyInduced_CCL_out200_noBG/knownResults.txt")

#interested in significance level, and whether a lot of CCLs have the motif:

#there are only 69 sequences so this is just an indicative analysis:
colnames(CCLncs_04Up_known_noBg)[6] <- "No. Target Sequences"

ggplot(CCLncs_04Up_known_noBg) + aes(x = -Log.P.value, y = `No. Target Sequences`) +
  geom_point()

#tracking down the top motifs:
#TGA6 - plant TF - bZIP motif - 

#some top motifs, are they DE?
fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE_2026filt.csv", header = T)
length(unique(fpkm_allGDE$EnsID))

#or at least expressed?
fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026filt.csv", header = T)
length(unique(fpkm_allG$EnsID))

