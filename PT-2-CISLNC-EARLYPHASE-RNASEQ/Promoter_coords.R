#promoter co-ordinates:

library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(rcompanion)
library(ggbeeswarm)

#import of all lncRNA-PCG pairs which are expressed in close 2D proximity or connected by the Zhao caSMC hiC loops to an expressed coding gene:
fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGClassUpdate.csv", header = T)

FPKM_CQV_OVERLAP_fpkm <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/FPKM_CQV_OVERLAP_fpkm.csv")
table(FPKM_CQV_OVERLAP_fpkm$IGV)#413 pass, 168 fail

#remove artefacts (push back to step 7?)
fpkm_allG_filt <- filter(fpkm_allG, grepl("chr", chr), !grepl("artefacts", GeneClassUpdate), !is.na(GeneClassUpdate))
#remove manual or Thresh4 fails
fpkm_allG_filt_manual <- filter(fpkm_allG_filt, 
                                !EnsID %in% filter(FPKM_CQV_OVERLAP_fpkm, IGV == "fail")$EnsID, #remove manual fails
)

fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE.csv", header = T)
fpkm_allGDE_filt <- filter(fpkm_allGDE, grepl("chr", chr), !grepl("artefacts", GeneClassUpdate), !is.na(GeneClassUpdate))
#remove manual fails
fpkm_allGDE <- filter(fpkm_allGDE_filt, EnsID %in% fpkm_allG_filt_manual$EnsID)

#clusters needed for later
#regulated within 4 hours:
fpkm_allGDE_Upwithin_4 <- filter(fpkm_allGDE, (Hour0_meanFPKM >=1 | Hour4_meanFPKM >=1) &
                                   (LogFC_0_4 >= log2(1.5) & preadj_0_4 <0.05))
fpkm_allGDE_Downwithin_4 <- filter(fpkm_allGDE, (Hour0_meanFPKM >=1 | Hour4_meanFPKM >=1) &
                                     (LogFC_0_4 < -log2(1.5) & preadj_0_4 <0.05))
#of remaining, within 8 hours:
fpkm_allGDE_Upwithin_8 <- filter(fpkm_allGDE, !EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                                 ((Hour0_meanFPKM >=1 | Hour8_meanFPKM >=1) &
                                    (LogFC_0_8 >= log2(1.5) & preadj_0_8 <0.05)
                                 )
                                 |
                                   ((Hour4_meanFPKM >=1 | Hour8_meanFPKM >=1) &
                                      (LogFC_4_8 >= log2(1.5) & preadj_4_8 <0.05)
                                   ))
fpkm_allGDE_Downwithin_8 <- filter(fpkm_allGDE, !EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                                   ((Hour0_meanFPKM >=1 | Hour8_meanFPKM >=1) &
                                      (LogFC_0_8 < -log2(1.5) & preadj_0_8 <0.05)
                                   )
                                   |
                                     ((Hour4_meanFPKM >=1 | Hour8_meanFPKM >=1) &
                                        (LogFC_4_8 < -log2(1.5) & preadj_4_8 <0.05)
                                     ))
#of remaining, within 24 hours:
fpkm_allGDE_Upwithin_24 <- filter(fpkm_allGDE, !EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                                  !EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID),
                                  ((Hour0_meanFPKM >=1 | Hour24_meanFPKM >=1) &
                                     (LogFC_0_24 >= log2(1.5) & preadj_0_24 <0.05)
                                  )
                                  |
                                    ((Hour4_meanFPKM >=1 | Hour24_meanFPKM >=1) &
                                       (LogFC_4_24 >= log2(1.5) & preadj_4_24 <0.05))
                                  |
                                    ((Hour8_meanFPKM >=1 | Hour24_meanFPKM >=1) &
                                       (LogFC_8_24 >= log2(1.5) & preadj_8_24 <0.05)
                                    ))
fpkm_allGDE_Downwithin_24 <- filter(fpkm_allGDE, !EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                                    !EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID),
                                    ((Hour0_meanFPKM >=1 | Hour24_meanFPKM >=1) &
                                       (LogFC_0_24 < -log2(1.5) & preadj_0_24 <0.05)
                                    )
                                    |
                                      ((Hour4_meanFPKM >=1 | Hour24_meanFPKM >=1) &
                                         (LogFC_4_24 < -log2(1.5) & preadj_4_24 <0.05))
                                    |
                                      ((Hour8_meanFPKM >=1 | Hour24_meanFPKM >=1) &
                                         (LogFC_8_24 < -log2(1.5) & preadj_8_24 <0.05)
                                      ))


#tidy up
fpkm_allG <- fpkm_allG_filt_manual
length(unique(fpkm_allG$EnsID))#10761
length(unique(fpkm_allGDE$EnsID))#4345

#annotate DE table with clusters
fpkm_allGDE$RegulationStart[fpkm_allGDE$EnsID %in% fpkm_allGDE_Upwithin_4$EnsID] <- "Induced <4hrs"
fpkm_allGDE$RegulationStart[fpkm_allGDE$EnsID %in% fpkm_allGDE_Downwithin_4$EnsID] <- "Repressed <4hrs"

fpkm_allGDE$RegulationStart[fpkm_allGDE$EnsID %in% fpkm_allGDE_Upwithin_8$EnsID] <- "Induced 4-8hrs"
fpkm_allGDE$RegulationStart[fpkm_allGDE$EnsID %in% fpkm_allGDE_Downwithin_8$EnsID] <- "Repressed 4-8hrs"

fpkm_allGDE$RegulationStart[fpkm_allGDE$EnsID %in% fpkm_allGDE_Upwithin_24$EnsID] <- "Induced 8-24hrs"
fpkm_allGDE$RegulationStart[fpkm_allGDE$EnsID %in% fpkm_allGDE_Downwithin_24$EnsID] <- "Repressed 8-24hrs"

#above lines seperates genes into distinct buckets
table(fpkm_allGDE$RegulationStart) #973 4hr induced, 559 8-24 hr repressed


#gene co-ordinates:
#table of all lncRNAs + CAGE sites if available + TSS based on CAGE if available:
Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime.csv", header = T)
Enhancer_lociII_DEsig_Enh <- Enhancer_lociII
length(unique(Enhancer_lociII_DEsig_Enh$EnsID))#selecting all lncs = 558 genes, if just enhancer = 77 (7/2021)
length(unique(Enhancer_lociII_DEsig_Enh$MSTRG_Tx_ID))#870 as expected

#get co-ords based on FANTOM TSS
Enhancer_lociII_DEsig_Enh$Enhancer_Coords <- paste(Enhancer_lociII_DEsig_Enh$chr, 
                                                   Enhancer_lociII_DEsig_Enh$Enh_Start, 
                                                   Enhancer_lociII_DEsig_Enh$Enh_Stop, sep = ",")

#table summary:
#info during matching 5p to CAGE cluster (CTSS)
#whether the outcome was a valid CAGE match (TIEScore >60, exon 1 not reduced <10%)
#enhancer designation (based on Genehancer and FANTOM annotation)
#differential expression during timecourse

#remove AS overlap artefacts:
Enhancer_lociII_DEsig_Enh <- filter(Enhancer_lociII_DEsig_Enh, EnsID %in% fpkm_allG$EnsID)

#option, keep only major isoforms for matching up to HiC:
#Enhancer_lociII_DEsig_Enhi <- filter(Enhancer_lociII_DEsig_Enh, MSTRG_Tx_ID %in% iso_pct_10$transcript_id)
#in isolated cases where this is missing, use the one in fpkm_allG:
#Enhancer_lociII_DEsig_Enh <- rbind(Enhancer_lociII_DEsig_Enhi, filter(Enhancer_lociII_DEsig_Enh, !EnsID %in% Enhancer_lociII_DEsig_Enhi$EnsID))

#some stats:
#428 lncs expressed
length(unique(Enhancer_lociII_DEsig_Enh$EnsID))
length(unique(Enhancer_lociII_DEsig_Enh$MSTRG_Tx_ID))#727 isoforms being assessed
#198 DE lncs over 24 hours:
length(unique(filter(Enhancer_lociII_DEsig_Enh, !is.na(DiffExprs))$EnsID))

#Get TSS from FANTOM and TSS from GENCODE in same column
Enhancer_lociII_DEsig_Enh$TSS_FANTOM_GENCODE <- Enhancer_lociII_DEsig_Enh$BestStart

#obtain TSS co-ords using these CAGE sites or just 5' limit from GENCODE/Stringtie transcripts for others
trial <- fpkm_allG
trial$Tx_start <- as.numeric(sapply(strsplit(sapply(strsplit(trial$Tx_Locus, ":"), "[[", 2), "-"), "[[", 1))
trial$Tx_stop <- as.numeric(gsub(" [+-]", "", sapply(strsplit(sapply(strsplit(trial$Tx_Locus, ":"), "[[", 2), "-"), "[[", 2)))

trial <- unique(trial[,c(2,7,61:62,10,49:50)])
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
length(unique(allGB$EnsID))#10761 apiece
length(unique(allGB$MSTRG_Tx_ID))#31058 Tx total 
#(multiple TSS per gene now)


#pairs post-HiC
AllLNC_AllPCG_2d3d <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_2d3d_Aug2025.csv")

length(unique(AllLNC_AllPCG_2d3d$EnsID))  #@250kp 381 lnc, @400kbp 391

CoRegPairs_04_48_24_extended <- filter(AllLNC_AllPCG_2d3d,
                                       #AllLNC_AllPCG_1,
                                       (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                     fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                       fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                        fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                       fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
)
#@250kbp 248 pairs
#@400kbp 364 pairs


#make granges for all gene promoters
allGB_GR <- makeGRangesFromDataFrame(allGB[,-7], 
                                     start.field = "Tx_start", 
                                     end.field = "Tx_stop", 
                                     seqnames.field = "chr", 
                                     strand.field = "str", keep.extra.columns = T)

#n.b. 2000 seems too large for HOMER, 500 in line with ISMARA (LISA differs in using ATAC/DNase to find open bits)
allGB_GR_promoters <- promoters(allGB_GR, upstream = 500, downstream = 0)

#flatten to aid next steps
allGB_GR_promoters_reduce <- reduce(allGB_GR_promoters)
length(allGB_GR_promoters_reduce)#14693 regions @2kbp, 16411 @500kbp

#prom ID for promoter co-ordinates:
allGB_GR_promoters_reduce$promID <- paste(as.character(seqnames(allGB_GR_promoters_reduce)),
                                          start(allGB_GR_promoters_reduce), end(allGB_GR_promoters_reduce), 
                                          allGB_GR_promoters_reduce$MSTRG_Tx_ID, sep = "_")

#add in gene ID per reduced prom:
remapindex <- findOverlaps(query = allGB_GR_promoters, subject = allGB_GR_promoters_reduce)

#return matched peaks from ReMap2022:
allPromRed_EnsID <- unique(data.frame("EnsID" = allGB_GR_promoters$EnsID[queryHits(remapindex)],
                                      "promID" = allGB_GR_promoters_reduce$promID[subjectHits(remapindex)]))

#one gene ID per promoter region:
trial <- split(allPromRed_EnsID, allPromRed_EnsID$promID)
triali <- lapply(trial, function(x){
  paste(x$EnsID, collapse = "-") # not ideal but other good seperators (_ or .) are taken already
})
names(triali)
allPromRed_EnsIDi <- data.frame("AllEnsID" = unlist(triali),
                                "promID" = names(triali))

#look at some dual promoter regions:
head(filter(allPromRed_EnsIDi, grepl("-", AllEnsID)), n = 50)

allGB_GR_promoters_reduce$AllEnsID <- allPromRed_EnsIDi$AllEnsID[match(allGB_GR_promoters_reduce$promID, allPromRed_EnsIDi$promID)]

#width summary:
summary(width(allGB_GR_promoters_reduce))
#mostly a v. standard 2.2kbp size or 500

#isolate a group of promoters
#easiest way is to make second col:
twoGeneProms <- filter(allPromRed_EnsIDi, grepl("-", AllEnsID))
twoGeneProms$EnsID1 <- sapply(strsplit(twoGeneProms$AllEnsID, "-"), "[[", 1)
twoGeneProms$EnsID2 <- sapply(strsplit(twoGeneProms$AllEnsID, "-"), "[[", 2)

oneGeneProms <- filter(allPromRed_EnsIDi, !grepl("-", AllEnsID))
oneGeneProms$EnsID1 <- oneGeneProms$AllEnsID
oneGeneProms$EnsID2 <- "noSecondGene"

allPromRed_EnsIDi <- rbind(oneGeneProms, twoGeneProms)

allGB_GR_promoters_reduce$EnsID1 <- allPromRed_EnsIDi$EnsID1[match(allGB_GR_promoters_reduce$promID, allPromRed_EnsIDi$promID)]
allGB_GR_promoters_reduce$EnsID2 <- allPromRed_EnsIDi$EnsID2[match(allGB_GR_promoters_reduce$promID, allPromRed_EnsIDi$promID)]


DEG_4UPproms <- allGB_GR_promoters_reduce[allGB_GR_promoters_reduce$EnsID1 %in% fpkm_allGDE_Upwithin_4$EnsID | 
                                            allGB_GR_promoters_reduce$EnsID2 %in% fpkm_allGDE_Upwithin_4$EnsID]
length(unique(DEG_4UPproms$promID))#1342 proms (1531 @500kbp)
length(unique(c(DEG_4UPproms$EnsID1, DEG_4UPproms$EnsID2)))#978 EnsIDs
#@500bp 1 additional
#@2kbp n.b. 5 additional genes, these are connected to a DEG promoter but are not themselves DE
unique(c(DEG_4UPproms$EnsID1, DEG_4UPproms$EnsID2))[!unique(c(DEG_4UPproms$EnsID1, DEG_4UPproms$EnsID2)) %in% fpkm_allGDE_Upwithin_4$EnsID]
#n.b.b. not bidirectional, strand specific prom reduction done

#bed files for e.g. AME/HOMER
DEG_4UPproms_bed <- data.frame(seqnames(DEG_4UPproms), 
                                start(DEG_4UPproms), 
                                end(DEG_4UPproms),
                               DEG_4UPproms$AllEnsID,
                                ".",
                                strand(DEG_4UPproms))

#write.table(DEG_4UPproms_bed, "DEG_4UPproms_bed.bed", quote = F, col.names = F, row.names = F, sep = "\t")

DEG_4UPproms_lncs <- DEG_4UPproms[DEG_4UPproms$EnsID1 %in% filter(fpkm_allG, grepl("Lnc|fide", GeneClassUpdate))$EnsID |
                                    DEG_4UPproms$EnsID2 %in% filter(fpkm_allG, grepl("Lnc|fide", GeneClassUpdate))$EnsID]

DEG_4UPproms_lncs_bed <- data.frame(seqnames(DEG_4UPproms_lncs), 
                               start(DEG_4UPproms_lncs), 
                               end(DEG_4UPproms_lncs),
                               DEG_4UPproms_lncs$AllEnsID,
                               ".",
                               strand(DEG_4UPproms_lncs))

#write.table(DEG_4UPproms_lncs_bed, "DEG_4UPproms_lncs_bed.bed", quote = F, col.names = F, row.names = F, sep = "\t")


#background proms
allGB_GR_promoters_reduce_bed <- data.frame(seqnames(allGB_GR_promoters_reduce), 
                               start(allGB_GR_promoters_reduce), 
                               end(allGB_GR_promoters_reduce),
                               allGB_GR_promoters_reduce$AllEnsID,
                               ".",
                               strand(allGB_GR_promoters_reduce))

#write.table(allGB_GR_promoters_reduce_bed, "allGB_GR_promoters_reduce_bed.bed", quote = F, col.names = F, row.names = F, sep = "\t")




#some verification of AME output
a <- 770
b <- 1342
c <- 23
d <- 1342

fisher.test(data.frame("hit" = c(a, b-a),
                       "Not"   = c(c, d-c)), alternative = "greater")$p * 2108

#lncRNAs
DEL_4UP <- filter(allGB, EnsID %in% filter(fpkm_allGDE_Upwithin_4, grepl("Lnc|fide", GeneClassUpdate))$EnsID)
length(unique(DEL_4UP$EnsID))#428 as expected
length(unique(DEL_4UP$MSTRG_Tx_ID))#735 as expected

DEL_4UP_GR <- makeGRangesFromDataFrame(DEL_4UP[,c(1:6)], 
                                          start.field = "Tx_start", 
                                          end.field = "Tx_stop", 
                                          seqnames.field = "chr", 
                                          strand.field = "str", keep.extra.columns = T)

earlyUp_promotersL <- promoters(DEL_4UP_GR)
earlyUp_promotersL <- reduce(earlyUp_promotersL)
earlyUp_promotersL <- data.frame(gsub("chr", "", seqnames(earlyUp_promotersL)), start(earlyUp_promotersL), end(earlyUp_promotersL))

#write.table(earlyUp_promotersL, "earlyUp_promotersL.bed", quote = F, col.names = F, row.names = F, sep = "\t")


#try compare lncs to this:
DEL_4UP_DEG_compare <- DEG_4UP_GR[!DEG_4UP_GR$MSTRG_Tx_ID %in% DEL_4UP_GR$MSTRG_Tx_ID]

earlyUp_promotersL_back <- promoters(DEL_4UP_DEG_compare)
earlyUp_promotersL_back <- reduce(earlyUp_promotersL_back)
earlyUp_promotersL_back <- data.frame(gsub("chr", "", seqnames(earlyUp_promotersL_back)), 
                                      start(earlyUp_promotersL_back), end(earlyUp_promotersL_back))

#write.table(earlyUp_promotersL_back, "earlyUp_promotersL_back.bed", quote = F, col.names = F, row.names = F, sep = "\t")


#or this, all expressed in 0-4hr:
expr_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >1 | Hour4_meanFPKM >1)$EnsID)

Expressed_4 <- filter(allGB, EnsID %in% expr_04)
length(unique(Expressed_4$EnsID))#428 as expected
length(unique(Expressed_4$MSTRG_Tx_ID))#735 as expected

Expressed_4_GR <- makeGRangesFromDataFrame(Expressed_4[,c(1:6)], 
                                       start.field = "Tx_start", 
                                       end.field = "Tx_stop", 
                                       seqnames.field = "chr", 
                                       strand.field = "str", keep.extra.columns = T)

early_promoters <- promoters(Expressed_4_GR)
early_promoters <- reduce(early_promoters)
early_promoters <- data.frame(gsub("chr", "", seqnames(early_promoters)), start(early_promoters), end(early_promoters))
early_promoters <- filter(early_promoters, !gsub..chr.......seqnames.early_promoters.. == "M")

write.table(early_promoters, "early_promoters.bed", quote = F, col.names = F, row.names = F, sep = "\t")
