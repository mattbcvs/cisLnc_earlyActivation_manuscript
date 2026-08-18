#### SMC biobank of Shu group, from contact with Charles Solomon ####
#repeating 11 may be needed if going from scratch

library(dplyr)
library(ggplot2)
library(GenomicRanges)

#import key tables
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


AllLNC_AllPCG_2d3d <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_2d3d_Aug2025.csv")

length(unique(AllLNC_AllPCG_2d3d$EnsID))  #391 lnc
length(unique(AllLNC_AllPCG_2d3d$EnsID.y))  #2081 lnc

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
#248 pairs
#364 @400kbp


#### All variants for all 2d/3d PCGs ####

#All variants needs to just become "lnc variants", first 2 columns of Charles table:
Shu_allVar25 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/PCGs_4_Charles_Summstats_hg38_160925.csv")

length(unique(Shu_allVar25$GeneName))#2009 of supplied list are eGene (match on ID might get more)

#variants from subsetting the tables for potential lnc targets:
targetVariants <- unique(Shu_allVar25[,1:3])
targetVariants$chr <- paste("chr", sapply(sapply(targetVariants$snps, strsplit, ":"), "[[", 1), sep = "")
targetVariants$coords <- sapply(sapply(targetVariants$snps, strsplit, ":"), "[[", 2)

SVSMC_pairedlnc <- unique(AllLNC_AllPCG_2d3d$EnsID)#391 lncs in the pairings


#isolate ranges for all lncs with a neighbouring PCG (2D/3D) - just genebody, no promoters:
allGB_LNCS <- filter(allGB, EnsID %in% filter(fpkm_allG, grepl("Lnc|fide", GeneClassUpdate))$EnsID)
length(unique(allGB_LNCS$EnsID))#428 as expected
length(unique(allGB_LNCS$MSTRG_Tx_ID))#735 as expected

allGB_LNCS_GR <- makeGRangesFromDataFrame(allGB_LNCS[,c(1:6)], 
                                          start.field = "Tx_start", 
                                          end.field = "Tx_stop", 
                                          seqnames.field = "chr", 
                                          strand.field = "str", keep.extra.columns = T)
PairedLncs_GR <- allGB_LNCS_GR[allGB_LNCS_GR$EnsID %in% SVSMC_pairedlnc]
PairedLncs_GR
#660 total lnc tx for those paired to a neighbour expressed PCG 

#check for lnc promoter/genebody overlapping variants:
targetVariants_GR <- makeGRangesFromDataFrame(targetVariants, start.field = "coords", end.field = "coords", keep.extra.columns = T)

Variantindex <- findOverlaps(query = targetVariants_GR, subject = PairedLncs_GR)
Variantoverlaps <- unique(data.frame("Variant" = targetVariants_GR$snps[queryHits(Variantindex)],
                                     "Tx" = PairedLncs_GR$MSTRG_Tx_ID[subjectHits(Variantindex)]))

length(unique(Variantoverlaps$Tx)) #554 tx
length(unique(Variantoverlaps$Variant)) #11747
table(table(Variantoverlaps$Tx)) #some tx with insane no. variants, probs long intron lncRNAs?


#variants in the lnc exon regions:
#dominant isoforms
#exon Granges (with gene id as extra column)
stringtie_gtf <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/PLAR_4timepoints/3PLAR_allgenv26_Timecourse_ff.gtf", 
                            header= F, stringsAsFactors=F)
table(stringtie_gtf$V3)
stringtie_gtf$MSTRG_Tx_ID <- sapply(stringtie_gtf$V9, function(x){
  y <- unlist(strsplit(x, "; "))
  y <- y[grepl("transcript", y)]
  return(y)
})
stringtie_gtf$MSTRG_Tx_ID <- gsub("transcript_id ", "", stringtie_gtf$MSTRG_Tx_ID)
#obtain the transcripts for paired lncs (all >1fpkm tx otherwise wouldn't be in this list)
stringtie_gtf_majorPairedLncs <- filter(stringtie_gtf, MSTRG_Tx_ID %in% PairedLncs_GR$MSTRG_Tx_ID)
stringtie_gtf_majorPairedLncs <- merge(stringtie_gtf_majorPairedLncs, fpkm_allG[,c(2,49)], by = "MSTRG_Tx_ID")
length(unique(stringtie_gtf_majorPairedLncs$EnsID))#391 genes
length(unique(stringtie_gtf_majorPairedLncs$MSTRG_Tx_ID))#660 transcripts

stringtie_gtf_majorPairedLncsExons <- filter(stringtie_gtf_majorPairedLncs, V3 == "exon")

PairedLncsExons_GR <- makeGRangesFromDataFrame(stringtie_gtf_majorPairedLncsExons[,c(2,5,6,1,11)], seqnames.field = "V1", 
                                               start.field = "V4", end.field = "V5", keep.extra.columns = T)

#overlap exons
VariantindexE <- findOverlaps(query = targetVariants_GR, subject = PairedLncsExons_GR)
VariantoverlapsE <- unique(data.frame("Variant" = targetVariants_GR$snps[queryHits(VariantindexE)],
                                      "Tx" = PairedLncsExons_GR$MSTRG_Tx_ID[subjectHits(VariantindexE)]))

length(unique(VariantoverlapsE$Tx)) #413 
length(unique(VariantoverlapsE$Variant)) #1181 
table(table(VariantoverlapsE$Tx)) #fewer but still some genes with insane no. variants


#additional variants in promoter region of lncs:
VariantindexP <- findOverlaps(query = targetVariants_GR, subject = promoters(PairedLncs_GR, upstream = 2000,
                                                                             downstream = 0))
VariantoverlapsP <- data.frame("Variant" = targetVariants_GR$snps[queryHits(VariantindexP)],
                               "Tx" = PairedLncs_GR$MSTRG_Tx_ID[subjectHits(VariantindexP)])
length(unique(VariantoverlapsP$Tx)) #485
length(unique(VariantoverlapsP$Variant)) #1314
table(table(VariantoverlapsP$Tx))


#overlap splice junctions (needs work)
PairedLncsExonLimits_GR <- makeGRangesFromDataFrame(data.frame("seqnames" = seqnames(PairedLncsExons_GR), 
                                                               "start" = c(start(PairedLncsExons_GR), end(PairedLncsExons_GR)), 
                                                               "end" = c(start(PairedLncsExons_GR), end(PairedLncsExons_GR))))
#add 50bp either way
PairedLncsSJ_GR <- flank(PairedLncsExonLimits_GR,  width = 50, both = T)

VariantindexSJ <- findOverlaps(query = targetVariants_GR, subject = PairedLncsSJ_GR)
VariantoverlapsSJ <- unique(data.frame("Variant" = targetVariants_GR$snps[queryHits(VariantindexSJ)],
                                       "Tx" = PairedLncsExons_GR$MSTRG_Tx_ID[subjectHits(VariantindexSJ)]))
length(unique(VariantoverlapsSJ$Tx)) #220
length(unique(VariantoverlapsSJ$Variant)) #428
table(table(VariantoverlapsSJ$Tx)) #fewer but one tx with many...


#overlap termination site (n.b. cannot filter by tx, gtf has no entries for non ENSG genes)
#for each exon:
Tx_3primePos <- sapply(split(PairedLncsExons_GR, PairedLncsExons_GR$MSTRG_Tx_ID), function(x){
  end(range(x))
})
Tx_3primeNeg <- sapply(split(PairedLncsExons_GR, PairedLncsExons_GR$MSTRG_Tx_ID), function(x){
  start(range(x))
})

Tx_3primePos <- Tx_3primePos[order(names(Tx_3primePos))]
Tx_3primeNeg <- Tx_3primeNeg[order(names(Tx_3primeNeg))]

trial <- unique(stringtie_gtf_majorPairedLncs[,c(1,2,8)])
trial <- trial[order(trial$MSTRG_Tx_ID),]
trial$TTS_Pos <- Tx_3primePos
trial$TTS_Neg <- Tx_3primeNeg
trial$TTS <- trial$TTS_Pos
trial$TTS[trial$V7 == "-"] <- trial$TTS_Neg[trial$V7 == "-"]
Tx_3prime <- trial

PairedLncsTTS_GR <- makeGRangesFromDataFrame(Tx_3prime, seqnames.field = "V1", strand.field = "V7",
                                             start.field = "TTS", end.field = "TTS", keep.extra.columns = T)

#add flank 2500bp downstream (area in which the HiC could connect)
#settled on 500bp, less chance of kicking out eQTLs with bad overlap
PairedLncsTTS_GR <- terminators(PairedLncsTTS_GR, downstream = 500, upstream = 0)

VariantindexTTS <- findOverlaps(query = targetVariants_GR, subject = PairedLncsTTS_GR)
VariantoverlapsTTS <- unique(data.frame("Variant" = targetVariants_GR$snps[queryHits(VariantindexTTS)],
                                        "Tx" = PairedLncsTTS_GR$MSTRG_Tx_ID[subjectHits(VariantindexTTS)]))
length(unique(VariantoverlapsTTS$Tx)) #503
length(unique(VariantoverlapsTTS$Variant)) #1893
table(table(VariantoverlapsTTS$Tx))

#remove Variant-gene pairs from exon or promoter to just get splice (and TES...)
#some variants may be splice for one gene, promoter/exon for another, so use pairs only to remove
#n.b. not a good label: refers to variant-lnc pairs not variant-eGene pairs
VariantoverlapsSJ$variantGene <- paste(VariantoverlapsSJ$Variant, VariantoverlapsSJ$Tx, sep = "-")
VariantoverlapsE$variantGene <- paste(VariantoverlapsE$Variant, VariantoverlapsE$Tx, sep = "-")
VariantoverlapsP$variantGene <- paste(VariantoverlapsP$Variant, VariantoverlapsP$Tx, sep = "-")
VariantoverlapsTTS$variantGene <- paste(VariantoverlapsTTS$Variant, VariantoverlapsTTS$Tx, sep = "-")


#remaining matches should be intron - and probs a bit non-specific:
Variantoverlaps$variantGene <- paste(Variantoverlaps$Variant, Variantoverlaps$Tx, sep = "-")
VariantoverlapsI <- filter(Variantoverlaps, !variantGene %in% c(VariantoverlapsP$variantGene, 
                                                                VariantoverlapsE$variantGene,
                                                                VariantoverlapsSJ$variantGene,
                                                                VariantoverlapsTTS$variantGene
))
length(unique(VariantoverlapsI$Tx)) #483
length(unique(VariantoverlapsI$Variant)) #10719
table(table(VariantoverlapsI$Tx))

#can subset the table by these lnc-overlapping variants:
#list of variants to take back to subset the gtex tables:
CisLncVariants <- unique(c(VariantoverlapsE$Variant, VariantoverlapsP$Variant, VariantoverlapsSJ$Variant, 
                           VariantoverlapsI$Variant, VariantoverlapsTTS$Variant))
length(unique(CisLncVariants))#14481

Shu_LncVar25 <- filter(Shu_allVar25, snps %in% CisLncVariants)
head(Shu_LncVar25)

#annotate variants
VariantoverlapsE$OverlapType <- "Exon"
VariantoverlapsP$OverlapType <- "Promoter"
VariantoverlapsSJ$OverlapType <- "Splice"
VariantoverlapsI$OverlapType <- "Intron"
VariantoverlapsTTS$OverlapType <- "TTS"

VariantoverlapsAll <- rbind(VariantoverlapsE, VariantoverlapsP, VariantoverlapsSJ, VariantoverlapsI, VariantoverlapsTTS)
length(unique(VariantoverlapsAll$Variant))#14481 as expected


#add info about multiple overlaps:
allGB_GR <- makeGRangesFromDataFrame(allGB[,c(1:6)], 
                                     start.field = "Tx_start", 
                                     end.field = "Tx_stop", 
                                     seqnames.field = "chr", 
                                     strand.field = "str", keep.extra.columns = T)
#31058 total tx to check

#check for lnc promoter/genebody/TTS overlapping variants:
#extend up from TSS by 2000kbp
allGB_GRpromoters <- promoters(allGB_GR, upstream = 2000)#2000 originally
#extend down from TTS by 2500kbp (HiC range)
allGB_GRterminators <- terminators(allGB_GR, downstream = 500)#2500 originally

#take union
allGB_GR <- GenomicRanges::punion(allGB_GR, allGB_GRpromoters)
allGB_GR <- GenomicRanges::punion(allGB_GR, allGB_GRterminators)

allGB_GR$EnsID <- allGB_GRpromoters$EnsID
allGB_GR$MSTRG_Tx_ID <- allGB_GRpromoters$MSTRG_Tx_ID

#focus only on  variants found in a lncRNA - are any of these overlapping another gene?
VariantoverlapsAll_GR <- targetVariants_GR[targetVariants_GR$snps %in% VariantoverlapsAll$Variant]

Variantindex <- findOverlaps(query = VariantoverlapsAll_GR, subject = allGB_GR)
Variantoverlaps <- unique(data.frame("Variant" = VariantoverlapsAll_GR$snps[queryHits(Variantindex)],
                                     "Tx" = allGB_GR$MSTRG_Tx_ID[subjectHits(Variantindex)]))

length(unique(Variantoverlaps$Tx)) #1226 tx
length(unique(Variantoverlaps$Variant)) #14481 variants, as expected

#add gene to a) variants overlapping lnc table:
trial <- merge(fpkm_allG[,c(2,49)], VariantoverlapsAll, by.x = "MSTRG_Tx_ID", by.y = "Tx")
#and to b) variants overlapping all genes table:
triali <- merge(fpkm_allG[,c(2,49)], Variantoverlaps, by.x = "MSTRG_Tx_ID", by.y = "Tx")

#variant-gene pairs in both table a) and b)
trial$variantGene <- paste(trial$Variant, trial$EnsID, sep = "-")
triali$variantGene <- paste(triali$Variant, triali$EnsID, sep = "-")

trial$variantTx <- paste(trial$Variant, trial$MSTRG_Tx_ID, sep = "-")
triali$variantTx <- paste(triali$Variant, triali$MSTRG_Tx_ID, sep = "-")

#ID variant-tx overlapping a lnc
#triali$lncMatch[triali$variantTx %in% trial$variantTx] <- "LncMatch"

#for each variant, check number of gene overlaps:
noGeneOverlaps <- as.data.frame(table(unique(triali[,c(2,3)])$Variant))

PCG_ID <- filter(fpkm_allG, grepl("coding|TF|CC", GeneClassUpdate))$EnsID
noPCGOverlaps <- as.data.frame(table(unique(filter(triali, EnsID %in% PCG_ID)[,c(2,3)])$Variant))
dim(noPCGOverlaps)

Other_ID <- filter(fpkm_allG, !grepl("coding|TF|CC", GeneClassUpdate))$EnsID
noOtherOverlaps <- as.data.frame(table(unique(filter(triali, EnsID %in% Other_ID)[,c(2,3)])$Variant))
dim(noOtherOverlaps)

trialii <- merge(noGeneOverlaps, noPCGOverlaps, by = "Var1", all.x = T)
trialii <- merge(trialii, noOtherOverlaps, by = "Var1")
colnames(trialii) <- c("Variant", "TotalOverlaps", "PCGOverlaps", "OtherOverlaps")
table(trialii$TotalOverlaps)
#@400kbp 62 overlap 3, 11718 overlap 1  (i.e. a lncRNA)

#check on IGV - all good (bear in mind promoter regions used too)
VariantoverlapsAll <- merge(VariantoverlapsAll, trialii, by = "Variant")

Shu_LncVar25_2 <- unique(merge(VariantoverlapsAll, Shu_LncVar25, by.x = "Variant", by.y = "snps"))
colnames(Shu_LncVar25_2)[c(2,8)] <- c("MSTRG_Tx_ID", "EnsID.y")
length(unique(Shu_LncVar25_2$Variant))
length(unique(VariantoverlapsAll$Variant))

length(unique(Shu_LncVar25_2$EnsID.y))#1538 lnc neighbours have an eQTL overlapping a lncRNA

head(Shu_LncVar25_2)

#add in lncRNA IDs
trial <- merge(fpkm_allG[,c(2,49)], Shu_LncVar25_2, by = "MSTRG_Tx_ID")
Shu_LncVar25_2 <- trial

#write.csv(Shu_LncVar25_2, "C:/users/mbennet5/Shu_LncVar25_shortTTS.csv", row.names = F)


#### import other bits for testing enrichment ####

#eqtl-egnenes that overlap lncs:
#Shu_LncVar <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Shu_LncVar25.csv")
Shu_LncVar <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Shu_LncVar25_shortTTS.csv")

#old for comparison
#Shu_LncVarOld <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Shu_LncVar.csv")

#no genes that each variant overlaps should be the same:
#colnames(Shu_LncVar)
#colnames(Shu_LncVarOld)
#trial <- unique(merge(Shu_LncVar[,c(3,6,7,8)], Shu_LncVarOld[,c(3,6,7,8)], by = "Variant"))

#trial$comp1 <- sapply(as.data.frame(t(trial[,c(2:4)])), paste, collapse = "-")
#trial$comp2 <- sapply(as.data.frame(t(trial[,c(5:7)])), paste, collapse = "-")

#failedOverlaps <- trial[!trial$comp1 == trial$comp2,]

#231 incorrect

#expression:
#From GEO for the paper - raw counts for the genes of interest in SVSMC - add full complement of genes here:
Shu_SVSMCpairedPCG_exprs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/PCGs_4_Charles_counts25.csv", header = F)

#need a colSum for all genes, done in R on eddie:
Shu_colSums <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CharlesTable_colSums.csv", header = T)
dim(Shu_colSums)
dim(Shu_SVSMCpairedPCG_exprs)

#gene lengths, would be better to have the RSEM...:
library(AnnotationHub)
library(ensembldb)
ah <- AnnotationHub()
query(ah, c("EnsDb", "v100", "Homo sapiens"))
txdb <- ah[["AH79689"]] #Ensembl 100

txsByGene=transcriptsBy(txdb,"gene")

lengthDataI <- lengthOf(txdb, of = "tx")
lengthDataI <- data.frame("tx_id" = names(lengthDataI), "tx_length" = lengthDataI)
txsByGene <- data.frame("tx_id" = unlist(txsByGene)$tx_id, "gene_id" = unlist(txsByGene)$gene_id)

trial <- merge(txsByGene, lengthDataI, by = "tx_id")
filter <- dplyr::filter
trial <- filter(trial, gene_id %in% Shu_SVSMCpairedPCG_exprs$V1)
trial <- split(trial, trial$gene_id)
lengthfiltI <- sapply(trial, function(x){
  mean(x$tx_length)
})

#matches to the table:
head(names(lengthfiltI))
head(Shu_SVSMCpairedPCG_exprs$V1)
sum(names(lengthfiltI) == Shu_SVSMCpairedPCG_exprs$V1)

#FPK - read counts by gene lengths
trial <- sapply(Shu_SVSMCpairedPCG_exprs[,-1], function(x){
  x/(lengthfiltI/1000)
})

#TPM - divide FPK by colSum/1000000
triali <- trial

for (i in 1:dim(trial)[2]){
  triali[,i] <- trial[,i]/(Shu_colSums[i,2]/1000000)
}

dim(triali)

#quick compare to a random RSEM sample, extract from isoform table:
#transcript_id  	  gene_id	            length	effective_length	expected_count	TPM	  FPKM	IsoPct
#ENST00000570567.1	ENSG00000006283.17	1255	  1011.41	          2.00	          0.03	0.08	100.00
#ENST00000005995.7	ENSG00000007038.10	1101    857.41	          0.00	          0.00	0.00	0.00
#ENST00000442282.1	ENSG00000011590.13	691	    448.06	          4.00	          0.11	0.35	52.56
#ENST00000218230.5	ENSG00000102109.8	  1050	  806.41	          26	            0.41	1.25	100

#sum of fpk counts for this sample: 78948650.1

#works:
(2/1.011)/(78948650.1/1000000)
(4/0.44806)/(78948650.1/1000000)
(26/0.806)/(78948650.1/1000000)

Shu_SVSMCpairedPCG_exprsTPM <- data.frame("EnsID_merge" = Shu_SVSMCpairedPCG_exprs$V1, triali)
Shu_SVSMCpairedPCG_exprs[1:5,1:5]
Shu_SVSMCpairedPCG_exprsTPM[1:5,1:5]
head(Shu_colSums)

lengthfiltI[Shu_SVSMCpairedPCG_exprs$V1[1]]

775/0.6397143/(30360137/1000000)

#may not remove anything other than v. low TPM
sum(rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))>1)/dim(Shu_SVSMCpairedPCG_exprsTPM)[1]#-5%
sum(rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))>5)/dim(Shu_SVSMCpairedPCG_exprsTPM)[1]#-18%
sum(rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))>10)/dim(Shu_SVSMCpairedPCG_exprsTPM)[1]#-32%
sum(rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))>20)/dim(Shu_SVSMCpairedPCG_exprsTPM)[1]#-52%
sum(rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))>30)/dim(Shu_SVSMCpairedPCG_exprsTPM)[1]#-65%
sum(rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))>50)/dim(Shu_SVSMCpairedPCG_exprsTPM)[1]#-77%

#identify eGenes at given p val
#original requested genes
Shu_allVar_p <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/PCGs_4_Charles_Summstats_hg38_160925.csv")
length(unique(Shu_allVar_p$gene)) #2009 eGenes from Shu/Charles

#note very diff p val distribution in this resource compared to GTEX
pThresh <- 0.05    # -0%
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < pThresh)
length(unique(Shu_allVar_pThresh$snps))/length(unique(Shu_allVar_p$snps))

pThresh <- 0.025    # -22%
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < pThresh)
length(unique(Shu_allVar_pThresh$snps))/length(unique(Shu_allVar_p$snps))

pThresh <- 0.01    # -40%
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < pThresh)
length(unique(Shu_allVar_pThresh$snps))/length(unique(Shu_allVar_p$snps))

pThresh <- 0.005   # -55%
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < pThresh)
length(unique(Shu_allVar_pThresh$snps))/length(unique(Shu_allVar_p$snps))

pThresh <- 0.001   # -62%
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < pThresh)
length(unique(Shu_allVar_pThresh$snps))/length(unique(Shu_allVar_p$snps))

pThresh <- 0.0001  # -75%
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < pThresh)
length(unique(Shu_allVar_pThresh$snps))/length(unique(Shu_allVar_p$snps))

pThresh <- 0.00001 # -75%
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < pThresh)
length(unique(Shu_allVar_pThresh$snps))/length(unique(Shu_allVar_p$snps))

pThresh <- 0.000001 # -85%
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < pThresh)
length(unique(Shu_allVar_pThresh$snps))/length(unique(Shu_allVar_p$snps))


#### start optimisation - is there a sweet spot where eQTLp predicts correct cisLnc pairs ####

#ready to run tests as done on GTEX:
Shu_fish_list <- list()

#allow finding pairs in each:
Shu_LncVar$pairs <- paste(Shu_LncVar$EnsID, Shu_LncVar$EnsID.y, sep = "-")

AllLNC_AllPCG_2d3d$pairs_merge <- gsub("\\.[0-9]*$", "", AllLNC_AllPCG_2d3d$pairs)
AllLNC_AllPCG_2d3d$EnsID_merge <- gsub("\\.[0-9]*", "", AllLNC_AllPCG_2d3d$EnsID.y)

#matching more or less to amount of variants removed by GTEX thresholds
pThresh_df <- data.frame("pThresh" = c(rep(0.05, 1),
                                       rep(0.025, 1),# -22%
                                       rep(0.01, 1), # -40%
                                       rep(0.001, 1),# -55%
                                       rep(0.0001, 1),# -62%
                                       rep(0.00001, 1))# -75%
                         )

#promoter:
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType == "Promoter", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
    )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_df_prom <- Shu_fish_list_df
Shu_fish_list_df_prom$eQTLoverlaps <- "Promoter"

#write.csv(Shu_fish_list_df_prom, "Shu_fish_list_df_prom.csv", row.names = F)

#somehow 2x less pairs found here at 1e-5, (15 now vs. 17 before)
#how would this happen? not excluding any lncs, PCGs, eGenes or eQTLs from before... just adding more
#no issue from including "fake" PCGs from PLAR (one found though, not issue as not DE)
#overlap code likely culprit: adding in the TTS likely means more overlaps have been found
#some argument to remove TTS, no sig from GTEX
#makes sense to keep for stringency, given the nature of lnc work

#modifications have improved the base unadjusted p value for these runs


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType == "Exon", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_df_exon <- Shu_fish_list_df
Shu_fish_list_df_exon$eQTLoverlaps <- "Exon"

#write.csv(Shu_fish_list_df_exon, "Shu_fish_list_df_exon.csv", row.names = F)


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType %in% c("Promoter", "Exon"), 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_df_promExon <- Shu_fish_list_df
Shu_fish_list_df_promExon$eQTLoverlaps <- "PromoterExon"

#write.csv(Shu_fish_list_df_promExon, "Shu_fish_list_df_promExon.csv", row.names = F)


#intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType == "Intron", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_df_intron <- Shu_fish_list_df
Shu_fish_list_df_intron$eQTLoverlaps <- "Intron"

#write.csv(Shu_fish_list_df_intron, "Shu_fish_list_df_intron.csv", row.names = F)


#locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              #OverlapType == "Intron", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_df_locus <- Shu_fish_list_df
Shu_fish_list_df_locus$eQTLoverlaps <- "Locus"

#write.csv(Shu_fish_list_df_locus, "Shu_fish_list_df_locus.csv", row.names = F)


#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType == "Splice", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_df_splice <- Shu_fish_list_df
Shu_fish_list_df_splice$eQTLoverlaps <- "Splice"

#write.csv(Shu_fish_list_df_splice, "Shu_fish_list_df_splice.csv", row.names = F)


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType == "TTS", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_df_TTS <- Shu_fish_list_df
Shu_fish_list_df_TTS$eQTLoverlaps <- "TTS"

write.csv(Shu_fish_list_df_TTS, "Shu_fish_list_df_TTS.csv", row.names = F)

#promoter only (strong), some signs of something in splice, nothing else even close


#### round ii optimisation - is there a sweet spot where eQTLp predicts correct SAME TIMEFRAME cisLnc pairs ####

CoRegPairs_04_48_24_extendedSame <- filter(AllLNC_AllPCG_2d3d,
                                           #AllLNC_AllPCG_1,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                         fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                          fpkm_allGDE_Downwithin_4$EnsID)) |
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                           fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                            fpkm_allGDE_Downwithin_8$EnsID)) |
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                           fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                             fpkm_allGDE_Downwithin_24$EnsID))
)

#ready to run tests as done on GTEX:
Shu_fish_list <- list()

#allow finding pairs in each:
Shu_LncVar$pairs <- paste(Shu_LncVar$EnsID, Shu_LncVar$EnsID.y, sep = "-")

AllLNC_AllPCG_2d3d$pairs_merge <- gsub("\\.[0-9]*$", "", AllLNC_AllPCG_2d3d$pairs)
AllLNC_AllPCG_2d3d$EnsID_merge <- gsub("\\.[0-9]*", "", AllLNC_AllPCG_2d3d$EnsID.y)

#matching more or less to amount of variants removed by GTEX thresholds
pThresh_df <- data.frame("pThresh" = c(rep(0.05, 1),
                                       rep(0.025, 1),# -22%
                                       rep(0.01, 1), # -40%
                                       rep(0.001, 1),# -55%
                                       rep(0.0001, 1),# -62%
                                       rep(0.00001, 1))# -75%
)

#promoter:
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType == "Promoter", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedSame$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedSame$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_same_dfprom <- Shu_fish_list_df
Shu_fish_list_same_dfprom$eQTLoverlaps <- "Promoter"

#write.csv(Shu_fish_list_same_dfprom, "Shu_fish_list_same_dfprom.csv", row.names = F)

#somehow 2x less pairs found here at 1e-5, (15 now vs. 17 before)
#how would this happen? not excluding any lncs, PCGs, eGenes or eQTLs from before... just adding more
#no issue from including "fake" PCGs from PLAR (one found though, not issue as not DE)
#overlap code likely culprit: adding in the TTS likely means more overlaps have been found
#some argument to remove TTS, no sig from GTEX
#makes sense to keep for stringency, given the nature of lnc work

#modifications have improved the base unadjusted p value for these runs


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType == "Exon", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedSame$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedSame$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_same_dfexon <- Shu_fish_list_df
Shu_fish_list_same_dfexon$eQTLoverlaps <- "Exon"

#write.csv(Shu_fish_list_same_dfexon, "Shu_fish_list_same_dfexon.csv", row.names = F)


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType %in% c("Promoter", "Exon"), 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedSame$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedSame$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_same_dfpromExon <- Shu_fish_list_df
Shu_fish_list_same_dfpromExon$eQTLoverlaps <- "PromoterExon"

#write.csv(Shu_fish_list_same_dfpromExon, "Shu_fish_list_same_dfpromExon.csv", row.names = F)


#intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType == "Intron", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedSame$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedSame$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_same_dfintron <- Shu_fish_list_df
Shu_fish_list_same_dfintron$eQTLoverlaps <- "Intron"

#write.csv(Shu_fish_list_same_dfintron, "Shu_fish_list_same_dfintron.csv", row.names = F)


#locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              #OverlapType == "Intron", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedSame$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedSame$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_same_dflocus <- Shu_fish_list_df
Shu_fish_list_same_dflocus$eQTLoverlaps <- "Locus"

#write.csv(Shu_fish_list_same_dflocus, "Shu_fish_list_same_dflocus.csv", row.names = F)


#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType == "Splice", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedSame$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedSame$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_same_dfsplice <- Shu_fish_list_df
Shu_fish_list_same_dfsplice$eQTLoverlaps <- "Splice"

#write.csv(Shu_fish_list_same_dfsplice, "Shu_fish_list_same_dfsplice.csv", row.names = F)


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType == "TTS", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedSame$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedSame$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_same_dfTTS <- Shu_fish_list_df
Shu_fish_list_same_dfTTS$eQTLoverlaps <- "TTS"

#write.csv(Shu_fish_list_same_dfTTS, "Shu_fish_list_same_dfTTS.csv", row.names = F)

#none sig


#### Higher TPM optimisation - is there a sweet spot where eQTLp predicts correct cisLnc pairs, using high TPM genes only? ####

#ready to run tests as done on GTEX:
Shu_fish_list <- list()

#allow finding pairs in each:
Shu_LncVar$pairs <- paste(Shu_LncVar$EnsID, Shu_LncVar$EnsID.y, sep = "-")

AllLNC_AllPCG_2d3d$pairs_merge <- gsub("\\.[0-9]*$", "", AllLNC_AllPCG_2d3d$pairs)
AllLNC_AllPCG_2d3d$EnsID_merge <- gsub("\\.[0-9]*", "", AllLNC_AllPCG_2d3d$EnsID.y)

#matching more or less to amount of variants removed by GTEX thresholds
pThresh_df <- data.frame("pThresh" = c(rep(0.05, 1),
                                       rep(0.025, 1),# -22%
                                       rep(0.01, 1), # -40%
                                       rep(0.001, 1),# -55%
                                       rep(0.0001, 1),# -62%
                                       rep(0.00001, 1))# -75%
)


#promoter:
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 50]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType == "Promoter", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_df_prom50 <- Shu_fish_list_df
Shu_fish_list_df_prom50$eQTLoverlaps <- "Promoter"

#worse result


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 50]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType == "Exon", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_df_exon <- Shu_fish_list_df
Shu_fish_list_df_exon$eQTLoverlaps <- "Exon"

#slightly more sig


#intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 50]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType == "Intron", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_df_intron <- Shu_fish_list_df
Shu_fish_list_df_intron$eQTLoverlaps <- "Intron"

#slight more signal


#locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 50]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              #OverlapType == "Intron", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_df_locus <- Shu_fish_list_df
Shu_fish_list_df_locus$eQTLoverlaps <- "Locus"

#slight more signal


#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 50]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType == "Splice", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_df_splice <- Shu_fish_list_df
Shu_fish_list_df_splice$eQTLoverlaps <- "Splice"

#worse


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 50]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType == "TTS", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_df_TTS <- Shu_fish_list_df
Shu_fish_list_df_TTS$eQTLoverlaps <- "TTS"


#no improvement at 30 or 50 TPM


#### Shorter dist optimisation - is there a sweet spot where eQTLp predicts correct cisLnc pairs, using high TPM genes only? ####

#ready to run tests as done on GTEX:
Shu_fish_list <- list()

#allow finding pairs in each:
Shu_LncVar$pairs <- paste(Shu_LncVar$EnsID, Shu_LncVar$EnsID.y, sep = "-")

AllLNC_AllPCG_2d3d$pairs_merge <- gsub("\\.[0-9]*$", "", AllLNC_AllPCG_2d3d$pairs)
AllLNC_AllPCG_2d3d$EnsID_merge <- gsub("\\.[0-9]*", "", AllLNC_AllPCG_2d3d$EnsID.y)

#matching more or less to amount of variants removed by GTEX thresholds
pThresh_df <- data.frame("pThresh" = c(rep(0.05, 1),
                                       rep(0.025, 1),# -22%
                                       rep(0.01, 1), # -40%
                                       rep(0.001, 1),# -55%
                                       rep(0.0001, 1),# -62%
                                       rep(0.00001, 1))# -75%
)


#promoter:
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 50]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d, DisLnc_PCG < 250, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < selectedp,
                                                              OverlapType == "Promoter", 
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  #all potential targets for potential cisLncs
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
  DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
  
  a <- dim(DELNC_DEPCG_1_eQTL)[1]
  b <- dim(DELNC_DEPCG_1)[1]
  c <- dim(DELNC_PCG_1_eQTL)[1]
  d <- dim(DELNC_PCG_1)[1]
  
  Shu_fish_list[[j]] <- c(
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
    fisher.test(data.frame("cisLnc" = c(a, b-a),
                           "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
    a,b,c,d,
    mean(statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d),
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d), 
         statmod::power.fisher.test(a/b, c/d, c, d))
  )
}

length(Shu_fish_list)
names(Shu_fish_list) <- paste("Run", pThresh_df[,1], sep = "_")

trial <- lapply(Shu_fish_list, as.data.frame)
trial <- lapply(trial, t)
trial <- lapply(trial, as.data.frame)

Shu_fish_list_df <- bind_rows(trial, .id = "Run")

Shu_fish_list_df$ThreshP <- as.numeric(sapply(sapply(Shu_fish_list_df$Run, strsplit, "_"), "[[", 2))
Shu_fish_list_df$p_adj <- p.adjust(Shu_fish_list_df$V1, method = "BH")

Shu_fish_list_df_promShort <- Shu_fish_list_df
Shu_fish_list_df_promShort$eQTLoverlaps <- "Promoter"

#worse result


#### visualise eQTL-supported loci - promoter eQTL matches, strength of p, number/position of eQTLs relative to lncRNA TSS etc ####

Shu_fish_list_df_prom <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Shu_SMC_Biobank/Shu_fish_list_df_prom.csv")

SMCBank_best <- filter(Shu_fish_list_df_prom, V5 >10)
SMCBank_best$p_adj <- p.adjust(SMCBank_best$V1, method = "BH")
SMCBank_best <- filter(SMCBank_best, V1 == min(V1))[1,]

Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < 0.0001)

#### make plot of all eQTL-linked pairs per test ####

#allow finding pairs in each:
Shu_LncVar$pairs <- paste(Shu_LncVar$EnsID, Shu_LncVar$EnsID.y, sep = "-")

AllLNC_AllPCG_2d3d$pairs_merge <- gsub("\\.[0-9]*$", "", AllLNC_AllPCG_2d3d$pairs)
AllLNC_AllPCG_2d3d$EnsID_merge <- gsub("\\.[0-9]*", "", AllLNC_AllPCG_2d3d$EnsID.y)

CoRegPairs_04_48_24_extended$pairs_merge <- gsub("\\.[0-9]*$", "", CoRegPairs_04_48_24_extended$pairs)

triali <- filter(AllLNC_AllPCG_2d3d, EnsID %in% CoRegPairs_04_48_24_extended$EnsID,
                 EnsID_merge %in% Shu_exprsG, 
                 EnsName.y %in% Shu_allVar_pThresh$GeneName)

triali$eQTL_validated_tissue <- "No"
triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < 0.0001,
                                                            OverlapType == "Promoter", 
                                                            TotalOverlaps == 1)$pairs] <- "Yes"
table(triali$eQTL_validated_tissue)

ShuBioBankSupported_AllVar <- filter(Shu_LncVar,  
                                     pairs %in% triali$pairs_merge,
                                     #pvalue < 0.00001,
                                     OverlapType == "Promoter", 
                                     TotalOverlaps == 1)

ShuBioBankSupported_AllVar <- unique(ShuBioBankSupported_AllVar[,c(2,9,3,13,14,19)])

#include the other testable lnc-PCG pairs:
trialii <- filter(triali, !pairs_merge %in% ShuBioBankSupported_AllVar$pairs)

trialii <- data.frame("EnsID" = trialii$EnsID,
                      "EnsID.y" = trialii$EnsID.y,
                      "Variant" = NA,
                      "pvalue" = NA,
                      "FDR" = NA,
                      "pairs"= trialii$pairs)

trialii <- rbind(ShuBioBankSupported_AllVar, trialii)

ShuBioBankSupported_AllVar <- trialii
length(unique(triali$pairs))
length(unique(ShuBioBankSupported_AllVar$pairs))

#order pairs by max eQTL sig.
trial <- sapply(split(ShuBioBankSupported_AllVar, ShuBioBankSupported_AllVar$pairs), function(x){
  min(x$pvalue)
})

bestEqtl <- data.frame("pairs" = names(trial), "eQTLbest" = trial)
bestEqtl <- bestEqtl[order(bestEqtl$pairs),]

ShuBioBankSupported_AllVar$PairType <- "Other CClncRNA-Neighbour pairs"
ShuBioBankSupported_AllVar$PairType[ShuBioBankSupported_AllVar$pairs %in% CoRegPairs_04_48_24_extended$pairs_merge] <- "CClncRNA-Target pairs"

ShuBioBankSupported_AllVar <- ShuBioBankSupported_AllVar[order(ShuBioBankSupported_AllVar$pairs),]
ShuBioBankSupported_AllVar$pairs <- as.factor(ShuBioBankSupported_AllVar$pairs)
ShuBioBankSupported_AllVar$pairs <- factor(ShuBioBankSupported_AllVar$pairs, 
                                           levels = levels(ShuBioBankSupported_AllVar$pairs)[order(bestEqtl$eQTLbest)])

#max eQTL p simpler to visualise:
bestEqtl$eQTLbest[is.na(bestEqtl$eQTLbest)] <- 1
ShuBioBankSupported_BestVar <- unique(merge(ShuBioBankSupported_AllVar[,-c(3:5)], bestEqtl, by = "pairs", all.x = T))

selectedp <- SMCBank_best$ThreshP

ggplot(ShuBioBankSupported_BestVar) + aes(y = -log10(eQTLbest), x = pairs, fill = PairType, color = PairType) +
  geom_bar(stat = "identity") +
  annotate("text", size = 6, x = SMCBank_best$V6*0.25, y = -log10(selectedp)*2.1, label = paste("eQTL p=", selectedp), fontface = 3, color = "grey60") +
  annotate("text", size = 6, x = SMCBank_best$V6*0.25, y = -log10(selectedp)*3.65, label = paste("Enrichment: p=", formatC(SMCBank_best$p_adj, format = "e", digits = 1),
                                                                                               "OR=", round(SMCBank_best$V2,1)), fontface = 3, color = "grey60") +
  geom_hline(yintercept = -log10(selectedp), linetype = "dashed") +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c(`CClncRNA-Target pairs` = "mediumorchid1",
                                `Other CClncRNA-Neighbour pairs` = "grey50")) +
  scale_fill_manual(values = c(`CClncRNA-Target pairs` = "mediumorchid1",
                               `Other CClncRNA-Neighbour pairs` = "grey50")) +
  coord_cartesian(ylim = c(0,75)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), legend.position = "none", panel.grid.major =  element_blank(),
        text = element_text(size=21)) +
  xlab(paste(SMCBank_best$V6, "CClncRNA-Neighbour pairs with\neGene in SMC BioBank", sep = " ")) +
  ylab("eQTL significance\n(maximum -log10p)")

#### other plots ####

triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
                 EnsName.y %in% Shu_allVar_pThresh$GeneName)

triali$eQTL_validated_tissue <- "No"
triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < 0.0001,
                                                            OverlapType == "Promoter", 
                                                            TotalOverlaps == 1)$pairs] <- "Yes"

DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")

ShuBioBankSupported <- DELNC_DEPCG_1_eQTL


#previous optimisations with TPM yielded diff results, 
#but not now:
#Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 50]
#Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < 0.05)

#triali <- filter(AllLNC_AllPCG_2d3d, EnsID_merge %in% Shu_exprsG, 
#                 EnsName.y %in% Shu_allVar_pThresh$GeneName)

#triali$eQTL_validated_tissue <- "No"
#triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar, pvalue < 0.05,
#                                                            OverlapType == "Promoter", 
#                                                            TotalOverlaps == 1
#                                                            )$pairs] <- "Yes"

#all potential targets for potential cisLncs
#DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
#DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")

#ShuBioBankSupported50 <- DELNC_DEPCG_1_eQTL
#only 2x lncRNAs from this are kept

#with a more nuanced approach, using more factors (e.g. eQTL no., eQTL distance, eQTL beta) would come more predictive value
#probs neg binomial too with shuffled pairs (FANTOM, which does smth similar with co-expression)


#number and significance of linked variants used to link up these pairs:
ShuBioBankSupported_Var <- filter(Shu_LncVar, pairs %in% ShuBioBankSupported$pairs_merge, 
       pvalue < 0.00001,
       OverlapType == "Promoter", 
       TotalOverlaps == 1)

ShuBioBankSupported_Var <- unique(ShuBioBankSupported_Var[,c(3,10,13,14,19)])

#generally more than one eQTL supports (12/15)
table(table(ShuBioBankSupported_Var$pairs))

sapply(split(ShuBioBankSupported_Var, ShuBioBankSupported_Var$pairs), function(x){
  summary(x$pvalue)
})


#n.b. SHOULD ALSO SHOW THE CO-REG WITHOUT GOOD EQTL?
#compare to all promoter-linked variants linked to all neighbours:
ShuBioBankSupported_AllVar <- filter(Shu_LncVar, 
                                     EnsID %in% ShuBioBankSupported$EnsID, 
                                     pairs %in% AllLNC_AllPCG_2d3d$pairs_merge,
                                  #pvalue < 0.00001,
                                  OverlapType == "Promoter", 
                                  TotalOverlaps == 1)

ShuBioBankSupported_AllVar <- unique(ShuBioBankSupported_AllVar[,c(2,9,3,10,13,14,19)])


#include any non-eGene expressed pairs? not used in the calculation but may be good to be complete:
trial <- filter(AllLNC_AllPCG_2d3d, EnsID %in% ShuBioBankSupported$EnsID, 
                !pairs_merge %in% ShuBioBankSupported_AllVar$pairs)

triali <- data.frame("EnsID" = trial$EnsID,
                    "EnsID.y" = trial$EnsID_merge,
                    "Variant" = NA,
                    "GeneName" = trial$EnsName.y,
                    "pvalue" = NA,
                    "FDR" = NA,
                    "pairs"= trial$pairs_merge)

trialii <- rbind(ShuBioBankSupported_AllVar, triali)

ShuBioBankSupported_AllVar <- trialii

#annotate pairs: DE/non-DE eQTL/non-eQTL/not eGene
CoRegPairs_04_48_24_extended$pairs_merge <- gsub("\\.[0-9]*$", "", CoRegPairs_04_48_24_extended$pairs)
#needed to find non-coreg eQTL linked:
ShuBioBankSupported_AllVar2 <- filter(Shu_LncVar, EnsID %in% ShuBioBankSupported$EnsID, 
                                     pvalue < 0.00001,
                                     OverlapType == "Promoter", 
                                     TotalOverlaps == 1)

ShuBioBankSupported_AllVar$PairType <- "Not CClncRNA Target"
#ShuBioBankSupported_AllVar$PairType[ShuBioBankSupported_AllVar$pairs %in% ShuBioBankSupported_AllVar2$pairs] <- "eQTL-linked but not CClncRNA Target"

ShuBioBankSupported_AllVar$PairType[ShuBioBankSupported_AllVar$pairs %in% CoRegPairs_04_48_24_extended$pairs_merge] <- "CClncRNA Target"
#ShuBioBankSupported_AllVar$PairType[ShuBioBankSupported_AllVar$pairs %in% CoRegPairs_04_48_24_extended$pairs_merge &
#                                      ShuBioBankSupported_AllVar$pairs %in% ShuBioBankSupported_Var$pairs] <- "eQTL-linked CClncRNA Target"
table(ShuBioBankSupported_AllVar$PairType)

Shu_BB_EQTL_LNCS <- unique(ShuBioBankSupported_AllVar$EnsID)

Shu_BB_EQTL_LNCS_PLOT <- list()

#max value for plot
summary(-log10(ShuBioBankSupported_AllVar$pvalue))

#try arrange in TSS order:
allGB$EnsID_merge <- gsub("\\.[0-9]*", "", allGB$EnsID)

trial <- allGB
trial$TSS <- trial$Tx_start
trial$TSS[trial$str == "-"] <- trial$Tx_stop[trial$str == "-"]

trial <- split(trial, trial$EnsID_merge)

trial <- sapply(trial, function(x){
  mean(x$TSS)
})

trial <- data.frame("EnsID_merge" = names(trial), "meanTSS" = trial)

trial <- merge(ShuBioBankSupported_AllVar, trial, by.x = "EnsID.y", by.y = "EnsID_merge")
trial <- merge(trial, unique(allGB[,c(3,8)]), by.x = "EnsID.y", by.y = "EnsID_merge")
ShuBioBankSupported_AllVar <- trial

#label closest to the lncRNA:
trial <- filter(AllLNC_AllPCG_2d3d, EnsID %in% ShuBioBankSupported$EnsID)
trial <- split(trial, trial$EnsID)

trial <- lapply(trial, function(x){
  x[x$DisLnc_PCG == min(x$DisLnc_PCG), c(13,17)]
})
trial <- bind_rows(trial)

ShuBioBankSupported_AllVar$closest2Lnc <- "Further"
ShuBioBankSupported_AllVar$closest2Lnc[ShuBioBankSupported_AllVar$pairs %in% trial$pairs] <- "Closest"

ShuBioBankSupported_AllVar$GeneName.y <- ShuBioBankSupported_AllVar$GeneName
ShuBioBankSupported_AllVar$GeneName.y[ShuBioBankSupported_AllVar$closest2Lnc == "Closest"] <- paste0(ShuBioBankSupported_AllVar$GeneName.y[ShuBioBankSupported_AllVar$closest2Lnc == "Closest"], "*")

#label lncRNA:
trial <- merge(ShuBioBankSupported_AllVar, unique(fpkm_allG[,c(2,5)]), by = "EnsID")

trial$LncName <- trial$EnsName
trial$LncName[is.na(trial$LncName)] <- trial$EnsID[is.na(trial$LncName)]

ShuBioBankSupported_AllVar <- trial

for (i in 1:length(Shu_BB_EQTL_LNCS)){
  
  selectedLocus <- filter(ShuBioBankSupported_AllVar, EnsID == Shu_BB_EQTL_LNCS[i])
  
  selectedLocus <- selectedLocus[order(selectedLocus$GeneName.y),]
  selectedLocus$GeneName.y <- as.factor(selectedLocus$GeneName.y)
  selectedLocus$GeneName.y <- factor(selectedLocus$GeneName.y, 
                                   levels = levels(selectedLocus$GeneName.y)[
                                     order(unique(selectedLocus[,c(12,9)])$meanTSS, decreasing = F)]
                                   )
  
  Shu_BB_EQTL_LNCS_PLOT[[i]] <- ggplot(selectedLocus) + 
    aes(x = GeneName.y, y = -log10(pvalue), color = PairType) +
    geom_hline(yintercept = -log10(1e-4), linetype = "dashed", colour = "grey60") +
    #geom_boxplot(outlier.shape = NA, width = 0.3) +
    #scale_y_continuous(limits = c(0,55)) +
    #scale_color_manual(values = c(`CClncRNA Target` = "mediumorchid1",
    #                              `eQTL-linked CClncRNA Target` = "mediumorchid4",
    #                              `Not CClncRNA Target` = "darkblue",
    #                              `eQTL-linked but not CClncRNA Target` = "steelblue")) +
    scale_color_manual(values = c(`CClncRNA Target` = "mediumorchid1",
                                  `Not CClncRNA Target` = "grey50")) +
    #coord_cartesian(ylim = c(0,20)) +
    geom_point(width = 0.1, size =  2.25,alpha = 0.75) +
    xlab("") +
    ylab("") +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(size=12.75)) +
    Seurat::RotatedAxis() +
    ggtitle(selectedLocus$LncName)
}

Shu_BB_EQTL_LNCS_PLOT[[1]]
Shu_BB_EQTL_LNCS_PLOT[[2]]
Shu_BB_EQTL_LNCS_PLOT[[3]]
Shu_BB_EQTL_LNCS_PLOT[[4]]
Shu_BB_EQTL_LNCS_PLOT[[5]]
Shu_BB_EQTL_LNCS_PLOT[[6]]
Shu_BB_EQTL_LNCS_PLOT[[7]]
Shu_BB_EQTL_LNCS_PLOT[[8]]
Shu_BB_EQTL_LNCS_PLOT[[9]]
Shu_BB_EQTL_LNCS_PLOT[[10]]#+theme(legend.position = "left")
Shu_BB_EQTL_LNCS_PLOT[[11]]#+theme(legend.position = "left")
Shu_BB_EQTL_LNCS_PLOT[[12]]#+theme(legend.position = "left")
Shu_BB_EQTL_LNCS_PLOT[[13]]#+theme(legend.position = "left")
Shu_BB_EQTL_LNCS_PLOT[[14]]

library(grid)
library(gridExtra)
grid.arrange(#top = textGrob("CClncRNA-targets supported by SMC Biobank eQTLs\nin lncRNA promoters\n", gp = gpar(fontface = 'bold', fontsize = 18)),
             bottom = textGrob("CClncRNA Neighbours", gp = gpar(fontface = 'bold', fontsize = 15)),
             left = textGrob("SMC biobank lncRNA promoter eQTL significance -log10(p)", gp = gpar(fontface = 'bold', fontsize = 15), rot = 90),
             grobs = Shu_BB_EQTL_LNCS_PLOT, ncol = 2, width = 10)
#per locus the top pvalue eQTLs match to a co-reg gene 

#modify 5,7,9,12 to have smaller x axis:
#Shu_BB_EQTL_LNCS_PLOT[[9]] <- Shu_BB_EQTL_LNCS_PLOT[[9]]+ theme(axis.text.x = element_text(size=7.5))

#### HMGA2-RP11221 were joined previously, would be nice to see the locus view of this and any other genes ####

#to ascertain if smth obvious has been missed
CoRegPairs_04_48_24_extended_T2_24 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CoRegPairs_04_48_24_extendedIII.csv")


#all promoter-linked variants linked to all neighbours:
Lnc2CheckEQTL <- "ENSG00000256268.1" #RP11-221N13.3, only 1 variant found, weak p val
Lnc2CheckEQTL <- "ENSG00000268858.2" #RP4-591, 0 variants found
Lnc2CheckEQTL <- "ENSG00000215039.6" #CD27-AS1, few eQTLs but best matches the co-reg
Lnc2CheckEQTL <- "MSTRG.12913" #MSTRG.12913, lots of strong eQTLs, preference for MTHFSD
Lnc2CheckEQTL <- "MSTRG.12914" #MSTRG.12914, one (weak) eQTL, preference for FOXL1
Lnc2CheckEQTL <- "MSTRG.12915" #MSTRG.12915, lots of (weak) eQTLs, preference for MTHFSD and other
Lnc2CheckEQTL <- "ENSG00000232949.1" #AC002480.4 no eQTLs
Lnc2CheckEQTL <- "ENSG00000232759.1" #AC002480.3 no eQTLs
Lnc2CheckEQTL <- "ENSG00000240476.1" #LINC00973, quite strong eQTLs linking the the co-reg gene
Lnc2CheckEQTL <- "MSTRG.24277" #LINC00973, quite strong eQTLs linking the the co-reg gene


CClncEQTL_AllVar <- filter(Shu_LncVar, EnsID %in% Lnc2CheckEQTL, 
                           #pvalue < 0.00001,
                           OverlapType == "Promoter", 
                           TotalOverlaps == 1)

CClncEQTL_AllVar <- unique(CClncEQTL_AllVar[,c(2,9,3,10,13,14,19)])


#include any non-eGene expressed pairs? not used in the calculation but may be good to be complete:
trial <- filter(AllLNC_AllPCG_2d3d, EnsID %in% Lnc2CheckEQTL, !pairs_merge %in% CClncEQTL_AllVar$pairs)

trial <- data.frame("EnsID" = trial$EnsID,
                    "EnsID.y" = trial$EnsID_merge,
                    "Variant" = NA,
                    "GeneName" = trial$EnsName.y,
                    "pvalue" = NA,
                    "FDR" = NA,
                    "pairs"= trial$pairs_merge)

trial <- rbind(CClncEQTL_AllVar, trial)

CClncEQTL_AllVar <- trial

#annotate pairs: DE/non-DE eQTL/non-eQTL/not eGene
CoRegPairs_04_48_24_extended$pairs_merge <- gsub("\\.[0-9]*$", "", CoRegPairs_04_48_24_extended$pairs)
#needed to find non-coreg eQTL linked:
CClncEQTL_AllVar2 <- filter(Shu_LncVar, EnsID %in% Lnc2CheckEQTL, 
                            pvalue < 0.00001,
                            OverlapType == "Promoter", 
                            TotalOverlaps == 1)

CClncEQTL_AllVar$PairType <- "Not CClncRNA Target"
#CClncEQTL_AllVar$PairType[CClncEQTL_AllVar$pairs %in% CClncEQTL_AllVar2$pairs] <- "eQTL-linked but not CClncRNA Target"

CClncEQTL_AllVar$PairType[CClncEQTL_AllVar$pairs %in% CoRegPairs_04_48_24_extended$pairs_merge] <- "CClncRNA Target"
#CClncEQTL_AllVar$PairType[CClncEQTL_AllVar$pairs %in% CoRegPairs_04_48_24_extended$pairs_merge &
#                                      CClncEQTL_AllVar$pairs %in% CClncEQTL_Var$pairs] <- "eQTL-linked CClncRNA Target"
table(CClncEQTL_AllVar$PairType)

#try arrange in TSS order:
allGB$EnsID_merge <- gsub("\\.[0-9]*", "", allGB$EnsID)

trial <- allGB
trial$TSS <- trial$Tx_start
trial$TSS[trial$str == "-"] <- trial$Tx_stop[trial$str == "-"]

trial <- split(trial, trial$EnsID_merge)

trial <- sapply(trial, function(x){
  mean(x$TSS)
})

trial <- data.frame("EnsID_merge" = names(trial), "meanTSS" = trial)

trial <- merge(CClncEQTL_AllVar, trial, by.x = "EnsID.y", by.y = "EnsID_merge")
trial <- merge(trial, unique(allGB[,c(3,8)]), by.x = "EnsID.y", by.y = "EnsID_merge")
CClncEQTL_AllVar <- trial

#label closest to the lncRNA:
trial <- filter(AllLNC_AllPCG_2d3d, EnsID %in% Lnc2CheckEQTL)
trial <- split(trial, trial$EnsID)

trial <- lapply(trial, function(x){
  x[x$DisLnc_PCG == min(x$DisLnc_PCG), c(13,17)]
})
trial <- bind_rows(trial)

CClncEQTL_AllVar$closest2Lnc <- "Further"
CClncEQTL_AllVar$closest2Lnc[CClncEQTL_AllVar$pairs %in% trial$pairs] <- "Closest"

CClncEQTL_AllVar$GeneName.y <- CClncEQTL_AllVar$GeneName
CClncEQTL_AllVar$GeneName.y[CClncEQTL_AllVar$closest2Lnc == "Closest"] <- paste0(CClncEQTL_AllVar$GeneName.y[CClncEQTL_AllVar$closest2Lnc == "Closest"], "*")

#label lncRNA:
trial <- merge(CClncEQTL_AllVar, unique(fpkm_allG[,c(2,5)]), by = "EnsID")

trial$LncName <- trial$EnsName
trial$LncName[is.na(trial$LncName)] <- trial$EnsID[is.na(trial$LncName)]

CClncEQTL_AllVar <- trial

selectedLocus <- CClncEQTL_AllVar

selectedLocus <- selectedLocus[order(selectedLocus$GeneName.y),]
selectedLocus$GeneName.y <- as.factor(selectedLocus$GeneName.y)
selectedLocus$GeneName.y <- factor(selectedLocus$GeneName.y, 
                                   levels = levels(selectedLocus$GeneName.y)[
                                     order(unique(selectedLocus[,c(12,9)])$meanTSS, decreasing = F)]
)

ggplot(selectedLocus) + 
  aes(x = GeneName.y, y = -log10(pvalue), color = PairType) +
  geom_hline(yintercept = -log10(1e-4), linetype = "dashed", colour = "grey60") +
  geom_boxplot(outlier.shape = NA, width = 0.3) +
  #scale_y_continuous(limits = c(0,55)) +
  #scale_color_manual(values = c(`CClncRNA Target` = "mediumorchid1",
  #                              `eQTL-linked CClncRNA Target` = "mediumorchid4",
  #                              `Not CClncRNA Target` = "darkblue",
  #                              `eQTL-linked but not CClncRNA Target` = "steelblue")) +
  scale_color_manual(values = c(`CClncRNA Target` = "mediumorchid1",
                                `Not CClncRNA Target` = "grey50")) +
  #coord_cartesian(ylim = c(0,20)) +
  geom_jitter(width = 0.1, alpha = 0.5, height = 0) +
  xlab("") +
  ylab("") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(size=7.5)) +
  Seurat::RotatedAxis() +
  ggtitle(selectedLocus$LncName)

wilcox.test(filter(selectedLocus, PairType == "Not CClncRNA Target")$pvalue, 
            filter(selectedLocus, PairType == "CClncRNA Target")$pvalue, alternative = "greater")

#### overlap identified eQTL-supported pairs with those from FANTOM analysis ####

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
length(unique(triali$pairs))#143 found both expressed in SVSMC

#part of the co-reg set:
triali <- filter(trial, pairs %in% CoRegPairs_04_48_24_extended$pairs, 
                 FANTpairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff |
                   pairs_noSuff %in% FANTOM_eQTL_pairs_cis$pairs_noSuff)
length(unique(triali$pairs))#18 found both expressed in SVSMC and in a co-reg pairing

#includes SNHG15-TBRG4

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