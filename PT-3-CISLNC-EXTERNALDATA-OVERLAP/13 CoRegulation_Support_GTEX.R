#Predicting cis-acting lncRNA via GTEXv8 eQTLs + testing some associations between lnc-PCG regulation pattern over time
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


#### GTEXv8 eQTLs subset to lnc variants ####

#targets to subset GTEX tables for 2081 PCGs near a lnc:
SVSMC_pairedPCG <- unique(AllLNC_AllPCG_2d3d$EnsID.y)
#write.csv(SVSMC_pairedPCG, "SVSMC_pairedPCG.csv", row.names = F)


#variants from subsetting the tables for potential lnc targets (in Eddie):
targetVariants <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/targetVariants",
                             header = F)
targetVariants$chr <- gsub("_.*", "", targetVariants$V1)
targetVariants$coords <- sapply(strsplit(targetVariants$V1, "_"), "[[", 2)

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
Variantoverlaps <- unique(data.frame("Variant" = targetVariants_GR$V1[queryHits(Variantindex)],
                                     "Tx" = PairedLncs_GR$MSTRG_Tx_ID[subjectHits(Variantindex)]))

length(unique(Variantoverlaps$Tx)) #561 tx
length(unique(Variantoverlaps$Variant)) #12179
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
VariantoverlapsE <- unique(data.frame("Variant" = targetVariants_GR$V1[queryHits(VariantindexE)],
                                      "Tx" = PairedLncsExons_GR$MSTRG_Tx_ID[subjectHits(VariantindexE)]))

length(unique(VariantoverlapsE$Tx)) #426 
length(unique(VariantoverlapsE$Variant)) #1470 
table(table(VariantoverlapsE$Tx)) #fewer but still some genes with insane no. variants


#additional variants in promoter region of lncs:
VariantindexP <- findOverlaps(query = targetVariants_GR, subject = promoters(PairedLncs_GR, downstream = 0))
VariantoverlapsP <- data.frame("Variant" = targetVariants_GR$V1[queryHits(VariantindexP)],
                               "Tx" = PairedLncs_GR$MSTRG_Tx_ID[subjectHits(VariantindexP)])
length(unique(VariantoverlapsP$Tx)) #496
length(unique(VariantoverlapsP$Variant)) #1549
table(table(VariantoverlapsP$Tx))


#overlap splice junctions (needs work)
PairedLncsExonLimits_GR <- makeGRangesFromDataFrame(data.frame("seqnames" = seqnames(PairedLncsExons_GR), 
                                             "start" = c(start(PairedLncsExons_GR), end(PairedLncsExons_GR)), 
                                             "end" = c(start(PairedLncsExons_GR), end(PairedLncsExons_GR))))
#add 50bp either way
PairedLncsSJ_GR <- flank(PairedLncsExonLimits_GR,  width = 50, both = T)

VariantindexSJ <- findOverlaps(query = targetVariants_GR, subject = PairedLncsSJ_GR)
VariantoverlapsSJ <- unique(data.frame("Variant" = targetVariants_GR$V1[queryHits(VariantindexSJ)],
                                       "Tx" = PairedLncsExons_GR$MSTRG_Tx_ID[subjectHits(VariantindexSJ)]))
length(unique(VariantoverlapsSJ$Tx)) #229
length(unique(VariantoverlapsSJ$Variant)) #496
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
#settled on 500bp, reduce chance of kicking out irrelevant overlaps
PairedLncsTTS_GR <- terminators(PairedLncsTTS_GR, downstream = 500, upstream = 0)

VariantindexTTS <- findOverlaps(query = targetVariants_GR, subject = PairedLncsTTS_GR)
VariantoverlapsTTS <- unique(data.frame("Variant" = targetVariants_GR$V1[queryHits(VariantindexTTS)],
                                       "Tx" = PairedLncsTTS_GR$MSTRG_Tx_ID[subjectHits(VariantindexTTS)]))
length(unique(VariantoverlapsTTS$Tx)) #510 (289 500bp)
length(unique(VariantoverlapsTTS$Variant)) #2255 (493 500bp)
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
length(unique(VariantoverlapsI$Tx)) #491
length(unique(VariantoverlapsI$Variant)) #10908
table(table(VariantoverlapsI$Tx))


#checking out individual eQTLs, 10 from each, see if in an incorrect place
#all have passed

#tag variants that overlap another genebody + promoter, some will be removed when treating overlapping/non-overlapping seperately

#<overlap all these variants to all transcripts>#

#label with no. other PCGs (exp/DE) overlapping, no. other lncs (exp/DE) overlapping

#list of variants to take back to subset the gtex tables:
CisLncVariants <- unique(c(VariantoverlapsE$Variant, VariantoverlapsP$Variant, VariantoverlapsSJ$Variant, 
                           VariantoverlapsI$Variant, VariantoverlapsTTS$Variant))
length(unique(CisLncVariants))#15352 (13956 with smaller TTS)
#write.csv(CisLncVariants, "variants.csv", row.names = F)

#### GTEXv8 eQTLs import ####

#pattern for eQTL files
eQTL_files <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/filtVar", pattern = '*signif_variant*', full.names = T)

eQTL_filesName <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/filtVar", pattern = "*signif_variant*")
eQTL_filesName <- sapply(strsplit(eQTL_filesName, "\\."), "[[", 1)
eQTL_filesName <- gsub("filtVar_", "", eQTL_filesName)

trial <- lapply(eQTL_files, function(x){
  read.delim(x, header = F)
})

for (i in 1:49){
  trial[[i]]$Tissue <- eQTL_filesName[i]
}

triali <- bind_rows(trial)

GTEX_4pairs <- triali
colnames(GTEX_4pairs) <- c("variant_id", "gene_id", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope",
                           "slope_se", "pval_nominal_threshold", "min_pval_nominal", "pval_beta", "tissueType")
length(unique(GTEX_4pairs$variant_id))#@250kbp 12966, @400kbp 15352


#worth reviewing the methods and the column info: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5776756/
#Within each tissue, cis-eQTLs were identified by linear regression, as implemented in FastQTL71, adjusting for 
#PEER factors, sex, genotyping platform, and three genotype-based principal components (PCs)

#We restricted our search to variants within 1 Mb of the TSS of each gene and, in the tissue of analysis, 
#minor allele frequencies ≥0.01 with the minor allele observed in at least 10 samples. 

#Nominal P values for each variant–gene pair were estimated using a two-tailed t-test. 

#The significance of the most highly associated variant per gene was determined from empirical P values, 
#extrapolated from a Beta distribution fitted to adaptive permutations with the setting –permute 1000 10000. 

#These empirical P values were subsequently corrected for multiple testing across genes using Storey’s q value method16. 
#To identify the list of all significant variant–gene pairs associated with eGenes, 
#variants with a nominal P value below the gene-level threshold were considered significant 
#and included in the final list of variant–gene pairs.


#pval_beta seems like the "q value" 

#a handful of the associations are weaker, not used to define eGenes: 
#"to obtain the list of eGenes, select the rows with 'qval' ≤ 0.05"

#all samples
unique(GTEX_4pairs$tissueType)
#vascular samples (artery only):
unique(GTEX_4pairs$tissueType)[c(4,5,6)]
#other muscly/mesenchymal samples
unique(GTEX_4pairs$tissueType)[c(1,2,4,5,6,21,28,29,34)]

#merge with lncRNAs
VariantoverlapsE$OverlapType <- "Exon"
VariantoverlapsP$OverlapType <- "Promoter"
VariantoverlapsSJ$OverlapType <- "Splice"
VariantoverlapsI$OverlapType <- "Intron"
VariantoverlapsTTS$OverlapType <- "TTS"

VariantoverlapsAll <- rbind(VariantoverlapsE, VariantoverlapsP, VariantoverlapsSJ, VariantoverlapsI, VariantoverlapsTTS)
length(unique(VariantoverlapsAll$Variant))#15352 as expected (13956 with TTS 500bp)

#categorise multiple overlaps - how many are to another expressed gene:
allGB_GR <- makeGRangesFromDataFrame(allGB[,c(1:6)], 
                                     start.field = "Tx_start", 
                                     end.field = "Tx_stop", 
                                     seqnames.field = "chr", 
                                     strand.field = "str", keep.extra.columns = T)
#31058 total tx to check

#check for lnc promoter/genebody/TTS overlapping variants:
#extend up from TSS by 2000kbp
allGB_GRpromoters <- promoters(allGB_GR, upstream = 2000)
#extend down from TTS by 2500kbp (HiC range) or 500bp (less stringent, chance of overlap)
allGB_GRterminators <- terminators(allGB_GR, downstream = 500)

#take union
allGB_GR <- GenomicRanges::punion(allGB_GR, allGB_GRpromoters)
allGB_GR <- GenomicRanges::punion(allGB_GR, allGB_GRterminators)

allGB_GR$EnsID <- allGB_GRpromoters$EnsID
allGB_GR$MSTRG_Tx_ID <- allGB_GRpromoters$MSTRG_Tx_ID

#focus only on  variants found in a lncRNA:
VariantoverlapsAll_GR <- targetVariants_GR[targetVariants_GR$V1 %in% VariantoverlapsAll$Variant]

Variantindex <- findOverlaps(query = VariantoverlapsAll_GR, subject = allGB_GR)
Variantoverlaps <- unique(data.frame("Variant" = VariantoverlapsAll_GR$V1[queryHits(Variantindex)],
                                     "Tx" = allGB_GR$MSTRG_Tx_ID[subjectHits(Variantindex)]))

length(unique(Variantoverlaps$Tx)) #1254 tx (1201)
length(unique(Variantoverlaps$Variant)) #15352 variants, as expected (13956)

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
#@250kbp 31 eQTLs overlap 3 genes ,10222 eQTLs overlap only 1 gene (i.e. a lncRNA)
#@400kbp 94 overlap 3, 11864 overlap 1  (i.e. a lncRNA)
#~w/500bp TTS, 37 overlap 3, 11034 overlap 1

#check on IGV - all good (bear in mind promoter regions used too)
VariantoverlapsAll <- merge(VariantoverlapsAll, trialii, by = "Variant")

GTEX_4pairs <- unique(merge(VariantoverlapsAll, GTEX_4pairs, by.x = "Variant", by.y = "variant_id"))
colnames(GTEX_4pairs)[c(2,8)] <- c("MSTRG_Tx_ID", "EnsID.y")
length(unique(GTEX_4pairs$Variant))
length(unique(VariantoverlapsAll$Variant))

length(unique(GTEX_4pairs$EnsID.y))#1285 lnc targets have an eQTL overlapping a lncRNA (1243)

head(GTEX_4pairs)

#now have obtained a final (massive) table of the SVSMC lnc-PCG pairs expressed in close genomic space which have a GTEX eQTL linkage
#as well as the position of the overlap and whether it overlaps multiple genes

#values used on GTEX browser are transformed e.g.
#-log10(filter(GTEX_4pairs, EnsID == "MSTRG.12913")$pval_nominal)

GTEX_4pairs$pval_nominal_l10 <- -log10(GTEX_4pairs$pval_nominal)

#there are other eQTLs which reach about 40 for FOXL1 but lots of eQTLs cover this lncRNA locus

#add in lncRNA IDs
trial <- merge(fpkm_allG[,c(2,49)], GTEX_4pairs, by = "MSTRG_Tx_ID")
GTEX_4pairs <- trial

#final table save:
colnames(GTEX_4pairs)
GTEX_4pairs$pairs <- paste(GTEX_4pairs$EnsID, GTEX_4pairs$EnsID.y, sep = "-")

#write.csv(GTEX_4pairs, "GTEX_4pairs.csv", row.names = F)

#CHECK IGV ON A FEW WITH MULTIPLE OVERLAPS AND IN TTS ETC THEN CONTINUE


#### does GTEX have any predictive power in finding cisLnc pairs? ####

#import table, lnc-PCG eQTL-eGenes:
GTEX_4pairs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/GTEX_4pairs_shortTSS.csv")
#GTEX_4pairs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/GTEX_4pairs.csv")

#test will be enrichment of eQTL-eGene pairs within CisLnc-target pairs vs. CisLnc-other pairs

#parameters to optimise will be 
#a) significance of the eQTL in the GTEX tissue  - measure of strength of variant association with eGene expression
#b) whether overlaps are in lncRNA promoter, exon, TTS, intron, promoter/exon or whole locus
#optional: 
#beta?
#expression level of the pairs in the GTEX tissue - simple measure of how much that gene is being used
#above benchmarks well in pilot

#done per tissue, to find the best predictive setting

#final statement will be 
#"for a given tissue, 
#eQTL-eGenes 
#(defined by eQTLs with p < x and eGenes with at least 1 eQTL <x and expression over y TPM) 
#were enriched/depleted/no effect amongst 
#CisLnc-target pairs vs. CisLnc-other gene pairs
#suggesting that eQTL-eGene pairings tend to select the same target pairs as the timing-based prediction pipeline"

#list of which of the lnc-paired PCGs are eGenes (including no eQTL overlap to a lncRNA)
eQTL_filesII <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX", pattern = "*signif_variant*", full.names = T)

#need the eGene variants BEFORE filtering for the lncRNAs to get this:
eQTL_filesNameII <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX", pattern = "*signif_variant*")
eQTL_filesNameII <- sapply(strsplit(eQTL_filesNameII, "\\."), "[[", 1)
eQTL_filesNameII <- gsub("filt_", "", eQTL_filesNameII)

trial <- lapply(eQTL_filesII, function(x){
  read.delim(x, header = F)
})

for (i in 1:49){
  trial[[i]]$Tissue <- eQTL_filesNameII[i]
}

for (i in 1:49){
  colnames(trial[[i]]) <- c("variant_id", "gene_id", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope",
                            "slope_se", "pval_nominal_threshold", "min_pval_nominal", "pval_beta", "tissueType")
}

#eGenes per tissue:
GTEX_4pairsAll <- trial
rm(trial)
#don't want to do that again:
#saveRDS(GTEX_4pairsAll, "GTEX_4pairsAll_eGenes.rds")

GTEX_4pairsAll <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/GTEX_4pairsAll_eGenes.rds")
names(GTEX_4pairsAll) <- eQTL_filesNameII

#can check here how many eGenes at a given threshold

#n.b. diff threshes required compared to Shu biobank? ideally keep the same
pThresh <- 0.05    # -0%
summary(sapply(GTEX_4pairsAll, function(x){
  length(unique(filter(x, pval_nominal < pThresh)$variant_id))/
    length(unique(x$variant_id))
}))
#pThresh <- 0.025    # -0%
#pThresh <- 0.01    # -0%
#pThresh <- 0.001   # -0%
pThresh <- 0.00005 # -5-15%
summary(sapply(GTEX_4pairsAll, function(x){
  length(unique(filter(x, pval_nominal < pThresh)$variant_id))/
    length(unique(x$variant_id))
}))
pThresh <- 0.00001 # -25-30%
summary(sapply(GTEX_4pairsAll, function(x){
  length(unique(filter(x, pval_nominal < pThresh)$variant_id))/
    length(unique(x$variant_id))
}))
pThresh <- 0.000001 # -40-45%
summary(sapply(GTEX_4pairsAll, function(x){
  length(unique(filter(x, pval_nominal < pThresh)$variant_id))/
    length(unique(x$variant_id))
}))
pThresh <- 0.0000001 # -50-55%
summary(sapply(GTEX_4pairsAll, function(x){
  length(unique(filter(x, pval_nominal < pThresh)$variant_id))/
    length(unique(x$variant_id))
}))
pThresh <- 0.00000001 # -60-70%
summary(sapply(GTEX_4pairsAll, function(x){
  length(unique(filter(x, pval_nominal < pThresh)$variant_id))/
    length(unique(x$variant_id))
}))

#import GTEX TPMs, for simple removal of any CClnc target that are not expressed
GTEX_exprs <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                         skip = 2)

#match column names up by removing the "_" and "." from each:
unique(GTEX_4pairs$tissueType)[gsub("-", "", gsub("_", "", unique(GTEX_4pairs$tissueType))) %in% 
                                 gsub("\\.", "", colnames(GTEX_exprs))]

#same order required for eQTL table, eGene tables to select correct tissue from seperate lists via a single number:
unique(GTEX_4pairs$tissueType)[2]
names(GTEX_4pairsAll)[2]
GTEX_4pairs <- GTEX_4pairs[order(GTEX_4pairs$tissueType),]
GTEX_4pairsAll <- GTEX_4pairsAll[order(names(GTEX_4pairsAll))]
#sorted:
unique(GTEX_4pairs$tissueType) == names(GTEX_4pairsAll)
unique(GTEX_4pairs$tissueType)[2]
names(GTEX_4pairsAll)[2]

#number of Cclnc neighbours that are eGenes per dataset:
exprsThreshTest <- list()
for (i in 1:49){
  #total number of Cclnc neighbours that are eGenes per dataset:
  AlleGenes <- unique(GTEX_4pairsAll[[i]]$gene_id)
  AtThresh <- GTEX_exprs$Name[GTEX_exprs[,2+i] >1]
  exprsThreshTest[[i]] <- sum(AlleGenes %in% AtThresh)/length(AlleGenes)
}
summary(unlist(exprsThreshTest)) #-5%


#### now set up/run optimisation ####

#selected p
pThresh_df <- data.frame("pThresh" = c(rep(0.05, 1),
                                       rep(0.00005, 1),# -5-15%
                                       rep(0.00001, 1), # -25-30%
                                       rep(0.000001, 1),# -40-45%
                                       rep(0.0000001, 1),# -50-55%
                                       rep(0.00000001, 1))# -60-70%
)

#store outputs here:
GTEX_eQTL_locusFish <- list()
p_test_list <- list()

#run once for each group of lncRNA overlap type
#whole locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                         #OverlapType == "Promoter", 
                                                         pval_nominal < selectedp,
                                                         TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs$tissueType)))[i], 
                                                        gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
    table(findExprs)
    triali <- filter(triali, EnsID.y %in% GTEX_exprs[findExprs, 1])
    
    triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
    
    #all potential targets for potential cisLncs
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_eQTL_locusFish[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_eQTL_locusFish[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_locusFish_dfi <- p_test_list_df
write.csv(GTEX_eQTL_locusFish_df, "GTEX_eQTL_locusFish_df.csv", row.names = F)
unique(filter(GTEX_eQTL_locusFish_df, p <0.05)$tissue)
#0x tissues (previously it was 3x but weak) (still 0 with TTS 500bp but improved)

#no significant tests before p adjust (at 250kbp there were 3x tissues)
#with HiC, same timeframe worked better


#promoter
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          OverlapType == "Promoter", 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs$tissueType)))[i], 
                                                        gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
    table(findExprs)
    triali <- filter(triali, EnsID.y %in% GTEX_exprs[findExprs, 1])
    
    triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
    
    #all potential targets for potential cisLncs
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_eQTL_locusFish[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_eQTL_locusFish[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_promoterFish_df <- p_test_list_df
write.csv(GTEX_eQTL_promoterFish_df, "GTEX_eQTL_promoterFish_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterFish_df, p <0.05)$tissue)
#8x tissues (previously it was 3x) (down to 6x with small TTS)

#seems like the hit no. stays similar whilst b and d increase (more for d though)


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          OverlapType == "Exon", 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs$tissueType)))[i], 
                                                        gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
    table(findExprs)
    triali <- filter(triali, EnsID.y %in% GTEX_exprs[findExprs, 1])
    
    triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
    
    #all potential targets for potential cisLncs
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_eQTL_locusFish[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_eQTL_locusFish[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_exonFish_df <- p_test_list_df
write.csv(GTEX_eQTL_exonFish_df, "GTEX_eQTL_exonFish_df.csv", row.names = F)
unique(filter(GTEX_eQTL_exonFish_df, p <0.05)$tissue)
#20x tissues (previously it was 27x with just 1 TPM more like ~11x)


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          OverlapType %in% c("Promoter", "Exon"), 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs$tissueType)))[i], 
                                                        gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
    table(findExprs)
    triali <- filter(triali, EnsID.y %in% GTEX_exprs[findExprs, 1])
    
    triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
    
    #all potential targets for potential cisLncs
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_eQTL_locusFish[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_eQTL_locusFish[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_promoterExonFish_df <- p_test_list_df
write.csv(GTEX_eQTL_promoterExonFish_df, "GTEX_eQTL_promoterExonFish_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterExonFish_df, p <0.05)$tissue)
#8x tissues (not done previously)
#unexpectedly weaker than exon alone, consider dropping


#Intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          OverlapType %in% c("Intron"), 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs$tissueType)))[i], 
                                                        gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
    table(findExprs)
    triali <- filter(triali, EnsID.y %in% GTEX_exprs[findExprs, 1])
    
    triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
    
    #all potential targets for potential cisLncs
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_eQTL_locusFish[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_eQTL_locusFish[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_IntronFish_df <- p_test_list_df
write.csv(GTEX_eQTL_IntronFish_df, "GTEX_eQTL_IntronFish_df.csv", row.names = F)
unique(filter(GTEX_eQTL_IntronFish_df, p <0.05)$tissue)
#3x tissue, weak (4x at TPM 1 previously, seems weaker) (1x if just using ttS 500BP )


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          OverlapType %in% c("TTS"), 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs$tissueType)))[i], 
                                                        gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
    table(findExprs)
    triali <- filter(triali, EnsID.y %in% GTEX_exprs[findExprs, 1])
    
    triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
    
    #all potential targets for potential cisLncs
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_eQTL_locusFish[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_eQTL_locusFish[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_TTSFish_df <- p_test_list_df
write.csv(GTEX_eQTL_TTSFish_df, "GTEX_eQTL_TTSFish_df.csv", row.names = F)
unique(filter(GTEX_eQTL_TTSFish_df, p <0.05)$tissue)
#1x tissue, weak, not done previously (2x TTS 500bp)


#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          OverlapType %in% c("Splice"), 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs$tissueType)))[i], 
                                                        gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
    table(findExprs)
    triali <- filter(triali, EnsID.y %in% GTEX_exprs[findExprs, 1])
    
    triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
    
    #all potential targets for potential cisLncs
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_eQTL_locusFish[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_eQTL_locusFish[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_SpliceFish_df <- p_test_list_df
write.csv(GTEX_eQTL_SpliceFish_df, "GTEX_eQTL_SpliceFish_df.csv", row.names = F)
unique(filter(GTEX_eQTL_SpliceFish_df, p <0.05)$tissue)
#0x tissue, weak, not done previously


#### now set up/run optimisation - same timeframe pairs ####

#selected p
#pThresh_df <- data.frame("pThresh" = c(rep(0.05, 1),
#                                       rep(0.00005, 1),
#                                       rep(0.00001, 1), 
#                                       rep(0.000001, 1),
#                                       rep(0.0000001, 1),
#                                       rep(0.00000001, 1))
#)

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

#store outputs here:
GTEX_eQTL_locusFish <- list()
p_test_list <- list()

#run once for each group of lncRNA overlap type
#whole locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          #OverlapType == "Promoter", 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs$tissueType)))[i], 
                                                        gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
    table(findExprs)
    triali <- filter(triali, EnsID.y %in% GTEX_exprs[findExprs, 1])
    
    triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
    
    #all potential targets for potential cisLncs
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedSame$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedSame$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_eQTL_locusFish[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_eQTL_locusFish[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_locusFish_same_df <- p_test_list_df
#write.csv(GTEX_eQTL_locusFish_same_df, "GTEX_eQTL_locusFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_locusFish_same_df, p <0.05)$tissue)
#0x tissues (previously it was 0x too)
#N.B. 7x tissues with TTS 500bp, big increase in sig relative to same/later (0x tissues)


#promoter
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          OverlapType == "Promoter", 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs$tissueType)))[i], 
                                                        gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
    table(findExprs)
    triali <- filter(triali, EnsID.y %in% GTEX_exprs[findExprs, 1])
    
    triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
    
    #all potential targets for potential cisLncs
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedSame$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedSame$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_eQTL_locusFish[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_eQTL_locusFish[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_promoterFish_same_df <- p_test_list_df
#write.csv(GTEX_eQTL_promoterFish_same_df, "GTEX_eQTL_promoterFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterFish_same_df, p <0.05)$tissue)
#5x tissues (previously it was 4x with 1 TPM)
#note more were found with same/later

#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          OverlapType == "Exon", 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs$tissueType)))[i], 
                                                        gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
    table(findExprs)
    triali <- filter(triali, EnsID.y %in% GTEX_exprs[findExprs, 1])
    
    triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
    
    #all potential targets for potential cisLncs
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedSame$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedSame$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_eQTL_locusFish[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_eQTL_locusFish[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_exonFish_same_df <- p_test_list_df
#write.csv(GTEX_eQTL_exonFish_same_df, "GTEX_eQTL_exonFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_exonFish_same_df, p <0.05)$tissue)
#13x tissues (previously it was 13 with 1 TPM)
#got 20x with same/later


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          OverlapType %in% c("Promoter", "Exon"), 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs$tissueType)))[i], 
                                                        gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
    table(findExprs)
    triali <- filter(triali, EnsID.y %in% GTEX_exprs[findExprs, 1])
    
    triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
    
    #all potential targets for potential cisLncs
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedSame$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedSame$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_eQTL_locusFish[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_eQTL_locusFish[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_promoterExonFish_same_df <- p_test_list_df
#write.csv(GTEX_eQTL_promoterExonFish_same_df, "GTEX_eQTL_promoterExonFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterExonFish_same_df, p <0.05)$tissue)
#6x tissues (not done previously)
#unexpectedly weaker than exon alone, consider dropping
#bit less then same/later (8x)


#Intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          OverlapType %in% c("Intron"), 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs$tissueType)))[i], 
                                                        gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
    table(findExprs)
    triali <- filter(triali, EnsID.y %in% GTEX_exprs[findExprs, 1])
    
    triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
    
    #all potential targets for potential cisLncs
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedSame$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedSame$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_eQTL_locusFish[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_eQTL_locusFish[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_IntronFish_same_df <- p_test_list_df
#write.csv(GTEX_eQTL_IntronFish_same_df, "GTEX_eQTL_IntronFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_IntronFish_same_df, p <0.05)$tissue)
#8x tissue, weak (10x at TPM 1 previously)
#previously 3x


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          OverlapType %in% c("TTS"), 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs$tissueType)))[i], 
                                                        gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
    table(findExprs)
    triali <- filter(triali, EnsID.y %in% GTEX_exprs[findExprs, 1])
    
    triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
    
    #all potential targets for potential cisLncs
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedSame$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedSame$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_eQTL_locusFish[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_eQTL_locusFish[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_TTSFish_same_df <- p_test_list_df
#write.csv(GTEX_eQTL_TTSFish_same_df, "GTEX_eQTL_TTSFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_TTSFish_same_df, p <0.05)$tissue)
#1x tissue, weak, not done previously


#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          OverlapType %in% c("Splice"), 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs$tissueType)))[i], 
                                                        gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
    table(findExprs)
    triali <- filter(triali, EnsID.y %in% GTEX_exprs[findExprs, 1])
    
    triali_eQTL <- filter(triali, eQTL_validated_tissue == "Yes")
    
    #all potential targets for potential cisLncs
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedSame$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedSame$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_eQTL_locusFish[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_eQTL_locusFish[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_SpliceFish_same_df <- p_test_list_df
#write.csv(GTEX_eQTL_SpliceFish_same_df, "GTEX_eQTL_SpliceFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_SpliceFish_same_df, p <0.05)$tissue)
#0x tissue, weak, not done previously


#### visualise eQTL-supported loci - promoter eQTL matches, strength of p, number/position of eQTLs relative to lncRNA TSS etc ####

#for tissues/overlap/pair types of interest

#promoter
GTEX_eQTL_promFish_df <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/GTEX_eQTL_promoterFish_df.csv")

AdSubProm_best <- filter(GTEX_eQTL_promFish_df, tissue == "Adipose_Subcutaneous", c>10)
AdSubProm_best$p_adj <- p.adjust(AdSubProm_best$p, method = "BH")
AdSubProm_best <- filter(AdSubProm_best, p == min(p))[1,]

AdVisc_best <- filter(GTEX_eQTL_promFish_df, tissue == "Adipose_Visceral_Omentum", c>10)
AdVisc_best$p_adj <- p.adjust(AdVisc_best$p, method = "BH")
AdVisc_best <- filter(AdVisc_best, p == min(p))[1,]

#exon-overlapping eQTLs that are from analysis of Skin, small intest, subcutaneous adipose, lung, breast, cerebellum: 
GTEX_eQTL_exonFish_df <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/GTEX_eQTL_exonFish_df.csv")

Skin_best <- filter(GTEX_eQTL_exonFish_df, tissue == "Skin_Not_Sun_Exposed_Suprapubic", c>10)
Skin_best$p_adj <- p.adjust(Skin_best$p, method = "BH")
Skin_best <- filter(Skin_best, p == min(p))[1,]

SmallInt_best <- filter(GTEX_eQTL_exonFish_df, tissue == "Small_Intestine_Terminal_Ileum", c>10)
SmallInt_best$p_adj <- p.adjust(SmallInt_best$p, method = "BH")
SmallInt_best <- filter(SmallInt_best, p == min(p))[1,]

AdSub_best <- filter(GTEX_eQTL_exonFish_df, tissue == "Adipose_Subcutaneous", c>10)
AdSub_best$p_adj <- p.adjust(AdSub_best$p, method = "BH")
AdSub_best <- filter(AdSub_best, p == min(p))[1,]

Lung_best <- filter(GTEX_eQTL_exonFish_df, tissue == "Lung", c>10)
Lung_best$p_adj <- p.adjust(Lung_best$p, method = "BH")
Lung_best <- filter(Lung_best, p == min(p))[1,]

Cereb_best <- filter(GTEX_eQTL_exonFish_df, tissue == "Brain_Cerebellum", c>10)
Cereb_best$p_adj <- p.adjust(Cereb_best$p, method = "BH")
Cereb_best <- filter(Cereb_best, p == min(p))[1,]

Prost_best <- filter(GTEX_eQTL_exonFish_df, tissue == "Prostate", c>10)
Prost_best$p_adj <- p.adjust(Prost_best$p, method = "BH")
Prost_best <- filter(Prost_best, p == min(p))[1,]

#TTS
GTEX_eQTL_TTSFish_df <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/GTEX_eQTL_TTSFish_df.csv")

Skin_best <- filter(GTEX_eQTL_TTSFish_df, tissue == "Skin_Not_Sun_Exposed_Suprapubic", c>10)
Skin_best$p_adj <- p.adjust(Skin_best$p, method = "BH")
Skin_best <- filter(Skin_best, p == min(p))[1,]


#same timeframe examples:
GTEX_eQTL_promFish_df <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/GTEX_eQTL_promoterFish_same_df.csv")

AdVisc_best <- filter(GTEX_eQTL_promFish_df, tissue == "Adipose_Visceral_Omentum", c>10)
AdVisc_best$p_adj <- p.adjust(AdVisc_best$p, method = "BH")
AdVisc_best <- filter(AdVisc_best, p == min(p))[1,]

GTEX_eQTL_exonFish_df <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/GTEX_eQTL_exonFish_same_df.csv")

ArtAo_best <- filter(GTEX_eQTL_exonFish_df, tissue == "Artery_Aorta", c>10)
ArtAo_best$p_adj <- p.adjust(ArtAo_best$p, method = "BH")
ArtAo_best <- filter(ArtAo_best, p == min(p))[1,]

GTEX_eQTL_intronFish_df <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/GTEX_eQTL_intronFish_same_df.csv")

Breast_best <- filter(GTEX_eQTL_intronFish_df, tissue == "Breast_Mammary_Tissue", c>10)
Breast_best$p_adj <- p.adjust(Breast_best$p, method = "BH")
Breast_best <- filter(Breast_best, p == min(p))[1,]

GTEX_eQTL_TTSFish_df <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/GTEX_eQTL_TTSFish_same_df.csv")

Nerve_best <- filter(GTEX_eQTL_TTSFish_df, tissue == "Nerve_Tibial", c>10)
Nerve_best$p_adj <- p.adjust(Nerve_best$p, method = "BH")
Nerve_best <- filter(Nerve_best, p == min(p))[1,]

ArtAo_best <- filter(GTEX_eQTL_TTSFish_df, tissue == "Artery_Aorta", c>10)
ArtAo_best$p_adj <- p.adjust(ArtAo_best$p, method = "BH")
ArtAo_best <- filter(ArtAo_best, p == min(p))[1,]


#input above selection:
tissue_best <- ArtAo_best
selectedOverlap <- "TTS"
selectedp <- as.numeric(sapply(strsplit(tissue_best$Run, "_"), "[[", 2))
Pairs2Test <- CoRegPairs_04_48_24_extendedSame

#
#### make plot of all eQTL-linked pairs per test ####

#plot all lnc-PCG pairs where PCG is a) expressed in tissue of interest and b) is an eGene in tissue of interest
findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", tissue_best$tissue)), 
                                                    gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
table(findExprs)

triali <- filter(AllLNC_AllPCG_2d3d, EnsID %in% Pairs2Test$EnsID,
                 EnsID.y %in% GTEX_exprs[findExprs, 1], 
                 EnsID.y %in% filter(GTEX_4pairsAll[[tissue_best$tissue]], 
                                     pval_nominal < selectedp)$gene_id)

triali$eQTL_validated_tissue <- "No"
triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(tissue_best$tissue, tissueType), 
                                                      OverlapType == selectedOverlap, 
                                                      pval_nominal < selectedp,
                                                      TotalOverlaps == 1)$pairs] <- "Yes"

tissueplot_AllVar <- filter(GTEX_4pairs, grepl(tissue_best$tissue, tissueType),
                                   pairs %in% triali$pairs,
                                   #pval_nominal < selectedp,
                                   OverlapType == selectedOverlap, 
                                   TotalOverlaps == 1)

tissueplot_AllVar <- unique(tissueplot_AllVar[,c(2,9,3,14,19,22)])

#include the other testable lnc-PCG pairs:
trialii <- filter(triali, !pairs %in% tissueplot_AllVar$pairs)

trialii <- data.frame("EnsID" = trialii$EnsID,
                    "EnsID.y" = trialii$EnsID.y,
                    "Variant" = NA,
                    "pval_nominal" = NA,
                    "pval_beta" = NA,
                    "pairs"= trialii$pairs)

trialii <- rbind(tissueplot_AllVar, trialii)

tissueplot_AllVar <- trialii
length(unique(triali$pairs))
length(unique(tissueplot_AllVar$pairs))

#order pairs by max eQTL sig.
trial <- sapply(split(tissueplot_AllVar, tissueplot_AllVar$pairs), function(x){
  min(x$pval_nominal)
})

bestEqtl <- data.frame("pairs" = names(trial), "eQTLbest" = trial)
bestEqtl <- bestEqtl[order(bestEqtl$pairs),]

tissueplot_AllVar$PairType <- "Other CClncRNA-Neighbour pairs"
tissueplot_AllVar$PairType[tissueplot_AllVar$pairs %in% Pairs2Test$pairs] <- "CClncRNA-Target pairs"

tissueplot_AllVar <- tissueplot_AllVar[order(tissueplot_AllVar$pairs),]
tissueplot_AllVar$pairs <- as.factor(tissueplot_AllVar$pairs)
tissueplot_AllVar$pairs <- factor(tissueplot_AllVar$pairs, levels = levels(tissueplot_AllVar$pairs)[order(bestEqtl$eQTLbest)])

#max eQTL p simpler to visualise:
bestEqtl$eQTLbest[is.na(bestEqtl$eQTLbest)] <- 1
tissueplot_BestVar <- unique(merge(tissueplot_AllVar[,-c(3:5)], bestEqtl, by = "pairs", all.x = T))

ggplot(tissueplot_BestVar) + aes(y = -log10(eQTLbest), x = pairs, fill = PairType, color = PairType) +
  geom_bar(stat = "identity") +
  annotate("text", size = 3, x = tissue_best$d*0.8, y = -log10(selectedp)*1.2, label = paste("eQTL p=", selectedp), fontface = 3, color = "grey60") +
  annotate("text", size = 3, x = tissue_best$d*0.76, y = -log10(selectedp)*1.45, label = paste("Enrichment: p=", formatC(tissue_best$p_adj, format = "e", digits = 1),
                                                                                          "OR=", round(tissue_best$OR,1)), fontface = 3, color = "grey60") +
  geom_hline(yintercept = -log10(selectedp), linetype = "dashed") +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c(`CClncRNA-Target pairs` = "mediumorchid1",
                                `Other CClncRNA-Neighbour pairs` = "grey50")) +
  scale_fill_manual(values = c(`CClncRNA-Target pairs` = "mediumorchid1",
                               `Other CClncRNA-Neighbour pairs` = "grey50")) +
  coord_cartesian(ylim = c(0,25)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), legend.position = "none", panel.grid.major =  element_blank()) +
  xlab(paste(tissue_best$d, "CClncRNA-Neighbour pairs with\neGene in", tissue_best$tissue, sep = " ")) +
  ylab("eQTL significance\n(maximum -log10p)")



#### other plots ####
tissueplot_AllVar$pval_nominal[is.na(tissueplot_AllVar$pval_nominal)] <- 1

ggplot(tissueplot_AllVar) + aes(y = -log10(pval_nominal), x = pairs, fill = PairType, color = PairType) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_point(size = 1.2) +
  annotate("text", x = tissue_best$d-50, y = -log10(selectedp)+3, label = paste("eQTLp=", selectedp), fontface = 3, color = "grey60") +
  geom_hline(yintercept = -log10(selectedp), linetype = "dashed") +
  scale_color_manual(values = c(`CClncRNA-Target pairs` = "mediumorchid1",
                                `Other CClncRNA-Neighbour pairs` = "grey50")) +
  scale_fill_manual(values = c(`CClncRNA-Target pairs` = "mediumorchid1",
                               `Other CClncRNA-Neighbour pairs` = "grey50")) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  xlab(paste(tissue_best$d, "CClncRNA-Neighbour pairs with\neGene in", tissue_best$tissue, sep = " ")) +
  ylab(paste("CCLncRNA ", selectedOverlap, "-overlapping\neQTLs -log10(p)", sep = ""))

#demonstrates concept, but also do we need all the non-eQTL pairs in the background?
ggplot(filter(tissueplot_AllVar, pval_nominal<1)) + aes(y = -log10(pval_nominal), x = pairs, fill = PairType, color = PairType) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_point() +
  annotate("text", x = tissue_best$c+15, y = -log10(selectedp)+4, label = paste("eQTLp=", selectedp), fontface = 3, color = "grey60") +
  geom_hline(yintercept = -log10(selectedp), linetype = "dashed") +
  scale_color_manual(values = c(`CClncRNA-Target pairs` = "mediumorchid1",
                                `Other CClncRNA-Neighbour pairs` = "grey50")) +
  scale_fill_manual(values = c(`CClncRNA-Target pairs` = "mediumorchid1",
                               `Other CClncRNA-Neighbour pairs` = "grey50")) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  xlab(paste(length(unique(filter(tissueplot_AllVar, pval_nominal<1)$pairs)), 
             "CClncRNA-Neighbour pairs with\neQTL-eGene in", tissue_best$tissue, sep = " ")) +
  ylab(paste("CCLncRNA ", selectedOverlap, "-overlapping\neQTLs -log10(p)", sep = ""))


#### run plot code to make of all confirmed cclnc locus per test ####
findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", tissue_best$tissue)), 
                                                    gsub("\\.", "", colnames(GTEX_exprs)))]] > 1
table(findExprs)

triali <- filter(AllLNC_AllPCG_2d3d, EnsID.y %in% GTEX_exprs[findExprs, 1], 
                 EnsID.y %in% filter(GTEX_4pairsAll[[tissue_best$tissue]], 
                                       pval_nominal < selectedp)$gene_id)

triali$eQTL_validated_tissue <- "No"
triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(tissue_best$tissue, tissueType), 
                                                      OverlapType == selectedOverlap, 
                                                      pval_nominal < selectedp,
                                                      TotalOverlaps == 1)$pairs] <- "Yes"

DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extended$pairs)
DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")

eQTLSupportedplot <- DELNC_DEPCG_1_eQTL


#number and significance of linked variants used to link up these pairs:
eQTLSupportedplot_Var <- filter(GTEX_4pairs, grepl(tissue_best$tissue, tissueType),
                                pairs %in% eQTLSupportedplot$pairs, 
                                  pval_nominal < selectedp,
                                  OverlapType == selectedOverlap, 
                                  TotalOverlaps == 1)

eQTLSupportedplot_Var <- unique(eQTLSupportedplot_Var[,c(3,9,14,19,22)])

#generally more than one eQTL supports (8/12)
table(table(eQTLSupportedplot_Var$pairs))

sapply(split(eQTLSupportedplot_Var, eQTLSupportedplot_Var$pairs), function(x){
  summary(x$pval_nominal)
})

#n.b. SHOULD ALSO SHOW THE CO-REG WITHOUT GOOD EQTL?
#compare to all variants linked to all neighbours:
eQTLSupportedplot_AllVar <- filter(GTEX_4pairs, grepl(tissue_best$tissue, tissueType),
                                   EnsID %in% eQTLSupportedplot$EnsID,
                                   pairs %in% AllLNC_AllPCG_2d3d$pairs,
                                   #pval_nominal < selectedp,
                                   OverlapType == selectedOverlap, 
                                   TotalOverlaps == 1)

eQTLSupportedplot_AllVar <- unique(eQTLSupportedplot_AllVar[,c(2,9,3,14,19,22)])


#include any non-eGene expressed pairs? not used in the calculation but may be good to visualise alongside to be complete:
#these are the pairs:
trial <- filter(AllLNC_AllPCG_2d3d, EnsID %in% eQTLSupportedplot$EnsID, 
                !pairs %in% eQTLSupportedplot_AllVar$pairs)

trial <- data.frame("EnsID" = trial$EnsID,
                    "EnsID.y" = trial$EnsID.y,
                    "Variant" = NA,
                    "pval_nominal" = NA,
                    "pval_beta" = NA,
                    "pairs"= trial$pairs)

trial <- rbind(eQTLSupportedplot_AllVar, trial)

eQTLSupportedplot_AllVar <- trial

#annotate pairs: DE/non-DE eQTL/non-eQTL/not eGene
#CoRegPairs_04_48_24_extended$pairs_merge <- gsub("\\.[0-9]*$", "", CoRegPairs_04_48_24_extended$pairs)

#needed to find non-coreg eQTL linked:
#ShuBioBankSupported_AllVar2 <- filter(Shu_LncVar, EnsID %in% ShuBioBankSupported$EnsID, 
#                                      pvalue < 0.00001,
#                                      OverlapType == "Promoter", 
#                                      TotalOverlaps == 1)

eQTLSupportedplot_AllVar$PairType <- "Not CClncRNA Target"
#ShuBioBankSupported_AllVar$PairType[ShuBioBankSupported_AllVar$pairs %in% ShuBioBankSupported_AllVar2$pairs] <- "eQTL-linked but not CClncRNA Target"

eQTLSupportedplot_AllVar$PairType[eQTLSupportedplot_AllVar$pairs %in% CoRegPairs_04_48_24_extended$pairs] <- "CClncRNA Target"
#ShuBioBankSupported_AllVar$PairType[ShuBioBankSupported_AllVar$pairs %in% CoRegPairs_04_48_24_extended$pairs_merge &
#                                      ShuBioBankSupported_AllVar$pairs %in% ShuBioBankSupported_Var$pairs] <- "eQTL-linked CClncRNA Target"
table(eQTLSupportedplot_AllVar$PairType)

EQTL_CCLNCS <- unique(eQTLSupportedplot_AllVar$EnsID)

EQTL_CCLNCS_PLOT <- list()

#max value for plot
summary(-log10(eQTLSupportedplot_AllVar$pval_nominal))

#try arrange in TSS order:
trial <- allGB
trial$TSS <- trial$Tx_start
trial$TSS[trial$str == "-"] <- trial$Tx_stop[trial$str == "-"]

trial <- split(trial, trial$EnsID)

trial <- sapply(trial, function(x){
  mean(x$TSS)
})

trial <- data.frame("EnsID" = names(trial), "meanTSS" = trial)

trial <- merge(eQTLSupportedplot_AllVar, trial, by.x = "EnsID.y", by.y = "EnsID")
trial <- merge(trial, unique(allGB[,c(3,2)]), by.x = "EnsID.y", by.y = "EnsID")
eQTLSupportedplot_AllVar <- trial

#label closest to the lncRNA:
trial <- filter(AllLNC_AllPCG_2d3d, EnsID %in% eQTLSupportedplot$EnsID)
trial <- split(trial, trial$EnsID)

trial <- lapply(trial, function(x){
  x[x$DisLnc_PCG == min(x$DisLnc_PCG), c(13,3)]
})
trial <- bind_rows(trial)

eQTLSupportedplot_AllVar$closest2Lnc <- "Further"
eQTLSupportedplot_AllVar$closest2Lnc[eQTLSupportedplot_AllVar$pairs %in% trial$pairs] <- "Closest"

#add gene names
eQTLSupportedplot_AllVar <- merge(eQTLSupportedplot_AllVar, unique(fpkm_allG[,c(2,5)]), by.x = "EnsID.y", by.y = "EnsID")

eQTLSupportedplot_AllVar$EnsName.y <- eQTLSupportedplot_AllVar$EnsName
eQTLSupportedplot_AllVar$EnsName.y[eQTLSupportedplot_AllVar$closest2Lnc == "Closest"] <- paste0(eQTLSupportedplot_AllVar$EnsName.y[eQTLSupportedplot_AllVar$closest2Lnc == "Closest"], "*")

#label lncRNA:
trial <- merge(eQTLSupportedplot_AllVar, unique(fpkm_allG[,c(2,5)]), by = "EnsID")
colnames(trial)[13] <- "EnsName"

trial$LncName <- trial$EnsName
trial$LncName[is.na(trial$LncName)] <- trial$EnsID[is.na(trial$LncName)]

eQTLSupportedplot_AllVar <- trial

for (i in 1:length(EQTL_CCLNCS)){
  
  selectedLocus <- filter(eQTLSupportedplot_AllVar, EnsID == EQTL_CCLNCS[i])
  
  selectedLocus <- selectedLocus[order(selectedLocus$EnsName.y),]
  selectedLocus$EnsName.y <- as.factor(selectedLocus$EnsName.y)
  selectedLocus$EnsName.y <- factor(selectedLocus$EnsName.y, 
                                     levels = levels(selectedLocus$EnsName.y)[
                                       order(unique(selectedLocus[,c(12,8)])$meanTSS, decreasing = F)]
                                    )
  
  EQTL_CCLNCS_PLOT[[i]] <- ggplot(selectedLocus) + 
    aes(x = EnsName.y, y = -log10(pval_nominal), color = PairType) +
    geom_hline(yintercept = -log10(selectedp), linetype = "dashed", colour = "grey60") +
    geom_boxplot(outlier.shape = NA, width = 0.3) +
    #scale_y_continuous(limits = c(0,55)) +
    #scale_color_manual(values = c(`CClncRNA Target` = "mediumorchid1",
    #                              `eQTL-linked CClncRNA Target` = "mediumorchid4",
    #                              `Not CClncRNA Target` = "darkblue",
    #                              `eQTL-linked but not CClncRNA Target` = "steelblue")) +
    scale_color_manual(values = c(`CClncRNA Target` = "mediumorchid1",
                                  `Not CClncRNA Target` = "grey50")) +
    #coord_cartesian(ylim = c(0,20)) +
    geom_point(alpha = 0.5) +
    xlab("") +
    ylab("") +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(size=7.5)) +
    Seurat::RotatedAxis() +
    ggtitle(selectedLocus$LncName)
}

#EQTL_CCLNCS_PLOT[[1]]
#EQTL_CCLNCS_PLOT[[2]]
#EQTL_CCLNCS_PLOT[[3]]
#EQTL_CCLNCS_PLOT[[4]]
#EQTL_CCLNCS_PLOT[[5]]
#EQTL_CCLNCS_PLOT[[6]]
#EQTL_CCLNCS_PLOT[[7]]
#EQTL_CCLNCS_PLOT[[8]] + theme(legend.position = "left")

#### make plot ####
library(grid)
library(gridExtra)
grid.arrange(top = textGrob("CClncRNA-targets supported by Skin Tissue eQTLs\nin lncRNA exons\n", gp = gpar(fontface = 'bold', fontsize = 18)),
             bottom = textGrob("", gp = gpar(fontface = 'bold', fontsize = 15)),
             left = textGrob("eQTL Significance -log10(p)", gp = gpar(fontface = 'bold', fontsize = 15), rot = 90),
             grobs = EQTL_CCLNCS_PLOT, ncol = 4, width = 10)
#per locus the top pvalue eQTLs match to a co-reg gene 

grid.arrange(top = textGrob("CClncRNA-targets supported by Intestine Tissue eQTLs\nin lncRNA exons\n", gp = gpar(fontface = 'bold', fontsize = 18)),
             bottom = textGrob("", gp = gpar(fontface = 'bold', fontsize = 15)),
             left = textGrob("eQTL Significance -log10(p)", gp = gpar(fontface = 'bold', fontsize = 15), rot = 90),
             grobs = EQTL_CCLNCS_PLOT, ncol = 4, width = 10)

grid.arrange(top = textGrob("CClncRNA-targets supported by Adipose Tissue eQTLs\nin lncRNA exons\n", gp = gpar(fontface = 'bold', fontsize = 18)),
             bottom = textGrob("", gp = gpar(fontface = 'bold', fontsize = 15)),
             left = textGrob("eQTL Significance -log10(p)", gp = gpar(fontface = 'bold', fontsize = 15), rot = 90),
             grobs = EQTL_CCLNCS_PLOT, ncol = 4, width = 10)

#modify 5,7,9,12 to have smaller x axis:
#Shu_BB_EQTL_LNCS_PLOT[[9]] <- Shu_BB_EQTL_LNCS_PLOT[[9]]+ theme(axis.text.x = element_text(size=7.5))


#### adding in eGene TPM (bit complex for CR, not proper-enough stats for genome?) ####

#trialling without expression at first (less p correction) but this chunk solves it:
#from GTEX online download repository, mean TPM:
GTEX_exprs <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                         skip = 2)

#match column names up by removing the "_" and "." from each:
unique(GTEX_4pairs$tissueType)[gsub("-", "", gsub("_", "", unique(GTEX_4pairs$tissueType))) %in% 
                                 gsub("\\.", "", colnames(GTEX_exprs))]

#same order required for eQTL table, eGene tables to select correct tissue via a number:
unique(GTEX_4pairs$tissueType)[2]
names(GTEX_4pairsAll)[2]


GTEX_4pairs <- GTEX_4pairs[order(GTEX_4pairs$tissueType),]
GTEX_4pairsAll <- GTEX_4pairsAll[order(names(GTEX_4pairsAll))]

#sorted:
unique(GTEX_4pairs$tissueType) == names(GTEX_4pairsAll)

#number of Cclnc neighbours that are eGenes per dataset:

exprsThreshTest <- list()
for (i in 1:49){
  #total number of Cclnc neighbours that are eGenes per dataset:
  AlleGenes <- unique(GTEX_4pairsAll[[i]]$gene_id)
  AtThresh <- GTEX_exprs$Name[GTEX_exprs[,2+i] >1]
  exprsThreshTest[[i]] <- sum(AlleGenes %in% AtThresh)/length(AlleGenes)
}
summary(unlist(exprsThreshTest)) #-5%
for (i in 1:49){
  #total number of Cclnc neighbours that are eGenes per dataset:
  AlleGenes <- unique(GTEX_4pairsAll[[i]]$gene_id)
  AtThresh <- GTEX_exprs$Name[GTEX_exprs[,2+i] >5]
  exprsThreshTest[[i]] <- sum(AlleGenes %in% AtThresh)/length(AlleGenes)
}
summary(unlist(exprsThreshTest))#-15%
for (i in 1:49){
  #total number of Cclnc neighbours that are eGenes per dataset:
  AlleGenes <- unique(GTEX_4pairsAll[[i]]$gene_id)
  AtThresh <- GTEX_exprs$Name[GTEX_exprs[,2+i] >10]
  exprsThreshTest[[i]] <- sum(AlleGenes %in% AtThresh)/length(AlleGenes)
}
summary(unlist(exprsThreshTest))#-30%
for (i in 1:49){
  #total number of Cclnc neighbours that are eGenes per dataset:
  AlleGenes <- unique(GTEX_4pairsAll[[i]]$gene_id)
  AtThresh <- GTEX_exprs$Name[GTEX_exprs[,2+i] >20]
  exprsThreshTest[[i]] <- sum(AlleGenes %in% AtThresh)/length(AlleGenes)
}
summary(unlist(exprsThreshTest))#-50%
for (i in 1:49){
  #total number of Cclnc neighbours that are eGenes per dataset:
  AlleGenes <- unique(GTEX_4pairsAll[[i]]$gene_id)
  AtThresh <- GTEX_exprs$Name[GTEX_exprs[,2+i] >30]
  exprsThreshTest[[i]] <- sum(AlleGenes %in% AtThresh)/length(AlleGenes)
}
summary(unlist(exprsThreshTest))#-65%
for (i in 1:49){
  #total number of Cclnc neighbours that are eGenes per dataset:
  AlleGenes <- unique(GTEX_4pairsAll[[i]]$gene_id)
  AtThresh <- GTEX_exprs$Name[GTEX_exprs[,2+i] >50]
  exprsThreshTest[[i]] <- sum(AlleGenes %in% AtThresh)/length(AlleGenes)
}
summary(unlist(exprsThreshTest))#-80%


