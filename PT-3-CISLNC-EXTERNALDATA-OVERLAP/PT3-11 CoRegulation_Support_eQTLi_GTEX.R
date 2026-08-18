library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(rcompanion)
library(ggbeeswarm)
library(rtracklayer)
library(GenomicRanges)
library(reshape2)

#for pairs formed in 2d, what % are eQTL linked amongst CisLnc neighbours, CisLnc targets 

#enriched?

#for eQTLs in which part of the lncRNA?

#within a certain distance only?

#exclude later timeframe pairings like with HiC? (if so then probs remove from rest of paper too i.e. fig2)

#rob young alt. idea was to do permutation and find which samples had best overlap to set of cis acting lncs


#### import and set-up tables without artefacts ####

fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026filt.csv", header = T)
length(unique(fpkm_allG$EnsID))#12740

fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE_2026filt.csv", header = T)
length(unique(fpkm_allGDE$EnsID))#5081

fpkm_allGDE_Upwithin_4 <- filter(fpkm_allGDE, RegulationStart == "Induced <4hrs")
fpkm_allGDE_Downwithin_4 <- filter(fpkm_allGDE, RegulationStart == "Repressed <4hrs")

fpkm_allGDE_Upwithin_8 <- filter(fpkm_allGDE, RegulationStart == "Induced 4-8hrs")
fpkm_allGDE_Downwithin_8 <- filter(fpkm_allGDE, RegulationStart == "Repressed 4-8hrs")

fpkm_allGDE_Upwithin_24 <- filter(fpkm_allGDE, RegulationStart == "Induced 8-24hrs")
fpkm_allGDE_Downwithin_24 <- filter(fpkm_allGDE, RegulationStart == "Repressed 8-24hrs")

fpkm_allGDE_within_4 <- rbind(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Downwithin_4)
fpkm_allGDE_within_8 <- rbind(fpkm_allGDE_Upwithin_8, fpkm_allGDE_Downwithin_8)
fpkm_allGDE_within_24 <- rbind(fpkm_allGDE_Upwithin_24, fpkm_allGDE_Downwithin_24)

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


#### import pairs ####

#all 2d pairs in 1mbp, plus any HiC up to 1Mbp
AllLNC_AllPCG_2d3d <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_2d3d_2026.csv")
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
281/8454 #1Mbp

#focus on pairs found with TSSs 1mbp
AllLNC_AllPCG_2d3d <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <1000)
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
266/8439

#some further filtering may be necessary:
#e.g. maximum range of GTEX eQTLs is 1mbp - so are longer range genes unfairly discriminated against?
#e.g. may be a minimum range of GTEX eQTLs (less likely)

#get eQTLs first


#### GTEXv8 eQTLs subset to lnc variants ####

#targets to subset GTEX tables for 5351 PCGs near a lnc:
SVSMC_pairedPCG <- unique(AllLNC_AllPCG_2d3d$EnsID.y)
#write.csv(SVSMC_pairedPCG, "SVSMC_pairedPCG_2026.csv", row.names = F)

#variants from subsetting the tables for potential lnc targets (in Eddie):
targetVariants <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/targetVariants_2026",
                             header = F)
targetVariants$chr <- gsub("_.*", "", targetVariants$V1)
targetVariants$coords <- sapply(strsplit(targetVariants$V1, "_"), "[[", 2)

SVSMC_pairedlnc <- unique(AllLNC_AllPCG_2d3d$EnsID)#579 lncs in the pairings

#isolate ranges for all lncs with a neighbouring PCG (2D/3D) - just genebody, no promoters:
allGB_LNCS <- filter(allGB, EnsID %in% filter(fpkm_allG, grepl("Lnc|fide", GeneClassUpdate))$EnsID)
length(unique(allGB_LNCS$EnsID))#597 as expected
length(unique(allGB_LNCS$MSTRG_Tx_ID))#1575 as expected

allGB_LNCS_GR <- makeGRangesFromDataFrame(allGB_LNCS[,c(1:6)], 
                                          start.field = "Tx_start", 
                                          end.field = "Tx_stop", 
                                          seqnames.field = "chr", 
                                          strand.field = "str", keep.extra.columns = T)
PairedLncs_GR <- allGB_LNCS_GR[allGB_LNCS_GR$EnsID %in% SVSMC_pairedlnc]
PairedLncs_GR
#1516 total lnc tx for those paired to a neighbour expressed PCG 

#check for lnc promoter/genebody overlapping variants:
targetVariants_GR <- makeGRangesFromDataFrame(targetVariants, start.field = "coords", end.field = "coords", keep.extra.columns = T)

Variantindex <- findOverlaps(query = targetVariants_GR, subject = PairedLncs_GR)
Variantoverlaps <- unique(data.frame("Variant" = targetVariants_GR$V1[queryHits(Variantindex)],
                                     "Tx" = PairedLncs_GR$MSTRG_Tx_ID[subjectHits(Variantindex)]))
length(unique(Variantoverlaps$Tx)) #1329 tx
length(unique(Variantoverlaps$Variant)) #23370
table(table(Variantoverlaps$Tx)) #some tx with insane no. variants, long intron lncRNAs? highly variable regions?


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
#obtain the transcripts for paired lncs (all >0.8fpkm tx otherwise wouldn't be in this list)
stringtie_gtf_majorPairedLncs <- filter(stringtie_gtf, MSTRG_Tx_ID %in% PairedLncs_GR$MSTRG_Tx_ID)
stringtie_gtf_majorPairedLncs <- merge(stringtie_gtf_majorPairedLncs, fpkm_allG[,c(2,47)], by = "MSTRG_Tx_ID")
length(unique(stringtie_gtf_majorPairedLncs$EnsID))#579 genes
length(unique(stringtie_gtf_majorPairedLncs$MSTRG_Tx_ID))#1516 transcripts

stringtie_gtf_majorPairedLncsExons <- filter(stringtie_gtf_majorPairedLncs, V3 == "exon")

PairedLncsExons_GR <- makeGRangesFromDataFrame(stringtie_gtf_majorPairedLncsExons[,c(2,5,6,1,11)], seqnames.field = "V1", 
                                               start.field = "V4", end.field = "V5", keep.extra.columns = T)

#overlap exons
VariantindexE <- findOverlaps(query = targetVariants_GR, subject = PairedLncsExons_GR)
VariantoverlapsE <- unique(data.frame("Variant" = targetVariants_GR$V1[queryHits(VariantindexE)],
                                      "Tx" = PairedLncsExons_GR$MSTRG_Tx_ID[subjectHits(VariantindexE)]))
length(unique(VariantoverlapsE$Tx)) #968 
length(unique(VariantoverlapsE$Variant)) #2849 
table(table(VariantoverlapsE$Tx)) #fewer but still some genes with insane no. variants

#additional variants in promoter region of lncs:
VariantindexP <- findOverlaps(query = targetVariants_GR, subject = promoters(PairedLncs_GR, downstream = 0))
VariantoverlapsP <- data.frame("Variant" = targetVariants_GR$V1[queryHits(VariantindexP)],
                               "Tx" = PairedLncs_GR$MSTRG_Tx_ID[subjectHits(VariantindexP)])
length(unique(VariantoverlapsP$Tx)) #1153
length(unique(VariantoverlapsP$Variant)) #2710
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
length(unique(VariantoverlapsSJ$Tx)) #530
length(unique(VariantoverlapsSJ$Variant)) #965
table(table(VariantoverlapsSJ$Tx)) #fewer 


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
#settled on 500bp, reduce chance of kicking out irrelevant overlaps for others
PairedLncsTTS_GR <- terminators(PairedLncsTTS_GR, downstream = 500, upstream = 0)

VariantindexTTS <- findOverlaps(query = targetVariants_GR, subject = PairedLncsTTS_GR)
VariantoverlapsTTS <- unique(data.frame("Variant" = targetVariants_GR$V1[queryHits(VariantindexTTS)],
                                        "Tx" = PairedLncsTTS_GR$MSTRG_Tx_ID[subjectHits(VariantindexTTS)]))
length(unique(VariantoverlapsTTS$Tx)) #694
length(unique(VariantoverlapsTTS$Variant)) #978
#n.b. some potential issue in above numbers (decreased unexpectedly when expanding number of variants, no idea why, could be a previous fault)
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
length(unique(VariantoverlapsI$Tx)) #1190
length(unique(VariantoverlapsI$Variant)) #21557
table(table(VariantoverlapsI$Tx))

#checking out individual eQTLs, 10 from each, see if in an incorrect place
#all have passed

#tag variants that overlap another genebody + promoter, some will be removed when treating overlapping/non-overlapping seperately

#<overlap all these variants to all transcripts>#

#label with no. other PCGs (exp/DE) overlapping, no. other lncs (exp/DE) overlapping

#list of variants to take back to subset the gtex tables:
CisLncVariants <- unique(c(VariantoverlapsE$Variant, VariantoverlapsP$Variant, VariantoverlapsSJ$Variant, 
                           VariantoverlapsI$Variant, VariantoverlapsTTS$Variant))
length(unique(CisLncVariants))#25966
#write.csv(CisLncVariants, "variants_2026.csv", row.names = F)


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
length(unique(GTEX_4pairs$variant_id))#25966 variants

#worth reviewing the methods and the column info: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5776756/
#Within each tissue, cis-eQTLs were identified by linear regression, as implemented in FastQTL[71], adjusting for 
#PEER factors, sex, genotyping platform, and three genotype-based principal components (PCs)

#We restricted our search to variants within 1 Mb of the TSS of each gene and, in the tissue of analysis, 
#minor allele frequencies ≥0.01 with the minor allele observed in at least 10 samples. 

#Nominal P values for each variant–gene pair were estimated using a two-tailed t-test. 

#The significance of the most highly associated variant per gene was determined from empirical P values, 
#extrapolated from a Beta distribution fitted to adaptive permutations with the setting –permute 1000 10000. 

#These empirical P values were subsequently corrected for multiple testing across genes using Storey’s q value method[16]. 
#To identify the list of all significant variant–gene pairs associated with eGenes, 
#variants with a nominal P value below the gene-level threshold were considered significant 
#and included in the final list of variant–gene pairs.

#pval_beta seems like the "q value" 

#a handful of the associations are weaker, not used to define eGenes: 
#"to obtain the list of eGenes, select the rows with 'qval' ≤ 0.05"

#all below the min pval nominal
dim(filter(GTEX_4pairs, min_pval_nominal >0.05))

#will handle p val filtering later

#all samples
unique(GTEX_4pairs$tissueType)
#vascular samples (artery only):
unique(GTEX_4pairs$tissueType)[c(4,5,6)]
#other muscly/mesenchymal samples
unique(GTEX_4pairs$tissueType)[c(1,2,4,5,6,21,28,29,34)]

#rob young idea was to do permutation and find which samples had best overlap to set of cis acting lncs

#merge with lncRNAs
VariantoverlapsE$OverlapType <- "Exon"
VariantoverlapsP$OverlapType <- "Promoter"
VariantoverlapsSJ$OverlapType <- "Splice"
VariantoverlapsI$OverlapType <- "Intron"
VariantoverlapsTTS$OverlapType <- "TTS"

VariantoverlapsAll <- rbind(VariantoverlapsE, VariantoverlapsP, VariantoverlapsSJ, VariantoverlapsI, VariantoverlapsTTS)
length(unique(VariantoverlapsAll$Tx))#1386 lncRNA transcripts have
length(unique(VariantoverlapsAll$Variant))#25966 variants overlapping them

#optional: filter tx more strongly (accept any drop outs)
fpkm_allG_hiXP <- filter(fpkm_allG, Tx_Max_Average > 0.8)
length(unique(fpkm_allG$EnsID))
length(unique(fpkm_allG_hiXP$EnsID))

allGB_hiXP <- filter(allGB, MSTRG_Tx_ID %in% fpkm_allG_hiXP$MSTRG_Tx_ID)

#categorise multiple overlaps - how many are to another expressed gene:
allGB_GR <- makeGRangesFromDataFrame(allGB_hiXP[,c(1:6)],
                                     #allGB[,c(1:6)], 
                                     start.field = "Tx_start", 
                                     end.field = "Tx_stop", 
                                     seqnames.field = "chr", 
                                     strand.field = "str", keep.extra.columns = T)
#42511 total tx to check (previously 31058)
#hard filter gets 30675

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

length(unique(Variantoverlaps$Tx)) #2384 transcript for the (1254 previously) (1579 with hard xp filter)
length(unique(Variantoverlaps$Variant)) #25966 variants overlapping a lncRNA (19701 with hard xp filter)

#add gene to a) variants overlapping lnc table:
trial <- merge(fpkm_allG[,c(2,47)], VariantoverlapsAll, by.x = "MSTRG_Tx_ID", by.y = "Tx")
#and to b) variants overlapping all genes table:
triali <- merge(fpkm_allG[,c(2,47)], Variantoverlaps, by.x = "MSTRG_Tx_ID", by.y = "Tx")

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
#w/2026: 178 overlap 3, 20603 overlap 1 (there are more lncRNAs now)
#hard xp: 128 overlap 3, 15067 overlap 1

#seems like a lot more overlaps... regardless of tx filtering (non-merged genes overlap or close prox...)

#check on IGV - all good (bear in mind promoter regions used too)
VariantoverlapsAll <- merge(VariantoverlapsAll, trialii, by = "Variant")

GTEX_4pairs <- unique(merge(VariantoverlapsAll, GTEX_4pairs, by.x = "Variant", by.y = "variant_id"))
colnames(GTEX_4pairs)[c(2,8)] <- c("MSTRG_Tx_ID", "EnsID.y")
length(unique(GTEX_4pairs$Variant))
length(unique(VariantoverlapsAll$Variant))

#add in lncRNA IDs
trial <- merge(fpkm_allG[,c(2,47)], GTEX_4pairs, by = "MSTRG_Tx_ID")
GTEX_4pairs <- trial

#remove genes with only low xp tx from here if doing this option
GTEX_4pairs_ <- filter(GTEX_4pairs, 
                       EnsID %in% fpkm_allG_hiXP$EnsID, 
                       EnsID.y %in% fpkm_allG_hiXP$EnsID)

length(unique(GTEX_4pairs$EnsID))#538 lnc targets have an eQTL overlapping a lncRNA (452 with hiXP)
length(unique(GTEX_4pairs$EnsID.y))#2466 lnc targets have an eQTL overlapping a lncRNA (2021 with hiXP)

GTEX_4pairs <- GTEX_4pairs_

#maximum and minimum eQTL-eGene distances from GTEX:

#consider removing lncRNA-PCGs not wholly contained within these max/mins?
summary(abs(GTEX_4pairs$tss_distance))

#final table save:
colnames(GTEX_4pairs)
GTEX_4pairs$pairs <- paste(GTEX_4pairs$EnsID, GTEX_4pairs$EnsID.y, sep = "-")

#write.csv(GTEX_4pairs, "GTEX_4pairs_2026_hiXP.csv", row.names = F)


#CHECK IGV ON A FEW WITH MULTIPLE OVERLAPS AND IN TTS ETC THEN CONTINUE

#now have obtained a final (massive) table of the SVSMC lnc-PCG pairs expressed in close genomic space which have a GTEX eQTL linkage
#as well as the position of the overlap and whether it overlaps multiple genes

#values used on GTEX browser are transformed e.g.
#-log10(filter(GTEX_4pairs, EnsID == "MSTRG.12913")$pval_nominal)

GTEX_4pairs$pval_nominal_l10 <- -log10(GTEX_4pairs$pval_nominal)

#there are other eQTLs which reach about 40 for FOXL1 but lots of eQTLs cover this lncRNA locus


#### does GTEX have any predictive power in finding cisLnc pairs? ####

#import table, lnc-PCG eQTL-eGenes:
GTEX_4pairs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/GTEX_4pairs_2026.csv")
#GTEX_4pairs_hiXP <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/GTEX_4pairs_2026_hiXP.csv")

length(unique(GTEX_4pairs$EnsID))#538 lnc targets have an eQTL overlapping a lncRNA (452 with hiXP)
length(unique(GTEX_4pairs$EnsID.y))#2466 lnc targets have an eQTL overlapping a lncRNA (2021 with hiXP)

#length(unique(GTEX_4pairs_hiXP$EnsID))#(452 with hiXP)
#length(unique(GTEX_4pairs_hiXP$EnsID.y))#(2021 with hiXP)

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
#(defined by eQTLs with p < x and eGenes with at least 1 eQTL <x and expression over 1 TPM) 
#were enriched/depleted/no effect amongst
#CisLnc-target pairs vs. CisLnc-other gene pairs
#suggesting that eQTL-eGene pairings tend to select the same target pairs as the timing-based prediction pipeline"

#might be wise to seperate same-timeframe and later-timeframe as with HiC

#list of which of the lnc-paired PCGs are eGenes (including no eQTL overlap to a lncRNA)
eQTL_filesII <- list.files("heavyFiles4R/", 
                           pattern = "*signif_variant*", full.names = T)

#need the eGene variants BEFORE filtering for the lncRNAs to get this:
eQTL_filesNameII <- list.files("heavyFiles4R/", pattern = "*signif_variant*")
eQTL_filesNameII <- sapply(strsplit(eQTL_filesNameII, "\\."), "[[", 1)
eQTL_filesNameII <- gsub("filt_", "", eQTL_filesNameII)

#trial <- lapply(eQTL_filesII, function(x){
#  read.delim(x, header = F)
#})

#for (i in 1:49){
#  trial[[i]]$Tissue <- eQTL_filesNameII[i]
#}

#for (i in 1:49){
#  colnames(trial[[i]]) <- c("variant_id", "gene_id", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope",
#                            "slope_se", "pval_nominal_threshold", "min_pval_nominal", "pval_beta", "tissueType")
#}

#eGenes per tissue:
#GTEX_4pairsAll <- trial
#rm(trial)
#don't want to do that again:
#saveRDS(GTEX_4pairsAll, "GTEX_4pairsAll_eGenes_2026.rds")

GTEX_4pairsAll <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/GTEX/GTEX_4pairsAll_eGenes_2026.rds")
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
pThresh <- 0.00005 # -4-12%
summary(sapply(GTEX_4pairsAll, function(x){
  length(unique(filter(x, pval_nominal < pThresh)$variant_id))/
    length(unique(x$variant_id))
}))
pThresh <- 0.00001 # -22-27%
summary(sapply(GTEX_4pairsAll, function(x){
  length(unique(filter(x, pval_nominal < pThresh)$variant_id))/
    length(unique(x$variant_id))
}))
pThresh <- 0.000001 # -40-43%
summary(sapply(GTEX_4pairsAll, function(x){
  length(unique(filter(x, pval_nominal < pThresh)$variant_id))/
    length(unique(x$variant_id))
}))
pThresh <- 0.0000001 # -50-56%
summary(sapply(GTEX_4pairsAll, function(x){
  length(unique(filter(x, pval_nominal < pThresh)$variant_id))/
    length(unique(x$variant_id))
}))
pThresh <- 0.00000001 # -57-65%
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
#GTEX_4pairs_hiXP <- GTEX_4pairs_hiXP[order(GTEX_4pairs_hiXP$tissueType),]
GTEX_4pairsAll <- GTEX_4pairsAll[order(names(GTEX_4pairsAll))]
#sorted:
unique(GTEX_4pairs$tissueType) == names(GTEX_4pairsAll)
unique(GTEX_4pairs$tissueType)[2]
#unique(GTEX_4pairs_hiXP$tissueType)[2]
names(GTEX_4pairsAll)[2]
colnames(GTEX_exprs)[2+2]

#number of CClnc neighbours that are eGenes per dataset:
exprsThreshTest <- list()
for (i in 1:49){
  #total number of Cclnc neighbours that are eGenes per dataset:
  AlleGenes <- unique(GTEX_4pairsAll[[i]]$gene_id)
  AtThresh <- GTEX_exprs$Name[GTEX_exprs[,2+i] >1]
  exprsThreshTest[[i]] <- sum(AlleGenes %in% AtThresh)/length(AlleGenes)
}
summary(unlist(exprsThreshTest)) #-5%

#lnc-PCG pairs not wholly within 1mbp have less chance to be picked up - remove
#e.g. a lncRNA with a TSS 1mbp away from a PCG, being transcribed away from the PCG, would be unfair to include
#because pretty much all of it's gene body is too far for an eQTL-eGene connection
AllLNC_AllPCG_2d3d <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_2d3d_2026.csv")
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
281/8454 #numbers for 1Mbp TSS dist + HiC

AllLNC_AllPCG_2d3d <- filter(AllLNC_AllPCG_2d3d, pair_range < 1000000)
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
260/8116 #numbers for within 1Mbp

#lnc_pcg_ranges <- list()
#for (i in 1:length(AllLNC_AllPCG_2d3d$pairs)){
#  lnc_pcg_ranges[[i]] <- range(allGB_LNCS_GR[allGB_LNCS_GR$EnsID %in% AllLNC_AllPCG_2d3d[i,2]], 
#                               allGB_PCGs_GR[allGB_PCGs_GR$EnsID == AllLNC_AllPCG_2d3d[i,1]], ignore.strand = T)
#}

#length(unique(AllLNC_AllPCG_1$EnsID))
#AllLNC_AllPCG_2d3d$pair_range <- sapply(lnc_pcg_ranges, width)
#AllLNC_AllPCG_2d3d$pair_range_bin <- NA
#AllLNC_AllPCG_2d3d$pair_range_bin[AllLNC_AllPCG_2d3d$pair_range > 1000000] <- ">1000kbp"
#AllLNC_AllPCG_2d3d$pair_range_bin[AllLNC_AllPCG_2d3d$pair_range < 1000000] <- "500-1000kbp"
#AllLNC_AllPCG_2d3d$pair_range_bin[AllLNC_AllPCG_2d3d$pair_range < 500000] <- "400-500kbp"
#AllLNC_AllPCG_2d3d$pair_range_bin[AllLNC_AllPCG_2d3d$pair_range < 400000] <- "300-400kbp"
#AllLNC_AllPCG_2d3d$pair_range_bin[AllLNC_AllPCG_2d3d$pair_range < 300000] <- "200-300kbp"
#AllLNC_AllPCG_2d3d$pair_range_bin[AllLNC_AllPCG_2d3d$pair_range < 200000] <- "100-200kbp"
#AllLNC_AllPCG_2d3d$pair_range_bin[AllLNC_AllPCG_2d3d$pair_range < 100000] <- "50-100kbp"
#AllLNC_AllPCG_2d3d$pair_range_bin[AllLNC_AllPCG_2d3d$pair_range < 50000] <- "25-50kbp"
#AllLNC_AllPCG_2d3d$pair_range_bin[AllLNC_AllPCG_2d3d$pair_range < 25000] <- "15-25kbp"
#AllLNC_AllPCG_2d3d$pair_range_bin[AllLNC_AllPCG_2d3d$pair_range < 15000] <- "<15kbp"
#table(AllLNC_AllPCG_2d3d$pair_range_bin)

#should save this to avoid doing every time
#write.csv(AllLNC_AllPCG_2d3d, "AllLNC_AllPCG_2d3d_2026.csv", row.names = F)


#### 1a same timeframe pairs ####

#selected p
pThresh_df <- data.frame("pThresh" = c(rep(0.05, 1),
                                       rep(0.00005, 1),
                                       rep(0.00001, 1), 
                                       rep(0.000001, 1),
                                       rep(0.0000001, 1),
                                       rep(0.00000001, 1)))

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
                                                                                                             fpkm_allGDE_Downwithin_24$EnsID)))
#456

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
#8x tissues, 0x previously

#p adjust:
#select best p per Tissue_pair_overlap:
trial <- split(GTEX_eQTL_locusFish_same_df, GTEX_eQTL_locusFish_same_df$tissue)

triali <- lapply(trial, function(x){
  y <- filter(x, c >10)
  y$V1_padj <- p.adjust(y$p, method = "BH")
  return(y[order(y$V1, decreasing = F),][1,])
})

triali <- bind_rows(triali)
unique(filter(triali, V1_padj <0.05)$tissue)
#3x tissues, 0x previously
#includes 2x SMC-enriched tissue


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
#4x tissues


#p adjust:
#select best p per Tissue_pair_overlap:
trial <- split(GTEX_eQTL_promoterFish_same_df, GTEX_eQTL_promoterFish_same_df$tissue)

triali <- lapply(trial, function(x){
  y <- filter(x, c >10)
  y$V1_padj <- p.adjust(y$p, method = "BH")
  return(y[order(y$V1, decreasing = F),][1,])
})

triali <- bind_rows(triali)
unique(filter(triali, V1_padj <0.05)$tissue)
#1x tissue only (brain)


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
#10x tissues

#p adjust:
#select best p per Tissue_pair_overlap:
trial <- split(GTEX_eQTL_exonFish_same_df, GTEX_eQTL_exonFish_same_df$tissue)

triali <- lapply(trial, function(x){
  y <- filter(x, c >10)
  y$V1_padj <- p.adjust(y$p, method = "BH")
  return(y[order(y$V1, decreasing = F),][1,])
})

triali <- bind_rows(triali)
unique(filter(triali, V1_padj <0.05)$tissue)
#2x brain tissues only


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
#9x tissues

#p adjust:
#select best p per Tissue_pair_overlap:
trial <- split(GTEX_eQTL_promoterExonFish_same_df, GTEX_eQTL_promoterExonFish_same_df$tissue)

triali <- lapply(trial, function(x){
  y <- filter(x, c >10)
  y$V1_padj <- p.adjust(y$p, method = "BH")
  return(y[order(y$V1, decreasing = F),][1,])
})

triali <- bind_rows(triali)
unique(filter(triali, V1_padj <0.05)$tissue)
#2x brain tissues only


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
#2x tissues


#p adjust:
#select best p per Tissue_pair_overlap:
trial <- split(GTEX_eQTL_IntronFish_same_df, GTEX_eQTL_IntronFish_same_df$tissue)

triali <- lapply(trial, function(x){
  y <- filter(x, c >10)
  y$V1_padj <- p.adjust(y$p, method = "BH")
  return(y[order(y$V1, decreasing = F),][1,])
})

triali <- bind_rows(triali)
unique(filter(triali, V1_padj <0.05)$tissue)
#0x tissues


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
#7x tissues

#p adjust:
#select best p per Tissue_pair_overlap:
trial <- split(GTEX_eQTL_TTSFish_same_df, GTEX_eQTL_TTSFish_same_df$tissue)

triali <- lapply(trial, function(x){
  y <- filter(x, c >10)
  y$V1_padj <- p.adjust(y$p, method = "BH")
  return(y[order(y$V1, decreasing = F),][1,])
})

triali <- bind_rows(triali)
unique(filter(triali, V1_padj <0.05)$tissue)
#1x tissues


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
#1x tissue

#p adjust:
#select best p per Tissue_pair_overlap:
trial <- split(GTEX_eQTL_SpliceFish_same_df, GTEX_eQTL_SpliceFish_same_df$tissue)

triali <- lapply(trial, function(x){
  y <- filter(x, c >10)
  y$V1_padj <- p.adjust(y$p, method = "BH")
  return(y[order(y$V1, decreasing = F),][1,])
})

triali <- bind_rows(triali)
unique(filter(triali, V1_padj <0.05)$tissue)
#0x tissues

#whole locus makes most sense, SMC-enriched tissue
#otherwise... brain tissue??

#some decent results
unique(filter(GTEX_eQTL_locusFish_same_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_promoterFish_same_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_exonFish_same_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_IntronFish_same_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_TTSFish_same_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_SpliceFish_same_df, p <0.05)$tissue)



#### 1b same timeframe pairs - hiXP transcripts only ####

#note that not recalculating no. genes overlapping the eQTL - so still considering low XP for that

#selected p
pThresh_df <- data.frame("pThresh" = c(rep(0.05, 1),
                                       rep(0.00005, 1),
                                       rep(0.00001, 1), 
                                       rep(0.000001, 1),
                                       rep(0.0000001, 1),
                                       rep(0.00000001, 1)))

AllLNC_AllPCG_2d3d_hi <- filter(AllLNC_AllPCG_2d3d, EnsID %in% fpkm_allG_hiXP$EnsID, EnsID.y %in% fpkm_allG_hiXP$EnsID)

CoRegPairs_04_48_24_extendedSame <- filter(AllLNC_AllPCG_2d3d_hi,
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
  
  for (i in 1:length(unique(GTEX_4pairs_hiXP$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_hi
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs_hiXP, grepl(unique(GTEX_4pairs_hiXP$tissueType)[i], tissueType), 
                                                          #OverlapType == "Promoter", 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs_hiXP$tissueType)))[i], 
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
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs_hiXP$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_locusFish_samehiXP_df <- p_test_list_df
#write.csv(GTEX_eQTL_locusFish_samehiXP_df, "GTEX_eQTL_locusFish_samehiXP_df.csv", row.names = F)
unique(filter(GTEX_eQTL_locusFish_samehiXP_df, p <0.05)$tissue)
#6x tissues, 

#p adjust:
#select best p per Tissue_pair_overlap:
trial <- split(GTEX_eQTL_locusFish_samehiXP_df, GTEX_eQTL_locusFish_samehiXP_df$tissue)

triali <- lapply(trial, function(x){
  y <- filter(x, c >10)
  y$V1_padj <- p.adjust(y$p, method = "BH")
  return(y[order(y$V1, decreasing = F),][1,])
})

triali <- bind_rows(triali)
unique(filter(triali, V1_padj <0.05)$tissue)
#0x tissues (3 with all)


#promoter
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs_hiXP$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_hi
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs_hiXP, grepl(unique(GTEX_4pairs_hiXP$tissueType)[i], tissueType), 
                                                          OverlapType == "Promoter", 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs_hiXP$tissueType)))[i], 
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
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs_hiXP$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_promoterFish_samehiXP_df <- p_test_list_df
#write.csv(GTEX_eQTL_promoterFish_samehiXP_df, "GTEX_eQTL_promoterFish_samehiXP_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterFish_samehiXP_df, p <0.05)$tissue)
#4x tissues

#p adjust:
#select best p per Tissue_pair_overlap:
trial <- split(GTEX_eQTL_promoterFish_samehiXP_df, GTEX_eQTL_promoterFish_samehiXP_df$tissue)

triali <- lapply(trial, function(x){
  y <- filter(x, c >10)
  y$V1_padj <- p.adjust(y$p, method = "BH")
  return(y[order(y$V1, decreasing = F),][1,])
})

triali <- bind_rows(triali)
unique(filter(triali, V1_padj <0.05)$tissue)
#1x tissue only (brain)


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs_hiXP$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_hi
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs_hiXP, grepl(unique(GTEX_4pairs_hiXP$tissueType)[i], tissueType), 
                                                          OverlapType == "Exon", 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs_hiXP$tissueType)))[i], 
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
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs_hiXP$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_exonFish_samehiXP_df <- p_test_list_df
#write.csv(GTEX_eQTL_exonFish_samehiXP_df, "GTEX_eQTL_exonFish_samehiXP_df.csv", row.names = F)
unique(filter(GTEX_eQTL_exonFish_samehiXP_df, p <0.05)$tissue)
#7x tissues

#p adjust:
#select best p per Tissue_pair_overlap:
trial <- split(GTEX_eQTL_exonFish_samehiXP_df, GTEX_eQTL_exonFish_samehiXP_df$tissue)

triali <- lapply(trial, function(x){
  y <- filter(x, c >10)
  y$V1_padj <- p.adjust(y$p, method = "BH")
  return(y[order(y$V1, decreasing = F),][1,])
})

triali <- bind_rows(triali)
unique(filter(triali, V1_padj <0.05)$tissue)
#1x brain tissues only

#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs_hiXP$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_hi
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs_hiXP, grepl(unique(GTEX_4pairs_hiXP$tissueType)[i], tissueType), 
                                                          OverlapType %in% c("Promoter", "Exon"), 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs_hiXP$tissueType)))[i], 
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
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs_hiXP$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_promoterExonFish_samehiXP_df <- p_test_list_df
#write.csv(GTEX_eQTL_promoterExonFish_samehiXP_df, "GTEX_eQTL_promoterExonFish_samehiXP_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterExonFish_samehiXP_df, p <0.05)$tissue)
#5x tissues

#p adjust:
#select best p per Tissue_pair_overlap:
trial <- split(GTEX_eQTL_promoterExonFish_samehiXP_df, GTEX_eQTL_promoterExonFish_samehiXP_df$tissue)

triali <- lapply(trial, function(x){
  y <- filter(x, c >10)
  y$V1_padj <- p.adjust(y$p, method = "BH")
  return(y[order(y$V1, decreasing = F),][1,])
})

triali <- bind_rows(triali)
unique(filter(triali, V1_padj <0.05)$tissue)
#1x brain tissues only


#Intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs_hiXP$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs_hiXP, grepl(unique(GTEX_4pairs_hiXP$tissueType)[i], tissueType), 
                                                          OverlapType %in% c("Intron"), 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs_hiXP$tissueType)))[i], 
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
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs_hiXP$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_IntronFish_samehiXP_df <- p_test_list_df
#write.csv(GTEX_eQTL_IntronFish_samehiXP_df, "GTEX_eQTL_IntronFish_samehiXP_df.csv", row.names = F)
unique(filter(GTEX_eQTL_IntronFish_samehiXP_df, p <0.05)$tissue)
#3x tissues

#p adjust:
#select best p per Tissue_pair_overlap:
trial <- split(GTEX_eQTL_IntronFish_samehiXP_df, GTEX_eQTL_IntronFish_samehiXP_df$tissue)

triali <- lapply(trial, function(x){
  y <- filter(x, c >10)
  y$V1_padj <- p.adjust(y$p, method = "BH")
  return(y[order(y$V1, decreasing = F),][1,])
})

triali <- bind_rows(triali)
unique(filter(triali, V1_padj <0.05)$tissue)
#0x tissues


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs_hiXP$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs_hiXP, grepl(unique(GTEX_4pairs_hiXP$tissueType)[i], tissueType), 
                                                          OverlapType %in% c("TTS"), 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs_hiXP$tissueType)))[i], 
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
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs_hiXP$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_TTSFish_samehiXP_df <- p_test_list_df
#write.csv(GTEX_eQTL_TTSFish_samehiXP_df, "GTEX_eQTL_TTSFish_samehiXP_df.csv", row.names = F)
unique(filter(GTEX_eQTL_TTSFish_samehiXP_df, p <0.05)$tissue)
#6x tissues

#p adjust:
#select best p per Tissue_pair_overlap:
trial <- split(GTEX_eQTL_TTSFish_samehiXP_df, GTEX_eQTL_TTSFish_samehiXP_df$tissue)

triali <- lapply(trial, function(x){
  y <- filter(x, c >10)
  y$V1_padj <- p.adjust(y$p, method = "BH")
  return(y[order(y$V1, decreasing = F),][1,])
})

triali <- bind_rows(triali)
unique(filter(triali, V1_padj <0.05)$tissue)
#2x tissues


#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs_hiXP$tissueType))){
    triali <- AllLNC_AllPCG_2d3d
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs_hiXP, grepl(unique(GTEX_4pairs_hiXP$tissueType)[i], tissueType), 
                                                          OverlapType %in% c("Splice"), 
                                                          pval_nominal < selectedp,
                                                          TotalOverlaps == 1)$pairs] <- "Yes"
    
    #remove anything not an eGene at the p thresh (i.e. impossible to be confirmed as CClnc target in this tissue):
    triali <- filter(triali, EnsID.y %in% filter(GTEX_4pairsAll[[i]], 
                                                 pval_nominal < selectedp)$gene_id)
    #remove anything un-expressed >1 TPM in this tissue:
    findExprs <- GTEX_exprs[,colnames(GTEX_exprs)[grepl(gsub("-", "", gsub("_", "", unique(GTEX_4pairs_hiXP$tissueType)))[i], 
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
  
  names(GTEX_eQTL_locusFish) <- unique(GTEX_4pairs_hiXP$tissueType)
  p_test_list[[j]] <- as.data.frame(t(bind_rows(GTEX_eQTL_locusFish, .id = "tissue")))
  colnames(p_test_list[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(p_test_list)
names(p_test_list) <- paste("Run", pThresh_df[,1], sep = "_")
p_test_list_df <- bind_rows(p_test_list, .id = "Run")
p_test_list_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(p_test_list_df))
GTEX_eQTL_SpliceFish_samehiXP_df <- p_test_list_df
#write.csv(GTEX_eQTL_SpliceFish_samehiXP_df, "GTEX_eQTL_SpliceFish_samehiXP_df.csv", row.names = F)
unique(filter(GTEX_eQTL_SpliceFish_samehiXP_df, p <0.05)$tissue)
#1x tissue

#p adjust:
#select best p per Tissue_pair_overlap:
trial <- split(GTEX_eQTL_SpliceFish_samehiXP_df, GTEX_eQTL_SpliceFish_samehiXP_df$tissue)

triali <- lapply(trial, function(x){
  y <- filter(x, c >10)
  y$V1_padj <- p.adjust(y$p, method = "BH")
  return(y[order(y$V1, decreasing = F),][1,])
})

triali <- bind_rows(triali)
unique(filter(triali, V1_padj <0.05)$tissue)
#0x tissues

#
unique(filter(GTEX_eQTL_locusFish_samehiXP_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_promoterFish_samehiXP_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_exonFish_samehiXP_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_IntronFish_samehiXP_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_TTSFish_samehiXP_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_SpliceFish_samehiXP_df, p <0.05)$tissue)


#### 1c 400kbp pairs, same timeframe pairs ####

#bit odd that totally different tissues to previous have been found
#a key diff is the genomic window being much wider @1mbp
#400kbp was the previous, so tried again here
#in theory, less noisy longer distance pairs

#selected p
pThresh_df <- data.frame("pThresh" = c(rep(0.05, 1),
                                       rep(0.00005, 1),
                                       rep(0.00001, 1), 
                                       rep(0.000001, 1),
                                       rep(0.0000001, 1),
                                       rep(0.00000001, 1)))

AllLNC_AllPCG_2d3d_400 <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <400)

CoRegPairs_04_48_24_extendedSame <- filter(AllLNC_AllPCG_2d3d_400,
                                           #AllLNC_AllPCG_1,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                         fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                          fpkm_allGDE_Downwithin_4$EnsID)) |
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                           fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                            fpkm_allGDE_Downwithin_8$EnsID)) |
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                           fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                             fpkm_allGDE_Downwithin_24$EnsID)))
#232

#store outputs here:
GTEX_eQTL_locusFish <- list()
p_test_list <- list()

#run once for each group of lncRNA overlap type
#whole locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400
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
GTEX_eQTL_locusFish_same400_df <- p_test_list_df
#write.csv(GTEX_eQTL_locusFish_same_df, "GTEX_eQTL_locusFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_locusFish_same400_df, p <0.05)$tissue)
#0x tissues (loss of 8x)


#promoter
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400
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
GTEX_eQTL_promoterFish_same400_df <- p_test_list_df
#write.csv(GTEX_eQTL_promoterFish_same_df, "GTEX_eQTL_promoterFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterFish_same400_df, p <0.05)$tissue)
#0x tissues (loss of 4x tissues)


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400
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
GTEX_eQTL_exonFish_same400_df <- p_test_list_df
#write.csv(GTEX_eQTL_exonFish_same_df, "GTEX_eQTL_exonFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_exonFish_same400_df, p <0.05)$tissue)
#1x tissues (loss of 10x)


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400
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
GTEX_eQTL_promoterExonFish_same400_df <- p_test_list_df
#write.csv(GTEX_eQTL_promoterExonFish_same_df, "GTEX_eQTL_promoterExonFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterExonFish_same400_df, p <0.05)$tissue)
#0x tissues (loss of 9x)


#Intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400
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
GTEX_eQTL_IntronFish_same400_df <- p_test_list_df
#write.csv(GTEX_eQTL_IntronFish_same_df, "GTEX_eQTL_IntronFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_IntronFish_same400_df, p <0.05)$tissue)
#0x tissues (loss of 2x)


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400
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
GTEX_eQTL_TTSFish_same400_df <- p_test_list_df
#write.csv(GTEX_eQTL_TTSFish_same_df, "GTEX_eQTL_TTSFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_TTSFish_same400_df, p <0.05)$tissue)
#0x tissues (loss of 7x)

#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400
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
GTEX_eQTL_SpliceFish_same400_df <- p_test_list_df
#write.csv(GTEX_eQTL_SpliceFish_same_df, "GTEX_eQTL_SpliceFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_SpliceFish_same400_df, p <0.05)$tissue)
#0x tissue

#big fail
unique(filter(GTEX_eQTL_locusFish_same400_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_promoterFish_same400_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_exonFish_same400_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_IntronFish_same400_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_TTSFish_same400_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_SpliceFish_same400_df, p <0.05)$tissue)


#### 1d 250kbp pairs, same timeframe pairs ####

#bit odd that totally different tissues to previous have been found
#a key diff is the genomic window being much wider @1mbp
#400kbp was the previous, so tried again here
#in theory, less noisy longer distance pairs

#selected p
pThresh_df <- data.frame("pThresh" = c(rep(0.05, 1),
                                       rep(0.00005, 1),
                                       rep(0.00001, 1), 
                                       rep(0.000001, 1),
                                       rep(0.0000001, 1),
                                       rep(0.00000001, 1)))

AllLNC_AllPCG_2d3d_250 <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <250)

CoRegPairs_04_48_24_extendedSame <- filter(AllLNC_AllPCG_2d3d_250,
                                           #AllLNC_AllPCG_1,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                         fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                          fpkm_allGDE_Downwithin_4$EnsID)) |
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                           fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                            fpkm_allGDE_Downwithin_8$EnsID)) |
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                           fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                             fpkm_allGDE_Downwithin_24$EnsID)))
#153

#store outputs here:
GTEX_eQTL_locusFish <- list()
p_test_list <- list()

#run once for each group of lncRNA overlap type
#whole locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_250
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
GTEX_eQTL_locusFish_same250_df <- p_test_list_df
#write.csv(GTEX_eQTL_locusFish_same_df, "GTEX_eQTL_locusFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_locusFish_same250_df, p <0.05)$tissue)
#8x tissues (but p much weaker)


#promoter
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_250
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
GTEX_eQTL_promoterFish_same250_df <- p_test_list_df
#write.csv(GTEX_eQTL_promoterFish_same_df, "GTEX_eQTL_promoterFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterFish_same250_df, p <0.05)$tissue)
#0x tissues


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_250
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
GTEX_eQTL_exonFish_same250_df <- p_test_list_df
#write.csv(GTEX_eQTL_exonFish_same_df, "GTEX_eQTL_exonFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_exonFish_same250_df, p <0.05)$tissue)
#2x tissues 


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_250
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
GTEX_eQTL_promoterExonFish_same250_df <- p_test_list_df
#write.csv(GTEX_eQTL_promoterExonFish_same_df, "GTEX_eQTL_promoterExonFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterExonFish_same250_df, p <0.05)$tissue)
#1x tissues (loss of 9x)


#Intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_250
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
GTEX_eQTL_IntronFish_same250_df <- p_test_list_df
#write.csv(GTEX_eQTL_IntronFish_same_df, "GTEX_eQTL_IntronFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_IntronFish_same250_df, p <0.05)$tissue)
#7x tissues (loss of 2x)


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_250
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
GTEX_eQTL_TTSFish_same250_df <- p_test_list_df
#write.csv(GTEX_eQTL_TTSFish_same_df, "GTEX_eQTL_TTSFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_TTSFish_same250_df, p <0.05)$tissue)
#0x tissues (loss of 7x)

#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_250
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
GTEX_eQTL_SpliceFish_same250_df <- p_test_list_df
#write.csv(GTEX_eQTL_SpliceFish_same_df, "GTEX_eQTL_SpliceFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_SpliceFish_same250_df, p <0.05)$tissue)
#0x tissue

#seems bit worse
unique(filter(GTEX_eQTL_locusFish_same250_df, p <0.05)$tissue) #8
unique(filter(GTEX_eQTL_promoterFish_same250_df, p <0.05)$tissue) #0
unique(filter(GTEX_eQTL_exonFish_same250_df, p <0.05)$tissue) #2
unique(filter(GTEX_eQTL_IntronFish_same250_df, p <0.05)$tissue) #7
unique(filter(GTEX_eQTL_TTSFish_same250_df, p <0.05)$tissue) #0
unique(filter(GTEX_eQTL_SpliceFish_same250_df, p <0.05)$tissue) #0

#compare to 1mbp
unique(filter(GTEX_eQTL_locusFish_same_df, p <0.05)$tissue) #8
unique(filter(GTEX_eQTL_promoterFish_same_df, p <0.05)$tissue) #4
unique(filter(GTEX_eQTL_exonFish_same_df, p <0.05)$tissue) #10
unique(filter(GTEX_eQTL_IntronFish_same_df, p <0.05)$tissue) #2
unique(filter(GTEX_eQTL_TTSFish_same_df, p <0.05)$tissue) #7
unique(filter(GTEX_eQTL_SpliceFish_same_df, p <0.05)$tissue) #1


#### 2a delayed timeframe pairs ####

#selected p
pThresh_df <- data.frame("pThresh" = c(rep(0.05, 1),
                                       rep(0.00005, 1),
                                       rep(0.00001, 1), 
                                       rep(0.000001, 1),
                                       rep(0.0000001, 1),
                                       rep(0.00000001, 1)))

#are eQTL-eGene connections predictive of these pairs:
CoRegPairs_04_48_24_extendedDelayed <- filter(AllLNC_AllPCG_2d3d,
                                           #AllLNC_AllPCG_1,
                                           (EnsID %in% fpkm_allGDE_within_4$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                          fpkm_allGDE_within_24$EnsID)) |
                                             (EnsID %in% fpkm_allGDE_within_8$EnsID & EnsID.y %in% fpkm_allGDE_within_24$EnsID))
#506 pairs

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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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
GTEX_eQTL_locusFish_del_df <- p_test_list_df
#write.csv(GTEX_eQTL_locusFish_del_df, "GTEX_eQTL_locusFish_del_df.csv", row.names = F)
unique(filter(GTEX_eQTL_locusFish_del_df, p <0.05)$tissue)
#0x tissues

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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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
GTEX_eQTL_promoterFish_del_df <- p_test_list_df
#write.csv(GTEX_eQTL_promoterFish_del_df, "GTEX_eQTL_promoterFish_del_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterFish_del_df, p <0.05)$tissue)
#0x tissues

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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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
GTEX_eQTL_exonFish_del_df <- p_test_list_df
#write.csv(GTEX_eQTL_exonFish_del_df, "GTEX_eQTL_exonFish_del_df.csv", row.names = F)
unique(filter(GTEX_eQTL_exonFish_del_df, p <0.05)$tissue)
#0x tissues

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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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
GTEX_eQTL_promoterExonFish_del_df <- p_test_list_df
#write.csv(GTEX_eQTL_promoterExonFish_del_df, "GTEX_eQTL_promoterExonFish_del_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterExonFish_del_df, p <0.05)$tissue)
#0x tissues

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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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
GTEX_eQTL_IntronFish_del_df <- p_test_list_df
#write.csv(GTEX_eQTL_IntronFish_del_df, "GTEX_eQTL_IntronFish_del_df.csv", row.names = F)
unique(filter(GTEX_eQTL_IntronFish_del_df, p <0.05)$tissue)
#0x tissues

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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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
GTEX_eQTL_TTSFish_del_df <- p_test_list_df
#write.csv(GTEX_eQTL_TTSFish_del_df, "GTEX_eQTL_TTSFish_del_df.csv", row.names = F)
unique(filter(GTEX_eQTL_TTSFish_del_df, p <0.05)$tissue)
#0x tissues

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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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
GTEX_eQTL_SpliceFish_del_df <- p_test_list_df
#write.csv(GTEX_eQTL_SpliceFish_del_df, "GTEX_eQTL_SpliceFish_del_df.csv", row.names = F)
unique(filter(GTEX_eQTL_SpliceFish_del_df, p <0.05)$tissue)
#0x tissue, weak, not done previously

#despite a ~25% bigger pool of 506 pairs, no enrichment of eQTL support

#big fail
unique(filter(GTEX_eQTL_locusFish_del_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_promoterFish_del_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_exonFish_del_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_IntronFish_del_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_TTSFish_del_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_SpliceFish_del_df, p <0.05)$tissue)


#### 2b delayed timeframe pairs 250 ####

#selected p
pThresh_df <- data.frame("pThresh" = c(rep(0.05, 1),
                                       rep(0.00005, 1),
                                       rep(0.00001, 1), 
                                       rep(0.000001, 1),
                                       rep(0.0000001, 1),
                                       rep(0.00000001, 1)))

#are eQTL-eGene connections predictive of these pairs:
CoRegPairs_04_48_24_extendedDelayed <- filter(AllLNC_AllPCG_2d3d_250,
                                              #AllLNC_AllPCG_1,
                                              (EnsID %in% fpkm_allGDE_within_4$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                      fpkm_allGDE_within_24$EnsID)) |
                                                (EnsID %in% fpkm_allGDE_within_8$EnsID & EnsID.y %in% fpkm_allGDE_within_24$EnsID))
#506 pairs

#store outputs here:
GTEX_eQTL_locusFish <- list()
p_test_list <- list()

#run once for each group of lncRNA overlap type
#whole locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_250
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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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
GTEX_eQTL_locusFish_del250_df <- p_test_list_df
#write.csv(GTEX_eQTL_locusFish_del250_df, "GTEX_eQTL_locusFish_del250_df.csv", row.names = F)
unique(filter(GTEX_eQTL_locusFish_del250_df, p <0.05)$tissue)
#0x tissues

#promoter
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_250
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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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
GTEX_eQTL_promoterFish_del250_df <- p_test_list_df
#write.csv(GTEX_eQTL_promoterFish_del250_df, "GTEX_eQTL_promoterFish_del250_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterFish_del250_df, p <0.05)$tissue)
#0x tissues

#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_250
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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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
GTEX_eQTL_exonFish_del250_df <- p_test_list_df
#write.csv(GTEX_eQTL_exonFish_del250_df, "GTEX_eQTL_exonFish_del250_df.csv", row.names = F)
unique(filter(GTEX_eQTL_exonFish_del250_df, p <0.05)$tissue)
#0x tissues

#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_250
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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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
GTEX_eQTL_promoterExonFish_del250_df <- p_test_list_df
#write.csv(GTEX_eQTL_promoterExonFish_del250_df, "GTEX_eQTL_promoterExonFish_del250_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterExonFish_del250_df, p <0.05)$tissue)
#0x tissues

#Intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_250
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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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
GTEX_eQTL_IntronFish_del250_df <- p_test_list_df
#write.csv(GTEX_eQTL_IntronFish_del250_df, "GTEX_eQTL_IntronFish_del250_df.csv", row.names = F)
unique(filter(GTEX_eQTL_IntronFish_del250_df, p <0.05)$tissue)
#0x tissues

#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_250
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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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
GTEX_eQTL_TTSFish_del250_df <- p_test_list_df
#write.csv(GTEX_eQTL_TTSFish_del250_df, "GTEX_eQTL_TTSFish_del250_df.csv", row.names = F)
unique(filter(GTEX_eQTL_TTSFish_del250_df, p <0.05)$tissue)
#0x tissues

#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_250
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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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
GTEX_eQTL_SpliceFish_del250_df <- p_test_list_df
#write.csv(GTEX_eQTL_SpliceFish_del250_df, "GTEX_eQTL_SpliceFish_del250_df.csv", row.names = F)
unique(filter(GTEX_eQTL_SpliceFish_del250_df, p <0.05)$tissue)
#0x tissue, weak, not done previously

#despite a ~25% bigger pool of 506 pairs, no enrichment of eQTL support

#big fail
unique(filter(GTEX_eQTL_locusFish_del250_df, p <0.05)$tissue) #1
unique(filter(GTEX_eQTL_promoterFish_del250_df, p <0.05)$tissue) #1
unique(filter(GTEX_eQTL_exonFish_del250_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_IntronFish_del250_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_TTSFish_del250_df, p <0.05)$tissue)
unique(filter(GTEX_eQTL_SpliceFish_del250_df, p <0.05)$tissue)



#### 3a 0-4hr induced same timeframe pairs ####

#selected p
pThresh_df <- data.frame("pThresh" = c(rep(0.05, 1),
                                       rep(0.00005, 1),
                                       rep(0.00001, 1), 
                                       rep(0.000001, 1),
                                       rep(0.0000001, 1),
                                       rep(0.00000001, 1)))

#filter to 0-4hr expressed genes:
fpkm_allG_04 <- filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)
AllLNC_AllPCG_2d3d_04 <- filter(AllLNC_AllPCG_2d3d, EnsID %in% fpkm_allG_04$EnsID, EnsID.y %in% fpkm_allG_04$EnsID)

dim(filter(AllLNC_AllPCG_2d3d_04, !loopMethod == "Neither"))
223/7222
length(unique(AllLNC_AllPCG_2d3d_04$EnsID.y))#4611 genes returned

CoRegPairs_04_48_24_extendedSame <- filter(AllLNC_AllPCG_2d3d_04,
                                           #AllLNC_AllPCG_1,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                        fpkm_allGDE_Downwithin_4$EnsID)))

#167 testable

#considering range? are co-regs just closer together and so more likely to get eQTLs?
#no great diffs in ranges of pairs
summary(AllLNC_AllPCG_2d3d_04$pair_range)
summary(CoRegPairs_04_48_24_extendedSame$pair_range)
boxplot(AllLNC_AllPCG_2d3d_04$pair_range, CoRegPairs_04_48_24_extendedSame$pair_range)

#slight more obvious diff in TSS distance of pairs
summary(AllLNC_AllPCG_2d3d_04$AbsDistLnc_PCG)
summary(CoRegPairs_04_48_24_extendedSame$AbsDistLnc_PCG)
boxplot(AllLNC_AllPCG_2d3d_04$AbsDistLnc_PCG, CoRegPairs_04_48_24_extendedSame$AbsDistLnc_PCG)

#so some of effect could be simply that the coregs are closer together? but seems v slight


#store outputs here:
GTEX_eQTL_locusFish <- list()
p_test_list <- list()

#run once for each group of lncRNA overlap type
#whole locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_04
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
GTEX_eQTL_locusFish_same_df_04up <- p_test_list_df
#write.csv(GTEX_eQTL_locusFish_same_df, "GTEX_eQTL_locusFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_locusFish_same_df_04up, p <0.05)$tissue)

#massive increase now, in ORs
#note reducing the amount of pair massively (e.g. 733 brain frontal, to 215 now)
#by focusing on upreg lncs rather than all lncs, early timeframe PCGs only
#good message in there
#tease out by repeating with other clusters - unlikely to get same effect

#promoter
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_04
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
GTEX_eQTL_promoterFish_same_df_04up <- p_test_list_df
#write.csv(GTEX_eQTL_promoterFish_same_df, "GTEX_eQTL_promoterFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterFish_same_df_04up, p <0.05)$tissue)

#again, big increases

#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_04
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
GTEX_eQTL_exonFish_same_df_04up <- p_test_list_df
#write.csv(GTEX_eQTL_exonFish_same_df, "GTEX_eQTL_exonFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_exonFish_same_df_04up, p <0.05)$tissue)
#10x tissues


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_04
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
GTEX_eQTL_promoterExonFish_same_df_04up <- p_test_list_df
#write.csv(GTEX_eQTL_promoterExonFish_same_df, "GTEX_eQTL_promoterExonFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterExonFish_same_df_04up, p <0.05)$tissue)
#9x tissues


#Intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_04
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
GTEX_eQTL_IntronFish_same_df_04up <- p_test_list_df
#write.csv(GTEX_eQTL_IntronFish_same_df, "GTEX_eQTL_IntronFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_IntronFish_same_df_04up, p <0.05)$tissue)
#2x tissues


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_04
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
GTEX_eQTL_TTSFish_same_df_04up <- p_test_list_df
#write.csv(GTEX_eQTL_TTSFish_same_df, "GTEX_eQTL_TTSFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_TTSFish_same_df_04up, p <0.05)$tissue)
#7x tissues


#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_04
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
GTEX_eQTL_SpliceFish_same_df_df_04up <- p_test_list_df
#write.csv(GTEX_eQTL_SpliceFish_same_df, "GTEX_eQTL_SpliceFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_SpliceFish_same_df_df_04up, p <0.05)$tissue)
#1x tissue

#big increases, significance with big ORs for:

#whole locus, promoter, exon
#intron, splice, TTS too but seem weaker


#### 4a 0-4hr repressed same timeframe pairs ####

#selected p
pThresh_df <- data.frame("pThresh" = c(rep(0.05, 1),
                                       rep(0.00005, 1),
                                       rep(0.00001, 1), 
                                       rep(0.000001, 1),
                                       rep(0.0000001, 1),
                                       rep(0.00000001, 1)))

#filter to 0-4hr expressed genes:
fpkm_allG_04 <- filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)
AllLNC_AllPCG_2d3d_04 <- filter(AllLNC_AllPCG_2d3d, EnsID %in% fpkm_allG_04$EnsID, EnsID.y %in% fpkm_allG_04$EnsID)

dim(filter(AllLNC_AllPCG_2d3d_04, !loopMethod == "Neither"))
223/7222
length(unique(AllLNC_AllPCG_2d3d_04$EnsID.y))#4611 genes returned

CoRegPairs_04_48_24_extendedSame <- filter(AllLNC_AllPCG_2d3d_04,
                                           #AllLNC_AllPCG_1,
                                           (EnsID %in% c(fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                        fpkm_allGDE_Downwithin_4$EnsID)))
#86 testable

#considering range? are co-regs just closer together and so more likely to get eQTLs?
#no great diffs in ranges of pairs
summary(AllLNC_AllPCG_2d3d_04$pair_range)
summary(CoRegPairs_04_48_24_extendedSame$pair_range)
boxplot(AllLNC_AllPCG_2d3d_04$pair_range, CoRegPairs_04_48_24_extendedSame$pair_range)

#slight more obvious diff in TSS distance of pairs
summary(AllLNC_AllPCG_2d3d_04$AbsDistLnc_PCG)
summary(CoRegPairs_04_48_24_extendedSame$AbsDistLnc_PCG)
boxplot(AllLNC_AllPCG_2d3d_04$AbsDistLnc_PCG, CoRegPairs_04_48_24_extendedSame$AbsDistLnc_PCG)

#so some of effect could be simply that the coregs are closer together? but seems v slight

#store outputs here:
GTEX_eQTL_locusFish <- list()
p_test_list <- list()

#run once for each group of lncRNA overlap type
#whole locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_04
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
GTEX_eQTL_locusFish_same_df_04down <- p_test_list_df
#write.csv(GTEX_eQTL_locusFish_same_df, "GTEX_eQTL_locusFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_locusFish_same_df_04down, p <0.05)$tissue)


#promoter
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_04
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
GTEX_eQTL_promoterFish_same_df_04down <- p_test_list_df
#write.csv(GTEX_eQTL_promoterFish_same_df, "GTEX_eQTL_promoterFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterFish_same_df_04down, p <0.05)$tissue)


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_04
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
GTEX_eQTL_exonFish_same_df_04down <- p_test_list_df
#write.csv(GTEX_eQTL_exonFish_same_df, "GTEX_eQTL_exonFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_exonFish_same_df_04down, p <0.05)$tissue)


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_04
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
GTEX_eQTL_promoterExonFish_same_df_04down <- p_test_list_df
#write.csv(GTEX_eQTL_promoterExonFish_same_df, "GTEX_eQTL_promoterExonFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterExonFish_same_df_04down, p <0.05)$tissue)


#Intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_04
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
GTEX_eQTL_IntronFish_same_df_04down <- p_test_list_df
#write.csv(GTEX_eQTL_IntronFish_same_df, "GTEX_eQTL_IntronFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_IntronFish_same_df_04down, p <0.05)$tissue)


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_04
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
GTEX_eQTL_TTSFish_same_df_04down <- p_test_list_df
#write.csv(GTEX_eQTL_TTSFish_same_df, "GTEX_eQTL_TTSFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_TTSFish_same_df_04down, p <0.05)$tissue)


#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_04
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
GTEX_eQTL_SpliceFish_same_df_04down <- p_test_list_df
#write.csv(GTEX_eQTL_SpliceFish_same_df, "GTEX_eQTL_SpliceFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_SpliceFish_same_df_04down, p <0.05)$tissue)


#v poor confirmation
#there are very few confirmed pairs full stop (nvm enrichment)


#### takeaway from above ####

#by far biggest jump in predictive value is by splitting the lncRNAs

#and focusing on induced

#there is no big improvement with eQTL p thresholding (but save all for reference)

#delayed looks useless (make a comparison to 0-4hr induced)

#distance change is interesting, but leave out, too complex for this approach


#### build a supplementary table ####

#more effective than running numerous times
#need to set it to run through the tissues
#for each of the lncRNA overlap types
#for each lncRNA cluster
GTEX_runs_parameters <- data.frame("Lnc_cluster" = c("Up4", "Up8", "Up24", "Down4", "Down8", "Down24", 
                                                     "Up4Del", "Down4Del", "All-same", "All-delayed"),
                                   "Background" = c("0-4", "0-8", "0-24", "0-4", "0-8", "0-24", 
                                                    "0-24", "0-24", "0-24", "0-24"),
                                   "Overlap_type" = c(rep("Promoter|Exon|Intron|Splice|TTS",10), 
                                                      rep("Promoter",10),
                                                      rep("Exon",10),
                                                      rep("Intron",10),
                                                      rep("Splice",10),
                                                      rep("TTS",10)))

#background genes
fpkm_allG_04 <- filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)
fpkm_allG_08 <- filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8 | Hour8_meanFPKM >0.8)

GTEX_runs_background <- list("0-4" = fpkm_allG_04, "0-8" = fpkm_allG_08, "0-24" = fpkm_allG)

#co-reg genes
GTEX_runs_selection <- list("Up4" = filter(AllLNC_AllPCG_2d3d,
                                            (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                         fpkm_allGDE_Downwithin_4$EnsID))),
                             "Down4" = filter(AllLNC_AllPCG_2d3d,
                                            (EnsID %in% c(fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                         fpkm_allGDE_Downwithin_4$EnsID))),
                             "Up8" = filter(AllLNC_AllPCG_2d3d,
                                            (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                         fpkm_allGDE_Downwithin_8$EnsID))),
                             "Down8" = filter(AllLNC_AllPCG_2d3d,
                                            (EnsID %in% c(fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                         fpkm_allGDE_Downwithin_8$EnsID))),
                             "Up24" = filter(AllLNC_AllPCG_2d3d,
                                            (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                         fpkm_allGDE_Downwithin_24$EnsID))),
                             "Down24" = filter(AllLNC_AllPCG_2d3d,
                                            (EnsID %in% c(fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                         fpkm_allGDE_Downwithin_24$EnsID))),
                            "Up4Del" = filter(AllLNC_AllPCG_2d3d,
                                              (EnsID %in% fpkm_allGDE_Upwithin_4$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                        fpkm_allGDE_within_24$EnsID))),
                            "Down4Del" = filter(AllLNC_AllPCG_2d3d,
                                                (EnsID %in% fpkm_allGDE_Downwithin_4$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                          fpkm_allGDE_within_24$EnsID))),
                            "All-same" = filter(AllLNC_AllPCG_2d3d,
                                                (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                              fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                               fpkm_allGDE_Downwithin_4$EnsID)) |
                                                  (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                                 fpkm_allGDE_Downwithin_8$EnsID)) |
                                                  (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                                  fpkm_allGDE_Downwithin_24$EnsID))),
                            "All-same/delayed"= filter(AllLNC_AllPCG_2d3d,
                                                       (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                     fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE$EnsID)) |
                                                         (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                       fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                                        fpkm_allGDE_within_24$EnsID)) |
                                                         (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                       fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                                         fpkm_allGDE_Downwithin_24$EnsID)))
                            )

GTEX_SuppTable_TissueRes <- list()
GTEX_SuppTable <- list()

selectedp <- 0.05

for (j in 1:length(GTEX_runs_parameters$Lnc_cluster)){
  
  selectedLncs <- GTEX_runs_selection[[ GTEX_runs_parameters[j,1] ]]
  selectedBackground <- GTEX_runs_background[[ GTEX_runs_parameters[j,2] ]]
  selectedOverlap <- GTEX_runs_parameters[j,3]
    
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    
    triali <- filter(AllLNC_AllPCG_2d3d, EnsID %in% selectedBackground$EnsID, EnsID.y %in% selectedBackground$EnsID)
    
    #option to test without 27167 and 27169:
    #triali <- filter(triali, grepl("MSTRG.27167|MSTRG.27169", EnsID))
    
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          grepl(selectedOverlap, OverlapType), 
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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% selectedLncs$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% selectedLncs$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_SuppTable_TissueRes[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_SuppTable_TissueRes[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_SuppTable_TissueRes) <- unique(GTEX_4pairs$tissueType)
  GTEX_SuppTable[[j]] <- as.data.frame(t(bind_rows(GTEX_SuppTable_TissueRes, .id = "tissue")))
  colnames(GTEX_SuppTable[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(GTEX_SuppTable)
GTEX_runs_parameters$Overlap_type_name <- GTEX_runs_parameters$Overlap_type
GTEX_runs_parameters$Overlap_type_name[1:8] <- "Locus"
names(GTEX_SuppTable) <- paste("Run", GTEX_runs_parameters[,1], GTEX_runs_parameters[,4], sep = "_")
GTEX_SuppTable_df <- bind_rows(GTEX_SuppTable, .id = "Run")
GTEX_SuppTable_df$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(GTEX_SuppTable_df))
rownames(GTEX_SuppTable_df) <- NULL

#write.csv(GTEX_SuppTable_df, "GTEX_SuppTable_df.csv", row.names = F)


#### build a supplementary table - 500kbp ####

#would this limit zscan locus?
dim(filter(AllLNC_AllPCG_2d3d, grepl("MSTRG.27167|MSTRG.27169", EnsID)))#28 pairs

AllLNC_AllPCG_2d3d_250 <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <250)
dim(filter(AllLNC_AllPCG_2d3d_250, grepl("MSTRG.27167|MSTRG.27169", EnsID)))#3 pairs only

AllLNC_AllPCG_2d3d_500 <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <500)
dim(filter(AllLNC_AllPCG_2d3d_500, grepl("MSTRG.27167|MSTRG.27169", EnsID)))#8 pairs

#more effective than running numerous times
#need to set it to run through the tissues
#for each of the lncRNA overlap types
#for each lncRNA cluster
GTEX_runs_parameters <- data.frame("Lnc_cluster" = c("Up4", "Up8", "Up24", "Down4", "Down8", "Down24", 
                                                     "Up4Del", "Down4Del", "All-same", "All-delayed"),
                                   "Background" = c("0-4", "0-8", "0-24", "0-4", "0-8", "0-24", 
                                                    "0-24", "0-24", "0-24", "0-24"),
                                   "Overlap_type" = c(rep("Promoter|Exon|Intron|Splice|TTS",10), 
                                                      rep("Promoter",10),
                                                      rep("Exon",10),
                                                      rep("Intron",10),
                                                      rep("Splice",10),
                                                      rep("TTS",10)))

#background genes
fpkm_allG_04 <- filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)
fpkm_allG_08 <- filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8 | Hour8_meanFPKM >0.8)

GTEX_runs_background <- list("0-4" = fpkm_allG_04, "0-8" = fpkm_allG_08, "0-24" = fpkm_allG)

#co-reg genes
GTEX_runs_selection <- list("Up4" = filter(AllLNC_AllPCG_2d3d_500,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                        fpkm_allGDE_Downwithin_4$EnsID))),
                            "Down4" = filter(AllLNC_AllPCG_2d3d_500,
                                             (EnsID %in% c(fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                            fpkm_allGDE_Downwithin_4$EnsID))),
                            "Up8" = filter(AllLNC_AllPCG_2d3d_500,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                        fpkm_allGDE_Downwithin_8$EnsID))),
                            "Down8" = filter(AllLNC_AllPCG_2d3d_500,
                                             (EnsID %in% c(fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                            fpkm_allGDE_Downwithin_8$EnsID))),
                            "Up24" = filter(AllLNC_AllPCG_2d3d_500,
                                            (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                          fpkm_allGDE_Downwithin_24$EnsID))),
                            "Down24" = filter(AllLNC_AllPCG_2d3d_500,
                                              (EnsID %in% c(fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                              fpkm_allGDE_Downwithin_24$EnsID))),
                            "Up4Del" = filter(AllLNC_AllPCG_2d3d_500,
                                              (EnsID %in% fpkm_allGDE_Upwithin_4$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                        fpkm_allGDE_within_24$EnsID))),
                            "Down4Del" = filter(AllLNC_AllPCG_2d3d_500,
                                                (EnsID %in% fpkm_allGDE_Downwithin_4$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                            fpkm_allGDE_within_24$EnsID))),
                            "All-same" = filter(AllLNC_AllPCG_2d3d,
                                                (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                              fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                               fpkm_allGDE_Downwithin_4$EnsID)) |
                                                  (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                                 fpkm_allGDE_Downwithin_8$EnsID)) |
                                                  (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                                  fpkm_allGDE_Downwithin_24$EnsID))),
                            "All-same/delayed"= filter(AllLNC_AllPCG_2d3d,
                                                       (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                     fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE$EnsID)) |
                                                         (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                       fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                                        fpkm_allGDE_within_24$EnsID)) |
                                                         (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                       fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                                         fpkm_allGDE_Downwithin_24$EnsID)))
                            )

GTEX_SuppTable_TissueRes <- list()
GTEX_SuppTable <- list()

selectedp <- 0.05

for (j in 1:length(GTEX_runs_parameters$Lnc_cluster)){
  
  selectedLncs <- GTEX_runs_selection[[ GTEX_runs_parameters[j,1] ]]
  selectedBackground <- GTEX_runs_background[[ GTEX_runs_parameters[j,2] ]]
  selectedOverlap <- GTEX_runs_parameters[j,3]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    
    triali <- filter(AllLNC_AllPCG_2d3d_500, EnsID %in% selectedBackground$EnsID, EnsID.y %in% selectedBackground$EnsID)
    
    #option to test without 27167 and 27169:
    #triali <- filter(triali, !grepl("MSTRG.27167|MSTRG.27169", EnsID))
    
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          grepl(selectedOverlap, OverlapType), 
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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% selectedLncs$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% selectedLncs$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_SuppTable_TissueRes[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_SuppTable_TissueRes[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_SuppTable_TissueRes) <- unique(GTEX_4pairs$tissueType)
  GTEX_SuppTable[[j]] <- as.data.frame(t(bind_rows(GTEX_SuppTable_TissueRes, .id = "tissue")))
  colnames(GTEX_SuppTable[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(GTEX_SuppTable)
GTEX_runs_parameters$Overlap_type_name <- GTEX_runs_parameters$Overlap_type
GTEX_runs_parameters$Overlap_type_name[1:8] <- "Locus"
names(GTEX_SuppTable) <- paste("Run", GTEX_runs_parameters[,1], GTEX_runs_parameters[,4], sep = "_")
GTEX_SuppTable_df500 <- bind_rows(GTEX_SuppTable, .id = "Run")
GTEX_SuppTable_df500$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(GTEX_SuppTable_df500))
rownames(GTEX_SuppTable_df500) <- NULL

#much poorer results, only 1 for the up in 4
#check esoph-gas:
#sig is depleted but not completely lost without the zscans

#write.csv(GTEX_SuppTable_df, "GTEX_SuppTable_df.csv", row.names = F)



#### build a supplementary table - 250kbp ####

#would this limit zscan locus?
dim(filter(AllLNC_AllPCG_2d3d, grepl("MSTRG.27167|MSTRG.27169", EnsID)))#28 pairs

AllLNC_AllPCG_2d3d_250 <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <250)
dim(filter(AllLNC_AllPCG_2d3d_250, grepl("MSTRG.27167|MSTRG.27169", EnsID)))#3 pairs only

AllLNC_AllPCG_2d3d_500 <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <500)
dim(filter(AllLNC_AllPCG_2d3d_500, grepl("MSTRG.27167|MSTRG.27169", EnsID)))#8 pairs

#try 250kbp, stringent
dim(filter(AllLNC_AllPCG_2d3d_250, !loopMethod == "Neither"))
185/2563

#more effective than running numerous times
#need to set it to run through the tissues
#for each of the lncRNA overlap types
#for each lncRNA cluster
GTEX_runs_parameters <- data.frame("Lnc_cluster" = c("Up4", "Up8", "Up24", "Down4", "Down8", "Down24", 
                                                     "Up4Del", "Down4Del", "All-same", "All-delayed"),
                                   "Background" = c("0-4", "0-8", "0-24", "0-4", "0-8", "0-24", 
                                                    "0-24", "0-24", "0-24", "0-24"),
                                   "Overlap_type" = c(rep("Promoter|Exon|Intron|Splice|TTS",10), 
                                                      rep("Promoter",10),
                                                      rep("Exon",10),
                                                      rep("Intron",10),
                                                      rep("Splice",10),
                                                      rep("TTS",10)))

#background genes
fpkm_allG_04 <- filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)
fpkm_allG_08 <- filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8 | Hour8_meanFPKM >0.8)

GTEX_runs_background <- list("0-4" = fpkm_allG_04, "0-8" = fpkm_allG_08, "0-24" = fpkm_allG)

#co-reg genes
GTEX_runs_selection <- list("Up4" = filter(AllLNC_AllPCG_2d3d_250,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                        fpkm_allGDE_Downwithin_4$EnsID))),
                            "Down4" = filter(AllLNC_AllPCG_2d3d_250,
                                             (EnsID %in% c(fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                            fpkm_allGDE_Downwithin_4$EnsID))),
                            "Up8" = filter(AllLNC_AllPCG_2d3d_250,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                        fpkm_allGDE_Downwithin_8$EnsID))),
                            "Down8" = filter(AllLNC_AllPCG_2d3d_250,
                                             (EnsID %in% c(fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                            fpkm_allGDE_Downwithin_8$EnsID))),
                            "Up24" = filter(AllLNC_AllPCG_2d3d_250,
                                            (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                          fpkm_allGDE_Downwithin_24$EnsID))),
                            "Down24" = filter(AllLNC_AllPCG_2d3d_250,
                                              (EnsID %in% c(fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                              fpkm_allGDE_Downwithin_24$EnsID))),
                            "Up4Del" = filter(AllLNC_AllPCG_2d3d_250,
                                              (EnsID %in% fpkm_allGDE_Upwithin_4$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                        fpkm_allGDE_within_24$EnsID))),
                            "Down4Del" = filter(AllLNC_AllPCG_2d3d_250,
                                                (EnsID %in% fpkm_allGDE_Downwithin_4$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                            fpkm_allGDE_within_24$EnsID))),
                            "All-same" = filter(AllLNC_AllPCG_2d3d_250,
                                                (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                              fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                               fpkm_allGDE_Downwithin_4$EnsID)) |
                                                  (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                                 fpkm_allGDE_Downwithin_8$EnsID)) |
                                                  (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                                  fpkm_allGDE_Downwithin_24$EnsID))),
                            "All-delayed"= filter(AllLNC_AllPCG_2d3d_250,
                                                       (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                     fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                                      fpkm_allGDE_within_24$EnsID)) |
                                                         (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                       fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_within_24$EnsID)))
                            )

GTEX_SuppTable_TissueRes <- list()
GTEX_SuppTable <- list()

selectedp <- 0.05

for (j in 1:length(GTEX_runs_parameters$Lnc_cluster)){
  
  selectedLncs <- GTEX_runs_selection[[ GTEX_runs_parameters[j,1] ]]
  selectedBackground <- GTEX_runs_background[[ GTEX_runs_parameters[j,2] ]]
  selectedOverlap <- GTEX_runs_parameters[j,3]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    
    triali <- filter(AllLNC_AllPCG_2d3d_250, EnsID %in% selectedBackground$EnsID, EnsID.y %in% selectedBackground$EnsID)
    
    #option to test without 27167 and 27169:
    #triali <- filter(triali, !grepl("MSTRG.27167|MSTRG.27169", EnsID))
    
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          grepl(selectedOverlap, OverlapType), 
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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% selectedLncs$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% selectedLncs$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_SuppTable_TissueRes[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_SuppTable_TissueRes[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_SuppTable_TissueRes) <- unique(GTEX_4pairs$tissueType)
  GTEX_SuppTable[[j]] <- as.data.frame(t(bind_rows(GTEX_SuppTable_TissueRes, .id = "tissue")))
  colnames(GTEX_SuppTable[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(GTEX_SuppTable)
GTEX_runs_parameters$Overlap_type_name <- GTEX_runs_parameters$Overlap_type
GTEX_runs_parameters$Overlap_type_name[1:10] <- "Locus"
names(GTEX_SuppTable) <- paste("Run", GTEX_runs_parameters[,1], GTEX_runs_parameters[,4], sep = "_")
GTEX_SuppTable_df250i <- bind_rows(GTEX_SuppTable, .id = "Run")
GTEX_SuppTable_df250$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(GTEX_SuppTable_df250))
rownames(GTEX_SuppTable_df250) <- NULL

#write.csv(GTEX_SuppTable_df250, "GTEX_SuppTable_df250.csv", row.names = F)


#
#### build a supplementary table - 250kbp with eQTL p = 0.00001 ####

#would this limit zscan locus?
dim(filter(AllLNC_AllPCG_2d3d, grepl("MSTRG.27167|MSTRG.27169", EnsID)))#28 pairs

AllLNC_AllPCG_2d3d_250 <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <250)
dim(filter(AllLNC_AllPCG_2d3d_250, grepl("MSTRG.27167|MSTRG.27169", EnsID)))#3 pairs only

AllLNC_AllPCG_2d3d_500 <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <500)
dim(filter(AllLNC_AllPCG_2d3d_500, grepl("MSTRG.27167|MSTRG.27169", EnsID)))#8 pairs

#more effective than running numerous times
#need to set it to run through the tissues
#for each of the lncRNA overlap types
#for each lncRNA cluster
GTEX_runs_parameters <- data.frame("Lnc_cluster" = c("Up4", "Up8", "Up24", "Down4", "Down8", "Down24", 
                                                     "Up4Del", "Down4Del", "All-same", "All-delayed"),
                                   "Background" = c("0-4", "0-8", "0-24", "0-4", "0-8", "0-24", 
                                                    "0-24", "0-24", "0-24", "0-24"),
                                   "Overlap_type" = c(rep("Promoter|Exon|Intron|Splice|TTS",10), 
                                                      rep("Promoter",10),
                                                      rep("Exon",10),
                                                      rep("Intron",10),
                                                      rep("Splice",10),
                                                      rep("TTS",10)))

#background genes
fpkm_allG_04 <- filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)
fpkm_allG_08 <- filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8 | Hour8_meanFPKM >0.8)

GTEX_runs_background <- list("0-4" = fpkm_allG_04, "0-8" = fpkm_allG_08, "0-24" = fpkm_allG)

#co-reg genes
GTEX_runs_selection <- list("Up4" = filter(AllLNC_AllPCG_2d3d_250,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                        fpkm_allGDE_Downwithin_4$EnsID))),
                            "Down4" = filter(AllLNC_AllPCG_2d3d_250,
                                             (EnsID %in% c(fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                            fpkm_allGDE_Downwithin_4$EnsID))),
                            "Up8" = filter(AllLNC_AllPCG_2d3d_250,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                        fpkm_allGDE_Downwithin_8$EnsID))),
                            "Down8" = filter(AllLNC_AllPCG_2d3d_250,
                                             (EnsID %in% c(fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                            fpkm_allGDE_Downwithin_8$EnsID))),
                            "Up24" = filter(AllLNC_AllPCG_2d3d_250,
                                            (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                          fpkm_allGDE_Downwithin_24$EnsID))),
                            "Down24" = filter(AllLNC_AllPCG_2d3d_250,
                                              (EnsID %in% c(fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                              fpkm_allGDE_Downwithin_24$EnsID))),
                            "Up4Del" = filter(AllLNC_AllPCG_2d3d_250,
                                              (EnsID %in% fpkm_allGDE_Upwithin_4$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                        fpkm_allGDE_within_24$EnsID))),
                            "Down4Del" = filter(AllLNC_AllPCG_2d3d_250,
                                                (EnsID %in% fpkm_allGDE_Downwithin_4$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                            fpkm_allGDE_within_24$EnsID))),
                            "All-same" = filter(AllLNC_AllPCG_2d3d,
                                                (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                              fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                               fpkm_allGDE_Downwithin_4$EnsID)) |
                                                  (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                                 fpkm_allGDE_Downwithin_8$EnsID)) |
                                                  (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                                  fpkm_allGDE_Downwithin_24$EnsID))),
                            "All-same/delayed"= filter(AllLNC_AllPCG_2d3d,
                                                       (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                     fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE$EnsID)) |
                                                         (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                       fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                                        fpkm_allGDE_within_24$EnsID)) |
                                                         (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                       fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                                         fpkm_allGDE_Downwithin_24$EnsID)))
)

GTEX_SuppTable_TissueRes <- list()
GTEX_SuppTable <- list()

selectedp <- 0.00001

for (j in 1:length(GTEX_runs_parameters$Lnc_cluster)){
  
  selectedLncs <- GTEX_runs_selection[[ GTEX_runs_parameters[j,1] ]]
  selectedBackground <- GTEX_runs_background[[ GTEX_runs_parameters[j,2] ]]
  selectedOverlap <- GTEX_runs_parameters[j,3]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    
    triali <- filter(AllLNC_AllPCG_2d3d_250, EnsID %in% selectedBackground$EnsID, EnsID.y %in% selectedBackground$EnsID)
    
    #option to test without 27167 and 27169:
    #triali <- filter(triali, !grepl("MSTRG.27167|MSTRG.27169", EnsID))
    
    triali$eQTL_validated_tissue <- "No"
    triali$eQTL_validated_tissue[triali$pairs %in% filter(GTEX_4pairs, grepl(unique(GTEX_4pairs$tissueType)[i], tissueType), 
                                                          grepl(selectedOverlap, OverlapType), 
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
    DELNC_DEPCG_1 <- filter(triali, pairs %in% selectedLncs$pairs)
    DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
    #other neighbours of cisLncs
    DELNC_PCG_1 <- filter(triali, EnsID %in% selectedLncs$EnsID)
    DELNC_PCG_1_eQTL <- filter(DELNC_PCG_1, eQTL_validated_tissue == "Yes")
    
    if (dim(DELNC_DEPCG_1_eQTL)[1] >0) {
      
      a <- dim(DELNC_DEPCG_1_eQTL)[1]
      b <- dim(DELNC_DEPCG_1)[1]
      c <- dim(DELNC_PCG_1_eQTL)[1]
      d <- dim(DELNC_PCG_1)[1]
      
      GTEX_SuppTable_TissueRes[[i]] <- c(
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
        fisher.test(data.frame("cisLnc" = c(a, b-a),
                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
        a,b,c,d)
    } else {
      GTEX_SuppTable_TissueRes[[i]] <- c(NA, NA, NA, NA, NA, NA)
    }
  }
  
  names(GTEX_SuppTable_TissueRes) <- unique(GTEX_4pairs$tissueType)
  GTEX_SuppTable[[j]] <- as.data.frame(t(bind_rows(GTEX_SuppTable_TissueRes, .id = "tissue")))
  colnames(GTEX_SuppTable[[j]]) <- c("p", "OR", "a", "b", "c", "d")
}

length(GTEX_SuppTable)
GTEX_runs_parameters$Overlap_type_name <- GTEX_runs_parameters$Overlap_type
GTEX_runs_parameters$Overlap_type_name[1:10] <- "Locus"
names(GTEX_SuppTable) <- paste("Run", GTEX_runs_parameters[,1], GTEX_runs_parameters[,4], sep = "_")
GTEX_SuppTable_df250_00001 <- bind_rows(GTEX_SuppTable, .id = "Run")
GTEX_SuppTable_df250_00001$tissue <- gsub("\\.\\.\\.[0-9]*", "", rownames(GTEX_SuppTable_df250_00001))
rownames(GTEX_SuppTable_df250_00001) <- NULL

#still poorer results, v. few below 0.01:
#promoter artery-aorta in updel:
#sig is depleted but not completely lost without the zscans

#but did identify 6x lncRNAs - not worth the faff! just use FANTOM in that case

#write.csv(GTEX_SuppTable_df, "GTEX_SuppTable_df.csv", row.names = F)



#### old code #### 0.00001
#### now set up/run optimisation - closer pairs, same/delayed timeframe pairs ####

#bit odd that totally different tissues to previous have been found
#a key diff is the genomic window being much wider @1mbp
#400kbp was the previous, so tried again here
#in theory, less noisy longer distance pairs
#another is use of both same/delayed together

#selected p
pThresh_df <- data.frame("pThresh" = c(rep(0.05, 1),
                                       rep(0.00005, 1),
                                       rep(0.00001, 1), 
                                       rep(0.000001, 1),
                                       rep(0.0000001, 1),
                                       rep(0.00000001, 1)))

AllLNC_AllPCG_2d3d_400 <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <400)

CoRegPairs_04_48_24_extended <- filter(AllLNC_AllPCG_2d3d_400,
                                           #AllLNC_AllPCG_1,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                         fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE$EnsID)) |
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                           fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                            fpkm_allGDE_within_24$EnsID)) |
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                           fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                             fpkm_allGDE_Downwithin_24$EnsID)))
#434

#store outputs here:
GTEX_eQTL_locusFish <- list()
p_test_list <- list()

#run once for each group of lncRNA overlap type
#whole locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400
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
GTEX_eQTL_locusFish_samedel400_df <- p_test_list_df
#write.csv(GTEX_eQTL_locusFish_same_df, "GTEX_eQTL_locusFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_locusFish_samedel400_df, p <0.05)$tissue)
#0x tissues 


#promoter
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400
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
GTEX_eQTL_promoterFish_samedel400_df <- p_test_list_df
#write.csv(GTEX_eQTL_promoterFish_same_df, "GTEX_eQTL_promoterFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterFish_samedel400_df, p <0.05)$tissue)
#1x tissues 


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400
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
GTEX_eQTL_exonFish_samedel400_df <- p_test_list_df
#write.csv(GTEX_eQTL_exonFish_same_df, "GTEX_eQTL_exonFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_exonFish_samedel400_df, p <0.05)$tissue)
#0x tissues 


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400
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
GTEX_eQTL_promoterExonFish_samedel400_df <- p_test_list_df
#write.csv(GTEX_eQTL_promoterExonFish_same_df, "GTEX_eQTL_promoterExonFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterExonFish_samedel400_df, p <0.05)$tissue)
#0x tissues (loss of 9x)


#Intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400
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
GTEX_eQTL_IntronFish_samedel400_df <- p_test_list_df
#write.csv(GTEX_eQTL_IntronFish_same_df, "GTEX_eQTL_IntronFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_IntronFish_samedel400_df, p <0.05)$tissue)
#0x tissues (loss of 2x)


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400
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
GTEX_eQTL_TTSFish_samedel400_df <- p_test_list_df
#write.csv(GTEX_eQTL_TTSFish_same_df, "GTEX_eQTL_TTSFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_TTSFish_samedel400_df, p <0.05)$tissue)
#0x tissues (loss of 7x)

#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400
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
GTEX_eQTL_SpliceFish_samedel400_df <- p_test_list_df
#write.csv(GTEX_eQTL_SpliceFish_same_df, "GTEX_eQTL_SpliceFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_SpliceFish_samedel400_df, p <0.05)$tissue)
#0x tissue



#### now set up/run optimisation - closer pairs, same/delayed timeframe pairs, well-expressed ####

#bit odd that totally different tissues to previous have been found
#a key diff is the genomic window being much wider @1mbp
#400kbp was the previous, so tried again here
#in theory, less noisy longer distance pairs
#another is use of both same/delayed together
#another is use of higher expressed genes

#selected p
pThresh_df <- data.frame("pThresh" = c(rep(0.05, 1),
                                       rep(0.00005, 1),
                                       rep(0.00001, 1), 
                                       rep(0.000001, 1),
                                       rep(0.0000001, 1),
                                       rep(0.00000001, 1)))

AllLNC_AllPCG_2d3d_400 <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <400)

#higher exp requirement:
AllLNC_AllPCG_2d3d_400_FPKM_1.5 <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <400, 
                                          EnsID %in% filter(fpkm_allG, fpkm_max_treatment>1.5)$EnsID, 
                                          EnsID.y %in% filter(fpkm_allG, fpkm_max_treatment>1.5)$EnsID)

CoRegPairs_04_48_24_extended <- filter(AllLNC_AllPCG_2d3d_400_FPKM_1.5,
                                       #AllLNC_AllPCG_1,
                                       (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                     fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE$EnsID)) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                       fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                        fpkm_allGDE_within_24$EnsID)) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                       fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                         fpkm_allGDE_Downwithin_24$EnsID)))
#321

#store outputs here:
GTEX_eQTL_locusFish <- list()
p_test_list <- list()

#run once for each group of lncRNA overlap type
#whole locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400_FPKM_1.5
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
GTEX_eQTL_locusFish_samedel400hiXP_df <- p_test_list_df
#write.csv(GTEX_eQTL_locusFish_same_df, "GTEX_eQTL_locusFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_locusFish_samedel400hiXP_df, p <0.05)$tissue)
#0x tissues 


#promoter
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400_FPKM_1.5
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
GTEX_eQTL_promoterFish_samedel400hiXP_df <- p_test_list_df
#write.csv(GTEX_eQTL_promoterFish_same_df, "GTEX_eQTL_promoterFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterFish_samedel400hiXP_df, p <0.05)$tissue)
#1x tissues 


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400_FPKM_1.5
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
GTEX_eQTL_exonFish_samedel400hiXP_df <- p_test_list_df
#write.csv(GTEX_eQTL_exonFish_same_df, "GTEX_eQTL_exonFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_exonFish_samedel400hiXP_df, p <0.05)$tissue)
#0x tissues 


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400_FPKM_1.5
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
GTEX_eQTL_promoterExonFish_samedel400hiXP_df <- p_test_list_df
#write.csv(GTEX_eQTL_promoterExonFish_same_df, "GTEX_eQTL_promoterExonFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_promoterExonFish_samedel400hiXP_df, p <0.05)$tissue)
#0x tissues (loss of 9x)


#Intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400_FPKM_1.5
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
GTEX_eQTL_IntronFish_samedel400hiXP_df <- p_test_list_df
#write.csv(GTEX_eQTL_IntronFish_same_df, "GTEX_eQTL_IntronFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_IntronFish_samedel400hiXP_df, p <0.05)$tissue)
#0x tissues 


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400_FPKM_1.5
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
GTEX_eQTL_TTSFish_samedel400hiXP_df <- p_test_list_df
#write.csv(GTEX_eQTL_TTSFish_same_df, "GTEX_eQTL_TTSFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_TTSFish_samedel400hiXP_df, p <0.05)$tissue)
#0x tissues 

#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  for (i in 1:length(unique(GTEX_4pairs$tissueType))){
    triali <- AllLNC_AllPCG_2d3d_400_FPKM_1.5
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
GTEX_eQTL_SpliceFish_samedel400hiXP_df <- p_test_list_df
#write.csv(GTEX_eQTL_SpliceFish_same_df, "GTEX_eQTL_SpliceFish_same_df.csv", row.names = F)
unique(filter(GTEX_eQTL_SpliceFish_samedel400hiXP_df, p <0.05)$tissue)
#0x tissue



