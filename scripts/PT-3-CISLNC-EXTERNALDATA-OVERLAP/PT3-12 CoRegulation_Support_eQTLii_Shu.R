#### SMC biobank of Shu group, from contact with Charles Solomon ####
#repeating 11 may be needed if going from scratch

library(dplyr)
library(ggplot2)
library(GenomicRanges)


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
length(unique(AllLNC_AllPCG_2d3d$EnsID.y))#5362

Shu_allVar26 <- read.csv("PCGs_4_Charles_Jan2026_Summstats_hg38_050126")
length(unique(Shu_allVar26$gene))#5180 IDs returned

length(unique(AllLNC_AllPCG_2d3d$EnsID.y))#5199 genes returned
sum(gsub("\\.[0-9]*", "", unique(AllLNC_AllPCG_2d3d$EnsID.y)) %in% Shu_allVar26$gene)#5021 (missing are the >1mbp)
5181/5362

#focus on pairs found with TSSs 1mbp
#AllLNC_AllPCG_2d3d <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <1000)
#dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
#266/8439

#some different filtering may be necessary:
#e.g. maximum range of GTEX eQTLs is 1mbp - likely the case with Shu too? so are longer range genes unfairly discriminated against?
#n.b. there may also be a minimum range of GTEX eQTLs (less likely)

AllLNC_AllPCG_2d3d <- filter(AllLNC_AllPCG_2d3d, pair_range <1000000)
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
260/8116

length(unique(AllLNC_AllPCG_2d3d$EnsID.y))#5199 genes returned
sum(gsub("\\.[0-9]*", "", unique(AllLNC_AllPCG_2d3d$EnsID.y)) %in% Shu_allVar26$gene)#5021 (missing are the >1mbp)
5021/5199

#missing small 3.4% amount of genes, presumed to not be eGenes


#### All variants for all 2d/3d PCGs ####

#All variants needs to just become "lnc variants", first 2 columns of Charles table:

#variants from subsetting the tables for potential lnc targets:
targetVariants <- unique(Shu_allVar26[,1:3])
targetVariants$chr <- paste("chr", sapply(sapply(targetVariants$snps, strsplit, ":"), "[[", 1), sep = "")
targetVariants$coords <- sapply(sapply(targetVariants$snps, strsplit, ":"), "[[", 2)

SVSMC_pairedlnc <- unique(AllLNC_AllPCG_2d3d$EnsID)#576 lncs in the pairings


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
#1508 total lnc tx for those paired to a neighbour expressed PCG 

#check for lnc promoter/genebody overlapping variants:
targetVariants_GR <- makeGRangesFromDataFrame(targetVariants, start.field = "coords", end.field = "coords", keep.extra.columns = T)

Variantindex <- findOverlaps(query = targetVariants_GR, subject = PairedLncs_GR)
Variantoverlaps <- unique(data.frame("Variant" = targetVariants_GR$snps[queryHits(Variantindex)],
                                     "Tx" = PairedLncs_GR$MSTRG_Tx_ID[subjectHits(Variantindex)]))

length(unique(Variantoverlaps$Tx)) #1301 tx
length(unique(Variantoverlaps$Variant)) #25462
table(table(Variantoverlaps$Tx)) #some tx with insane no. variants, probs long intron lncRNAs? or highly mutated regions


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
stringtie_gtf_majorPairedLncs <- merge(stringtie_gtf_majorPairedLncs, fpkm_allG[,c(2,47)], by = "MSTRG_Tx_ID")
length(unique(stringtie_gtf_majorPairedLncs$EnsID))#576 genes
length(unique(stringtie_gtf_majorPairedLncs$MSTRG_Tx_ID))#1508 transcripts

stringtie_gtf_majorPairedLncsExons <- filter(stringtie_gtf_majorPairedLncs, V3 == "exon")

PairedLncsExons_GR <- makeGRangesFromDataFrame(stringtie_gtf_majorPairedLncsExons[,c(2,5,6,1,11)], seqnames.field = "V1", 
                                               start.field = "V4", end.field = "V5", keep.extra.columns = T)

#overlap exons
VariantindexE <- findOverlaps(query = targetVariants_GR, subject = PairedLncsExons_GR)
VariantoverlapsE <- unique(data.frame("Variant" = targetVariants_GR$snps[queryHits(VariantindexE)],
                                      "Tx" = PairedLncsExons_GR$MSTRG_Tx_ID[subjectHits(VariantindexE)]))

length(unique(VariantoverlapsE$Tx)) #968 
length(unique(VariantoverlapsE$Variant)) #2707 
table(table(VariantoverlapsE$Tx)) #fewer but still some genes with insane no. variants


#additional variants in promoter region of lncs:
VariantindexP <- findOverlaps(query = targetVariants_GR, subject = promoters(PairedLncs_GR, upstream = 2000,
                                                                             downstream = 0))
VariantoverlapsP <- data.frame("Variant" = targetVariants_GR$snps[queryHits(VariantindexP)],
                               "Tx" = PairedLncs_GR$MSTRG_Tx_ID[subjectHits(VariantindexP)])
length(unique(VariantoverlapsP$Tx)) #1162
length(unique(VariantoverlapsP$Variant)) #2478
table(table(VariantoverlapsP$Tx))


#overlap splice junctions (needs work)
PairedLncsExonLimits_GR <- makeGRangesFromDataFrame(data.frame("seqnames" = seqnames(PairedLncsExons_GR), 
                                                               "start" = c(start(PairedLncsExons_GR), end(PairedLncsExons_GR)), 
                                                               "end" = c(start(PairedLncsExons_GR), end(PairedLncsExons_GR))))
#add 50bp either way
PairedLncsSJ_GR <- flank(PairedLncsExonLimits_GR, width = 50, both = T)

VariantindexSJ <- findOverlaps(query = targetVariants_GR, subject = PairedLncsSJ_GR)
VariantoverlapsSJ <- unique(data.frame("Variant" = targetVariants_GR$snps[queryHits(VariantindexSJ)],
                                       "Tx" = PairedLncsExons_GR$MSTRG_Tx_ID[subjectHits(VariantindexSJ)]))
length(unique(VariantoverlapsSJ$Tx)) #493
length(unique(VariantoverlapsSJ$Variant)) #883
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
#settled on 500bp, less chance of kicking out eQTLs with bad overlap
PairedLncsTTS_GR <- terminators(PairedLncsTTS_GR, downstream = 500, upstream = 0)

VariantindexTTS <- findOverlaps(query = targetVariants_GR, subject = PairedLncsTTS_GR)
VariantoverlapsTTS <- unique(data.frame("Variant" = targetVariants_GR$snps[queryHits(VariantindexTTS)],
                                        "Tx" = PairedLncsTTS_GR$MSTRG_Tx_ID[subjectHits(VariantindexTTS)]))
length(unique(VariantoverlapsTTS$Tx)) #689
length(unique(VariantoverlapsTTS$Variant)) #858
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
length(unique(VariantoverlapsI$Tx)) #1153
length(unique(VariantoverlapsI$Variant)) #23724
table(table(VariantoverlapsI$Tx))

#can subset the table by these lnc-overlapping variants:
#list of variants to take back to subset the gtex tables:
CisLncVariants <- unique(c(VariantoverlapsE$Variant, VariantoverlapsP$Variant, VariantoverlapsSJ$Variant, 
                           VariantoverlapsI$Variant, VariantoverlapsTTS$Variant))
length(unique(CisLncVariants))#27861

Shu_LncVar26 <- filter(Shu_allVar26, snps %in% CisLncVariants)
head(Shu_LncVar26)

#annotate variants
VariantoverlapsE$OverlapType <- "Exon"
VariantoverlapsP$OverlapType <- "Promoter"
VariantoverlapsSJ$OverlapType <- "Splice"
VariantoverlapsI$OverlapType <- "Intron"
VariantoverlapsTTS$OverlapType <- "TTS"

VariantoverlapsAll <- rbind(VariantoverlapsE, VariantoverlapsP, VariantoverlapsSJ, VariantoverlapsI, VariantoverlapsTTS)
length(unique(VariantoverlapsAll$Variant))#27861 as expected

#add info about multiple overlaps:
allGB_GR <- makeGRangesFromDataFrame(allGB[,c(1:6)], 
                                     start.field = "Tx_start", 
                                     end.field = "Tx_stop", 
                                     seqnames.field = "chr", 
                                     strand.field = "str", keep.extra.columns = T)
#42511 total tx to check

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

length(unique(Variantoverlaps$Tx)) #2335 tx
length(unique(Variantoverlaps$Variant)) #27861 variants, as expected

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
#@1mbp 188 overlap 3, 23260 overlap 1  (i.e. a lncRNA)

#check on IGV - all good (bear in mind promoter regions used too)
VariantoverlapsAll <- merge(VariantoverlapsAll, trialii, by = "Variant")

Shu_LncVar26_2 <- unique(merge(VariantoverlapsAll, Shu_LncVar26, by.x = "Variant", by.y = "snps"))
colnames(Shu_LncVar26_2)[c(2,8)] <- c("MSTRG_Tx_ID", "EnsID.y")
length(unique(Shu_LncVar26_2$Variant))
length(unique(VariantoverlapsAll$Variant))

length(unique(Shu_LncVar26_2$EnsID.y))#3501 lnc neighbours have an eQTL overlapping a lncRNA

head(Shu_LncVar26_2)

#add in lncRNA IDs
trial <- merge(fpkm_allG[,c(2,47)], Shu_LncVar26_2, by = "MSTRG_Tx_ID")
Shu_LncVar26_2 <- trial

#write.csv(Shu_LncVar26_2, "Shu_LncVar26.csv", row.names = F)


#### import other bits for testing enrichment ####

#eqtl-egnenes that overlap lncs:
Shu_LncVar <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Shu_LncVar26.csv")
length(unique(Shu_LncVar$EnsID.y))#3501 eGenes with an eQTL overlapping a lncRNA

#expression:
#From GEO for the paper - raw counts for the genes of interest in SVSMC - download here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189300
#then subset in Eddie, import lncRNA-neighbouring complement of genes here:
Shu_SVSMCpairedPCG_exprs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/PCGs_4_Charles_counts26.csv", header = F)

length((Shu_SVSMCpairedPCG_exprs$V1))#5343
sum(Shu_SVSMCpairedPCG_exprs$V1 %in% gsub("\\.[0-9]*", "", Shu_LncVar$EnsID.y))

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

2204/1.256/(30360137/1000000)

#may not remove anything other than v. low TPM
sum(rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))>1)/dim(Shu_SVSMCpairedPCG_exprsTPM)[1]#-5%
sum(rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))>5)/dim(Shu_SVSMCpairedPCG_exprsTPM)[1]#-18%
sum(rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))>10)/dim(Shu_SVSMCpairedPCG_exprsTPM)[1]#-32%
sum(rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))>20)/dim(Shu_SVSMCpairedPCG_exprsTPM)[1]#-52%
sum(rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))>30)/dim(Shu_SVSMCpairedPCG_exprsTPM)[1]#-65%
sum(rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))>50)/dim(Shu_SVSMCpairedPCG_exprsTPM)[1]#-77%

#check mki67 - are SMCs activated?:
Shu_mki67 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/charles_mki67.csv", header = F)

#estimate FPKM, gene length from timecourse: 12314.66
Shu_mki67_fpk <- as.numeric(Shu_mki67[,-1])/12314.66
Shu_mki67_fpkm <- Shu_mki67_fpk/(Shu_colSums[,2]/1000000)
summary(Shu_mki67_fpkm)
#barely expressed

#others? CXCL8:
summary(as.numeric(filter(Shu_SVSMCpairedPCG_exprsTPM, EnsID_merge %in% c("ENSG00000169429"))))
#IL6
summary(as.numeric(filter(Shu_SVSMCpairedPCG_exprsTPM, EnsID_merge %in% c("ENSG00000136244"))))
#they look pretty activated, low level relative to IP for sure tho


#identify eGenes at given p val
#original requested genes
Shu_allVar_p <- read.csv("PCGs_4_Charles_Jan2026_Summstats_hg38_050126.csv")
length(unique(Shu_allVar_p$gene)) #5181 eGenes from Shu/Charles

#note very diff p val distribution in this resource compared to GTEX
pThresh <- 0.05    # -0%
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < pThresh)
length(unique(Shu_allVar_pThresh$snps))/length(unique(Shu_allVar_p$snps))

pThresh <- 0.025    # -22%
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < pThresh)
length(unique(Shu_allVar_pThresh$snps))/length(unique(Shu_allVar_p$snps))

pThresh <- 0.01    # -35%
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < pThresh)
length(unique(Shu_allVar_pThresh$snps))/length(unique(Shu_allVar_p$snps))

pThresh <- 0.005   # -45%
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < pThresh)
length(unique(Shu_allVar_pThresh$snps))/length(unique(Shu_allVar_p$snps))

pThresh <- 0.001   # -55%
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < pThresh)
length(unique(Shu_allVar_pThresh$snps))/length(unique(Shu_allVar_p$snps))

pThresh <- 0.0001  # -65%
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < pThresh)
length(unique(Shu_allVar_pThresh$snps))/length(unique(Shu_allVar_p$snps))

pThresh <- 0.00001 # -70%
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < pThresh)
length(unique(Shu_allVar_pThresh$snps))/length(unique(Shu_allVar_p$snps))

pThresh <- 0.000001 # -75%
Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < pThresh)
length(unique(Shu_allVar_pThresh$snps))/length(unique(Shu_allVar_p$snps))

#needs more thresholding than GTEX did to reduce

#note also the way more samples (like 3x the amount) may contribute to this effect (more power, higher p)

#note also that thresholding did not add much value in GTEX


#### 1a SAME TIMEFRAME cisLnc pairs -first optimisation ####

#double check numbers
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
260/8116
length(unique(AllLNC_AllPCG_2d3d$EnsID.y))#5199 genes returned


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
#456 testable

#ready to run tests as dpair_range_bin#ready to run tests as done on GTEX:
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

#0x significant tests


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

#some sig from:
#a) promExon
#b) whole locus
#c) TTS (at eQTLp=0.025)

#positive: enrichment of eQTL support amongst co-regs

#negative: not significant in isolation, adding in intron improves, so not prom/exon alone, weak OR

#the odds ratios are not great, 1.3 max

#considering range? are co-regs just closer together and so more likely to get eQTLs?
#no great diffs in ranges of pairs
summary(AllLNC_AllPCG_2d3d$pair_range)
summary(CoRegPairs_04_48_24_extendedSame$pair_range)
boxplot(AllLNC_AllPCG_2d3d$pair_range, CoRegPairs_04_48_24_extendedSame$pair_range)

#slight more obvious diff in TSS distance of pairs
summary(AllLNC_AllPCG_2d3d$AbsDistLnc_PCG)
summary(CoRegPairs_04_48_24_extendedSame$AbsDistLnc_PCG)
boxplot(AllLNC_AllPCG_2d3d$AbsDistLnc_PCG, CoRegPairs_04_48_24_extendedSame$AbsDistLnc_PCG)

#so some of effect could be simply that the coregs are closer together? but seems v slight


#### 1b SAME TIMEFRAME cisLnc pairs -second optimisation - pair dist ####

#improve by focusing on narrower range? expect more true cisReg connections

#double check numbers
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
260/8116
length(unique(AllLNC_AllPCG_2d3d$EnsID.y))#5199 genes returned

#250kbp also seemed of interest in GTEX, still some sig
AllLNC_AllPCG_2d3d_250 <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <250)
dim(filter(AllLNC_AllPCG_2d3d_250, !loopMethod == "Neither"))
185/2563
length(unique(AllLNC_AllPCG_2d3d_250$EnsID.y))#2152 genes returned

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

#no range diff
summary(AllLNC_AllPCG_2d3d_250$pair_range)
summary(CoRegPairs_04_48_24_extendedSame$pair_range)
boxplot(AllLNC_AllPCG_2d3d_250$pair_range, CoRegPairs_04_48_24_extendedSame$pair_range)

#no strong TSS diff
summary(AllLNC_AllPCG_2d3d_250$AbsDistLnc_PCG)
summary(CoRegPairs_04_48_24_extendedSame$AbsDistLnc_PCG)
boxplot(AllLNC_AllPCG_2d3d_250$AbsDistLnc_PCG, CoRegPairs_04_48_24_extendedSame$AbsDistLnc_PCG)

#promoter:
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_250, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfprom_250 <- Shu_fish_list_df
Shu_fish_list_same_dfprom_250$eQTLoverlaps <- "Promoter"

#write.csv(Shu_fish_list_same_dfprom, "Shu_fish_list_same_dfprom.csv", row.names = F)

#0x significant tests (again)


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_250, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfexon_250 <- Shu_fish_list_df
Shu_fish_list_same_dfexon_250$eQTLoverlaps <- "Exon"

#write.csv(Shu_fish_list_same_dfexon, "Shu_fish_list_same_dfexon.csv", row.names = F)

#no good

#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_250, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfpromExon_250 <- Shu_fish_list_df
Shu_fish_list_same_dfpromExon_250$eQTLoverlaps <- "PromoterExon"

#write.csv(Shu_fish_list_same_dfpromExon, "Shu_fish_list_same_dfpromExon.csv", row.names = F)

#no good 0x

#intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_250, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfintron_250 <- Shu_fish_list_df
Shu_fish_list_same_dfintron_250$eQTLoverlaps <- "Intron"

#write.csv(Shu_fish_list_same_dfintron, "Shu_fish_list_same_dfintron.csv", row.names = F)

#no good


#locus
for (j in 1:length(pThresh_df$pThresh)){
    
    selectedp <- pThresh_df[j,1]
    
    Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
    Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
    
    triali <- filter(AllLNC_AllPCG_2d3d_250, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dflocus_250 <- Shu_fish_list_df
Shu_fish_list_same_dflocus_250$eQTLoverlaps <- "Locus"

#write.csv(Shu_fish_list_same_dflocus, "Shu_fish_list_same_dflocus.csv", row.names = F)

#no good


#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_250, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfsplice_250 <- Shu_fish_list_df
Shu_fish_list_same_dfsplice_250$eQTLoverlaps <- "Splice"

#write.csv(Shu_fish_list_same_dfsplice, "Shu_fish_list_same_dfsplice.csv", row.names = F)

#no good


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_250, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfTTS_250 <- Shu_fish_list_df
Shu_fish_list_same_dfTTS_250$eQTLoverlaps <- "TTS"

#write.csv(Shu_fish_list_same_dfTTS, "Shu_fish_list_same_dfTTS.csv", row.names = F)

#no good

#in 250kbp no significance


#### 1c SAME TIMEFRAME cisLnc pairs -third optimisation - no overlap criteria ####

#has removing the overlaps helped?

#double check numbers
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
260/8116
length(unique(AllLNC_AllPCG_2d3d$EnsID.y))#5199 genes returned

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

#ready to run tests as done on GTEX:


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
                                                              #TotalOverlaps == 1
                                                              )$pairs] <- "Yes"
  
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

Shu_fish_list_same_dfprom_withOverlaps <- Shu_fish_list_df
Shu_fish_list_same_dfprom_withOverlaps$eQTLoverlaps <- "Promoter"

#write.csv(Shu_fish_list_same_dfprom, "Shu_fish_list_same_dfprom.csv", row.names = F)

#baseline sig, improvement


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
                                                              #TotalOverlaps == 1
                                                              )$pairs] <- "Yes"
  
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

Shu_fish_list_same_dfexon_withOverlaps <- Shu_fish_list_df
Shu_fish_list_same_dfexon_withOverlaps$eQTLoverlaps <- "Exon"

#write.csv(Shu_fish_list_same_dfexon, "Shu_fish_list_same_dfexon.csv", row.names = F)

#quite close


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
                                                              #TotalOverlaps == 1
                                                              )$pairs] <- "Yes"
  
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

Shu_fish_list_same_dfpromExon_withOverlaps <- Shu_fish_list_df
Shu_fish_list_same_dfpromExon_withOverlaps$eQTLoverlaps <- "PromoterExon"

#write.csv(Shu_fish_list_same_dfpromExon, "Shu_fish_list_same_dfpromExon.csv", row.names = F)

#still baseline sig, improves prom alone slightly

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
                                                              #TotalOverlaps == 1
                                                              )$pairs] <- "Yes"
  
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

Shu_fish_list_same_dfintron_withOverlaps <- Shu_fish_list_df
Shu_fish_list_same_dfintron_withOverlaps$eQTLoverlaps <- "Intron"

#write.csv(Shu_fish_list_same_dfintron, "Shu_fish_list_same_dfintron.csv", row.names = F)

#close


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
                                                              #TotalOverlaps == 1
                                                              )$pairs] <- "Yes"
  
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

Shu_fish_list_same_dflocus_withOverlaps <- Shu_fish_list_df
Shu_fish_list_same_dflocus_withOverlaps$eQTLoverlaps <- "Locus"

#write.csv(Shu_fish_list_same_dflocus, "Shu_fish_list_same_dflocus.csv", row.names = F)

#sig, bit worse off than promExon

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
                                                              #TotalOverlaps == 1
                                                              )$pairs] <- "Yes"
  
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

Shu_fish_list_same_dfsplice_withOverlaps <- Shu_fish_list_df
Shu_fish_list_same_dfsplice_withOverlaps$eQTLoverlaps <- "Splice"

#write.csv(Shu_fish_list_same_dfsplice, "Shu_fish_list_same_dfsplice.csv", row.names = F)

#nope


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
                                                              #TotalOverlaps == 1
                                                              )$pairs] <- "Yes"
  
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

Shu_fish_list_same_dfTTS_withOverlaps <- Shu_fish_list_df
Shu_fish_list_same_dfTTS_withOverlaps$eQTLoverlaps <- "TTS"

#write.csv(Shu_fish_list_same_dfTTS, "Shu_fish_list_same_dfTTS.csv", row.names = F)

#nope

#without overlap considerations some sig from:
#a) prom
#b) promExon
#c) whole locus

#positive: back to previous expectations, prom/exon do better, improvement

#negative: still fairly weak ORs, including potential for e.g. bidirectional promoters

#do not recommend






#### 2a DELAYED TIMEFRAME cisLnc pairs -first optimisation ####

#double check numbers
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
260/8116
length(unique(AllLNC_AllPCG_2d3d$EnsID.y))#5199 genes returned

fpkm_allGDE_within_4 <- rbind(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Downwithin_4)
fpkm_allGDE_within_8 <- rbind(fpkm_allGDE_Upwithin_8, fpkm_allGDE_Downwithin_8)
fpkm_allGDE_within_24 <- rbind(fpkm_allGDE_Upwithin_24, fpkm_allGDE_Downwithin_24)

CoRegPairs_04_48_24_extendedDelayed <- filter(AllLNC_AllPCG_2d3d,
                                              #AllLNC_AllPCG_1,
                                              (EnsID %in% fpkm_allGDE_within_4$EnsID & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                      fpkm_allGDE_within_24$EnsID)) |
                                                (EnsID %in% fpkm_allGDE_within_8$EnsID & EnsID.y %in% fpkm_allGDE_within_24$EnsID))
#506 pairs testable

#ready to run tests as dpair_range_bin#ready to run tests as done on GTEX:
Shu_fish_list <- list()

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
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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

Shu_fish_list_same_dfprom_del <- Shu_fish_list_df
Shu_fish_list_same_dfprom_del$eQTLoverlaps <- "Promoter"

#write.csv(Shu_fish_list_same_dfprom, "Shu_fish_list_same_dfprom.csv", row.names = F)

#0x significant tests


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
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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

Shu_fish_list_same_dfexon_del <- Shu_fish_list_df
Shu_fish_list_same_dfexon_del$eQTLoverlaps <- "Exon"

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
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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

Shu_fish_list_same_dfpromExon_del <- Shu_fish_list_df
Shu_fish_list_same_dfpromExon_del$eQTLoverlaps <- "PromoterExon"

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
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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

Shu_fish_list_same_dfintron_del <- Shu_fish_list_df
Shu_fish_list_same_dfintron_del$eQTLoverlaps <- "Intron"

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
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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

Shu_fish_list_same_dflocus_del <- Shu_fish_list_df
Shu_fish_list_same_dflocus_del$eQTLoverlaps <- "Locus"

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
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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

Shu_fish_list_same_dfsplice_del <- Shu_fish_list_df
Shu_fish_list_same_dfsplice_del$eQTLoverlaps <- "Splice"

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
  DELNC_DEPCG_1 <- filter(triali, pairs %in% CoRegPairs_04_48_24_extendedDelayed$pairs)
  DELNC_DEPCG_1_eQTL <- filter(DELNC_DEPCG_1, eQTL_validated_tissue == "Yes")
  #other neighbours of cisLncs
  DELNC_PCG_1 <- filter(triali, EnsID %in% CoRegPairs_04_48_24_extendedDelayed$EnsID)
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

Shu_fish_list_same_dfTTS_del <- Shu_fish_list_df
Shu_fish_list_same_dfTTS_del$eQTLoverlaps <- "TTS"

#write.csv(Shu_fish_list_same_dfTTS, "Shu_fish_list_same_dfTTS.csv", row.names = F)

#no sig found, not even a hint of sig


#### 3a EARLY SAME TIMEFRAME cisLnc pairs -first optimisation ####

#double check numbers
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
260/8116
length(unique(AllLNC_AllPCG_2d3d$EnsID.y))#5199 genes returned

#filter to 0-4hr expressed genes:
fpkm_allG_04 <- filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)
AllLNC_AllPCG_2d3d_04 <- filter(AllLNC_AllPCG_2d3d, EnsID %in% fpkm_allG_04$EnsID, EnsID.y %in% fpkm_allG_04$EnsID)

dim(filter(AllLNC_AllPCG_2d3d_04, !loopMethod == "Neither"))
223/7222
length(unique(AllLNC_AllPCG_2d3d_04$EnsID.y))#4611 genes returned

CoRegPairs_04_48_24_extendedSame <- filter(AllLNC_AllPCG_2d3d_04,
                                           #AllLNC_AllPCG_1,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                         fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                          fpkm_allGDE_Downwithin_4$EnsID)))
)
#253 testable

#ready to run tests as dpair_range_bin#ready to run tests as done on GTEX:
Shu_fish_list <- list()

#allow finding pairs in each:
Shu_LncVar$pairs <- paste(Shu_LncVar$EnsID, Shu_LncVar$EnsID.y, sep = "-")

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
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfprom04 <- Shu_fish_list_df
Shu_fish_list_same_dfprom04$eQTLoverlaps <- "Promoter"

#write.csv(Shu_fish_list_same_dfprom, "Shu_fish_list_same_dfprom.csv", row.names = F)

#0x significant tests


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfexon04 <- Shu_fish_list_df
Shu_fish_list_same_dfexon04$eQTLoverlaps <- "Exon"

#write.csv(Shu_fish_list_same_dfexon, "Shu_fish_list_same_dfexon.csv", row.names = F)


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfpromExon04 <- Shu_fish_list_df
Shu_fish_list_same_dfpromExon04$eQTLoverlaps <- "PromoterExon"

#write.csv(Shu_fish_list_same_dfpromExon, "Shu_fish_list_same_dfpromExon.csv", row.names = F)


#intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfintron04 <- Shu_fish_list_df
Shu_fish_list_same_dfintron04$eQTLoverlaps <- "Intron"

#write.csv(Shu_fish_list_same_dfintron, "Shu_fish_list_same_dfintron.csv", row.names = F)


#locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dflocus04 <- Shu_fish_list_df
Shu_fish_list_same_dflocus04$eQTLoverlaps <- "Locus"

#write.csv(Shu_fish_list_same_dflocus, "Shu_fish_list_same_dflocus.csv", row.names = F)


#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfsplice04 <- Shu_fish_list_df
Shu_fish_list_same_dfsplice04$eQTLoverlaps <- "Splice"

#write.csv(Shu_fish_list_same_dfsplice, "Shu_fish_list_same_dfsplice.csv", row.names = F)


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfTTS04 <- Shu_fish_list_df
Shu_fish_list_same_dfTTS04$eQTLoverlaps <- "TTS"

#write.csv(Shu_fish_list_same_dfTTS, "Shu_fish_list_same_dfTTS.csv", row.names = F)

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

#summary:
#some sig from whole locus only, OR reaches 1.45

#positive: ORs do increase a bit, one OR of 1.45 (still modest compared to previous)

#negative: less significance, distance appears to be a factor



#### 3b EARLY SAME TIMEFRAME cisLnc pairs -second optimisation - pair dist ####

#double check numbers
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
260/8116
length(unique(AllLNC_AllPCG_2d3d$EnsID.y))#5199 genes returned

#filter to 0-4hr expressed genes:
fpkm_allG_04 <- filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)
AllLNC_AllPCG_2d3d_04 <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <250, 
                                EnsID %in% fpkm_allG_04$EnsID, EnsID.y %in% fpkm_allG_04$EnsID)

dim(filter(AllLNC_AllPCG_2d3d_04, !loopMethod == "Neither"))
159/2298
length(unique(AllLNC_AllPCG_2d3d_04$EnsID.y))#1921 genes returned

CoRegPairs_04_48_24_extendedSame <- filter(AllLNC_AllPCG_2d3d_04,
                                           #AllLNC_AllPCG_1,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                         fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                          fpkm_allGDE_Downwithin_4$EnsID)))
)
#94 testable

#ready to run tests as dpair_range_bin#ready to run tests as done on GTEX:
Shu_fish_list <- list()

#allow finding pairs in each:
Shu_LncVar$pairs <- paste(Shu_LncVar$EnsID, Shu_LncVar$EnsID.y, sep = "-")

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
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfprom04_250 <- Shu_fish_list_df
Shu_fish_list_same_dfprom04_250$eQTLoverlaps <- "Promoter"

#write.csv(Shu_fish_list_same_dfprom, "Shu_fish_list_same_dfprom.csv", row.names = F)

#0x significant tests


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfexon04_250 <- Shu_fish_list_df
Shu_fish_list_same_dfexon04_250$eQTLoverlaps <- "Exon"

#write.csv(Shu_fish_list_same_dfexon, "Shu_fish_list_same_dfexon.csv", row.names = F)


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfpromExon04_250 <- Shu_fish_list_df
Shu_fish_list_same_dfpromExon04_250$eQTLoverlaps <- "PromoterExon"

#write.csv(Shu_fish_list_same_dfpromExon, "Shu_fish_list_same_dfpromExon.csv", row.names = F)


#intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfintron04_250 <- Shu_fish_list_df
Shu_fish_list_same_dfintron04_250$eQTLoverlaps <- "Intron"

#write.csv(Shu_fish_list_same_dfintron, "Shu_fish_list_same_dfintron.csv", row.names = F)


#locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dflocus04_250 <- Shu_fish_list_df
Shu_fish_list_same_dflocus04_250$eQTLoverlaps <- "Locus"

#write.csv(Shu_fish_list_same_dflocus, "Shu_fish_list_same_dflocus.csv", row.names = F)


#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfsplice04_250 <- Shu_fish_list_df
Shu_fish_list_same_dfsplice04_250$eQTLoverlaps <- "Splice"

#write.csv(Shu_fish_list_same_dfsplice, "Shu_fish_list_same_dfsplice.csv", row.names = F)


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfTTS04_250 <- Shu_fish_list_df
Shu_fish_list_same_dfTTS04_250$eQTLoverlaps <- "TTS"

#write.csv(Shu_fish_list_same_dfTTS, "Shu_fish_list_same_dfTTS.csv", row.names = F)

#considering range? are co-regs just closer together and so more likely to get eQTLs?
#no great diffs in ranges of pairs
summary(AllLNC_AllPCG_2d3d_04$pair_range)
summary(CoRegPairs_04_48_24_extendedSame$pair_range)
boxplot(AllLNC_AllPCG_2d3d_04$pair_range, CoRegPairs_04_48_24_extendedSame$pair_range)

#slight more obvious diff in TSS distance of pairs
summary(AllLNC_AllPCG_2d3d_04$AbsDistLnc_PCG)
summary(CoRegPairs_04_48_24_extendedSame$AbsDistLnc_PCG)
boxplot(AllLNC_AllPCG_2d3d_04$AbsDistLnc_PCG, CoRegPairs_04_48_24_extendedSame$AbsDistLnc_PCG)

#no distance change

#summary:
#no sig across the board


#### 3c EARLY SAME TIMEFRAME cisLnc pairs -third optimisation - hi XP####

#double check numbers
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
260/8116
length(unique(AllLNC_AllPCG_2d3d$EnsID.y))#5199 genes returned

#filter to 0-4hr expressed genes of hiXP:
fpkm_allG_04 <- filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)
fpkm_allG_04_hi <- filter(fpkm_allG, Hour0_meanFPKM >10 | Hour4_meanFPKM >10)

AllLNC_AllPCG_2d3d_04 <- filter(AllLNC_AllPCG_2d3d, EnsID %in% fpkm_allG_04$EnsID, EnsID.y %in% fpkm_allG_04_hi$EnsID)

dim(filter(AllLNC_AllPCG_2d3d_04, !loopMethod == "Neither"))
161/5074
length(unique(AllLNC_AllPCG_2d3d_04$EnsID.y))#3243 genes returned

CoRegPairs_04_48_24_extendedSame <- filter(AllLNC_AllPCG_2d3d_04,
                                           #AllLNC_AllPCG_1,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                         fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                          fpkm_allGDE_Downwithin_4$EnsID)))
#180 testable

#ready to run tests as dpair_range_bin#ready to run tests as done on GTEX:
Shu_fish_list <- list()

#allow finding pairs in each:
Shu_LncVar$pairs <- paste(Shu_LncVar$EnsID, Shu_LncVar$EnsID.y, sep = "-")

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
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 10]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfprom04_hiXP <- Shu_fish_list_df
Shu_fish_list_same_dfprom04_hiXP$eQTLoverlaps <- "Promoter"

#write.csv(Shu_fish_list_same_dfprom, "Shu_fish_list_same_dfprom.csv", row.names = F)

#0x significant tests


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 10]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfexon04_hiXP <- Shu_fish_list_df
Shu_fish_list_same_dfexon04_hiXP$eQTLoverlaps <- "Exon"

#write.csv(Shu_fish_list_same_dfexon, "Shu_fish_list_same_dfexon.csv", row.names = F)


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 10]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfpromExon04_hiXP <- Shu_fish_list_df
Shu_fish_list_same_dfpromExon04_hiXP$eQTLoverlaps <- "PromoterExon"

#write.csv(Shu_fish_list_same_dfpromExon, "Shu_fish_list_same_dfpromExon.csv", row.names = F)


#intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 10]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfintron04_hiXP <- Shu_fish_list_df
Shu_fish_list_same_dfintron04_hiXP$eQTLoverlaps <- "Intron"

#write.csv(Shu_fish_list_same_dfintron, "Shu_fish_list_same_dfintron.csv", row.names = F)


#locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 10]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dflocus04_hiXP <- Shu_fish_list_df
Shu_fish_list_same_dflocus04_hiXP$eQTLoverlaps <- "Locus"

#write.csv(Shu_fish_list_same_dflocus, "Shu_fish_list_same_dflocus.csv", row.names = F)


#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 10]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfsplice04_hiXP <- Shu_fish_list_df
Shu_fish_list_same_dfsplice04_hiXP$eQTLoverlaps <- "Splice"

#write.csv(Shu_fish_list_same_dfsplice, "Shu_fish_list_same_dfsplice.csv", row.names = F)


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 10]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfTTS04_hiXP <- Shu_fish_list_df
Shu_fish_list_same_dfTTS04_hiXP$eQTLoverlaps <- "TTS"

#write.csv(Shu_fish_list_same_dfTTS, "Shu_fish_list_same_dfTTS.csv", row.names = F)

#considering range? are co-regs just closer together and so more likely to get eQTLs?
#some diffs in ranges of pairs
summary(AllLNC_AllPCG_2d3d_04$pair_range)
summary(CoRegPairs_04_48_24_extendedSame$pair_range)
boxplot(AllLNC_AllPCG_2d3d_04$pair_range, CoRegPairs_04_48_24_extendedSame$pair_range)

# diff in TSS distance of pairs
summary(AllLNC_AllPCG_2d3d_04$AbsDistLnc_PCG)
summary(CoRegPairs_04_48_24_extendedSame$AbsDistLnc_PCG)
boxplot(AllLNC_AllPCG_2d3d_04$AbsDistLnc_PCG, CoRegPairs_04_48_24_extendedSame$AbsDistLnc_PCG)

#so some of effect could be simply that the coregs are closer together? but seems v slight

#summary:
#no sig at all from higher XP of targets @5 TPM or @20TPM



#### 3d EARLY SAME TIMEFRAME cisLnc pairs -fourth optimisation - co-induceds####

#double check numbers
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
260/8116
length(unique(AllLNC_AllPCG_2d3d$EnsID.y))#5199 genes returned

#filter to 0-4hr expressed genes:
fpkm_allG_04 <- filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)
AllLNC_AllPCG_2d3d_04 <- filter(AllLNC_AllPCG_2d3d, EnsID %in% fpkm_allG_04$EnsID, EnsID.y %in% fpkm_allG_04$EnsID)

dim(filter(AllLNC_AllPCG_2d3d_04, !loopMethod == "Neither"))
223/7222
length(unique(AllLNC_AllPCG_2d3d_04$EnsID.y))#4611 genes returned

CoRegPairs_04_48_24_extendedSame <- filter(AllLNC_AllPCG_2d3d_04,
                                           #AllLNC_AllPCG_1,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID)))
)
#119 testable

Shu_fish_list <- list()

#allow finding pairs in each:
Shu_LncVar$pairs <- paste(Shu_LncVar$EnsID, Shu_LncVar$EnsID.y, sep = "-")

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
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfprom04_coI <- Shu_fish_list_df
Shu_fish_list_same_dfprom04_coI$eQTLoverlaps <- "Promoter"

#write.csv(Shu_fish_list_same_dfprom, "Shu_fish_list_same_dfprom.csv", row.names = F)

#0x significant tests


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfexon04_coI <- Shu_fish_list_df
Shu_fish_list_same_dfexon04_coI$eQTLoverlaps <- "Exon"

#write.csv(Shu_fish_list_same_dfexon, "Shu_fish_list_same_dfexon.csv", row.names = F)


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfpromExon04_coI <- Shu_fish_list_df
Shu_fish_list_same_dfpromExon04_coI$eQTLoverlaps <- "PromoterExon"

#write.csv(Shu_fish_list_same_dfpromExon, "Shu_fish_list_same_dfpromExon.csv", row.names = F)


#intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfintron04_coI <- Shu_fish_list_df
Shu_fish_list_same_dfintron04_coI$eQTLoverlaps <- "Intron"

#write.csv(Shu_fish_list_same_dfintron, "Shu_fish_list_same_dfintron.csv", row.names = F)


#locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dflocus04_coI <- Shu_fish_list_df
Shu_fish_list_same_dflocus04_coI$eQTLoverlaps <- "Locus"

#write.csv(Shu_fish_list_same_dflocus, "Shu_fish_list_same_dflocus.csv", row.names = F)


#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfsplice04_coI <- Shu_fish_list_df
Shu_fish_list_same_dfsplice04_coI$eQTLoverlaps <- "Splice"

#write.csv(Shu_fish_list_same_dfsplice, "Shu_fish_list_same_dfsplice.csv", row.names = F)


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfTTS04_coI <- Shu_fish_list_df
Shu_fish_list_same_dfTTS04_coI$eQTLoverlaps <- "TTS"

#write.csv(Shu_fish_list_same_dfTTS, "Shu_fish_list_same_dfTTS.csv", row.names = F)

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

#summary:
#some sig from whole locus only, splice, TTS, with p adjusting

#positive: ORs do increase a bit, one OR of 1.45 (still modest compared to previous)

#negative: distance appears to be a factor, requires p adjusting




#### 3e EARLY SAME TIMEFRAME cisLnc pairs -fifth optimisation - concordants ####

#double check numbers
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
260/8116
length(unique(AllLNC_AllPCG_2d3d$EnsID.y))#5199 genes returned

#filter to 0-4hr expressed genes:
fpkm_allG_04 <- filter(fpkm_allG, Hour0_meanFPKM >0.8 | Hour4_meanFPKM >0.8)
AllLNC_AllPCG_2d3d_04 <- filter(AllLNC_AllPCG_2d3d, EnsID %in% fpkm_allG_04$EnsID, EnsID.y %in% fpkm_allG_04$EnsID)

dim(filter(AllLNC_AllPCG_2d3d_04, !loopMethod == "Neither"))
223/7222
length(unique(AllLNC_AllPCG_2d3d_04$EnsID.y))#4611 genes returned

CoRegPairs_04_48_24_extendedSame <- filter(AllLNC_AllPCG_2d3d_04,
                                           #AllLNC_AllPCG_1,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID)) |
                                             (EnsID %in% c(fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Downwithin_4$EnsID))
                                           )

#163 testable

Shu_fish_list <- list()

#allow finding pairs in each:
Shu_LncVar$pairs <- paste(Shu_LncVar$EnsID, Shu_LncVar$EnsID.y, sep = "-")

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
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfprom04_conc <- Shu_fish_list_df
Shu_fish_list_same_dfprom04_conc$eQTLoverlaps <- "Promoter"

#write.csv(Shu_fish_list_same_dfprom, "Shu_fish_list_same_dfprom.csv", row.names = F)

#0x significant tests


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfexon04_conc <- Shu_fish_list_df
Shu_fish_list_same_dfexon04_conc$eQTLoverlaps <- "Exon"

#write.csv(Shu_fish_list_same_dfexon, "Shu_fish_list_same_dfexon.csv", row.names = F)


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfpromExon04_conc <- Shu_fish_list_df
Shu_fish_list_same_dfpromExon04_conc$eQTLoverlaps <- "PromoterExon"

#write.csv(Shu_fish_list_same_dfpromExon, "Shu_fish_list_same_dfpromExon.csv", row.names = F)


#intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfintron04_conc <- Shu_fish_list_df
Shu_fish_list_same_dfintron04_conc$eQTLoverlaps <- "Intron"

#write.csv(Shu_fish_list_same_dfintron, "Shu_fish_list_same_dfintron.csv", row.names = F)


#locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dflocus04_conc <- Shu_fish_list_df
Shu_fish_list_same_dflocus04_conc$eQTLoverlaps <- "Locus"

#write.csv(Shu_fish_list_same_dflocus, "Shu_fish_list_same_dflocus.csv", row.names = F)


#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfsplice04_conc <- Shu_fish_list_df
Shu_fish_list_same_dfsplice04_conc$eQTLoverlaps <- "Splice"

#write.csv(Shu_fish_list_same_dfsplice, "Shu_fish_list_same_dfsplice.csv", row.names = F)


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfTTS04_conc <- Shu_fish_list_df
Shu_fish_list_same_dfTTS04_conc$eQTLoverlaps <- "TTS"

#write.csv(Shu_fish_list_same_dfTTS, "Shu_fish_list_same_dfTTS.csv", row.names = F)

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

#summary:
#some sig from whole locus only




#### 4a INDUCED EARLY SAME TIMEFRAME cisLnc pairs -first optimisation ####

#double check numbers
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
260/8116
length(unique(AllLNC_AllPCG_2d3d$EnsID.y))#5199 genes returned

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
)
#167 testable

#ready to run tests as dpair_range_bin#ready to run tests as done on GTEX:
Shu_fish_list <- list()

#allow finding pairs in each:
Shu_LncVar$pairs <- paste(Shu_LncVar$EnsID, Shu_LncVar$EnsID.y, sep = "-")

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
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfprom04up <- Shu_fish_list_df
Shu_fish_list_same_dfprom04up$eQTLoverlaps <- "Promoter"

#write.csv(Shu_fish_list_same_dfprom, "Shu_fish_list_same_dfprom.csv", row.names = F)

#0x significant tests


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfexon04up <- Shu_fish_list_df
Shu_fish_list_same_dfexon04up$eQTLoverlaps <- "Exon"

#write.csv(Shu_fish_list_same_dfexon, "Shu_fish_list_same_dfexon.csv", row.names = F)


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfpromExon04up <- Shu_fish_list_df
Shu_fish_list_same_dfpromExon04up$eQTLoverlaps <- "PromoterExon"

#write.csv(Shu_fish_list_same_dfpromExon, "Shu_fish_list_same_dfpromExon.csv", row.names = F)


#intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfintron04up <- Shu_fish_list_df
Shu_fish_list_same_dfintron04up$eQTLoverlaps <- "Intron"

#write.csv(Shu_fish_list_same_dfintron, "Shu_fish_list_same_dfintron.csv", row.names = F)


#locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dflocus04up <- Shu_fish_list_df
Shu_fish_list_same_dflocus04up$eQTLoverlaps <- "Locus"

#write.csv(Shu_fish_list_same_dflocus, "Shu_fish_list_same_dflocus.csv", row.names = F)


#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfsplice04up <- Shu_fish_list_df
Shu_fish_list_same_dfsplice04up$eQTLoverlaps <- "Splice"

#write.csv(Shu_fish_list_same_dfsplice, "Shu_fish_list_same_dfsplice.csv", row.names = F)


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfTTS04up <- Shu_fish_list_df
Shu_fish_list_same_dfTTS04up$eQTLoverlaps <- "TTS"

#write.csv(Shu_fish_list_same_dfTTS, "Shu_fish_list_same_dfTTS.csv", row.names = F)

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

#summary:
#some sig from whole locus only, OR reaches 1.45

#positive: ORs do increase a bit, one OR of 1.45 (still modest compared to previous)

#negative: less significance, distance appears to be a factor




#### 5a REPRESSED EARLY SAME TIMEFRAME cisLnc pairs -first optimisation ####

#double check numbers
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
260/8116
length(unique(AllLNC_AllPCG_2d3d$EnsID.y))#5199 genes returned

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

#ready to run tests as dpair_range_bin#ready to run tests as done on GTEX:
Shu_fish_list <- list()

#allow finding pairs in each:
Shu_LncVar$pairs <- paste(Shu_LncVar$EnsID, Shu_LncVar$EnsID.y, sep = "-")

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
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfprom04down <- Shu_fish_list_df
Shu_fish_list_same_dfprom04down$eQTLoverlaps <- "Promoter"

#write.csv(Shu_fish_list_same_dfprom, "Shu_fish_list_same_dfprom.csv", row.names = F)

#0x significant tests


#exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfexon04down <- Shu_fish_list_df
Shu_fish_list_same_dfexon04down$eQTLoverlaps <- "Exon"

#write.csv(Shu_fish_list_same_dfexon, "Shu_fish_list_same_dfexon.csv", row.names = F)


#promoter + exon
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfpromExon04down <- Shu_fish_list_df
Shu_fish_list_same_dfpromExon04down$eQTLoverlaps <- "PromoterExon"

#write.csv(Shu_fish_list_same_dfpromExon, "Shu_fish_list_same_dfpromExon.csv", row.names = F)


#intron
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfintron04down <- Shu_fish_list_df
Shu_fish_list_same_dfintron04down$eQTLoverlaps <- "Intron"

#write.csv(Shu_fish_list_same_dfintron, "Shu_fish_list_same_dfintron.csv", row.names = F)


#locus
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dflocus04down <- Shu_fish_list_df
Shu_fish_list_same_dflocus04down$eQTLoverlaps <- "Locus"

#write.csv(Shu_fish_list_same_dflocus, "Shu_fish_list_same_dflocus.csv", row.names = F)


#Splice
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfsplice04down <- Shu_fish_list_df
Shu_fish_list_same_dfsplice04down$eQTLoverlaps <- "Splice"

#write.csv(Shu_fish_list_same_dfsplice, "Shu_fish_list_same_dfsplice.csv", row.names = F)


#TTS
for (j in 1:length(pThresh_df$pThresh)){
  
  selectedp <- pThresh_df[j,1]
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(AllLNC_AllPCG_2d3d_04, EnsID_merge %in% Shu_exprsG, 
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

Shu_fish_list_same_dfTTS04down <- Shu_fish_list_df
Shu_fish_list_same_dfTTS04down$eQTLoverlaps <- "TTS"

#write.csv(Shu_fish_list_same_dfTTS, "Shu_fish_list_same_dfTTS.csv", row.names = F)



#### takeaways from above ####

#unexpectedly, splitting to individual lnc groups doesn't improve as much as it did in GTEX

#nor do other factors like increasing eQTLp shortening the distance or increasing the pcg xp

#or going to specific types like co-induced/concordant

#seems like GTEX data has better predictive value, but on less pairs, 

#whilst the biobank gives a lot of pairs, with low OR - be wary of including...


#### build a supplementary table ####

AllLNC_AllPCG_2d3d <- filter(AllLNC_AllPCG_2d3d, pair_range <1000000)
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
260/8116

#more effective than running numerous times
#for each of the lncRNA overlap types
#for each lncRNA cluster
#can use the same table as before
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

Shu_LncVar$pairs <- paste(Shu_LncVar$EnsID, Shu_LncVar$EnsID.y, sep = "-")

AllLNC_AllPCG_2d3d$pairs_merge <- gsub("\\.[0-9]*$", "", AllLNC_AllPCG_2d3d$pairs)
AllLNC_AllPCG_2d3d$EnsID_merge <- gsub("\\.[0-9]*", "", AllLNC_AllPCG_2d3d$EnsID.y)

Shu_SuppTable <- list()

selectedp <- 0.05

for (j in 1:length(GTEX_runs_parameters$Lnc_cluster)){
  
  selectedLncs <- GTEX_runs_selection[[ GTEX_runs_parameters[j,1] ]]
  selectedBackground <- GTEX_runs_background[[ GTEX_runs_parameters[j,2] ]]
  selectedOverlap <- GTEX_runs_parameters[j,3]
  
  #select pairs to test
  triali <- filter(AllLNC_AllPCG_2d3d, EnsID %in% selectedBackground$EnsID, EnsID.y %in% selectedBackground$EnsID)
  
  #triali <- filter(triali, !grepl("MSTRG.27167|MSTRG.27169", EnsID))
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar,
                                                        grepl(selectedOverlap, OverlapType), 
                                                        pvalue < selectedp,
                                                        TotalOverlaps == 1)$pairs] <- "Yes"
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(triali, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
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
    
    Shu_SuppTable[[j]] <- c(
      fisher.test(data.frame("cisLnc" = c(a, b-a),
                             "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
      fisher.test(data.frame("cisLnc" = c(a, b-a),
                             "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
      a,b,c,d)
  } else {
    Shu_SuppTable[[j]] <- c(NA, NA, NA, NA, NA, NA)
  }
  
}

length(Shu_SuppTable)
GTEX_runs_parameters$Overlap_type_name <- GTEX_runs_parameters$Overlap_type
GTEX_runs_parameters$Overlap_type_name[1:10] <- "Locus"
names(Shu_SuppTable) <- paste("Run", GTEX_runs_parameters[,c(1)], GTEX_runs_parameters[,c(4)], sep = "_")
Shu_SuppTable_df <- as.data.frame(t(bind_rows(Shu_SuppTable, .id = "Run")))
Shu_SuppTable_df$Run <- row.names(Shu_SuppTable_df)
colnames(Shu_SuppTable_df) <- c("p", "OR", "a", "b", "c", "d", "Run")
Shu_SuppTable_df$tissue <- "SMC_biobank"
Shu_SuppTable_df <- Shu_SuppTable_df[,c(7,1:6,8)]

#write.csv(Shu_SuppTable_df, "Shu_SuppTable_df.csv", row.names = F)


#### build a supplementary table -250kbp ####

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
#for each of the lncRNA overlap types
#for each lncRNA cluster
#can use the same table as before
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
                            "All-same/delayed"= filter(AllLNC_AllPCG_2d3d_250,
                                                       (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                     fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE$EnsID)) |
                                                         (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                       fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_within_8$EnsID, 
                                                                                                                        fpkm_allGDE_within_24$EnsID)) |
                                                         (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                       fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                                         fpkm_allGDE_Downwithin_24$EnsID)))
)

Shu_LncVar$pairs <- paste(Shu_LncVar$EnsID, Shu_LncVar$EnsID.y, sep = "-")

AllLNC_AllPCG_2d3d_250$pairs_merge <- gsub("\\.[0-9]*$", "", AllLNC_AllPCG_2d3d_250$pairs)
AllLNC_AllPCG_2d3d_250$EnsID_merge <- gsub("\\.[0-9]*", "", AllLNC_AllPCG_2d3d_250$EnsID.y)

Shu_SuppTable <- list()

for (j in 1:length(GTEX_runs_parameters$Lnc_cluster)){
  
  selectedLncs <- GTEX_runs_selection[[ GTEX_runs_parameters[j,1] ]]
  selectedBackground <- GTEX_runs_background[[ GTEX_runs_parameters[j,2] ]]
  selectedOverlap <- GTEX_runs_parameters[j,3]
  
  #select pairs to test
  triali <- filter(AllLNC_AllPCG_2d3d_250, EnsID %in% selectedBackground$EnsID, EnsID.y %in% selectedBackground$EnsID)
  
  #triali <- filter(triali, !grepl("MSTRG.27167|MSTRG.27169", EnsID))
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar,
                                                              grepl(selectedOverlap, OverlapType), 
                                                              pvalue < selectedp,
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < selectedp)
  
  triali <- filter(triali, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
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
    
    Shu_SuppTable[[j]] <- c(
      fisher.test(data.frame("cisLnc" = c(a, b-a),
                             "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
      fisher.test(data.frame("cisLnc" = c(a, b-a),
                             "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
      a,b,c,d)
  } else {
    Shu_SuppTable[[j]] <- c(NA, NA, NA, NA, NA, NA)
  }
  
}

length(Shu_SuppTable)
GTEX_runs_parameters$Overlap_type_name <- GTEX_runs_parameters$Overlap_type
GTEX_runs_parameters$Overlap_type_name[1:10] <- "Locus"
names(Shu_SuppTable) <- paste("Run", GTEX_runs_parameters[,c(1)], GTEX_runs_parameters[,c(4)], sep = "_")
Shu_SuppTable_df_250 <- as.data.frame(t(bind_rows(Shu_SuppTable, .id = "Run")))
Shu_SuppTable_df_250$Run <- row.names(Shu_SuppTable_df_250)
colnames(Shu_SuppTable_df_250) <- c("p", "OR", "a", "b", "c", "d", "Run")
Shu_SuppTable_df_250$tissue <- "SMC_biobank"
Shu_SuppTable_df_250 <- Shu_SuppTable_df_250[,c(7,1:6,8)]

#write.csv(Shu_SuppTable_df, "Shu_SuppTable_df.csv", row.names = F)

#### build a supplementary table -250kbp p = 0.025####

#would this limit zscan locus?
dim(filter(AllLNC_AllPCG_2d3d, grepl("MSTRG.27167|MSTRG.27169", EnsID)))#28 pairs

AllLNC_AllPCG_2d3d_250 <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <250)
dim(filter(AllLNC_AllPCG_2d3d_250, grepl("MSTRG.27167|MSTRG.27169", EnsID)))#3 pairs only

AllLNC_AllPCG_2d3d_500 <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <500)
dim(filter(AllLNC_AllPCG_2d3d_500, grepl("MSTRG.27167|MSTRG.27169", EnsID)))#8 pairs

#try 250kbp, stringent
dim(filter(AllLNC_AllPCG_2d3d_250, !loopMethod == "Neither"))
185/2563

#extra power in this SMC biobank, extra eQTLs, extra need for thresholding?
#try p = 0.025

#more effective than running numerous times
#for each of the lncRNA overlap types
#for each lncRNA cluster
#can use the same table as before
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

Shu_LncVar$pairs <- paste(Shu_LncVar$EnsID, Shu_LncVar$EnsID.y, sep = "-")

AllLNC_AllPCG_2d3d_250$pairs_merge <- gsub("\\.[0-9]*$", "", AllLNC_AllPCG_2d3d_250$pairs)
AllLNC_AllPCG_2d3d_250$EnsID_merge <- gsub("\\.[0-9]*", "", AllLNC_AllPCG_2d3d_250$EnsID.y)

Shu_SuppTable <- list()

for (j in 1:length(GTEX_runs_parameters$Lnc_cluster)){
  
  selectedLncs <- GTEX_runs_selection[[ GTEX_runs_parameters[j,1] ]]
  selectedBackground <- GTEX_runs_background[[ GTEX_runs_parameters[j,2] ]]
  selectedOverlap <- GTEX_runs_parameters[j,3]
  
  #select pairs to test
  triali <- filter(AllLNC_AllPCG_2d3d_250, EnsID %in% selectedBackground$EnsID, EnsID.y %in% selectedBackground$EnsID)
  
  #triali <- filter(triali, !grepl("MSTRG.27167|MSTRG.27169", EnsID))
  
  triali$eQTL_validated_tissue <- "No"
  triali$eQTL_validated_tissue[triali$pairs_merge %in% filter(Shu_LncVar,
                                                              grepl(selectedOverlap, OverlapType), 
                                                              pvalue < 0.025,
                                                              TotalOverlaps == 1)$pairs] <- "Yes"
  
  Shu_exprsG <- Shu_SVSMCpairedPCG_exprsTPM$EnsID_merge[rowMeans(as.matrix(Shu_SVSMCpairedPCG_exprsTPM[,-1]))> 1]
  Shu_allVar_pThresh <- filter(Shu_allVar_p, pvalue < 0.025)
  
  triali <- filter(triali, EnsID_merge %in% Shu_exprsG, 
                   EnsName.y %in% Shu_allVar_pThresh$GeneName)
  
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
    
    Shu_SuppTable[[j]] <- c(
      fisher.test(data.frame("cisLnc" = c(a, b-a),
                             "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p.value,
      fisher.test(data.frame("cisLnc" = c(a, b-a),
                             "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate,
      a,b,c,d)
  } else {
    Shu_SuppTable[[j]] <- c(NA, NA, NA, NA, NA, NA)
  }
  
}

length(Shu_SuppTable)
GTEX_runs_parameters$Overlap_type_name <- GTEX_runs_parameters$Overlap_type
GTEX_runs_parameters$Overlap_type_name[1:10] <- "Locus"
names(Shu_SuppTable) <- paste("Run", GTEX_runs_parameters[,c(1)], GTEX_runs_parameters[,c(4)], sep = "_")
Shu_SuppTable_df_250_p025 <- as.data.frame(t(bind_rows(Shu_SuppTable, .id = "Run")))
Shu_SuppTable_df_250_p025$Run <- row.names(Shu_SuppTable_df_250_p025)
colnames(Shu_SuppTable_df_250_p025) <- c("p", "OR", "a", "b", "c", "d", "Run")
Shu_SuppTable_df_250_p025$tissue <- "SMC_biobank"
Shu_SuppTable_df_250_p025 <- Shu_SuppTable_df_250_p025[,c(7,1:6,8)]

#write.csv(Shu_SuppTable_df_250_p025, "Shu_SuppTable_df_250_p025.csv", row.names = F)
