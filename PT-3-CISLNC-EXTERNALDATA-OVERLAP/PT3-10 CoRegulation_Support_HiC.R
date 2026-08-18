#Predicting cis-acting lncRNA via the Zhao HiC set + testing some associations between lnc-PCG regulation pattern over time
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(rcompanion)
library(ggbeeswarm)
library(rtracklayer)
library(GenomicRanges)
library(reshape2)

#### import and set-up tables without artefacts ####

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

#above lines seperates genes into distinct buckets:
table(fpkm_allGDE$RegulationStart) 

#from neighbour associations, all lncRNA-PCGs within a 250kbp window (2d/naive approach)
AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_info_2026.csv", header = T)

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


#### addtl. pairs from caSMC loops - convert Zhao et al. to hg38 ####

#Quanyi used 5kb resolution to get loops from FitHiC

#import of FitHiC loops across pooled samples
caSMC_Zhao_Loops1 <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Candidate gene selection/epigenetics/enhancers_superenhancers/CASMC ATAC HiC Zhao 2020/GSM4212721_FitHiC_allloops.txt",
                                header = F)
#colnames from https://github.com/ay-lab/fithic#output
colnames(caSMC_Zhao_Loops1) <- c("chr1",	"fragmentMid1",	"chr2",	"fragmentMid2",	"contactCount",	"p-value", "q-value")

#
chain <- import.chain("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/hg19ToHg38.over.chain")
trial <- makeGRangesFromDataFrame(caSMC_Zhao_Loops1, start.field = "fragmentMid1", end.field = "fragmentMid1", seqnames.field = "chr1", keep.extra.columns = T)
trial$locID_hg19 <- paste(caSMC_Zhao_Loops1$chr1, caSMC_Zhao_Loops1$fragmentMid1, caSMC_Zhao_Loops1$fragmentMid2, sep = "_")
#liftOver fragment peak 1
triali <- data.frame(unlist(liftOver(trial, chain)))[,-c(3:5)]
colnames(triali)[1:2] <- c("chr1", "fragmentMid1")
trialii <- makeGRangesFromDataFrame(triali, start.field = "fragmentMid2", end.field = "fragmentMid2", 
                                    seqnames.field = "chr2", keep.extra.columns = T)
#liftOver fragment peak 2
trialiii <- data.frame(unlist(liftOver(trialii, chain)))[,-c(3:5)]
colnames(trialiii)[1:2] <- c("chr2", "fragmentMid2")

#changes in distance with liftOver? hg38 dist:
trialiii$dist_hg38 <- (trialiii$fragmentMid2 - trialiii$fragmentMid1)/1000

caSMC_Zhao_Loops1$locID_hg19 <- paste(caSMC_Zhao_Loops1$chr1, 
                                      caSMC_Zhao_Loops1$fragmentMid1, caSMC_Zhao_Loops1$fragmentMid2, sep = "_")

trialiv <- unique(merge(trialiii, caSMC_Zhao_Loops1, by = "locID_hg19"))
trialiv$dist_hg19 <- (trialiv$fragmentMid2.y - trialiv$fragmentMid1.y)/1000

#mostly fine though with some outliers
ggplot(trialiv) + aes(x = dist_hg38, y = dist_hg19) +
  geom_point(alpha = 0.2)

summary(trialiv$dist_hg19) #330 median
summary(trialiv$dist_hg38) #325 median

#distance changes are rendered invalid - key input used in FitHiC confidence assignment:
#hence these are the valid lifted over loops:
trialiv <- filter(trialiv, dist_hg38 == dist_hg19)

caSMC_Zhao_Loops1_filt <- filter(trialiii, locID_hg19 %in% trialiv$locID_hg19)

#import of HiCCUPS loops - alternative process/tool - across pooled samples
caSMC_Zhao_Loops2 <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/BioinfGroupResources/Candidate gene selection/epigenetics/enhancers_superenhancers/CASMC ATAC HiC Zhao 2020/GSM4212721_HiCCUPS_allloop.bedpe",
                                header = T)
caSMC_Zhao_Loops2 <- caSMC_Zhao_Loops2[-1,]
caSMC_Zhao_Loops2 <- caSMC_Zhao_Loops2[,-c(7:10)]

colnames(caSMC_Zhao_Loops2) <- c("chromosome1", "x1",    "x2",    "chromosome2",    "y1",    "y2",    "color",    "observed",
                                 "expected_bottom_left",    "expected_donut",    "expected_horizontal",    "expected_vertical",    
                                 "fdr_bottom_left",    "fdr_donut",    "fdr_horizontal",    "fdr_vertical",    
                                 "number_collapsed",    "centroid1",    "centroid2",    "radius")

caSMC_Zhao_Loops2$chromosome1 <- paste("chr", caSMC_Zhao_Loops2$chromosome1, sep = "")
caSMC_Zhao_Loops2$chromosome2 <- paste("chr", caSMC_Zhao_Loops2$chromosome2, sep = "")

#
#chain <- import.chain("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/hg19ToHg38.over.chain")
trial <- makeGRangesFromDataFrame(caSMC_Zhao_Loops2, start.field = "centroid1", end.field = "centroid1", 
                                  seqnames.field = "chromosome1", keep.extra.columns = T)
trial$locID_hg19 <- paste(caSMC_Zhao_Loops2$chromosome1, caSMC_Zhao_Loops2$centroid1, caSMC_Zhao_Loops2$centroid2, sep = "_")
#liftOver fragment peak 1
triali <- data.frame(unlist(liftOver(trial, chain)))[,-c(3:5)]
colnames(triali)[1:2] <- c("chr1", "fragmentMid1")

trialii <- makeGRangesFromDataFrame(triali, start.field = "centroid2", end.field = "centroid2", 
                                    seqnames.field = "chromosome2", keep.extra.columns = T)
#liftOver fragment peak 2
trialiii <- data.frame(unlist(liftOver(trialii, chain)))[,-c(3:5)]
colnames(trialiii)[1:2] <- c("chr2", "fragmentMid2")

#changes in distance with liftOver? hg38 dist:
trialiii$dist_hg38 <- (trialiii$fragmentMid2 - trialiii$fragmentMid1)/1000

caSMC_Zhao_Loops2$locID_hg19 <- paste(caSMC_Zhao_Loops2$chromosome1, 
                                      caSMC_Zhao_Loops2$centroid1, caSMC_Zhao_Loops2$centroid2, sep = "_")

trialiv <- unique(merge(trialiii, caSMC_Zhao_Loops2, by = "locID_hg19"))
trialiv$dist_hg19 <- (trialiv$centroid2 - trialiv$centroid1)/1000

#mostly fine though with some outliers
ggplot(trialiv) + aes(x = dist_hg38, y = dist_hg19) +
  geom_point(alpha = 0.2)

summary(trialiv$dist_hg19) #295 median
summary(trialiv$dist_hg38) #292.5 median

#distance changes are rendered invalid - key input used in FitHiC confidence assignment:
trialiv <- filter(trialiv, dist_hg38 == dist_hg19)

caSMC_Zhao_Loops2_filt <- filter(trialiii, locID_hg19 %in% trialiv$locID_hg19)

#much less HiCCUPs (~11k) than FitHiC (~100k)


#### identify lncRNA contact sites ####

#process to get the overlap regions

#get all tx for the gene

#flatten into a single region with one 5' and 3' extremity

#if considering involving some promoter analysis later - remember their may be multiple sites within these regions


#now convert to ranges to be overlapped
#isolate loops, where one or other fragment overlaps an expressed lncRNA:
allGB_LNCS <- filter(allGB, EnsID %in% filter(fpkm_allG, grepl("Lnc|fide", GeneClassUpdate))$EnsID)
length(unique(allGB_LNCS$EnsID))#597 as expected

allGB_LNCS_GR <- makeGRangesFromDataFrame(allGB_LNCS[,-7], 
                                          start.field = "Tx_start", 
                                          end.field = "Tx_stop", 
                                          seqnames.field = "chr", 
                                          strand.field = "str", keep.extra.columns = T)

#extend on TSS by 2000kbp
allGB_LNCS_GRpromoters <- promoters(allGB_LNCS_GR, upstream = 2000)

#take union
allGB_LNCS_GR <- GenomicRanges::punion(allGB_LNCS_GR, allGB_LNCS_GRpromoters)
allGB_LNCS_GR$EnsID <- allGB_LNCS_GRpromoters$EnsID
allGB_LNCS_GR$MSTRG_Tx_ID <- allGB_LNCS_GRpromoters$MSTRG_Tx_ID

#FitHiC frag1 to lncRNAs:
trial <- makeGRangesFromDataFrame(caSMC_Zhao_Loops1_filt, 
                                  start.field = "fragmentMid1", 
                                  end.field = "fragmentMid1", seqnames.field = "chr1", keep.extra.columns = T)
caSMC_Zhao_Loops1_filt_f1 <- flank(trial, width = 2500, both = T)

trial <- makeGRangesFromDataFrame(caSMC_Zhao_Loops1_filt, 
                                  start.field = "fragmentMid2", 
                                  end.field = "fragmentMid2", seqnames.field = "chr2", keep.extra.columns = T)
caSMC_Zhao_Loops1_filt_f2 <- flank(trial, width = 2500, both = T)

#query lncNA granges with FitHiC Frag1:
Loop1Frag1index <- findOverlaps(query = caSMC_Zhao_Loops1_filt_f1, subject = allGB_LNCS_GR)
Loop1Frag1overlaps <- unique(data.frame("Loop1Frag1" = caSMC_Zhao_Loops1_filt_f1$locID_hg19[queryHits(Loop1Frag1index)],
                                        "LncRNA" = allGB_LNCS_GR$MSTRG_Tx_ID[subjectHits(Loop1Frag1index)]))
Loop1Frag1overlaps <- merge(Loop1Frag1overlaps, unique(fpkm_allG[,c(2,47)]), by.x = "LncRNA", by.y = "MSTRG_Tx_ID")
head(Loop1Frag1overlaps)
caSMC_Zhao_Loops1_filt_f1[caSMC_Zhao_Loops1_filt_f1$locID_hg19 == "chr1_1627500_1682500"]
allGB_LNCS_GR[allGB_LNCS_GR$MSTRG_Tx_ID == "MSTRG.124.1"]
length(unique(Loop1Frag1overlaps$Loop1Frag1))#992 contacts found within a gene in
length(unique(Loop1Frag1overlaps$EnsID))#205 lncRNAs found with a HiC contact in

#FitHiC frag2 to lncRNAs:
Loop1Frag2index <- findOverlaps(query = caSMC_Zhao_Loops1_filt_f2, subject = allGB_LNCS_GR)
Loop1Frag2overlaps <- unique(data.frame("Loop1Frag2" = caSMC_Zhao_Loops1_filt_f2$locID_hg19[queryHits(Loop1Frag2index)],
                                        "LncRNA" = allGB_LNCS_GR$MSTRG_Tx_ID[subjectHits(Loop1Frag2index)]))
Loop1Frag2overlaps <- merge(Loop1Frag2overlaps, unique(fpkm_allG[,c(2,47)]), by.x = "LncRNA", by.y = "MSTRG_Tx_ID")
head(Loop1Frag2overlaps)
caSMC_Zhao_Loops1_filt_f2[caSMC_Zhao_Loops1_filt_f2$locID_hg19 == "chr1_16117500_16177500"]
allGB_LNCS_GR[allGB_LNCS_GR$EnsID == "ENSG00000179743.4"]
length(unique(Loop1Frag2overlaps$Loop1Frag2))#1124 loops found contacting a total of:
length(unique(Loop1Frag2overlaps$EnsID))#202 lncRNAs found with a HiC contact in Frag2

#loops with a lncRNA in one or other fragment:
caSMC_Zhao_Loops1_filt_lncs <- filter(caSMC_Zhao_Loops1_filt, 
                                      locID_hg19 %in% Loop1Frag1overlaps$Loop1Frag1 | locID_hg19 %in% Loop1Frag2overlaps$Loop1Frag2)

#HiCCUP frag1 to lncRNAs:
trial <- makeGRangesFromDataFrame(caSMC_Zhao_Loops2_filt, 
                                  start.field = "fragmentMid1", 
                                  end.field = "fragmentMid1", seqnames.field = "chr1", keep.extra.columns = T)
caSMC_Zhao_Loops2_filt_f1 <- flank(trial, width = 2500, both = T)

trial <- makeGRangesFromDataFrame(caSMC_Zhao_Loops2_filt, 
                                  start.field = "fragmentMid2", 
                                  end.field = "fragmentMid2", seqnames.field = "chr2", keep.extra.columns = T)
caSMC_Zhao_Loops2_filt_f2 <- flank(trial, width = 2500, both = T)

#query lncNA granges with HiCCUP Frag1:
Loop2Frag1index <- findOverlaps(query = caSMC_Zhao_Loops2_filt_f1, subject = allGB_LNCS_GR)
Loop2Frag1overlaps <- unique(data.frame("Loop2Frag1" = caSMC_Zhao_Loops2_filt_f1$locID_hg19[queryHits(Loop2Frag1index)],
                                        "LncRNA" = allGB_LNCS_GR$MSTRG_Tx_ID[subjectHits(Loop2Frag1index)]))
Loop2Frag1overlaps <- merge(Loop2Frag1overlaps, unique(fpkm_allG[,c(2,47)]), by.x = "LncRNA", by.y = "MSTRG_Tx_ID")
head(Loop2Frag1overlaps)
caSMC_Zhao_Loops2_filt_f1[caSMC_Zhao_Loops2_filt_f1$locID_hg19 == "chr10_81590000_81645000"]
allGB_LNCS_GR[allGB_LNCS_GR$EnsID == "ENSG00000225484.6"]
length(unique(Loop2Frag1overlaps$Loop2Frag1))#141 loops contacting inside a genebody
length(unique(Loop2Frag1overlaps$EnsID))#96 lncRNAs found with a HiC contact in


#HiCCUP frag2 to lncRNAs:
Loop2Frag2index <- findOverlaps(query = caSMC_Zhao_Loops2_filt_f2, subject = allGB_LNCS_GR)
Loop2Frag2overlaps <- unique(data.frame("Loop2Frag2" = caSMC_Zhao_Loops2_filt_f2$locID_hg19[queryHits(Loop2Frag2index)],
                                        "LncRNA" = allGB_LNCS_GR$MSTRG_Tx_ID[subjectHits(Loop2Frag2index)]))
Loop2Frag2overlaps <- merge(Loop2Frag2overlaps, unique(fpkm_allG[,c(2,47)]), by.x = "LncRNA", by.y = "MSTRG_Tx_ID")
head(Loop2Frag2overlaps)
caSMC_Zhao_Loops2_filt_f2[caSMC_Zhao_Loops2_filt_f2$locID_hg19 == "chr10_31075000_31615000"]
allGB_LNCS_GR[allGB_LNCS_GR$EnsID == "ENSG00000237036.4"]
length(unique(Loop2Frag2overlaps$Loop2Frag2))#123 HiC found with a genebody in
length(unique(Loop2Frag2overlaps$EnsID))#82 lncRNAs found with a HiC contact in Frag2

#loops with a lncRNA in one or other fragment:
caSMC_Zhao_Loops2_filt_lncs <- filter(caSMC_Zhao_Loops2_filt, 
                                      locID_hg19 %in% Loop2Frag1overlaps$Loop2Frag1 | locID_hg19 %in% Loop2Frag2overlaps$Loop2Frag2)

#expecting similar lncRNAs from FitHiC and HiCCUPs - but maybe more in FitHiC
FitHiC_LNCS <- unique(c(Loop1Frag1overlaps$EnsID, Loop1Frag2overlaps$EnsID))
HiCCUP_LNCS <- unique(c(Loop2Frag1overlaps$EnsID, Loop2Frag2overlaps$EnsID))

#library("ggVennDiagram")
ggVennDiagram::ggVennDiagram(list("FitHiC_lncs" = FitHiC_LNCS, "HiCCUP_lncs" = HiCCUP_LNCS), label_alpha = 0)
#136 shared

#clear that lncRNAs contacted by both approaches have lower p values

#HiCCUP finds nearly all FitHiC lncRNAs, FitHiC finds more, suspect FitHiC is permissive


#### HiC validation 1: position and strength of contacts relative to lncRNA gene ####

#where are contacts? lnc promoter? genebody? are stronger found at promoter? 
#non-random association would be promising
#call promoter 0% and TES 100% and find fragment centre within:

#add lncRNA coords to the matched loops:
Loop1Frag1overlaps <- merge(Loop1Frag1overlaps, allGB_LNCS[,c(1,3:6)], by.x = "LncRNA", by.y = "MSTRG_Tx_ID")
#add frag1 centroid:
caSMC_Zhao_Loops1_filt_f1$centre <- as.data.frame(caSMC_Zhao_Loops1_filt_f1)[2] + 2500
Loop1Frag1overlaps <- merge(Loop1Frag1overlaps, as.data.frame(caSMC_Zhao_Loops1_filt_f1)[,c(8,10,11,13)], by.x = "Loop1Frag1", by.y = "locID_hg19")

#lnc locus length
Loop1Frag1overlaps$lncLocLength <- Loop1Frag1overlaps$Tx_stop - Loop1Frag1overlaps$Tx_start
#dist from contact point to stop point
Loop1Frag1overlaps$lncLocContact <- Loop1Frag1overlaps$Tx_stop - Loop1Frag1overlaps$start
#as percentage of lnc total
Loop1Frag1overlaps$lncLocContact_perc <- Loop1Frag1overlaps$lncLocContact/Loop1Frag1overlaps$lncLocLength*100

#for +ve, if contact lower than start then call perc -10% - so they stick out clearly on any graph
Loop1Frag1overlaps$lncLocContact_perc[Loop1Frag1overlaps$str == "+" & Loop1Frag1overlaps$start < Loop1Frag1overlaps$Tx_start] <- -10
#for -ve, if contact is higher than stop, then -10%
Loop1Frag1overlaps$lncLocContact_perc[Loop1Frag1overlaps$str == "-" & Loop1Frag1overlaps$start > Loop1Frag1overlaps$Tx_stop] <- -10

#for +ve, subtract % from 100 so pos and neg are measured from TSS
Loop1Frag1overlaps$lncLocContact_perc[Loop1Frag1overlaps$str == "+"] <- 100- Loop1Frag1overlaps$lncLocContact_perc[Loop1Frag1overlaps$str == "+"]

#for +ve, if contact higher than end, call +110%
Loop1Frag1overlaps$lncLocContact_perc[Loop1Frag1overlaps$str == "+" & Loop1Frag1overlaps$start > Loop1Frag1overlaps$Tx_stop] <- 110
#for -ve, if contact is lower than start, then 110%
Loop1Frag1overlaps$lncLocContact_perc[Loop1Frag1overlaps$str == "-" & Loop1Frag1overlaps$start < Loop1Frag1overlaps$Tx_start] <- 110

#for FitHiC frags, the contacts seem to be biased to TES but also TSS
ggplot(Loop1Frag1overlaps) + aes(x = lncLocContact_perc) +
  geom_histogram()
#described previously? expecting same from frag2 FitHiC


#add lncRNA coords to the matched loops:
Loop1Frag2overlaps <- merge(Loop1Frag2overlaps, allGB_LNCS[,c(1,3:6)], by.x = "LncRNA", by.y = "MSTRG_Tx_ID")
#add frag1 centroid:
caSMC_Zhao_Loops1_filt_f2$centre <- as.data.frame(caSMC_Zhao_Loops1_filt_f2)[2] + 2500
Loop1Frag2overlaps <- merge(Loop1Frag2overlaps, as.data.frame(caSMC_Zhao_Loops1_filt_f2)[,c(8,10,11,13)], by.x = "Loop1Frag2", by.y = "locID_hg19")

#lnc locus length
Loop1Frag2overlaps$lncLocLength <- Loop1Frag2overlaps$Tx_stop - Loop1Frag2overlaps$Tx_start
#dist from contact point to stop point
Loop1Frag2overlaps$lncLocContact <- Loop1Frag2overlaps$Tx_stop - Loop1Frag2overlaps$start
#as percentage of lnc total
Loop1Frag2overlaps$lncLocContact_perc <- Loop1Frag2overlaps$lncLocContact/Loop1Frag2overlaps$lncLocLength*100
#for +ve, if contact lower than start then call perc -10% - so they stick out clearly on any graph
Loop1Frag2overlaps$lncLocContact_perc[Loop1Frag2overlaps$str == "+" & Loop1Frag2overlaps$start < Loop1Frag2overlaps$Tx_start] <- -10
#for -ve, if contact is higher than stop, then -10%
Loop1Frag2overlaps$lncLocContact_perc[Loop1Frag2overlaps$str == "-" & Loop1Frag2overlaps$start > Loop1Frag2overlaps$Tx_stop] <- -10

#for +ve, subtract % from 100 so pos and neg are treated same
Loop1Frag2overlaps$lncLocContact_perc[Loop1Frag2overlaps$str == "+"] <- 100- Loop1Frag2overlaps$lncLocContact_perc[Loop1Frag2overlaps$str == "+"]

#for +ve, if contact higher than end, call +110%
Loop1Frag2overlaps$lncLocContact_perc[Loop1Frag2overlaps$str == "+" & Loop1Frag2overlaps$start > Loop1Frag2overlaps$Tx_stop] <- 110
#for -ve, if contact is lower than start, then 110%
Loop1Frag2overlaps$lncLocContact_perc[Loop1Frag2overlaps$str == "-" & Loop1Frag2overlaps$start < Loop1Frag2overlaps$Tx_start] <- 110

ggplot(Loop1Frag2overlaps) + aes(x = lncLocContact_perc) +
  geom_histogram()
#promoter bias lessened, TES bias still enriched


#HiCCUPs
#add lncRNA coords to the matched loops:
Loop2Frag1overlaps <- merge(Loop2Frag1overlaps, allGB_LNCS[,c(1,3:6)], by.x = "LncRNA", by.y = "MSTRG_Tx_ID")
#add frag1 centroid:
caSMC_Zhao_Loops2_filt_f1$centre <- as.data.frame(caSMC_Zhao_Loops2_filt_f1)[2] + 2500
Loop2Frag1overlaps <- merge(Loop2Frag1overlaps, as.data.frame(caSMC_Zhao_Loops2_filt_f1)[,c(13,18:21,24:26)], by.x = "Loop2Frag1", by.y = "locID_hg19")

#lnc locus length
Loop2Frag1overlaps$lncLocLength <- Loop2Frag1overlaps$Tx_stop - Loop2Frag1overlaps$Tx_start
#dist from contact point to stop point
Loop2Frag1overlaps$lncLocContact <- Loop2Frag1overlaps$Tx_stop - Loop2Frag1overlaps$start
#as percentage of lnc total
Loop2Frag1overlaps$lncLocContact_perc <- Loop2Frag1overlaps$lncLocContact/Loop2Frag1overlaps$lncLocLength * 100
#for +ve, if contact lower than start then call perc -10% - so they stick out clearly on any graph
Loop2Frag1overlaps$lncLocContact_perc[Loop2Frag1overlaps$str == "+" & Loop2Frag1overlaps$start < Loop2Frag1overlaps$Tx_start] <- - 10
#for -ve, if contact is higher than stop, then -10%
Loop2Frag1overlaps$lncLocContact_perc[Loop2Frag1overlaps$str == "-" & Loop2Frag1overlaps$start > Loop2Frag1overlaps$Tx_stop] <- - 10

#for +ve, subtract % from 100 so pos and neg are treated same
Loop2Frag1overlaps$lncLocContact_perc[Loop2Frag1overlaps$str == "+"] <- 100 - Loop2Frag1overlaps$lncLocContact_perc[Loop2Frag1overlaps$str == "+"]

#for +ve, if contact higher than end, call +110%
Loop2Frag1overlaps$lncLocContact_perc[Loop2Frag1overlaps$str == "+" & Loop2Frag1overlaps$start > Loop2Frag1overlaps$Tx_stop] <- 110
#for -ve, if contact is lower than start, then 110%
Loop2Frag1overlaps$lncLocContact_perc[Loop2Frag1overlaps$str == "-" & Loop2Frag1overlaps$start < Loop2Frag1overlaps$Tx_start] <- 110

#with HiCCUPs the trend to TTS or TES is much stronger, all gene body contacts are deprioritised
ggplot(Loop2Frag1overlaps) + aes(x = lncLocContact_perc) +
  geom_histogram()


#HiCCUPs Frag2
Loop2Frag2overlaps <- merge(Loop2Frag2overlaps, allGB_LNCS[,c(1,3:6)], by.x = "LncRNA", by.y = "MSTRG_Tx_ID")
#add frag1 centroid:
caSMC_Zhao_Loops2_filt_f2$centre <- as.data.frame(caSMC_Zhao_Loops2_filt_f2)[2] + 2500
Loop2Frag2overlaps <- merge(Loop2Frag2overlaps, as.data.frame(caSMC_Zhao_Loops2_filt_f2)[,c(13,18:21,24:26)], by.x = "Loop2Frag2", by.y = "locID_hg19")

#lnc locus length
Loop2Frag2overlaps$lncLocLength <- Loop2Frag2overlaps$Tx_stop - Loop2Frag2overlaps$Tx_start
#dist from contact point to stop point
Loop2Frag2overlaps$lncLocContact <- Loop2Frag2overlaps$Tx_stop - Loop2Frag2overlaps$start
#as percentage of lnc total
Loop2Frag2overlaps$lncLocContact_perc <- Loop2Frag2overlaps$lncLocContact/Loop2Frag2overlaps$lncLocLength * 100
#for +ve, if contact lower than start then call perc -10% - so they stick out clearly on any graph
Loop2Frag2overlaps$lncLocContact_perc[Loop2Frag2overlaps$str == "+" & Loop2Frag2overlaps$start < Loop2Frag2overlaps$Tx_start] <- - 10
#for -ve, if contact is higher than stop, then -10%
Loop2Frag2overlaps$lncLocContact_perc[Loop2Frag2overlaps$str == "-" & Loop2Frag2overlaps$start > Loop2Frag2overlaps$Tx_stop] <- - 10

#for +ve, subtract % from 100 so pos and neg are treated same
Loop2Frag2overlaps$lncLocContact_perc[Loop2Frag2overlaps$str == "+"] <- 100 - Loop2Frag2overlaps$lncLocContact_perc[Loop2Frag2overlaps$str == "+"]

#for +ve, if contact higher than end, call +110%
Loop2Frag2overlaps$lncLocContact_perc[Loop2Frag2overlaps$str == "+" & Loop2Frag2overlaps$start > Loop2Frag2overlaps$Tx_stop] <- 110
#for -ve, if contact is lower than start, then 110%
Loop2Frag2overlaps$lncLocContact_perc[Loop2Frag2overlaps$str == "-" & Loop2Frag2overlaps$start < Loop2Frag2overlaps$Tx_start] <- 110

ggplot(Loop2Frag2overlaps) + aes(x = lncLocContact_perc) +
  geom_histogram()


#### HiC find lnc-PCG connected pairs ####

#of lnc-overlapping loops, isolate loops where, non-lncRNA overlap fragment overlaps an expressed PCG
allGB_PCGs <- filter(allGB, EnsID %in% filter(fpkm_allG, EnsType == "protein_coding", grepl("TF|CC|coding", GeneClassUpdate))$EnsID)

allGB_PCGs_GR <- makeGRangesFromDataFrame(allGB_PCGs[,-7], 
                                          start.field = "Tx_start", 
                                          end.field = "Tx_stop", 
                                          seqnames.field = "chr", 
                                          strand.field = "str", keep.extra.columns = T)

#extend on TSS by 2000kbp
allGB_PCGs_GRpromoters <- promoters(allGB_PCGs_GR, upstream = 2000)

#take union
allGB_PCGs_GR <- GenomicRanges::punion(allGB_PCGs_GR, allGB_PCGs_GRpromoters)
allGB_PCGs_GR$EnsID <- allGB_PCGs_GRpromoters$EnsID
allGB_PCGs_GR$MSTRG_Tx_ID <- allGB_PCGs_GRpromoters$MSTRG_Tx_ID

#isolate lncRNA-contacted loops
FitHiC_LNCS_loops <- unique(c(Loop1Frag1overlaps$Loop1Frag1, Loop1Frag2overlaps$Loop1Frag2))
caSMC_Zhao_Loops1_filt_lncs <- filter(caSMC_Zhao_Loops1_filt, locID_hg19 %in% FitHiC_LNCS_loops)

#for lncs contacted by Frag1, need to check Frag2 for expressed PCGs
Loop1Frag1Lncs_Frag2 <- filter(caSMC_Zhao_Loops1_filt_lncs, locID_hg19 %in% Loop1Frag1overlaps$Loop1Frag1)[,c(1:2,8)]
Loop1Frag1Lncs_Frag2 <- makeGRangesFromDataFrame(Loop1Frag1Lncs_Frag2, 
                                                 start.field = "fragmentMid2", 
                                                 end.field = "fragmentMid2", 
                                                 seqnames.field = "chr2", keep.extra.columns = T)
Loop1Frag1Lncs_Frag2 <- flank(Loop1Frag1Lncs_Frag2, width = 2500, both = T)


#can now overlap loops with lncRNA at Frag1 with all PCGs that overlap the other fragment:
Loop1Frag1Lncs_Frag2index <- findOverlaps(query = Loop1Frag1Lncs_Frag2, subject = allGB_PCGs_GR)
Loop1Frag1Lncs_Frag2overlaps <- unique(data.frame("Loop1Frag1Lncs_Frag2" = Loop1Frag1Lncs_Frag2$locID_hg19[queryHits(Loop1Frag1Lncs_Frag2index)],
                                                  "PCG" = allGB_PCGs_GR$MSTRG_Tx_ID[subjectHits(Loop1Frag1Lncs_Frag2index)]))
Loop1Frag1Lncs_Frag2overlaps <- merge(Loop1Frag1Lncs_Frag2overlaps, unique(fpkm_allG[,c(2,47)]), by.x = "PCG", by.y = "MSTRG_Tx_ID")
head(Loop1Frag1Lncs_Frag2overlaps)
Loop1Frag1Lncs_Frag2[Loop1Frag1Lncs_Frag2$locID_hg19 == "chr15_40987500_41182500"]
allGB_PCGs_GR[allGB_PCGs_GR$MSTRG_Tx_ID == "ENST00000220509.9"]

length(unique(Loop1Frag1Lncs_Frag2overlaps$Loop1Frag1Lncs_Frag2))#366 loops that contact a lnc via frag1 have a PCG at the other end
length(unique(Loop1Frag1overlaps$Loop1Frag1))#there were 992 of these
366/992 #37%
length(unique(Loop1Frag1Lncs_Frag2overlaps$EnsID))#125 PCGs at other end


#for lncs contacted by Frag2, need to check Frag1 for expressed PCGs
Loop1Frag2Lncs_Frag1 <- filter(caSMC_Zhao_Loops1_filt_lncs, locID_hg19 %in% Loop1Frag2overlaps$Loop1Frag2)[,c(3:4,8)]
Loop1Frag2Lncs_Frag1 <- makeGRangesFromDataFrame(Loop1Frag2Lncs_Frag1, 
                                                 start.field = "fragmentMid1", 
                                                 end.field = "fragmentMid1", 
                                                 seqnames.field = "chr1", keep.extra.columns = T)
Loop1Frag2Lncs_Frag1 <- flank(Loop1Frag2Lncs_Frag1, width = 2500, both = T)

#can now overlap loops with lncRNA at Frag2 with all PCGs that overlap the other fragment:
Loop1Frag2Lncs_Frag1index <- findOverlaps(query = Loop1Frag2Lncs_Frag1, subject = allGB_PCGs_GR)
Loop1Frag2Lncs_Frag1overlaps <- unique(data.frame("Loop1Frag2Lncs_Frag1" = Loop1Frag2Lncs_Frag1$locID_hg19[queryHits(Loop1Frag2Lncs_Frag1index)],
                                                  "PCG" = allGB_PCGs_GR$MSTRG_Tx_ID[subjectHits(Loop1Frag2Lncs_Frag1index)]))
Loop1Frag2Lncs_Frag1overlaps <- merge(Loop1Frag2Lncs_Frag1overlaps, unique(fpkm_allG[,c(2,47)]), by.x = "PCG", by.y = "MSTRG_Tx_ID")
head(Loop1Frag2Lncs_Frag1overlaps)
Loop1Frag2Lncs_Frag1[Loop1Frag2Lncs_Frag1$locID_hg19 == "chr16_3112500_3177500"]
allGB_PCGs_GR[allGB_PCGs_GR$MSTRG_Tx_ID == "ENST00000008180.13"]

length(unique(Loop1Frag2Lncs_Frag1overlaps$Loop1Frag2Lncs_Frag1))#363 loops that contact a lnc via frag2 have a PCG at the other end
length(unique(Loop1Frag2overlaps$Loop1Frag2))#there were 1124 of these
363/1124 #32% - seems a higher rate of PCG matching for lncs with a f2 contact...
length(unique(Loop1Frag2Lncs_Frag1overlaps$EnsID))#123 PCGs at other end

#for FitHiC, define lncRNA-PCG contact loci
#this table has loop which contact a lnc at frag1:
head(Loop1Frag1overlaps)
#this table has which of these loops contact a PCG at frag2:
head(Loop1Frag1Lncs_Frag2overlaps)
#sanity check, all latter found in former:
dim(Loop1Frag1Lncs_Frag2overlaps)
sum(Loop1Frag1Lncs_Frag2overlaps$Loop1Frag1Lncs_Frag2 %in% Loop1Frag1overlaps$Loop1Frag1)

#merge on the 2 columns containing the hg19-defined locus ID
Loop1Frag1_LNC_PCG <- merge(Loop1Frag1overlaps, Loop1Frag1Lncs_Frag2overlaps, by.x = "Loop1Frag1", by.y = "Loop1Frag1Lncs_Frag2", 
                            all.x = T)

#equivalent for Frag2 lnc loops:
Loop1Frag2_LNC_PCG <- merge(Loop1Frag2overlaps, Loop1Frag2Lncs_Frag1overlaps, by.x = "Loop1Frag2", by.y = "Loop1Frag2Lncs_Frag1", 
                            all.x = T)


#HiCCUPs
HiCCUPs_LNCS_loops <- unique(c(Loop2Frag1overlaps$Loop2Frag1, Loop2Frag2overlaps$Loop2Frag2))

caSMC_Zhao_Loops2_filt_lncs <- filter(caSMC_Zhao_Loops2_filt, locID_hg19 %in% HiCCUPs_LNCS_loops)

#for lncs contacted by Frag1, need to check Frag2 for expressed PCGs
Loop2Frag1Lncs_Frag2 <- filter(caSMC_Zhao_Loops2_filt_lncs, locID_hg19 %in% Loop2Frag1overlaps$Loop2Frag1)[,c(1:2,21)]
Loop2Frag1Lncs_Frag2 <- makeGRangesFromDataFrame(Loop2Frag1Lncs_Frag2, 
                                                 start.field = "fragmentMid2", 
                                                 end.field = "fragmentMid2", 
                                                 seqnames.field = "chr2", keep.extra.columns = T)
Loop2Frag1Lncs_Frag2 <- flank(Loop2Frag1Lncs_Frag2, width = 2500, both = T)


#can now overlap loops with lncRNA at Frag1 with all PCGs that overlap the other fragment:
Loop2Frag1Lncs_Frag2index <- findOverlaps(query = Loop2Frag1Lncs_Frag2, subject = allGB_PCGs_GR)
Loop2Frag1Lncs_Frag2overlaps <- unique(data.frame("Loop2Frag1Lncs_Frag2" = Loop2Frag1Lncs_Frag2$locID_hg19[queryHits(Loop2Frag1Lncs_Frag2index)],
                                                  "PCG" = allGB_PCGs_GR$MSTRG_Tx_ID[subjectHits(Loop2Frag1Lncs_Frag2index)]))
Loop2Frag1Lncs_Frag2overlaps <- merge(Loop2Frag1Lncs_Frag2overlaps, unique(fpkm_allG[,c(2,47)]), by.x = "PCG", by.y = "MSTRG_Tx_ID")
head(Loop2Frag1Lncs_Frag2overlaps)
Loop2Frag1Lncs_Frag2[Loop2Frag1Lncs_Frag2$locID_hg19 == "chr4_129495833_129737500"]
allGB_PCGs_GR[allGB_PCGs_GR$MSTRG_Tx_ID == "ENST00000226319.10"]

length(unique(Loop2Frag1Lncs_Frag2overlaps$Loop2Frag1Lncs_Frag2))#44 loops that contact a lnc via frag1 have a PCG at the other end
length(unique(Loop2Frag1overlaps$Loop2Frag1))#there were 141 of these
44/141 #31%
length(unique(Loop2Frag1Lncs_Frag2overlaps$EnsID))#44 PCGs at other end


#for lncs contacted by Frag2, need to check Frag1 for expressed PCGs
Loop2Frag2Lncs_Frag1 <- filter(caSMC_Zhao_Loops2_filt_lncs, locID_hg19 %in% Loop2Frag2overlaps$Loop2Frag2)[,c(3:4,21)]
Loop2Frag2Lncs_Frag1 <- makeGRangesFromDataFrame(Loop2Frag2Lncs_Frag1, 
                                                 start.field = "fragmentMid1", 
                                                 end.field = "fragmentMid1", 
                                                 seqnames.field = "chr1", keep.extra.columns = T)
Loop2Frag2Lncs_Frag1 <- flank(Loop2Frag2Lncs_Frag1, width = 2500, both = T)

#can now overlap loops with lncRNA at Frag2 with all PCGs that overlap the other fragment:
Loop2Frag2Lncs_Frag1index <- findOverlaps(query = Loop2Frag2Lncs_Frag1, subject = allGB_PCGs_GR)
Loop2Frag2Lncs_Frag1overlaps <- unique(data.frame("Loop2Frag2Lncs_Frag1" = Loop2Frag2Lncs_Frag1$locID_hg19[queryHits(Loop2Frag2Lncs_Frag1index)],
                                                  "PCG" = allGB_PCGs_GR$MSTRG_Tx_ID[subjectHits(Loop2Frag2Lncs_Frag1index)]))
Loop2Frag2Lncs_Frag1overlaps <- merge(Loop2Frag2Lncs_Frag1overlaps, unique(fpkm_allG[,c(2,47)]), by.x = "PCG", by.y = "MSTRG_Tx_ID")
head(Loop2Frag2Lncs_Frag1overlaps)
Loop2Frag2Lncs_Frag1[Loop2Frag2Lncs_Frag1$locID_hg19 == "chr6_52375000_52525000"]
allGB_PCGs_GR[allGB_PCGs_GR$MSTRG_Tx_ID == "ENST00000182527.3"]

length(unique(Loop2Frag2Lncs_Frag1overlaps$Loop2Frag2Lncs_Frag1))#42 loops that contact a lnc via frag2 have a PCG at the other end
length(unique(Loop2Frag2overlaps$Loop2Frag2))#there were 123 of these
42/123 #34.1% - again seems a higher rate of PCG matching for lncs with a f2 contact...
length(unique(Loop2Frag2Lncs_Frag1overlaps$EnsID))#40 PCGs at other end

#for HiCCUPs, define lncRNA-PCG contact loci
#this table has loop which contact a lnc at frag1:
head(Loop2Frag1overlaps)
#this table has which of these loops contact a PCG at frag2:
head(Loop2Frag1Lncs_Frag2overlaps)
#sanity check, all latter found in former:
dim(Loop2Frag1Lncs_Frag2overlaps)
sum(Loop2Frag1Lncs_Frag2overlaps$Loop2Frag1Lncs_Frag2 %in% Loop2Frag1overlaps$Loop2Frag1)

#merge on the 2 columns containing the hg19 locus ID
Loop2Frag1_LNC_PCG <- merge(Loop2Frag1overlaps, Loop2Frag1Lncs_Frag2overlaps, by.x = "Loop2Frag1", by.y = "Loop2Frag1Lncs_Frag2", 
                            all.x = T)

#equivalent for Frag2 lnc loops:
Loop2Frag2_LNC_PCG <- merge(Loop2Frag2overlaps, Loop2Frag2Lncs_Frag1overlaps, by.x = "Loop2Frag2", by.y = "Loop2Frag2Lncs_Frag1", 
                            all.x = T)


#same/unique pairs FitHiC/HiCCUPs:
Loop1Frag1_LNC_PCG$pairs <- paste(Loop1Frag1_LNC_PCG$EnsID.x, Loop1Frag1_LNC_PCG$EnsID.y, sep = "-")
Loop1Frag2_LNC_PCG$pairs <- paste(Loop1Frag2_LNC_PCG$EnsID.x, Loop1Frag2_LNC_PCG$EnsID.y, sep = "-")
Loop2Frag1_LNC_PCG$pairs <- paste(Loop2Frag1_LNC_PCG$EnsID.x, Loop2Frag1_LNC_PCG$EnsID.y, sep = "-")
Loop2Frag2_LNC_PCG$pairs <- paste(Loop2Frag2_LNC_PCG$EnsID.x, Loop2Frag2_LNC_PCG$EnsID.y, sep = "-")

Loop1Frag1_LNC_PCG$pairs[is.na(Loop1Frag1_LNC_PCG$EnsID.y)] <- NA
Loop1Frag2_LNC_PCG$pairs[is.na(Loop1Frag2_LNC_PCG$EnsID.y)] <- NA
Loop2Frag1_LNC_PCG$pairs[is.na(Loop2Frag1_LNC_PCG$EnsID.y)] <- NA
Loop2Frag2_LNC_PCG$pairs[is.na(Loop2Frag2_LNC_PCG$EnsID.y)] <- NA

FitHiC_LNC_PCG_pairs <- unique(c(Loop1Frag1_LNC_PCG$pairs, Loop1Frag2_LNC_PCG$pairs))
HiCCUP_LNC_PCG_pairs <- unique(c(Loop2Frag1_LNC_PCG$pairs, Loop2Frag2_LNC_PCG$pairs))

#overlap of 41
sum(HiCCUP_LNC_PCG_pairs %in% FitHiC_LNC_PCG_pairs)
69/88 #78% of HiCCUPs found in HiC
69/263 #26% of FitHiC

#overlap is higher p FitHiC?
FHChi1 <- filter(Loop1Frag1_LNC_PCG, q.value < 0.0005, !is.na(pairs))$pairs
FHChi2 <- filter(Loop1Frag2_LNC_PCG, q.value < 0.0005, !is.na(pairs))$pairs

FitHiC_LNC_PCG_pairshi <- unique(c(FHChi1, FHChi2))

sum(HiCCUP_LNC_PCG_pairs %in% FitHiC_LNC_PCG_pairshi)
53/88 #60%
53/263 #20%

#lots of loops formed outside of the 2D limits:
HiCCUP_LNC_PCG_pairs[!HiCCUP_LNC_PCG_pairs %in% AllLNC_AllPCG_1$pairs]
FitHiC_LNC_PCG_pairs[!FitHiC_LNC_PCG_pairs %in% AllLNC_AllPCG_1$pairs]
FitHiC_LNC_PCG_pairshi[!FitHiC_LNC_PCG_pairshi %in% AllLNC_AllPCG_1$pairs]


#### Combine all pairs and annotate ####

#2D neighbours
#all 2d pairs:
AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_info_2026.csv", header = T)
#AllLNC_AllPCG_1 <- filter(AllLNC_AllPCG_1, DisLnc_PCG <600)
colnames(AllLNC_AllPCG_1)
#3D neighbours
#combine FitHiC, HiCCUPs
colnames(Loop1Frag1_LNC_PCG)[1] <- "loopID"
Loop1Frag1_LNC_PCG$loopFragLnc <- "Frag1"
colnames(Loop1Frag2_LNC_PCG)[1] <- "loopID"
Loop1Frag2_LNC_PCG$loopFragLnc <- "Frag2"
Loop1_LNC_PCG <- rbind(Loop1Frag1_LNC_PCG, Loop1Frag2_LNC_PCG)

colnames(Loop2Frag1_LNC_PCG)[1] <- "loopID"
Loop2Frag1_LNC_PCG$loopFragLnc <- "Frag1"
colnames(Loop2Frag2_LNC_PCG)[1] <- "loopID"
Loop2Frag2_LNC_PCG$loopFragLnc <- "Frag2"
Loop2_LNC_PCG <- rbind(Loop2Frag1_LNC_PCG, Loop2Frag2_LNC_PCG)
#take FDR single value from HiCCUPs:
Loop2_LNC_PCG$worstFDR <- Biobase::rowMax(as.matrix(Loop2_LNC_PCG[,c(9:12)]))

#combine FitHiC, HiCCUPs, just take minimum info
colnames(Loop1_LNC_PCG)
colnames(Loop2_LNC_PCG)

trial <- Loop1_LNC_PCG[,c(3,15,16)]
triali <- Loop2_LNC_PCG[,c(3,19,20)]

Loop_LNC_PCG <- unique(rbind(trial, triali))
#not considering non-connected lncs (for now)
Loop_LNC_PCG <- filter(Loop_LNC_PCG, !is.na(pairs))

Loop_LNC_PCG$loopMethod <- "FitHiC"
Loop_LNC_PCG$loopMethod[Loop_LNC_PCG$pairs %in% Loop2_LNC_PCG$pairs] <- "HiCCUPs"
Loop_LNC_PCG$loopMethod[Loop_LNC_PCG$pairs %in% Loop2_LNC_PCG$pairs & Loop_LNC_PCG$pairs %in% Loop1_LNC_PCG$pairs] <- "Both"
table(Loop_LNC_PCG$loopMethod)
length(unique(Loop_LNC_PCG$pairs))#281 gene pairs

#key info, IDs, loopID
trial <- AllLNC_AllPCG_1[,c(2,1,7)]
trial$loopMethod <- "Neither"

triali <- Loop_LNC_PCG
colnames(triali) <- colnames(trial)

#remove 2D pairs covered in 3D already:
trial <- trial[!trial$pairs %in% triali$pairs,]

AllLNC_AllPCG_2d3d <- unique(rbind(trial, triali))
#8453 pairs, mostly unlooped
table(AllLNC_AllPCG_2d3d$loopMethod)

#additional info, name, up/down, timeframe, biotype, expression per timeframe etc:
colnames(fpkm_allG)
trial <- unique(merge(AllLNC_AllPCG_2d3d, fpkm_allG[,c(2,3,58)], by = "EnsID"))
trial <- unique(merge(trial, fpkm_allG[,c(2,3,58)], by.x = "EnsID.y", by.y = "EnsID"))
colnames(fpkm_allGDE)
trial <- unique(merge(trial, fpkm_allGDE[,c(1,46)], by = "EnsID",  all.x = T))
trial <- unique(merge(trial, fpkm_allGDE[,c(1,46)], by.x = "EnsID.y", by.y = "EnsID",  all.x = T))

colnames(trial)[9:10] <- c("Lnc_Cluster", "PCG_Cluster")

trial$Lnc_Timeframe <- gsub(".* ", "", trial$Lnc_Cluster)
trial$PCG_Timeframe <- gsub(".* ", "", trial$PCG_Cluster)

#plot distance between
#add in TSS and distance between, allTSS:
allGB$TSS_FANTOM_GENCODE <- allGB$Tx_start
allGB$TSS_FANTOM_GENCODE[allGB$str == "-"] <- allGB$Tx_stop[allGB$str == "-"]

#n.b. multiple TSS per gene
triali <- unique(merge(trial, allGB[,c(2,8)], by.x = "EnsID", by.y = "EnsID"))
triali <- unique(merge(triali, allGB[,c(2,8)], by.x = "EnsID.y", by.y = "EnsID"))
triali$AbsDistLnc_PCG <- abs(triali$TSS_FANTOM_GENCODE.x - triali$TSS_FANTOM_GENCODE.y)/1000
triali$DistLnc_PCG <- (triali$TSS_FANTOM_GENCODE.x - triali$TSS_FANTOM_GENCODE.y)/1000

#shortest distance per pair:
trialii <- split(triali, triali$pairs)

trialii[[27]]
trialiii <- lapply(trialii, function(x){
  #remove the TSS cols
  unique(x[x$AbsDistLnc_PCG == min((x$AbsDistLnc_PCG)),-c(13:14)])
})

trial <- bind_rows(trialiii)

#lnc-PCG relationship, need strand
colnames(fpkm_allG)
trial <- unique(merge(trial, fpkm_allG[,c(2,8)], by = "EnsID"))
trial <- unique(merge(trial, fpkm_allG[,c(2,8)], by.x = "EnsID.y", by.y = "EnsID"))

#need to import overlaps
FPKM_CQV_OVERLAP_fpkm <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SameoverlapsG_lncs.csv")
#need the ensID col
triali <- unique(merge(FPKM_CQV_OVERLAP_fpkm, fpkm_allG[,1:2], by.x = "MSTRG_ID.y", by.y = "MSTRG_ID"))
triali$pairs <- paste(triali$EnsID.x, triali$EnsID, sep = "-")

trial$LNC_PCG_Type <- NA
#same strand overlaps
trial$LNC_PCG_Type[trial$pairs %in% triali$pairs & trial$str.x == trial$str.y] <- "Sense Overlap"
#opposite strand overlaps
trial$LNC_PCG_Type[trial$pairs %in% triali$pairs & !trial$str.x == trial$str.y] <- "Antisense Overlap"
#intergenic overlaps
trial$LNC_PCG_Type[!trial$pairs %in% triali$pairs] <- "Intergenic"
#where the lncRNA is enhancer annotated
#trial$LNC_PCG_Type[!trial$pairs %in% triali$pairs & trial$EnsID %in% filter(fpkm_allGDE,  grepl("ELnc", GeneClassUpdate))$EnsID] <- "Intergenic Enhancer"
#no overlap but the transcription appears divergent:
trial$LNC_PCG_Type[!trial$pairs %in% triali$pairs & !trial$str.x == trial$str.y & trial$AbsDistLnc_PCG <3] <- "Bidirectional"
table(trial$LNC_PCG_Type)#187 AS, 150 bidir, 5098 intergenic, 1 sense overlap

colnames(trial)[17] <- "LncRNA-PCG Relationship"

AllLNC_AllPCG_2d3d <- trial

length(unique(AllLNC_AllPCG_2d3d$EnsID)) #579
length(unique(AllLNC_AllPCG_2d3d$EnsID.y)) #5362

#write.csv(AllLNC_AllPCG_2d3d, "AllLNC_AllPCG_2d3d_2026.csv", row.names = F)


#### HiC validation 2 - set-up: pick out CClncRNA targets vs. other neighbours, compare to DEL-DEP ####

#all 2d pairs in 1mbp, plus any HiC up to 1Mbp
AllLNC_AllPCG_2d3d <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_2d3d_2026.csv")
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
281/8454 #1Mbp
281/5438 #previous 600kbp object

#for pairs formed in 2d, what % are HiC linked amongst CisLnc neighbours, CisLnc targets 
#(3d are obviously 100% HiC so excluded)
#similarly pairs with both loci in 15kbp cannot be considered - all 0% HiC:
#examine loop distances in FitHiC/HiCCUps:
#summary(caSMC_Zhao_Loops1$fragmentMid2 - caSMC_Zhao_Loops1$fragmentMid1)
#summary(caSMC_Zhao_Loops2$centroid2 - caSMC_Zhao_Loops2$centroid1)
#lowest possible distance is 15kbp confirmed

#GR for lncs
allGB_LNCS <- filter(allGB, EnsID %in% filter(fpkm_allG, grepl("Lnc|fide", GeneClassUpdate))$EnsID)
length(unique(allGB_LNCS$EnsID))#597 as expected

allGB_LNCS_GR <- makeGRangesFromDataFrame(allGB_LNCS[,-7], 
                                          start.field = "Tx_start", 
                                          end.field = "Tx_stop", 
                                          seqnames.field = "chr", 
                                          strand.field = "str", keep.extra.columns = T)

#extend on TSS by 2000kbp
allGB_LNCS_GRpromoters <- promoters(allGB_LNCS_GR, upstream = 2000)

#take union
allGB_LNCS_GR <- GenomicRanges::punion(allGB_LNCS_GR, allGB_LNCS_GRpromoters)
allGB_LNCS_GR$EnsID <- allGB_LNCS_GRpromoters$EnsID
allGB_LNCS_GR$MSTRG_Tx_ID <- allGB_LNCS_GRpromoters$MSTRG_Tx_ID

#GR for PCGs
allGB_PCGs <- filter(allGB, EnsID %in% filter(fpkm_allG, EnsType == "protein_coding", grepl("TF|CC|coding", GeneClassUpdate))$EnsID)

allGB_PCGs_GR <- makeGRangesFromDataFrame(allGB_PCGs[,-7], 
                                          start.field = "Tx_start", 
                                          end.field = "Tx_stop", 
                                          seqnames.field = "chr", 
                                          strand.field = "str", keep.extra.columns = T)

#extend on TSS by 2000kbp
allGB_PCGs_GRpromoters <- promoters(allGB_PCGs_GR, upstream = 2000)

#take union
allGB_PCGs_GR <- GenomicRanges::punion(allGB_PCGs_GR, allGB_PCGs_GRpromoters)
allGB_PCGs_GR$EnsID <- allGB_PCGs_GRpromoters$EnsID
allGB_PCGs_GR$MSTRG_Tx_ID <- allGB_PCGs_GRpromoters$MSTRG_Tx_ID

#import 2d pairs:
AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_info_2026.csv", header = T)

#remove any which are too close for HiC (expecting small minority)
#flattens all tx (+ promoter regions) in all pairs into a single range
lnc_pcg_ranges <- list()
for (i in 1:length(AllLNC_AllPCG_1$pairs)){
  lnc_pcg_ranges[[i]] <- range(allGB_LNCS_GR[allGB_LNCS_GR$EnsID %in% AllLNC_AllPCG_1[i,2]], 
                               allGB_PCGs_GR[allGB_PCGs_GR$EnsID == AllLNC_AllPCG_1[i,1]], ignore.strand = T)
}

length(unique(AllLNC_AllPCG_1$EnsID))
AllLNC_AllPCG_1$pair_range <- sapply(lnc_pcg_ranges, width)
AllLNC_AllPCG_1$pair_range_bin <- NA
AllLNC_AllPCG_1$pair_range_bin[AllLNC_AllPCG_1$pair_range > 500000] <- ">500kbp"
AllLNC_AllPCG_1$pair_range_bin[AllLNC_AllPCG_1$pair_range < 500000] <- "400-500kbp"
AllLNC_AllPCG_1$pair_range_bin[AllLNC_AllPCG_1$pair_range < 400000] <- "300-400kbp"
AllLNC_AllPCG_1$pair_range_bin[AllLNC_AllPCG_1$pair_range < 300000] <- "200-300kbp"
AllLNC_AllPCG_1$pair_range_bin[AllLNC_AllPCG_1$pair_range < 200000] <- "100-200kbp"
AllLNC_AllPCG_1$pair_range_bin[AllLNC_AllPCG_1$pair_range < 100000] <- "50-100kbp"
AllLNC_AllPCG_1$pair_range_bin[AllLNC_AllPCG_1$pair_range < 50000] <- "25-50kbp"
AllLNC_AllPCG_1$pair_range_bin[AllLNC_AllPCG_1$pair_range < 25000] <- "15-25kbp"
AllLNC_AllPCG_1$pair_range_bin[AllLNC_AllPCG_1$pair_range < 15000] <- "<15kbp"

table(AllLNC_AllPCG_1$pair_range_bin)
AllLNC_AllPCG_1_ <- filter(AllLNC_AllPCG_1, AbsDistLnc_PCG <1000) #can come back later and this line to do more filtering if needed

AllLNC_AllPCG_HiCTest <- filter(AllLNC_AllPCG_1_, pair_range >=15000)

AllLNC_AllPCG_HiCTest_naive <- filter(AllLNC_AllPCG_HiCTest, pairs %in% AllLNC_AllPCG_1_$pairs)
#total pairs
dim(filter(AllLNC_AllPCG_HiCTest_naive))
#with HiC
dim(filter(AllLNC_AllPCG_HiCTest_naive, pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs))
266/8400 #1mbp

#with closer pairs, a higher % of HiC loop support is found

AllLNC_AllPCG_HiCTest <- filter(AllLNC_AllPCG_1_, pair_range >=15000)

#2d/naive version - for paper - take only the pairs found in the set distance
#(should be no diff from this line)
AllLNC_AllPCG_HiCTest_naive <- filter(AllLNC_AllPCG_HiCTest, pairs %in% AllLNC_AllPCG_1_$pairs)


#
#### HiC validation 2 - all timeframes together ####

#considering all lnc-PCG neighours within 800kbp, are cclnc pairs enriched with hic?

#total coreg pairs (in the test)
CoRegPairs_04_48_24_extended_naive <- filter(AllLNC_AllPCG_1_,
                                             #AllLNC_AllPCG_1,
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                           fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                               (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                             fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                              fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                               (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                             fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
)
#1000kbp 1001
#800kbp 821
#600kbp 633 pairs
#500kbp 525 pairs
#250kbp 282 pairs

#enrichment irrespective of time
#for a given cclnc, is hic support found more readily amongst it's co-reg neighbours over it's other neighbours

#total pairs in the test that are co-regulated
b <- dim(filter(AllLNC_AllPCG_HiCTest_naive, pairs %in% CoRegPairs_04_48_24_extended_naive$pairs))[1]
#and those with HiC
a <- dim(filter(AllLNC_AllPCG_HiCTest_naive, pairs %in% CoRegPairs_04_48_24_extended_naive$pairs,
                pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs))[1]
44/996#1000

#enrichment of HiC connections amongst lnc co-reg pairs?

#all cclnc neighbour pairs:
d <- dim(filter(AllLNC_AllPCG_HiCTest_naive, EnsID %in% CoRegPairs_04_48_24_extended_naive$EnsID))[1]
#with HiC
c <- dim(filter(AllLNC_AllPCG_HiCTest_naive, EnsID %in% CoRegPairs_04_48_24_extended_naive$EnsID,
                pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs))[1]
98/2844#1000
fisher.test(data.frame("cisLnc" = c(a, b-a),
                       "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")


#useful at this stage, to have some control pairs
#enrichment of HiC connections amongst nonDE-lncs to any DE PCG
#all nonDE lncs with a DE PCG
bi <- dim(filter(AllLNC_AllPCG_HiCTest_naive, 
                 !EnsID %in% fpkm_allGDE$EnsID, EnsID.y %in% fpkm_allGDE$EnsID))[1]
#with HiC support
ai <- dim(filter(AllLNC_AllPCG_HiCTest_naive, 
                 !EnsID %in% fpkm_allGDE$EnsID, EnsID.y %in% fpkm_allGDE$EnsID,
                 pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs))[1]
65/2110#1000

di <- dim(filter(AllLNC_AllPCG_HiCTest_naive, !EnsID %in% fpkm_allGDE$EnsID))[1] #1380 non-DE lnc-PCG pairs
#with HiC
ci <- dim(filter(AllLNC_AllPCG_HiCTest_naive, !EnsID %in% fpkm_allGDE$EnsID,
                 pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs))[1]
151/5371#1000

fisher.test(data.frame("cisLnc" = c(ai, bi-ai),
                       "Not"   = c(ci-ai, di-ci-(bi-ai))), alternative = "greater")

#there is no observable tendency for closely found lncs/DE PCGs to have HiC contact regardless of whether the lnc is DE

#same timeframe equivalent:
CoRegPairs_04_48_24_extended_naiveSame <- filter(AllLNC_AllPCG_1,
                                                 (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                               fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                                fpkm_allGDE_Downwithin_4$EnsID)) |
                                                   (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                 fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID)) |
                                                   (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                 fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
                                                 )

#enrichment irrespective of time
bii <- dim(filter(AllLNC_AllPCG_HiCTest_naive, pairs %in% CoRegPairs_04_48_24_extended_naiveSame$pairs))[1]
#with HiC
aii <- dim(filter(AllLNC_AllPCG_HiCTest_naive, pairs %in% CoRegPairs_04_48_24_extended_naiveSame$pairs,
                  pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs))[1]
28/466#1000

#all cclnc neighbour pairs:
dii <- dim(filter(AllLNC_AllPCG_HiCTest_naive, EnsID %in% CoRegPairs_04_48_24_extended_naiveSame$EnsID))[1]
##cclnc target pairs
cii <- dim(filter(AllLNC_AllPCG_HiCTest_naive, EnsID %in% CoRegPairs_04_48_24_extended_naiveSame$EnsID,
                  pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs))[1]
86/2689#1000

fisher.test(data.frame("cisLnc" = c(aii, bii-aii),
                       "Not"   = c(cii-aii, dii-cii-(bii-aii))), alternative = "greater")


#### HiC validation 2 - early timeframe plots ####

#per timeframe, per induced/repressed analysis
DEG_cluster <- list(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Upwithin_8, fpkm_allGDE_Upwithin_24,
                    fpkm_allGDE_Downwithin_4, fpkm_allGDE_Downwithin_8, fpkm_allGDE_Downwithin_24,
                    rbind(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Downwithin_4),
                    rbind(fpkm_allGDE_Upwithin_8, fpkm_allGDE_Downwithin_8),
                    rbind(fpkm_allGDE_Upwithin_24, fpkm_allGDE_Downwithin_24)
                    )

#total coreg pairs (in the test)
CoRegPairs_04_48_24_extended_naive <- filter(AllLNC_AllPCG_1_,
                                             #AllLNC_AllPCG_1,
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                           fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                               (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                             fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                              fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                               (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                             fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
)

#same + delayed timeframe pairings:
LoopFish_SameDelayedPairs <- list()

for (i in 1:length(DEG_cluster)){
  #for co-reg lncs per cluster, all neighbour pairs capable of forming a HiC pair
  CCLNC_AllPCG_HiCTest <- filter(AllLNC_AllPCG_HiCTest, 
                                 #all lncs with a DE PCG near
                                 EnsID %in% CoRegPairs_04_48_24_extended_naive$EnsID,
                                 #from a given cluster
                                 EnsID %in% DEG_cluster[[i]]$EnsID)
  #of these, how many are co-reg, same/later timeframe
  CCLNC_targetPCG_HiCTest <- filter(CCLNC_AllPCG_HiCTest, 
                                    pairs %in% CoRegPairs_04_48_24_extended_naive$pairs)
  
  #per selection and background, identify the pairs which are HiC looped, count total pairs too
  a <- length(unique(filter(CCLNC_targetPCG_HiCTest, pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs)$pairs))
  b <- length(unique(CCLNC_targetPCG_HiCTest$pairs))
  c <- length(unique(filter(CCLNC_AllPCG_HiCTest, pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs)$pairs))
  d <- length(unique(CCLNC_AllPCG_HiCTest$pairs))
  
  LoopFish_SameDelayedPairs[[i]] <- data.frame(a,b,c,d,
                                               a/b,
                                               c/d,
                                               fisher.test(data.frame("cisLnc" = c(a, b-a),
                                                                      "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p,
                                               fisher.test(data.frame("cisLnc" = c(a, b-a),
                                                                      "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate)
  colnames(LoopFish_SameDelayedPairs[[i]]) <- c("a", "b", "c", "d", "SelectPercLooped", "BackPerLooped", "p", "OR")
}

names(LoopFish_SameDelayedPairs) <- c("Up4", "Up8", "Up24", "Down4", "Down8", "Down24", "Either4", "Either8", "Either24")

LoopFish_SameDelayedPairs_df_1000 <- bind_rows(LoopFish_SameDelayedPairs, .id = "LncTiming")


#same timeframe equivalent:
CoRegPairs_04_48_24_extended_naiveSame <- filter(AllLNC_AllPCG_1_,
                                                 (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                               fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                                fpkm_allGDE_Downwithin_4$EnsID)) |
                                                   (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                 fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID)) |
                                                   (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                 fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
)

LoopFish_SamePairs <- list()

for (i in 1:length(DEG_cluster)){
  #for co-reg lncs per cluster, all neighbour pairs capable of forming a HiC pair
  CCLNC_AllPCG_HiCTest <- filter(AllLNC_AllPCG_HiCTest, 
                                 EnsID %in% CoRegPairs_04_48_24_extended_naiveSame$EnsID, 
                                 EnsID %in% DEG_cluster[[i]]$EnsID)
  #just the CClncRNA-target pairs now
  CCLNC_targetPCG_HiCTest <- filter(CCLNC_AllPCG_HiCTest, 
                                    pairs %in% CoRegPairs_04_48_24_extended_naiveSame$pairs)
  
  a <- length(unique(filter(CCLNC_targetPCG_HiCTest, pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs)$pairs))
  b <- length(unique(CCLNC_targetPCG_HiCTest$pairs))
  c <- length(unique(filter(CCLNC_AllPCG_HiCTest, pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs)$pairs))
  d <- length(unique(CCLNC_AllPCG_HiCTest$pairs))
  
  LoopFish_SamePairs[[i]] <- data.frame(a,b,c,d,
                                        a/b,
                                        c/d,
                                        fisher.test(data.frame("cisLnc" = c(a, b-a),
                                                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p,
                                        fisher.test(data.frame("cisLnc" = c(a, b-a),
                                                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate)
  colnames(LoopFish_SamePairs[[i]]) <- c("a", "b", "c", "d", "SelectPercLooped", "BackPerLooped", "p", "OR")
}

names(LoopFish_SamePairs) <- c("Up4", "Up8", "Up24", "Down4", "Down8", "Down24", "Either4", "Either8", "Either24")

LoopFish_SamePairs_df_1000 <- bind_rows(LoopFish_SamePairs, .id = "LncTiming")


#collect a non-DE lnc comparison for 0-4hr to display alongside (later timepoint can be displayed without e.g. in supplement)
fpkm_allG_04 <- filter(fpkm_allG, (Hour0_meanFPKM>0.8 | Hour4_meanFPKM>0.8), !EnsID %in% fpkm_allGDE$EnsID)
fpkm_allGDE_within_4 <- rbind(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Downwithin_4)

#non-DE lncs pairs
CCLNC_AllPCG_HiCTest <- filter(AllLNC_AllPCG_HiCTest, EnsID %in% fpkm_allG_04$EnsID)

#no. non-DEL neighbour pairs where the neighbour is 0-4hr DEs
CCLNC_targetPCG_HiCTest <- filter(CCLNC_AllPCG_HiCTest, EnsID.y %in% fpkm_allGDE_within_4$EnsID)

a <- length(unique(filter(CCLNC_targetPCG_HiCTest, pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs)$pairs))
b <- length(unique(CCLNC_targetPCG_HiCTest$pairs))
c <- length(unique(filter(CCLNC_AllPCG_HiCTest, pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs)$pairs))
d <- length(unique(CCLNC_AllPCG_HiCTest$pairs))

data.frame(a/b,
           c/d,
           "p" = fisher.test(data.frame("cisLnc" = c(a, b-a),
                                        "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p,
           "or" =fisher.test(data.frame("cisLnc" = c(a, b-a),
                                        "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate)

nonDElnc_hicP <- fisher.test(data.frame("cisLnc" = c(a, b-a),
                                        "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p

#plot for hic support in early, 0-4hr induced/repressed/nonDE
DEL_PCG_type <- rbind("Induced" = LoopFish_SamePairs_df_250[1,c(2,4,6,7)],
                      "Repressed" = LoopFish_SamePairs_df_250[4,c(2,4,6,7)],
                      "Non-DE" = c(a,c,a/b,c/d)
                      #"Regulated" = LoopFish_SamePairs_df[7,c(2,4,6,7)]
                      )

colnames(DEL_PCG_type) <- c("a", "c", "0-4hr DE neighbour", "Any neighbour")
DEL_PCG_type$Type <- rownames(DEL_PCG_type)

library(reshape2)
DEL_PCG_typei <- melt(DEL_PCG_type[,3:5])
DEL_PCG_typei <- cbind(DEL_PCG_typei, melt(DEL_PCG_type[,1:2]))
DEL_PCG_type <- DEL_PCG_typei

DEL_PCG_type$value <- DEL_PCG_type$value*100

colnames(DEL_PCG_type)[c(2,5)] <- c("Paired with", "noPairs")

DEL_PCG_type$Type <- factor(DEL_PCG_type$Type)
DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,3,1)])

DEL_PCG_type$`Paired with` <- factor(DEL_PCG_type$`Paired with`)
DEL_PCG_type$`Paired with` <- factor(DEL_PCG_type$`Paired with`, levels = levels(DEL_PCG_type$`Paired with`)[c(2,1)])

ggplot(DEL_PCG_type) + aes(y = Type, x = value, fill = `Paired with`, label = noPairs) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "grey60") +
  scale_fill_manual(values = c(`0-4hr DE neighbour` = "mediumorchid", `Any neighbour` = "grey60")) +
  geom_text(aes(y = Type, x = value+1.5, fill = `Paired with`, label = noPairs),
            position = position_dodge(width = 0.7), color = "black", size =5.5) +
  ylab("") +
  xlab("% HiC Connected\nPairs") +
  scale_x_continuous(breaks = seq(0,15,5), limits = c(0,18)) +
  theme_minimal() +
  #ggtitle("LncRNAs with a DE PCG neighbour", subtitle = "(same timeframe)") +
  theme(text = element_text(size=24))# + Seurat::RotatedAxis()
p.adjust(c(LoopFish_SamePairs_df_1000$p[c(1,4)], nonDElnc_hicP), method = "bonferroni")
p.adjust(c(LoopFish_SamePairs_df_250$p[c(1,4)], nonDElnc_hicP), method = "bonferroni")


#### Option to build a supp table with relevant tests ####

#250kbp pairs
#1000kbp pairs
#by temporal cluster (omit "either" ?)
#with all fisher's details (abcd, %, OR, p)
#same
#same/delayed

trial <- rbind(LoopFish_SamePairs_df_250[1:6,],
               LoopFish_SamePairs_df_1000[1:6,],
               LoopFish_SameDelayedPairs_df_250[1:6,],
               LoopFish_SameDelayedPairs_df_1000[1:6,])
trial$Range <- as.factor(c(rep("250kbp",6), rep("1mbp",6), rep("250kbp",6), rep("1mbp",6)))
trial$Lnc_PCG_timeframe <- c(rep("Same timeframe",12), rep("Same/delayed timeframe",12))
HiCSuppTable <- trial

HiCSuppTable$LncTiming <- gsub("Up", "Induced ", HiCSuppTable$LncTiming)
HiCSuppTable$LncTiming <- gsub("Down", "Repressed ", HiCSuppTable$LncTiming)
HiCSuppTable$LncTiming <- gsub("4", "4hr ", HiCSuppTable$LncTiming)
HiCSuppTable$LncTiming <- gsub("8", "8hr ", HiCSuppTable$LncTiming)

HiCSuppTable$LncTiming <- as.factor(HiCSuppTable$LncTiming)
HiCSuppTable$LncTiming <- factor(HiCSuppTable$LncTiming, levels(HiCSuppTable$LncTiming)[c(2,5,3,6,1,4)])

#add annotaiton of sig (bonf corrected x3 as above)
HiCSuppTable$a_annotated <- HiCSuppTable$a
HiCSuppTable$a_annotated[c(1,7)] <- paste0(HiCSuppTable$a_annotated[c(1,7)], "**")
HiCSuppTable$a_annotated[c(13,19)] <- paste0(HiCSuppTable$a_annotated[c(13,19)], "*")

#simple 1mbp figure:
ggplot(filter(HiCSuppTable, Lnc_PCG_timeframe == "Same timeframe", Range == "1mbp")) + aes(x = LncTiming, y = a) +
  geom_bar(stat= "identity", position = position_dodge(width=0.8), fill = "mediumorchid", color = "grey40") +
  ylab("No. HiC Connected\nCClncRNA-Targets\n(Same Timeframe)") +
  xlab("") +
  #geom_text(aes(x = LncTiming, y = a+0.6, label =a),
  #          position = position_dodge(width=0.8), size =4.5) +
  theme_minimal() + 
  scale_y_continuous(limits = c(0,15)) +
  theme(text = element_text(size=20), legend.position = "none") + Seurat::RotatedAxis()


#show lack of pairs, and testing power in later timeframes relative to first:
ggplot(filter(HiCSuppTable, Lnc_PCG_timeframe == "Same timeframe")) + aes(x = LncTiming, y = a, fill = Range) +
  geom_bar(stat= "identity", position = position_dodge(width=0.8), color = "grey40") +
  ylab("No. HiC Connected\nCClncRNA-Targets") +
  xlab("LncRNA Regulation\nTimeframe") +
  geom_text(aes(x = LncTiming, y = a+0.6, label =a_annotated ),
            position = position_dodge(width=0.8), size =4.5) +
  theme_minimal() + 
  theme(text = element_text(size=20), legend.position = "none") + Seurat::RotatedAxis()

#generally higher % in cclnc-t pairs
#significance and good OR achieved for both 250/1mbp:
ggplot(filter(HiCSuppTable, Lnc_PCG_timeframe == "Same timeframe")) + aes(x = LncTiming, y = OR, fill = Range) +
  geom_bar(stat= "identity", position = position_dodge(width=0.9), color = "grey40") +
  ylab("Odds Ratio - HiC in\nCClncRNA-Targets") +
  xlab("LncRNA Regulation\nTimeframe") +
  geom_hline(yintercept = 1, color = "grey70", linetype = "dashed") +
  geom_text(aes(x = LncTiming, y = OR+0.5, label =a_annotated),position = position_dodge(width=0.9)) +
  theme_minimal() + 
  theme(text = element_text(size=20))#, legend.position = "none") + Seurat::RotatedAxis()

#show weakening of effect with same/delayed
ggplot(filter(HiCSuppTable, Range == "1mbp")) + aes(x = LncTiming, y = OR, fill = Lnc_PCG_timeframe) +
  geom_bar(stat= "identity", position = position_dodge(width=0.9), color = "grey40") +
  ylab("Odds Ratio - HiC in\nCClncRNA-Targets") +
  xlab("LncRNA Regulation\nTimeframe") +
  geom_hline(yintercept = 1, color = "grey70", linetype = "dashed") +
  geom_text(aes(x = LncTiming, y = OR+0.25, label =a_annotated),position = position_dodge(width=0.9)) +
  theme_minimal() + 
  theme(text = element_text(size=20), legend.position = "none") + Seurat::RotatedAxis()

ggplot(filter(HiCSuppTable, Range == "250kbp")) + aes(x = LncTiming, y = OR, fill = Lnc_PCG_timeframe) +
  geom_bar(stat= "identity", position = position_dodge(width=0.9), color = "grey40") +
  ylab("Odds Ratio - HiC in\nCClncRNA-Targets") +
  xlab("LncRNA Regulation\nTimeframe") +
  geom_hline(yintercept = 1, color = "grey70", linetype = "dashed") +
  geom_text(aes(x = LncTiming, y = OR+0.25, label =a_annotated),position = position_dodge(width=0.9)) +
  theme_minimal()  + 
  theme(text = element_text(size=20))#, legend.position = "none") + Seurat::RotatedAxis()

#could make the argument that same/delayed also has good OR and significance
#but significance IS reduced a lot to *, despite having more pairs available


#### delayed pairs only ####

#total coreg pairs (in the test)
CoRegPairs_04_48_24_extended_naiveDel <- filter(AllLNC_AllPCG_1_,
                                             (EnsID %in% fpkm_allGDE_within_4$EnsID & 
                                                EnsID.y %in% fpkm_allGDE$EnsID & 
                                                !EnsID.y %in% fpkm_allGDE_within_4$EnsID) |
                                               (EnsID %in% fpkm_allGDE_within_8$EnsID & 
                                                  EnsID.y %in% fpkm_allGDE_within_24$EnsID))
#531 pairs

#same + delayed timeframe pairings:
LoopFish_DelayedPairs <- list()

for (i in 1:length(DEG_cluster)){
  #for co-reg lncs per cluster, all neighbour pairs capable of forming a HiC pair
  CCLNC_AllPCG_HiCTest <- filter(AllLNC_AllPCG_HiCTest, 
                                 #all lncs with a DE PCG near
                                 EnsID %in% CoRegPairs_04_48_24_extended_naiveDel$EnsID,
                                 #from a given cluster
                                 EnsID %in% DEG_cluster[[i]]$EnsID)
  #of these, how many are co-reg, same/later timeframe
  CCLNC_targetPCG_HiCTest <- filter(CCLNC_AllPCG_HiCTest, 
                                    pairs %in% CoRegPairs_04_48_24_extended_naiveDel$pairs)
  
  #per selection and background, identify the pairs which are HiC looped, count total pairs too
  a <- length(unique(filter(CCLNC_targetPCG_HiCTest, pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs)$pairs))
  b <- length(unique(CCLNC_targetPCG_HiCTest$pairs))
  c <- length(unique(filter(CCLNC_AllPCG_HiCTest, pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs)$pairs))
  d <- length(unique(CCLNC_AllPCG_HiCTest$pairs))
  
  LoopFish_DelayedPairs[[i]] <- data.frame(a,b,c,d,
                                               a/b,
                                               c/d,
                                               fisher.test(data.frame("cisLnc" = c(a, b-a),
                                                                      "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p,
                                               fisher.test(data.frame("cisLnc" = c(a, b-a),
                                                                      "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate)
  colnames(LoopFish_DelayedPairs[[i]]) <- c("a", "b", "c", "d", "SelectPercLooped", "BackPerLooped", "p", "OR")
}

names(LoopFish_DelayedPairs) <- c("Up4", "Up8", "Up24", "Down4", "Down8", "Down24", "Either4", "Either8", "Either24")

LoopFish_DelPairs_df_1000 <- bind_rows(LoopFish_DelayedPairs, .id = "LncTiming")

#up and down 4 figure:
DEL_PCG_type <- rbind("Induced" = LoopFish_DelPairs_df_1000[1,c(2,4,6,7)],
                      "Repressed" = LoopFish_DelPairs_df_1000[4,c(2,4,6,7)])
                      #"Regulated" = LoopFish_SamePairs_df[7,c(2,4,6,7)]

colnames(DEL_PCG_type) <- c("a", "c", "Later DE neighbour", "Any neighbour")
DEL_PCG_type$Type <- rownames(DEL_PCG_type)

library(reshape2)
DEL_PCG_typei <- melt(DEL_PCG_type[,3:5])
DEL_PCG_typei <- cbind(DEL_PCG_typei, melt(DEL_PCG_type[,1:2]))
DEL_PCG_type <- DEL_PCG_typei

DEL_PCG_type$value <- DEL_PCG_type$value*100

colnames(DEL_PCG_type)[c(2,5)] <- c("Paired with", "noPairs")

DEL_PCG_type$Type <- factor(DEL_PCG_type$Type)
DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,1)])

DEL_PCG_type$`Paired with` <- factor(DEL_PCG_type$`Paired with`)
DEL_PCG_type$`Paired with` <- factor(DEL_PCG_type$`Paired with`, levels = levels(DEL_PCG_type$`Paired with`)[c(2,1)])

ggplot(DEL_PCG_type) + aes(y = Type, x = value, fill = `Paired with`) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "grey60") +
  scale_fill_manual(values = c(`Later DE neighbour` = "mediumorchid", `Any neighbour` = "grey60")) +
  #geom_text(aes(y = Type, x = value+0.45, label = noPairs),
  #          position = position_dodge(width = 0.7), color = "black", size =5.5) +
  ylab("") +
  xlab("% HiC Connected\nPairs") +
  scale_x_continuous(breaks = seq(0,4,2), limits = c(0,4)) +
  theme_minimal() +
  #ggtitle("LncRNAs with a DE PCG neighbour", subtitle = "(same timeframe)") +
  theme(text = element_text(size=24))# + Seurat::RotatedAxis()

#### return to HiC validation 1: position of contacts relative to co-regulated cisLnc (requires build up code to be run) - merge across fragments/loops #### 

#all lncs contacted by frag1 in HiCCUPs or FitHiC
#add in info on fragment no.:
#add in column for FitHiC or HiCCUPs p val
Loop1Frag1_LNC_PCG$Loop1Frag1 <- paste0("Frag1:", Loop1Frag1_LNC_PCG$Loop1Frag1)
Loop1Frag1_LNC_PCG$LoopTool <- "FitHiC"
Loop1Frag1_LNC_PCG$FitHiC_q <- Loop1Frag1_LNC_PCG$q.value
Loop1Frag1_LNC_PCG$HiCCUPs_bestFDR <- NA

Loop1Frag2_LNC_PCG$Loop1Frag2 <- paste0("Frag2:", Loop1Frag2_LNC_PCG$Loop1Frag2)
Loop1Frag2_LNC_PCG$LoopTool <- "FitHiC"
Loop1Frag2_LNC_PCG$FitHiC_q <- Loop1Frag2_LNC_PCG$q.value
Loop1Frag2_LNC_PCG$HiCCUPs_bestFDR <- NA

Loop2Frag1_LNC_PCG$Loop2Frag1 <- paste0("Frag1:", Loop2Frag1_LNC_PCG$Loop1Frag1)
Loop2Frag1_LNC_PCG$LoopTool <- "HiCCUPs"
Loop2Frag1_LNC_PCG$HiCCUPs_bestFDR <- Biobase::rowMin(as.matrix(Loop2Frag1_LNC_PCG[,c(9:12)]))
Loop2Frag1_LNC_PCG$FitHiC_q <- NA

Loop2Frag2_LNC_PCG$Loop2Frag2 <- paste0("Frag2:", Loop2Frag2_LNC_PCG$Loop1Frag2)
Loop2Frag2_LNC_PCG$LoopTool <- "HiCCUPs"
Loop2Frag2_LNC_PCG$HiCCUPs_bestFDR <- Biobase::rowMin(as.matrix(Loop2Frag2_LNC_PCG[,c(9:12)]))
Loop2Frag2_LNC_PCG$FitHiC_q <- NA

#combine all together:
colnames(Loop1Frag1_LNC_PCG)[1] <- "Loop_ID"
colnames(Loop1Frag2_LNC_PCG)[1] <- "Loop_ID"
colnames(Loop2Frag1_LNC_PCG)[1] <- "Loop_ID"
colnames(Loop2Frag2_LNC_PCG)[1] <- "Loop_ID"

colnames(Loop1Frag1_LNC_PCG)
colnames(Loop1Frag2_LNC_PCG)
colnames(Loop2Frag1_LNC_PCG)
colnames(Loop2Frag2_LNC_PCG)

#both methods, Frag1, Frag2
Both_Both_LNC_PCG_contactPoints <- rbind(Loop1Frag1_LNC_PCG[,c(1,3,7,13,16:19)], Loop1Frag2_LNC_PCG[,c(1,3,7,13,16:19)],
                                         Loop2Frag1_LNC_PCG[,c(1,3,7,17,20:21,23,22)], Loop2Frag2_LNC_PCG[,c(1,3,7,17,20:21,23,22)])
Both_Both_LNC_PCG_contactPoints <- filter(Both_Both_LNC_PCG_contactPoints, !is.na(pairs))

#match:
length(unique(Both_Both_LNC_PCG_contactPoints$pairs))
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))

#all HiC lncRNA contact point for all HiC-connected lnc-PCGs:
ggplot(unique(filter(Both_Both_LNC_PCG_contactPoints)[,c(2,4)])) + 
  aes(x = lncLocContact_perc) + #), fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme_minimal() + 
  theme(legend.position = "none") +
  xlab("All LncRNAs HiC Contact Point\n(% of transcript length)") +
  ylab("Freq.")+
  scale_x_continuous(breaks = seq(0,100,20)) #+scale_y_continuous(limits = c(0,97))
length(unique(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$EnsID))
98/200
#49% of hic contacted lncs have a TTS contact

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin <- NA
Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<0] <- "TSS"
Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>100] <- "TTS"

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>0 & 
                                                         Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<=10] <- "0-10%"

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>10 & 
                                                         Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<=20] <- "10-20%"

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>20 & 
                                                         Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<=30] <- "20-30%"

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>30 & 
                                                         Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<=40] <- "30-40%"

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>40 & 
                                                         Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<=50] <- "40-50%"

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>50 & 
                                                         Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<=60] <- "50-60%"

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>60 & 
                                                         Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<=70] <- "60-70%"

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>70 & 
                                                         Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<=80] <- "70-80%"

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>80 & 
                                                         Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<=90] <- "80-90%"

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>90 & 
                                                         Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<=100] <- "90-100%"
table(Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin)

ggplot(Both_Both_LNC_PCG_contactPoints) + 
  aes(x = lncLocContact_perc_bin, y = -log10(FitHiC_q)) + #), fill = EnsID.x) +
  geom_boxplot() + theme_minimal() + 
  theme(legend.position = "none") +
#  coord_cartesian(ylim = c(0,10)) +
  xlab("All LncRNAs HiC Contact Point\n(% of transcript length)") +
  ylab("FitHiC_q") + Seurat::RotatedAxis()

ggplot(Both_Both_LNC_PCG_contactPoints) + 
  aes(x = lncLocContact_perc_bin, y = -log10(HiCCUPs_bestFDR)) + #), fill = EnsID.x) +
  geom_boxplot() + theme_minimal() + 
  geom_point(alpha = 0.05) +
  theme(legend.position = "none") +
  #coord_cartesian(ylim = c(0,10)) +
  xlab("All LncRNAs HiC Contact Point\n(% of transcript length)") +
  ylab("HiCCUPs_bestFDR") + Seurat::RotatedAxis()

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin2 <- NA
Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin2[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<0] <- "TSS"
Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin2[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>100] <- "TTS"

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin2[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>0 & 
                                                         Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<=20] <- "0-20%"

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin2[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>20 & 
                                                         Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<=40] <- "20-40%"

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin2[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>40 & 
                                                         Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<=60] <- "40-60%"

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin2[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>60 & 
                                                         Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<=80] <- "60-80%"

Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin2[Both_Both_LNC_PCG_contactPoints$lncLocContact_perc>80 & 
                                                         Both_Both_LNC_PCG_contactPoints$lncLocContact_perc<=100] <- "80-100%"
table(Both_Both_LNC_PCG_contactPoints$lncLocContact_perc_bin2)

ggplot(Both_Both_LNC_PCG_contactPoints) + 
  aes(x = lncLocContact_perc_bin2, y = -log10(FitHiC_q)) + #), fill = EnsID.x) +
  geom_boxplot() + theme_minimal() + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0,10)) +
  xlab("All LncRNAs HiC Contact Point\n(% of transcript length)") +
  ylab("FitHiC_q") + Seurat::RotatedAxis()

ggplot(Both_Both_LNC_PCG_contactPoints) + 
  aes(x = lncLocContact_perc_bin2, y = -log10(HiCCUPs_bestFDR)) + #), fill = EnsID.x) +
  geom_boxplot() + theme_minimal() + 
  geom_point(alpha = 0.05) +
  theme(legend.position = "none") +
  #coord_cartesian(ylim = c(0,10)) +
  xlab("All LncRNAs HiC Contact Point\n(% of transcript length)") +
  ylab("HiCCUPs_bestFDR") + Seurat::RotatedAxis()

#some slight sign of increased contact strength towards 3' end of lncRNA, and decreased towards 5'

#same timeframe CClncRNAs:
EarlyInducedLncs_HiC_confirmed <- filter(AllLNC_AllPCG_2d3d, Lnc_Cluster == "Induced <4hrs", PCG_Timeframe == "<4hrs", !loopMethod == "Neither")
EarlyRepressedLncs_HiC_confirmed <- filter(AllLNC_AllPCG_2d3d, Lnc_Cluster == "Repressed <4hrs", PCG_Timeframe == "<4hrs", !loopMethod == "Neither")
EarlyLncs_HiC_confirmed <- filter(AllLNC_AllPCG_2d3d, Lnc_Timeframe == "<4hrs", PCG_Timeframe == "<4hrs", !loopMethod == "Neither")
unique(EarlyLncs_HiC_confirmed$EnsID)

#early same timeframe cclncs
ggplot(unique(filter(Both_Both_LNC_PCG_contactPoints, 
                     pairs %in% EarlyLncs_HiC_confirmed$pairs)[,c(2,4)])) + 
  aes(x = lncLocContact_perc) +#, fill = EnsID.x) +
  geom_histogram(color = "black", fill = "olivedrab3", bins = 22) + theme_minimal() +
  xlab("0-4hr CCLncRNAs HiC Contact Point\n(% of transcript length)") +
  ylab("Freq.") +
  scale_y_continuous(limits = c(0,13)) +
  scale_x_continuous(breaks = seq(0,100,20))# + Seurat::RotatedAxis()
unique(EarlyLncs_HiC_confirmed$EnsID)
11/17 #TTS contacted

#early induced same timeframe cclncs
ggplot(unique(filter(Both_Both_LNC_PCG_contactPoints, 
                     pairs %in% EarlyInducedLncs_HiC_confirmed$pairs)[,c(2,4)])) + 
  aes(x = lncLocContact_perc) +#, fill = EnsID.x) +
  geom_histogram(color = "black", fill = "olivedrab3",bins = 22) + theme_minimal() +
  xlab("All LncRNAs HiC Contact Point\n(% of transcript length)") +
  ylab("Freq.") +
  scale_y_continuous(limits = c(0,8))+
  scale_x_continuous(breaks = seq(0,100,20))
unique(EarlyInducedLncs_HiC_confirmed$EnsID)
7/12 #TTS contacted

#early repressed same timeframe cclncs
ggplot(unique(filter(Both_Both_LNC_PCG_contactPoints, 
                     pairs %in% EarlyRepressedLncs_HiC_confirmed$pairs)[,c(1,3)])) + 
  aes(x = lncLocContact_perc) +#, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme_minimal() +
  xlab("All LncRNAs HiC Contact Point\n(% of transcript length)") +
  ylab("Freq.") +
  scale_y_continuous(limits = c(0,8))
unique(EarlyRepressedLncs_HiC_confirmed$EnsID)
4/5 #TTS contacted

#are TTS contact points stronger than others?
Both_Both_LNC_PCG_contactPoints_earlyCC <- filter(Both_Both_LNC_PCG_contactPoints, pairs %in% EarlyLncs_HiC_confirmed$pairs)

ggplot(Both_Both_LNC_PCG_contactPoints_earlyCC) + 
  aes(x = lncLocContact_perc_bin2, y = -log10(FitHiC_q)) + #), fill = EnsID.x) +
  geom_boxplot() + theme_minimal() + 
  geom_point(alpha = 0.05) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0,20)) +
  xlab("All LncRNAs HiC Contact Point\n(% of transcript length)") +
  ylab("FitHiC_q") + Seurat::RotatedAxis()

ggplot(Both_Both_LNC_PCG_contactPoints_earlyCC) + 
  aes(x = lncLocContact_perc_bin2, y = -log10(HiCCUPs_bestFDR)) + #), fill = EnsID.x) +
  geom_boxplot() + theme_minimal() + 
  geom_point(alpha = 0.05) +
  theme(legend.position = "none") +
  #coord_cartesian(ylim = c(0,10)) +
  xlab("All LncRNAs HiC Contact Point\n(% of transcript length)") +
  ylab("HiCCUPs_bestFDR") + Seurat::RotatedAxis()

#summary: HiC contact points to lncRNAs are non-random
#particularly lnc-PCG connections, and cclnc-pcg connections which occur at TSS/TES primarily
#some signs of increased strength towards 3'

#see old code for frag1/2 diffs:
#Frag2 and Frag1 have differing biases, Frag2 more TES biased and less TSS biased - same pattern in FitHiC and HiCCUPs
#unclear why this would be, some bio process related to +ve/-ve strand? 
#more likely an artefact during data processing - to investigate in other HiC data/lnc overlaps

#notable for early cclncs where HiC has v. strong support of the pairings....
#the TES bias is kept, for Frag1 and particularly Frag2, whilst TSS bias falls away in Frag1
#way more frag2 contacts for early cclncs


#
unique(filter(Both_Both_LNC_PCG_contactPoints, lncLocContact_perc == 110,
              EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
              pairs %in% CoRegPairs_04_48_24_extended$pairs)$EnsID.x)


#which direction are the lncs heading?
ggplot(unique(filter(Both_Frag1_LNC_PCG_contactPoints, pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(1,2,3)])) + 
  aes(x = lncLocContact_perc, fill = str) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  scale_fill_manual(values = c("+" = "red", "-" = "blue")) +
  xlab("LncRNA Contact Point (% of transcript length)")

ggplot(unique(filter(Both_Frag1_LNC_PCG_contactPoints, EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                     pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(1,2,3)])) + 
  aes(x = lncLocContact_perc, fill = str) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  scale_fill_manual(values = c("+" = "red", "-" = "blue")) +
  xlab("LncRNA Contact Point (% of transcript length)")

#For Frag1: 
#the TSS contacted lncs are being transcribed from just outside of the loop edge and moving away from the loop (-ve direction)
#the TES contacted lncs are being transcribed from far away from the loop edge and move towards it (+ve direction)


ggplot(unique(filter(Both_Frag2_LNC_PCG_contactPoints, pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(1,2,3)])) + 
  aes(x = lncLocContact_perc, fill = str) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  scale_fill_manual(values = c("+" = "red", "-" = "blue")) +
  xlab("LncRNA Contact Point (% of transcript length)")

ggplot(unique(filter(Both_Frag2_LNC_PCG_contactPoints, EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                     pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(1,2,3)])) + 
  aes(x = lncLocContact_perc, fill = str) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  scale_fill_manual(values = c("+" = "red", "-" = "blue")) +
  xlab("LncRNA Contact Point (% of transcript length)")

#For Frag2: 
#TSS contacted lncRNAs in a minority
#the TES contacted lncs are being transcribed from inside the loop and finish at the outer edge


#some quite odd patterning, it's unclear why TES enriched at F2 only and not at F1
#lncRNAs give these loops directionality...? loops tend to follow the direction of their lncRNA? both +ve, both contained within the loop

#expecting some kind of CTCF link too... TBC (there will be an RBP/TF thing at somepoint)

#REPEAT FOR EMPHASIS
#most likely explanation of frag2 bias is artefactual
#more likely an artefact during data processing - to investigate in other HiC data/lnc overlaps...


#### old code ####
#plot for hic support in early, 0-4hr induced/repressed/nonDE
#collect a non-DE lnc comparison for 0-4hr to display alongside (later timepoint can be displayed without e.g. in supplement)
CCLNC_AllPCG_HiCTest <- filter(AllLNC_AllPCG_HiCTest,
                               EnsID %in% filter(fpkm_allG_04, !EnsID %in% fpkm_allGDE$EnsID)$EnsID)
#no. non-DEL neighbour pairs where the neighbour is 0-4hr  DEs
CCLNC_targetPCG_HiCTest <- filter(CCLNC_AllPCG_HiCTest, 
                                  EnsID.y %in% fpkm_allGDE$EnsID)

a <- length(unique(filter(CCLNC_targetPCG_HiCTest, pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs)$pairs))
b <- length(unique(CCLNC_targetPCG_HiCTest$pairs))
c <- length(unique(filter(CCLNC_AllPCG_HiCTest, pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs)$pairs))
d <- length(unique(CCLNC_AllPCG_HiCTest$pairs))

data.frame(a/b,
           c/d,
           "p" = fisher.test(data.frame("cisLnc" = c(a, b-a),
                                        "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p,
           "or" =fisher.test(data.frame("cisLnc" = c(a, b-a),
                                        "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate)

nonDElncSameLater_hicP <- fisher.test(data.frame("cisLnc" = c(a, b-a),
                                                 "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p

DEL_PCG_type <- rbind("Induced" = LoopFish_SameDelayedPairs_df_400[1,c(2,4,6,7)],
                      "Repressed" = LoopFish_SameDelayedPairs_df_400[4,c(2,4,6,7)],
                      "Non-DE" = c(a,c,a/b,c/d)
                      #"Regulated" = LoopFish_SamePairs_df[7,c(2,4,6,7)]
)
colnames(DEL_PCG_type) <- c("a", "c", "Any DE neighbour", "Any neighbour")
DEL_PCG_type$Type <- rownames(DEL_PCG_type)

library(reshape2)
DEL_PCG_typei <- melt(DEL_PCG_type[,3:5])
DEL_PCG_typei <- cbind(DEL_PCG_typei, melt(DEL_PCG_type[,1:2]))
DEL_PCG_type <- DEL_PCG_typei

DEL_PCG_type$value <- DEL_PCG_type$value*100

colnames(DEL_PCG_type)[c(2,5)] <- c("Paired with", "noPairs")

DEL_PCG_type$Type <- factor(DEL_PCG_type$Type)
DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,3,1)])

DEL_PCG_type$`Paired with` <- factor(DEL_PCG_type$`Paired with`)
DEL_PCG_type$`Paired with` <- factor(DEL_PCG_type$`Paired with`, levels = levels(DEL_PCG_type$`Paired with`)[c(2,1)])

ggplot(DEL_PCG_type) + aes(y = Type, x = value, fill = `Paired with`, label = noPairs) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "grey60") +
  scale_fill_manual(values = c(`Any DE neighbour` = "mediumorchid", `Any neighbour` = "grey60")) +
  geom_label(position = position_dodge(width = 0.7), color = "white", size =5.5) +
  ylab("") +
  xlab("% HiC Connected") +
  scale_x_continuous(breaks = seq(0,15,5), limits = c(0,12.5)) +
  theme_minimal() +
  #ggtitle("LncRNAs with a DE PCG neighbour", subtitle = "(same timeframe)") +
  theme(text = element_text(size=24))# + Seurat::RotatedAxis()
p.adjust(c(LoopFish_SameDelayedPairs_df_400$p[c(1,4)],nonDElncSameLater_hicP), method = "bonferroni")


#### HiC validation 2 - later timeframe plots for supp ####

#plot for hic support in 4-8hr induced/repressed/nonDE
DEL_PCG_type <- rbind("Induced" = LoopFish_SamePairs_df[2,c(2,4,6,7)],
                      "Repressed" = LoopFish_SamePairs_df[5,c(2,4,6,7)]
                      #"Regulated" = LoopFish_SamePairs_df[7,c(2,4,6,7)]
)

colnames(DEL_PCG_type) <- c("a", "c", "Same timeframe neighbour", "Any neighbour")
DEL_PCG_type$Type <- rownames(DEL_PCG_type)

library(reshape2)
DEL_PCG_typei <- melt(DEL_PCG_type[,3:5])
DEL_PCG_typei <- cbind(DEL_PCG_typei, melt(DEL_PCG_type[,1:2]))
DEL_PCG_type <- DEL_PCG_typei

DEL_PCG_type$value <- DEL_PCG_type$value*100

colnames(DEL_PCG_type)[c(2,5)] <- c("Paired with", "noPairs")

DEL_PCG_type$Type <- factor(DEL_PCG_type$Type)
DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,3,1)])

DEL_PCG_type$`Paired with` <- factor(DEL_PCG_type$`Paired with`)
DEL_PCG_type$`Paired with` <- factor(DEL_PCG_type$`Paired with`, levels = levels(DEL_PCG_type$`Paired with`)[c(2,1)])

ggplot(DEL_PCG_type) + aes(y = Type, x = value, fill = `Paired with`) + #, label = noPairs) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "grey60") +
  scale_fill_manual(values = c(`Same timeframe neighbour` = "mediumorchid", `Any neighbour` = "grey60")) +
 # geom_label(position = position_dodge(width = 0.7), color = "white", size =5) +
  ylab("") +
  xlab("% HiC Connected") +
  scale_x_continuous(breaks = seq(0,15,5), limits = c(0,10)) +
  theme_minimal() +
  #ggtitle("LncRNAs with a DE PCG neighbour", subtitle = "(same timeframe)") +
  theme(text = element_text(size=24))# + Seurat::RotatedAxis()
p.adjust(c(LoopFish_SamePairs_df$p[c(2,5)]), method = "bonferroni")


#plot for hic support in 8-24hr induced/repressed/nonDE
DEL_PCG_type <- rbind("Induced" = LoopFish_SamePairs_df[3,c(2,4,6,7)],
                      "Repressed" = LoopFish_SamePairs_df[6,c(2,4,6,7)]
                      #"Regulated" = LoopFish_SamePairs_df[7,c(2,4,6,7)]
)

colnames(DEL_PCG_type) <- c("a", "c", "8-24hr DE neighbour", "Any neighbour")
DEL_PCG_type$Type <- rownames(DEL_PCG_type)

library(reshape2)
DEL_PCG_typei <- melt(DEL_PCG_type[,3:5])
DEL_PCG_typei <- cbind(DEL_PCG_typei, melt(DEL_PCG_type[,1:2]))
DEL_PCG_type <- DEL_PCG_typei

DEL_PCG_type$value <- DEL_PCG_type$value*100

colnames(DEL_PCG_type)[c(2,5)] <- c("Paired with", "noPairs")

DEL_PCG_type$Type <- factor(DEL_PCG_type$Type)
DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,3,1)])

DEL_PCG_type$`Paired with` <- factor(DEL_PCG_type$`Paired with`)
DEL_PCG_type$`Paired with` <- factor(DEL_PCG_type$`Paired with`, levels = levels(DEL_PCG_type$`Paired with`)[c(2,1)])

ggplot(DEL_PCG_type) + aes(y = Type, x = value, fill = `Paired with`) +#, label = noPairs) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "grey60") +
  scale_fill_manual(values = c(`8-24hr DE neighbour` = "mediumorchid", `Any neighbour` = "grey60")) +
  #geom_label(position = position_dodge(width = 0.7), color = "white", size =5) +
  ylab("") +
  xlab("% HiC Connected") +
  scale_x_continuous(breaks = seq(0,15,5), limits = c(0,10)) +
  theme_minimal() +
  #ggtitle("LncRNAs with a DE PCG neighbour", subtitle = "(same timeframe)") +
  theme(text = element_text(size=24))# + Seurat::RotatedAxis()
p.adjust(c(LoopFish_SamePairs_df$p[c(3,6)]), method = "bonferroni")


#plot for hic support in 4-8hr induced/repressed/nonDE
DEL_PCG_type <- rbind("Induced" = LoopFish_SameDelayedPairs_df_400[2,c(2,4,6,7)],
                      "Repressed" = LoopFish_SameDelayedPairs_df_400[5,c(2,4,6,7)]
                      #"Regulated" = LoopFish_SamePairs_df[7,c(2,4,6,7)]
                      )

colnames(DEL_PCG_type) <- c("a", "c", "Same/later timeframe neighbour", "Any neighbour")
DEL_PCG_type$Type <- rownames(DEL_PCG_type)

library(reshape2)
DEL_PCG_typei <- melt(DEL_PCG_type[,3:5])
DEL_PCG_typei <- cbind(DEL_PCG_typei, melt(DEL_PCG_type[,1:2]))
DEL_PCG_type <- DEL_PCG_typei

DEL_PCG_type$value <- DEL_PCG_type$value*100

colnames(DEL_PCG_type)[c(2,5)] <- c("Paired with", "noPairs")

DEL_PCG_type$Type <- factor(DEL_PCG_type$Type)
DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,3,1)])

DEL_PCG_type$`Paired with` <- factor(DEL_PCG_type$`Paired with`)
DEL_PCG_type$`Paired with` <- factor(DEL_PCG_type$`Paired with`, levels = levels(DEL_PCG_type$`Paired with`)[c(2,1)])

ggplot(DEL_PCG_type) + aes(y = Type, x = value, fill = `Paired with`, label = noPairs) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "grey60") +
  scale_fill_manual(values = c(`Same/later timeframe neighbour` = "mediumorchid", `Any neighbour` = "grey60")) +
  geom_label(position = position_dodge(width = 0.7), color = "white", size =5) +
  ylab("") +
  xlab("% HiC Connected") +
  scale_x_continuous(breaks = seq(0,15,10), limits = c(0,15)) +
  theme_minimal() +
  #ggtitle("LncRNAs with a DE PCG neighbour", subtitle = "(same timeframe)") +
  theme(text = element_text(size=24))# + Seurat::RotatedAxis()
p.adjust(c(LoopFish_SameDelayedPairs_df_400$p[c(2,5)]), method = "bonferroni")


#plot for hic support in 4-8hr induced/repressed/nonDE
DEL_PCG_type <- rbind("Induced" = LoopFish_SameDelayedPairs_df_400[3,c(2,4,6,7)],
                      "Repressed" = LoopFish_SameDelayedPairs_df_400[6,c(2,4,6,7)]
                      #"Regulated" = LoopFish_SamePairs_df[7,c(2,4,6,7)]
)

colnames(DEL_PCG_type) <- c("a", "c", "Same/later timeframe neighbour", "Any neighbour")
DEL_PCG_type$Type <- rownames(DEL_PCG_type)

library(reshape2)
DEL_PCG_typei <- melt(DEL_PCG_type[,3:5])
DEL_PCG_typei <- cbind(DEL_PCG_typei, melt(DEL_PCG_type[,1:2]))
DEL_PCG_type <- DEL_PCG_typei

DEL_PCG_type$value <- DEL_PCG_type$value*100

colnames(DEL_PCG_type)[c(2,5)] <- c("Paired with", "noPairs")

DEL_PCG_type$Type <- factor(DEL_PCG_type$Type)
DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,3,1)])

DEL_PCG_type$`Paired with` <- factor(DEL_PCG_type$`Paired with`)
DEL_PCG_type$`Paired with` <- factor(DEL_PCG_type$`Paired with`, levels = levels(DEL_PCG_type$`Paired with`)[c(2,1)])

ggplot(DEL_PCG_type) + aes(y = Type, x = value, fill = `Paired with`, label = noPairs) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "grey60") +
  scale_fill_manual(values = c(`Same/later timeframe neighbour` = "mediumorchid", `Any neighbour` = "grey60")) +
  geom_label(position = position_dodge(width = 0.7), color = "white", size =5) +
  ylab("") +
  xlab("% HiC Connected") +
  scale_x_continuous(breaks = seq(0,15,10), limits = c(0,15)) +
  theme_minimal() +
  #ggtitle("LncRNAs with a DE PCG neighbour", subtitle = "(same timeframe)") +
  theme(text = element_text(size=24))# + Seurat::RotatedAxis()
p.adjust(c(LoopFish_SameDelayedPairs_df_400$p[c(2,5)]), method = "bonferroni")


#### overall numbers of hic connected genes etc ####

AllLNC_AllPCG_2d3d <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_2d3d_Aug2025.csv")
dim(filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither"))
179/2552

#total coreg pairs (in the test)
CoRegPairs_04_48_24_extended <- filter(AllLNC_AllPCG_2d3d,
                                       #AllLNC_AllPCG_1,
                                       (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                     fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                       fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                        fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                       fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
)#364 @400kbp

#same timeframe equivalent:
CoRegPairs_04_48_24_extendedSame <- filter(AllLNC_AllPCG_2d3d,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                         fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                          fpkm_allGDE_Downwithin_4$EnsID)) |
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                           fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID)) |
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                           fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
)#200 @400kbp

unique(filter(CoRegPairs_04_48_24_extended)$EnsID) #115 total

unique(filter(CoRegPairs_04_48_24_extended, !loopMethod == "Neither")$EnsID) #32 HiC looped

unique(filter(CoRegPairs_04_48_24_extended, EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID))$EnsID) #47 up
unique(filter(CoRegPairs_04_48_24_extended, EnsID %in% c(fpkm_allGDE_Downwithin_4$EnsID))$EnsID) #27 down
unique(filter(CoRegPairs_04_48_24_extended, EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID))$EnsID) #74 total

unique(filter(CoRegPairs_04_48_24_extended, !loopMethod == "Neither", 
              EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID))$EnsID) #24 HiC looped
unique(filter(CoRegPairs_04_48_24_extended, !loopMethod == "Neither", 
              EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID))$EnsID) #24 HiC looped
unique(filter(CoRegPairs_04_48_24_extended, !loopMethod == "Neither", 
              EnsID %in% c(fpkm_allGDE_Downwithin_4$EnsID))$EnsID) #24 HiC looped


#### return to HiC validation 1: position of contacts relative to co-regulated cisLnc (requires build up code to be run) - fragment/loop checks #### 

#is the HiC contact position non-random? helps to indicate biological function and validate the connections are not some technical thing

#For FitHiC
length(unique(Loop1Frag1_LNC_PCG$EnsID.x))#137 lncs are contacted via frag1
length(unique(Loop1Frag2_LNC_PCG$EnsID.x))#136 lncs are contacted via frag2
length(unique(c(Loop1Frag1_LNC_PCG$EnsID.x, 
                Loop1Frag2_LNC_PCG$EnsID.x)))#217 lncs are hic-contacted

length(unique(filter(Loop1Frag1_LNC_PCG, !is.na(PCG))$EnsID.x))#63 lncs loop-connected via frag1 to a PCG 
length(unique(filter(Loop1Frag2_LNC_PCG, !is.na(PCG))$EnsID.x))#76 lncs loop-connected via frag2 to a PCG
length(unique(c(filter(Loop1Frag1_LNC_PCG, !is.na(PCG))$EnsID.x, 
                filter(Loop1Frag2_LNC_PCG, !is.na(PCG))$EnsID.x)))#128 lncs are hic-contacted to a pcg

#loops where lncRNA overlaps a FitHiC Frag
#using histograms of 22 bins: 20 bins across the 100% then 2 for TTS and TES
#frag1 contacts (always the most upstream on the sense DNA strand)
ggplot(Loop1Frag1_LNC_PCG) + aes(x = lncLocContact_perc) +
  geom_histogram(bins = 22)
#frag2 contacts (always the most downstream on the sense DNA strand)
ggplot(Loop1Frag2_LNC_PCG) + aes(x = lncLocContact_perc) +
  geom_histogram(bins = 22)
#lnc loops at frag1 biased to behind TTS and after TES, frag2 v. highly biased to the TES
#why would frags show differences? frag1 is always upstream of frag2 on sense strand
#analysis quirk (most likely)? or some reason why one strand is favoured?

#show contact points per lnc:
head(unique(filter(Loop1Frag1_LNC_PCG)[,c(3,13)]))
ggplot(unique(filter(Loop1Frag1_LNC_PCG)[,c(3,13)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + #theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")
ggplot(unique(filter(Loop1Frag2_LNC_PCG)[,c(3,13)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + #theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")


#subset to lnc-PCG contact loops
ggplot(filter(Loop1Frag1_LNC_PCG, !is.na(PCG))) + aes(x = lncLocContact_perc) +
  geom_histogram(bins = 22)
ggplot(filter(Loop1Frag2_LNC_PCG, !is.na(PCG))) + aes(x = lncLocContact_perc) +
  geom_histogram(bins = 22)
#biases to extremities become slightly stronger?

#show contact points per lnc:
ggplot(unique(filter(Loop1Frag1_LNC_PCG, !is.na(PCG))[,c(3,13)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")
ggplot(unique(filter(Loop1Frag2_LNC_PCG, !is.na(PCG))[,c(3,13)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")


#subset to cclnc-target contact loops
ggplot(filter(Loop1Frag1_LNC_PCG, pairs %in% CoRegPairs_04_48_24_extended$pairs)) + aes(x = lncLocContact_perc) +
  geom_histogram(bins = 22)
ggplot(filter(Loop1Frag2_LNC_PCG, pairs %in% CoRegPairs_04_48_24_extended$pairs)) + aes(x = lncLocContact_perc) +
  geom_histogram(bins = 22)
#biases to extremities maintained, particularly for frag2 

#show contact points per lnc:
ggplot(unique(filter(Loop1Frag1_LNC_PCG, pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(3,13)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")
ggplot(unique(filter(Loop1Frag2_LNC_PCG, pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(3,13)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")


#early cclnc-target loops
ggplot(filter(Loop1Frag1_LNC_PCG, EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
              pairs %in% CoRegPairs_04_48_24_extended$pairs)) + aes(x = lncLocContact_perc) +
  geom_histogram(bins = 22)
ggplot(filter(Loop1Frag2_LNC_PCG, EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
              pairs %in% CoRegPairs_04_48_24_extended$pairs)) + aes(x = lncLocContact_perc) +
  geom_histogram(bins = 22)
#biases to the TES extremity maintained, frag1 no longer has much TSS bias

#show contact points per lnc:
ggplot(unique(filter(Loop1Frag1_LNC_PCG, EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                     pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(3,13)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")
ggplot(unique(filter(Loop1Frag2_LNC_PCG, EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                     pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(3,13)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")



#HiCCUPs loops where lncRNA overlaps a FitHiC Frag
#using histograms of 22 bins: 20 bins across the 100% then 2 for TTS and TES
ggplot(Loop2Frag1_LNC_PCG) + aes(x = lncLocContact_perc) +
  geom_histogram(bins = 22)
ggplot(Loop2Frag2_LNC_PCG) + aes(x = lncLocContact_perc) +
  geom_histogram(bins = 22)
#lnc loops at frag1 biased to behind TTS and after TES, frag2 v. highly biased to the TES
#why would frags show differences? frag1 is always upstream of frag2 on sense strand
#now in 2x analysis methods
#might be something upstream, in the HiC library generation itself?

#show contact points per lnc:
head(unique(filter(Loop2Frag1_LNC_PCG)[,c(3,17)]))
ggplot(unique(filter(Loop2Frag1_LNC_PCG)[,c(3,17)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")
ggplot(unique(filter(Loop2Frag2_LNC_PCG)[,c(3,17)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")


#subset to lnc-PCG contact loops
ggplot(filter(Loop2Frag1_LNC_PCG, !is.na(PCG))) + aes(x = lncLocContact_perc) +
  geom_histogram(bins = 22)
ggplot(filter(Loop2Frag2_LNC_PCG, !is.na(PCG))) + aes(x = lncLocContact_perc) +
  geom_histogram(bins = 22)
#biases to extremities become stronger, presumably less noise, more biological function relevant loops

#show contact points per lnc:
ggplot(unique(filter(Loop2Frag1_LNC_PCG, !is.na(PCG))[,c(3,17)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")
ggplot(unique(filter(Loop2Frag2_LNC_PCG, !is.na(PCG))[,c(3,17)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")


#subset to cclnc-target contact loops
ggplot(filter(Loop2Frag1_LNC_PCG, pairs %in% CoRegPairs_04_48_24_extended$pairs)) + aes(x = lncLocContact_perc) +
  geom_histogram(bins = 22)
ggplot(filter(Loop2Frag2_LNC_PCG, pairs %in% CoRegPairs_04_48_24_extended$pairs)) + aes(x = lncLocContact_perc) +
  geom_histogram(bins = 22)
#biases to extremities maintained, particularly for frag2

#show contact points per lnc:
ggplot(unique(filter(Loop2Frag1_LNC_PCG, pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(3,17)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")
ggplot(unique(filter(Loop2Frag2_LNC_PCG, pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(3,17)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")


#early cclnc-target loops
ggplot(filter(Loop2Frag1_LNC_PCG, EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
              pairs %in% CoRegPairs_04_48_24_extended$pairs)) + aes(x = lncLocContact_perc) +
  geom_histogram(bins = 22)
ggplot(filter(Loop2Frag2_LNC_PCG, EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
              pairs %in% CoRegPairs_04_48_24_extended$pairs)) + aes(x = lncLocContact_perc) +
  geom_histogram(bins = 22)
#biases to the TES extremity maintained, frag1 no longer has much TSS bias

#show contact points per lnc:
ggplot(unique(filter(Loop2Frag1_LNC_PCG, EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                     pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(3,17)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")
ggplot(unique(filter(Loop2Frag2_LNC_PCG, EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                     pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(3,17)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")


#### return to HiC validation 1: position of contacts relative to co-regulated cisLnc (requires build up code to be run) - merge across fragments/loops #### 

#all lncs contacted by frag1 in HiCCUPs or FitHiC
Both_Frag1_LNC_PCG_contactPoints <- unique(rbind(Loop1Frag1_LNC_PCG[,c(3,7,13,16)], Loop2Frag1_LNC_PCG[,c(3,7,17,20)]))

#all cclncs
ggplot(unique(filter(Both_Frag1_LNC_PCG_contactPoints,
                     pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(1,3)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")

#early cclncs
ggplot(unique(filter(Both_Frag1_LNC_PCG_contactPoints, EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                     pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(1,3)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")
#frag1 


#all lncs contacted by frag2 in HiCCUPs or FitHiC
Both_Frag2_LNC_PCG_contactPoints <- unique(rbind(Loop1Frag2_LNC_PCG[,c(3,7,13,16)], Loop2Frag2_LNC_PCG[,c(3,7,17,20)]))

#all cclncs
ggplot(unique(filter(Both_Frag2_LNC_PCG_contactPoints,
                     pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(1,3)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")
#early cclncs
ggplot(unique(filter(Both_Frag2_LNC_PCG_contactPoints, EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                     pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(1,3)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")
#frag2 bias is way more noticeable to TES than frag1


#both Frag1, Frag2
Both_Both_LNC_PCG_contactPoints <- rbind(Both_Frag1_LNC_PCG_contactPoints, Both_Frag2_LNC_PCG_contactPoints)
#all cclncs
ggplot(unique(filter(Both_Both_LNC_PCG_contactPoints,
                     pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(1,3)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("LncRNA Contact Point (% of transcript length)")
#early cclncs
ggplot(unique(filter(Both_Both_LNC_PCG_contactPoints, EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                     pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(1,3)])) + 
  aes(x = lncLocContact_perc, fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  xlab("Early CCLncRNA HiC Contact Point (% of transcript length)")
#early cclncs
ggplot(unique(filter(Both_Both_LNC_PCG_contactPoints, EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                     pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(1,3)])) + 
  aes(x = lncLocContact_perc) + #), fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme_minimal() + 
  theme(legend.position = "none") +
  xlab("Early CCLncRNA HiC Contact Point\n(% of transcript length)") +
  ylab("Freq.")
#all lncs as comparison
ggplot(unique(filter(Both_Both_LNC_PCG_contactPoints,
                     pairs %in% AllLNC_AllPCG_2d3d$pairs)[,c(1,3)])) + 
  aes(x = lncLocContact_perc) + #), fill = EnsID.x) +
  geom_histogram(color = "black", bins = 22) + theme_minimal() + 
  theme(legend.position = "none") +
  xlab("All LncRNAs HiC Contact Point\n(% of transcript length)") +
  ylab("Freq.")
#some enrichment seems likely - keep as picture?


#summary so far: HiC contact points to lncRNAs are non-random
#particularly lnc-PCG connections, and cclnc-pcg connections which occur at TSS/TES primarily
#Frag2 and Frag1 have differing biases, Frag2 more TES biased and less TSS biased - same pattern in FitHiC and HiCCUPs
#unclear why this would be, some bio process related to +ve/-ve strand? 
#more likely an artefact during data processing - to investigate in other HiC data/lnc overlaps

#notable for early cclncs where HiC has v. strong support of the pairings....
#the TES bias is kept, for Frag1 and particularly Frag2, whilst TSS bias falls away in Frag1
#way more frag2 contacts for early cclncs


#
unique(filter(Both_Both_LNC_PCG_contactPoints, lncLocContact_perc == 110,
              EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
              pairs %in% CoRegPairs_04_48_24_extended$pairs)$EnsID.x)


#which direction are the lncs heading?
ggplot(unique(filter(Both_Frag1_LNC_PCG_contactPoints, pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(1,2,3)])) + 
  aes(x = lncLocContact_perc, fill = str) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  scale_fill_manual(values = c("+" = "red", "-" = "blue")) +
  xlab("LncRNA Contact Point (% of transcript length)")

ggplot(unique(filter(Both_Frag1_LNC_PCG_contactPoints, EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                     pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(1,2,3)])) + 
  aes(x = lncLocContact_perc, fill = str) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  scale_fill_manual(values = c("+" = "red", "-" = "blue")) +
  xlab("LncRNA Contact Point (% of transcript length)")

#For Frag1: 
#the TSS contacted lncs are being transcribed from just outside of the loop edge and moving away from the loop (-ve direction)
#the TES contacted lncs are being transcribed from far away from the loop edge and move towards it (+ve direction)


ggplot(unique(filter(Both_Frag2_LNC_PCG_contactPoints, pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(1,2,3)])) + 
  aes(x = lncLocContact_perc, fill = str) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  scale_fill_manual(values = c("+" = "red", "-" = "blue")) +
  xlab("LncRNA Contact Point (% of transcript length)")

ggplot(unique(filter(Both_Frag2_LNC_PCG_contactPoints, EnsID.x %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                     pairs %in% CoRegPairs_04_48_24_extended$pairs)[,c(1,2,3)])) + 
  aes(x = lncLocContact_perc, fill = str) +
  geom_histogram(color = "black", bins = 22) + theme(legend.position = "none") +
  scale_fill_manual(values = c("+" = "red", "-" = "blue")) +
  xlab("LncRNA Contact Point (% of transcript length)")

#For Frag2: 
#TSS contacted lncRNAs in a minority
#the TES contacted lncs are being transcribed from inside the loop and finish at the outer edge


#some quite odd patterning, it's unclear why TES enriched at F2 only and not at F1
#lncRNAs give these loops directionality...? loops tend to follow the direction of their lncRNA? both +ve, both contained within the loop

#expecting some kind of CTCF link too... TBC (there will be an RBP/TF thing at somepoint)

#REPEAT FOR EMPHASIS
#most likely explanation of frag2 bias is artefactual
#more likely an artefact during data processing - to investigate in other HiC data/lnc overlaps...

#######################################
#### not reviewed beyond here yet ####

#delayed, less precise? 
#CoRegPairs_04_48_24_extended_delayed <- filter(AllLNC_AllPCG_2d3d,
#                                               #AllLNC_AllPCG_1,
#                                               (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
#                                                             fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
#                                                                                                              fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
#                                                 (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
#                                                               fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)))

#same timeframe pairings:
#LoopFish_DelayPairs <- list()

#for (i in 1:length(DEG_cluster)){
#just induced
#  CCLNC_AllPCG_HiCTest <- filter(AllLNC_AllPCG_HiCTest, EnsID %in% CoRegPairs_04_48_24_extended_delayed$EnsID, EnsID %in% DEG_cluster[[i]]$EnsID)
#  #just the CClncRNA-target pairs now
#  CCLNC_targetPCG_HiCTest <- filter(CCLNC_AllPCG_HiCTest, pairs %in% CoRegPairs_04_48_24_extended_delayed$pairs, EnsID %in% DEG_cluster[[i]]$EnsID)
#
#  a <- length(unique(filter(CCLNC_targetPCG_HiCTest, pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs)$pairs))
#  b <- length(unique(CCLNC_targetPCG_HiCTest$pairs))
#  c <- length(unique(filter(CCLNC_AllPCG_HiCTest, pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs)$pairs))
#  d <- length(unique(CCLNC_AllPCG_HiCTest$pairs))
#  
#  LoopFish_DelayPairs[[i]] <- data.frame(a,b,c,d,
#                                        a/b,
#                                        c/d,
#                                        fisher.test(data.frame("cisLnc" = c(a, b-a),
#                                                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$p,
#                                        fisher.test(data.frame("cisLnc" = c(a, b-a),
#                                                               "Not"   = c(c-a, d-c-(b-a))), alternative = "greater")$estimate)
#  colnames(LoopFish_DelayPairs[[i]]) <- c("a", "b", "c", "d", "SelectPercLooped", "BackPerLooped", "p", "OR")
#}

#names(LoopFish_DelayPairs) <- c("Up4", "Up8", "Up24", "Down4", "Down8", "Down24", "Either4", "Either8", "Either24")
#LoopFish_DelayPairs_df <- bind_rows(LoopFish_DelayPairs, .id = "LncTiming")


#NonDELnc-DEPCG pairs - 377 pairs
#NonDE_DE <- filter(AllLNC_AllPCG_2d3d, !EnsID %in% fpkm_allGDE$EnsID, EnsID.y %in% fpkm_allGDE$EnsID)
#NonDE_AllPCG_HiCTest <- filter(AllLNC_AllPCG_HiCTest, EnsID %in% NonDE_DE$EnsID)
#NonDE_DEPCG_HiCTest <- filter(NonDE_AllPCG_HiCTest, pairs %in% NonDE_DE$pairs)
#background
#length(unique(NonDE_AllPCG_HiCTest$pairs))#855 neighbours at loci of nonDE lncs to DE PCGs
#length(unique(filter(NonDE_AllPCG_HiCTest, pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs)$pairs)) #44

#selection:
#length(unique(NonDE_DEPCG_HiCTest$pairs))#359
#length(unique(filter(NonDE_DEPCG_HiCTest, pairs %in% filter(AllLNC_AllPCG_2d3d, !loopMethod == "Neither")$pairs)$pairs)) #17
#44/855
#17/359

#aiv <- 17
#biv <- 359
#civ <- 44
#div <- 855

#c(fisher.test(data.frame("cisLnc" = c(aiv, biv-aiv),
#                         "Not"   = c(civ-aiv, div-civ-(biv-aiv))), alternative = "greater")$p,
#  fisher.test(data.frame("cisLnc" = c(aiv, biv-aiv),
#                         "Not"   = c(civ-aiv, div-civ-(biv-aiv))), alternative = "greater")$estimate)
#very non-significant
#shows that effect specifically links DE PCGs to DE lncs not non-DE lncs

#######################################
#### plotting ####

#lncRNA - PCG pairs where the PCG was co-regulated in the same timeframe, a later timeframe or both.
#we tested if the number of HiC connected pairs amongst these sets
#was more than would be expected compared to the number of HiC connected pairs amongst all the PCG neighbours (DE or non-DE) 
#found for each set of lncRNAs
#this showed x

#we also repeated the test for pairs where the PCG was co-regulated in a prior timeframe to the lncRNA showing y

#finally we checked a set of nonDE lncRNAs paired to nonDE PCGs finding z

#this confirmed an overlap between HiC collected from SMCs and particular sets of co-regulated lnc-PCG pairs from the SVSMC

#particularly those where the lncRNA changes at the same time or prior to the PCG

#this suggests loci where the lncRNA could play an initiating function in the change in the PCG 
#are more likely to be enriched with chromatin contacts than others
#and so have additional evidence supporting a cis-acting functional role
#timing of lnc and target regulation may be an effective criteria to use when identifying cis-acting lncs

#early regulated/induced/repressed
DEL_PCG_type <- rbind("Induced" = LoopFish_SamePairs_df[1,c(2,4,6,7)],
                      "Repressed" = LoopFish_SamePairs_df[4,c(2,4,6,7)],
                      "Regulated" = LoopFish_SamePairs_df[7,c(2,4,6,7)])
colnames(DEL_PCG_type) <- c("a", "c", "CCLncRNA targets\n(same timeframe)", "All neighbours")
DEL_PCG_type$Type <- rownames(DEL_PCG_type)

library(reshape2)
DEL_PCG_typei <- melt(DEL_PCG_type[,3:5])
DEL_PCG_typei <- cbind(DEL_PCG_typei, melt(DEL_PCG_type[,1:2]))
DEL_PCG_type <- DEL_PCG_typei

DEL_PCG_type$Type <- as.factor(DEL_PCG_type$Type)
#DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,3,1)])
DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,1,3)])

DEL_PCG_type$value <- DEL_PCG_type$value*100

colnames(DEL_PCG_type)[c(2,5)] <- c("LncRNA-PCG Pairings", "noPairs")

ggplot(DEL_PCG_type) + aes(y = Type, x = value, fill = `LncRNA-PCG Pairings`) + #, label = noPairs) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "grey60") +
  scale_fill_manual(values = c(`CCLncRNA targets\n(same timeframe)` = "mediumorchid", `All neighbours` = "grey60")) +
  #geom_label(position = position_dodge(width = 0.7), color = "white", size =5) +
  #ylab("\nEarly CClncRNAs\n(with target in same timeframe)") +
  xlab("% HiC Connected") +
  scale_x_continuous(breaks = seq(0,15,5), limits = c(0,15)) +
  theme_minimal() +
  #ggtitle("LncRNAs with a DE PCG neighbour", subtitle = "(same timeframe)") +
  theme(text = element_text(size=24))# + Seurat::RotatedAxis()

p.adjust(LoopFish_SamePairs_df$p[c(1,4,7)], method = "bonferroni")

#medium induced/repressed
DEL_PCG_type <- rbind("Induced" = LoopFish_SamePairs_df[2,c(2,4,6,7)],
                      "Repressed" = LoopFish_SamePairs_df[5,c(2,4,6,7)],
                      "Regulated" = LoopFish_SamePairs_df[8,c(2,4,6,7)])
colnames(DEL_PCG_type) <- c("a", "c", "Co-regulated neighbours\n(same timeframe)", "All neighbours")
DEL_PCG_type$Type <- rownames(DEL_PCG_type)

library(reshape2)
DEL_PCG_typei <- melt(DEL_PCG_type[,3:5])
DEL_PCG_typei <- cbind(DEL_PCG_typei, melt(DEL_PCG_type[,1:2]))
DEL_PCG_type <- DEL_PCG_typei

DEL_PCG_type$Type <- as.factor(DEL_PCG_type$Type)
#DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,3,1)])
DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,1,3)])

DEL_PCG_type$value <- DEL_PCG_type$value*100

colnames(DEL_PCG_type)[c(2,5)] <- c("LncRNA-PCG Pairings", "noPairs")

ggplot(DEL_PCG_type) + aes(x = Type, y = value, fill = `LncRNA-PCG Pairings`, label = noPairs) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "grey60") +
  scale_fill_manual(values = c(`Co-regulated neighbours\n(same timeframe)` = "mediumorchid", `All neighbours` = "grey60")) +
  geom_label(position = position_dodge(width = 0.7), color = "white") +
  xlab("\n4-8hr CClncRNAs\n(with target in same timeframe)") +
  ylab("% HiC Connected") +
  scale_y_continuous(breaks = seq(0,15,5), limits = c(0,15)) +
  theme_minimal() +
  #ggtitle("LncRNAs with a DE PCG neighbour", subtitle = "(same timeframe)") +
  theme(text = element_text(size=16)) + Seurat::RotatedAxis()

#all ns already obvz

#late induced/repressed
DEL_PCG_type <- rbind("Induced" = LoopFish_SamePairs_df[3,c(2,4,6,7)],
                      "Repressed" = LoopFish_SamePairs_df[6,c(2,4,6,7)],
                      "Regulated" = LoopFish_SamePairs_df[9,c(2,4,6,7)])
colnames(DEL_PCG_type) <- c("a", "c", "Co-regulated neighbours\n(same timeframe)", "All neighbours")
DEL_PCG_type$Type <- rownames(DEL_PCG_type)

library(reshape2)
DEL_PCG_typei <- melt(DEL_PCG_type[,3:5])
DEL_PCG_typei <- cbind(DEL_PCG_typei, melt(DEL_PCG_type[,1:2]))
DEL_PCG_type <- DEL_PCG_typei

DEL_PCG_type$Type <- as.factor(DEL_PCG_type$Type)
#DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,3,1)])
DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,1,3)])

DEL_PCG_type$value <- DEL_PCG_type$value*100

colnames(DEL_PCG_type)[c(2,5)] <- c("LncRNA-PCG Pairings", "noPairs")

ggplot(DEL_PCG_type) + aes(x = Type, y = value, fill = `LncRNA-PCG Pairings`, label = noPairs) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "grey60") +
  scale_fill_manual(values = c(`Co-regulated neighbours\n(same timeframe)` = "mediumorchid", `All neighbours` = "grey60")) +
  #geom_label(position = position_dodge(width = 0.7), color = "white") +
  xlab("\n8-24hr CClncRNAs\n(with target in same timeframe target)") +
  ylab("% HiC Connected") +
  scale_y_continuous(breaks = seq(0,15,5), limits = c(0,15)) +
  theme_minimal() +
  #ggtitle("LncRNAs with a DE PCG neighbour", subtitle = "(same timeframe)") +
  theme(text = element_text(size=16)) + Seurat::RotatedAxis()

#all ns obvs


## same/later targets
#early
DEL_PCG_type <- rbind("Induced" = LoopFish_SameDelayedPairs_df[1,c(2,4,6,7)],
                      "Repressed" = LoopFish_SameDelayedPairs_df[4,c(2,4,6,7)],
                      "Regulated" = LoopFish_SameDelayedPairs_df[7,c(2,4,6,7)])
colnames(DEL_PCG_type) <- c("a", "c", "CCLncRNA targets\n(same or later timeframe)", "All neighbours")
DEL_PCG_type$Type <- rownames(DEL_PCG_type)

library(reshape2)
DEL_PCG_typei <- melt(DEL_PCG_type[,3:5])
DEL_PCG_typei <- cbind(DEL_PCG_typei, melt(DEL_PCG_type[,1:2]))
DEL_PCG_type <- DEL_PCG_typei

DEL_PCG_type$Type <- as.factor(DEL_PCG_type$Type)
#DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,3,1)])
DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,1,3)])

DEL_PCG_type$value <- DEL_PCG_type$value*100

colnames(DEL_PCG_type)[c(2,5)] <- c("LncRNA-PCG Pairings", "noPairs")

ggplot(DEL_PCG_type) + aes(x = Type, y = value, fill = `LncRNA-PCG Pairings`, label = noPairs) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "grey60") +
  scale_fill_manual(values = c(`CCLncRNA targets\n(same or later timeframe)` = "orchid2", `All neighbours` = "grey60")) +
  geom_label(position = position_dodge(width = 0.7), color = "white", size = 5) +
  xlab("\nEarly CCLncRNAs\n(with target in same/later timeframe)") +
  ylab("% HiC Connected") +
  scale_y_continuous(breaks = seq(0,15,5), limits = c(0,15)) +
  theme_minimal() +
  #ggtitle("LncRNAs with a DE PCG neighbour", subtitle = "(same timeframe)") +
  theme(text = element_text(size=24))# + Seurat::RotatedAxis()

p.adjust(LoopFish_SameDelayedPairs_df$p[c(1,4,7)], method = "bonferroni")


#medium
DEL_PCG_type <- rbind("Induced" = LoopFish_SameDelayedPairs_df[2,c(2,4,6,7)],
                      "Repressed" = LoopFish_SameDelayedPairs_df[5,c(2,4,6,7)],
                      "Regulated" = LoopFish_SameDelayedPairs_df[8,c(2,4,6,7)])
colnames(DEL_PCG_type) <- c("a", "c", "Co-regulated neighbours\n(same/later timeframe)", "All neighbours")
DEL_PCG_type$Type <- rownames(DEL_PCG_type)

library(reshape2)
DEL_PCG_typei <- melt(DEL_PCG_type[,3:5])
DEL_PCG_typei <- cbind(DEL_PCG_typei, melt(DEL_PCG_type[,1:2]))
DEL_PCG_type <- DEL_PCG_typei

DEL_PCG_type$Type <- as.factor(DEL_PCG_type$Type)
#DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,3,1)])
DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,1,3)])

DEL_PCG_type$value <- DEL_PCG_type$value*100

colnames(DEL_PCG_type)[c(2,5)] <- c("LncRNA-PCG Pairings", "noPairs")

ggplot(DEL_PCG_type) + aes(x = Type, y = value, fill = `LncRNA-PCG Pairings`, label = noPairs) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "grey60") +
  scale_fill_manual(values = c(`Co-regulated neighbours\n(same/later timeframe)` = "orchid1", `All neighbours` = "grey60")) +
  #geom_label(position = position_dodge(width = 0.7), color = "white") +
  xlab("\n4-8hr CClncRNAs\n(same/later timeframe target)") +
  ylab("% HiC Connected") +
  scale_y_continuous(breaks = seq(0,20,5), limits = c(0,20)) +
  theme_minimal() +
  #ggtitle("LncRNAs with a DE PCG neighbour", subtitle = "(same timeframe)") +
  theme(text = element_text(size=16)) + Seurat::RotatedAxis()

p.adjust(LoopFish_SameDelayedPairs_df$p[c(2,5,8)], method = "bonferroni")


#later
DEL_PCG_type <- rbind("Induced" = LoopFish_SameDelayedPairs_df[3,c(2,4,6,7)],
                      "Repressed" = LoopFish_SameDelayedPairs_df[6,c(2,4,6,7)],
                      "Regulated" = LoopFish_SameDelayedPairs_df[9,c(2,4,6,7)])
colnames(DEL_PCG_type) <- c("a", "c", "Co-regulated neighbours", "All neighbours")
DEL_PCG_type$Type <- rownames(DEL_PCG_type)

library(reshape2)
DEL_PCG_typei <- melt(DEL_PCG_type[,3:5])
DEL_PCG_typei <- cbind(DEL_PCG_typei, melt(DEL_PCG_type[,1:2]))
DEL_PCG_type <- DEL_PCG_typei

DEL_PCG_type$Type <- as.factor(DEL_PCG_type$Type)
#DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,3,1)])
DEL_PCG_type$Type <- factor(DEL_PCG_type$Type, levels = levels(DEL_PCG_type$Type)[c(2,1,3)])

DEL_PCG_type$value <- DEL_PCG_type$value*100

colnames(DEL_PCG_type)[c(2,5)] <- c("LncRNA-PCG Pairings", "noPairs")

ggplot(DEL_PCG_type) + aes(x = Type, y = value, fill = `LncRNA-PCG Pairings`, label = noPairs) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "grey60") +
  scale_fill_manual(values = c(`Co-regulated neighbours` = "orchid1", `All neighbours` = "grey60")) +
  geom_label(position = position_dodge(width = 0.7), color = "white") +
  xlab("\n8-24hr CClncRNAs\n(with target in same/later timeframe)") +
  ylab("% HiC Connected") +
  #scale_y_continuous(breaks = seq(0,15,5), limits = c(0,15)) +
  theme_minimal() +
  #ggtitle("LncRNAs with a DE PCG neighbour", subtitle = "(same timeframe)") +
  theme(text = element_text(size=16)) + Seurat::RotatedAxis()

p.adjust(LoopFish_SameDelayedPairs_df$p[c(2,5,8)], method = "bonferroni")


#### (save for after eQTL) HiC validation 3: recovery of known cis lncs in the HiC pairs ####

ControlCisLncs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/ControlCisLncs.csv")[,-1]

#merge:
fpkm_allG$EnsID_merge <- gsub("\\.[0-9]*", "", fpkm_allG$EnsID)

#27x lncs with strong signs of cis-activity from perturbation expressed somewhere in the dataset:
dim(unique(filter(fpkm_allG, EnsID_merge %in% ControlCisLncs$Ens_ID)[,c(2,5)]))

#9x are also in the caSMC FitHiC loops connecting lncs and PCGs: 
trial <- unique(filter(fpkm_allG, EnsID_merge %in% ControlCisLncs$Ens_ID)$EnsID[
  filter(fpkm_allG, EnsID_merge %in% ControlCisLncs$Ens_ID)$EnsID %in% sapply(strsplit(FitHiC_LNC_PCG_pairs, "-"), "[[", 1)])
9/27

#unique(filter(fpkm_allG, EnsName %in% CisLncExpected | EnsID %in% CisLncExpected)$EnsID[
#  filter(fpkm_allG, EnsName %in% CisLncExpected | EnsID %in% CisLncExpected)$EnsID %in% sapply(strsplit(FitHiC_LNC_PCG_pairshi, "-"), "[[", 1)])
#6/25 #many in the strong set, including AC002480.4 which is immune response linked

#any additional from HiCCUPs?
triali <- unique(filter(fpkm_allG, EnsID_merge %in% ControlCisLncs$Ens_ID)$EnsID[
  filter(fpkm_allG, EnsID_merge %in% ControlCisLncs$Ens_ID)$EnsID %in% sapply(strsplit(HiCCUP_LNC_PCG_pairs, "-"), "[[", 1)])
3/27

#all HiCCUPs found in FitHiC, total is 9 of 25 found
unique(trial, triali)
9/27

#what about all lncs:
expressedLncs <- unique(filter(fpkm_allG, grepl("fide|Lnc", GeneClassUpdate))$EnsID)

length(unique(expressedLncs[expressedLncs %in% sapply(strsplit(FitHiC_LNC_PCG_pairs, "-"), "[[", 1)|
                              expressedLncs %in% sapply(strsplit(HiCCUP_LNC_PCG_pairs, "-"), "[[", 1)]))

138/428
#v similar rate of HiC link to a PCG in cislnc and all lncs


#also MSTRG.12914-FOXL1, INKILN-CXCL8, RP11-HMGA2, DCBLD2(found in former work, putative link to GeneHancer)

#manual check for expected genes:
# foxl1 86,576,368-86,582,160
# 12913 86,681,935-86,696,468
# 12914 86,761,216-86,764,032

filter(caSMC_Zhao_Loops1_filt, 
       chr1 == "chr16",
       fragmentMid2 < 86764032, #less than 12914 R-hand extremity
       fragmentMid1 > 86576368) #more than foxl1 L-hand extremity



#### GO/KEGG/REACTOME for targets of 0-4hr lncs ####

library(clusterProfiler)
library(org.Hs.eg.db)

#first background, anything setting them apart from early up DEGs generally?
EarlyCoInduced <- filter(AllLNC_AllPCG_2d3d,
                         (EnsID %in% fpkm_allGDE_Upwithin_4$EnsID & 
                            EnsID.y %in% fpkm_allGDE_Upwithin_4$EnsID))

EarlyCoInduced$EnsID_merge.y <- gsub("\\.[0-9]*", "", EarlyCoInduced$EnsID.y)

fpkm_PCGDE_Upwithin_4 <- filter(fpkm_allGDE_Upwithin_4, grepl("protein_coding", EnsType))
fpkm_PCGDE_Upwithin_4$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE_Upwithin_4$EnsID)

#note high p thresh - 0 terms
CoUpEarly_UpEarly_GO <- enrichGO(gene          = unique(EarlyCoInduced$EnsID_merge.y),
                                 universe      = unique(fpkm_PCGDE_Upwithin_4$EnsID_merge.y),
                                 keyType       = "ENSEMBL",
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "all",    
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.1,
                                 qvalueCutoff  = 0.1,
                                 readable      = TRUE)

fpkm_PCGDE <- filter(fpkm_allGDE, grepl("protein_coding", EnsType))
fpkm_PCGDE$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE$EnsID)

#second background, anything setting apart from DEGs generally
CoUpEarly_DEGs_GO <- enrichGO(gene          = unique(EarlyCoInduced$EnsID_merge.y),
                                 universe      = unique(fpkm_PCGDE$EnsID_merge.y),
                                 keyType       = "ENSEMBL",
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "all",    
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.1,
                                 qvalueCutoff  = 0.1,
                                 readable      = TRUE)

#mostly cytokine/chemokine, zinc related, one or two loci only... not good to use
CoUpEarly_DEGs_GO_df <- as.data.frame(CoUpEarly_DEGs_GO)
#6x epigenetic gene regulation lower down (>0.05) (HMGA2/TET3/DNMT1)


#third background, anything setting apart from EGs generally
fpkm_PCG <- filter(fpkm_allG, grepl("protein_coding", EnsType))
fpkm_PCG$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCG$EnsID)

CoUpEarly_EGs_GO <- enrichGO(gene          = unique(EarlyCoInduced$EnsID_merge.y),
                              universe      = unique(fpkm_PCG$EnsID_merge.y),
                              keyType       = "ENSEMBL",
                              OrgDb         = org.Hs.eg.db,
                              ont           = "all",    
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.1,
                              qvalueCutoff  = 0.1,
                              readable      = TRUE)

CoUpEarly_EGs_GO_df <- as.data.frame(CoUpEarly_EGs_GO)
#as above, but good sig now, incl. on the epigenetic thing (3x genes only)


#first background, anything setting them apart from early Down DEGs generally?
EarlyCoRepressed <- filter(AllLNC_AllPCG_2d3d,
                         (EnsID %in% fpkm_allGDE_Downwithin_4$EnsID & 
                            EnsID.y %in% fpkm_allGDE_Downwithin_4$EnsID))

EarlyCoRepressed$EnsID_merge.y <- gsub("\\.[0-9]*", "", EarlyCoRepressed$EnsID.y)

fpkm_PCGDE_Downwithin_4 <- filter(fpkm_allGDE_Downwithin_4, grepl("protein_coding", EnsType))
fpkm_PCGDE_Downwithin_4$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE_Downwithin_4$EnsID)

#note high p thresh - 0 terms
CoDownEarly_DownEarly_GO <- enrichGO(gene          = unique(EarlyCoRepressed$EnsID_merge.y),
                                     universe      = unique(fpkm_PCGDE_Downwithin_4$EnsID_merge.y),
                                     keyType       = "ENSEMBL",
                                     OrgDb         = org.Hs.eg.db,
                                     ont           = "all",    
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.1,
                                     qvalueCutoff  = 0.1,
                                     readable      = TRUE)

#TFs, but coming from 2x loci mainly, collective lncRNA targeting? not rly despite good sig
CoDownEarly_DownEarly_GO_df <- as.data.frame(CoDownEarly_DownEarly_GO)

fpkm_PCGDE <- filter(fpkm_allGDE, grepl("protein_coding", EnsType))
fpkm_PCGDE$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE$EnsID)

#second background, anything setting apart from DEGs generally
CoDownEarly_DEGs_GO <- enrichGO(gene          = unique(EarlyCoRepressed$EnsID_merge.y),
                                universe      = unique(fpkm_PCGDE$EnsID_merge.y),
                                keyType       = "ENSEMBL",
                                OrgDb         = org.Hs.eg.db,
                                ont           = "all",    
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.1,
                                qvalueCutoff  = 0.1,
                                readable      = TRUE)

#picking up on the HOX genes again, developmental effect collectively being targeted? not rly despite good sig
CoDownEarly_DEGs_GO_df <- as.data.frame(CoDownEarly_DEGs_GO)


#third background, anything setting apart from EGs generally
fpkm_PCG <- filter(fpkm_allG, grepl("protein_coding", EnsType))
fpkm_PCG$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCG$EnsID)

CoDownEarly_EGs_GO <- enrichGO(gene          = unique(EarlyCoRepressed$EnsID_merge.y),
                               universe      = unique(fpkm_PCG$EnsID_merge.y),
                               keyType       = "ENSEMBL",
                               OrgDb         = org.Hs.eg.db,
                               ont           = "all",    
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.1,
                               qvalueCutoff  = 0.1,
                               readable      = TRUE)

CoDownEarly_EGs_GO_df <- as.data.frame(CoDownEarly_EGs_GO)
#as above


#joint testing of concordant early targets

ConcEarly_EGs_GO <- enrichGO(gene          = unique(c(unique(EarlyCoRepressed$EnsID_merge.y),
                                               unique(EarlyCoInduced$EnsID_merge.y))),
                               universe      = unique(fpkm_PCG$EnsID_merge.y),
                               keyType       = "ENSEMBL",
                               OrgDb         = org.Hs.eg.db,
                               ont           = "all",    
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.1,
                               qvalueCutoff  = 0.1,
                               readable      = TRUE)

ConcEarly_EGs_GO_df <- as.data.frame(ConcEarly_EGs_GO)

#for enrichR
EarlyCoInduced_con <- bitr(unique(EarlyCoInduced$EnsID_merge.y), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
write.csv(unique(EarlyCoInduced_con$ENTREZID), "EarlyCoInduced_entrez.csv", row.names = F)

EarlyCoRepressed_con <- bitr(unique(EarlyCoRepressed$EnsID_merge.y), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
write.csv(unique(EarlyCoRepressed_con$ENTREZID), "EarlyCoRepressed_entrez.csv", row.names = F)
#nothing seperately or together


#now all early targets of early lncs
#first background, anything setting them apart from early up DEGs generally?
EarlyTargets <- filter(AllLNC_AllPCG_2d3d,
                         (EnsID %in% fpkm_allGDE_within_4$EnsID & 
                            EnsID.y %in% fpkm_allGDE_within_4$EnsID))

EarlyTargets$EnsID_merge.y <- gsub("\\.[0-9]*", "", EarlyTargets$EnsID.y)

fpkm_PCGDE <- filter(fpkm_allGDE, grepl("protein_coding", EnsType))
fpkm_PCGDE$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCGDE$EnsID)

#note high p thresh - 0 terms
EarlyTargets_GO <- enrichGO(gene          = unique(EarlyTargets$EnsID_merge.y),
                                 universe      = unique(fpkm_PCGDE$EnsID_merge.y),
                                 keyType       = "ENSEMBL",
                                 OrgDb         = org.Hs.eg.db,
                                 ont           = "all",    
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.1,
                                 qvalueCutoff  = 0.1,
                                 readable      = TRUE)


#no additional info
EarlyTargets_GO_df <- as.data.frame(EarlyTargets_GO)

#2nd background, anything setting apart from EGs generally
fpkm_PCG <- filter(fpkm_allG, grepl("protein_coding", EnsType))
fpkm_PCG$EnsID_merge.y <- gsub("\\.[0-9]*", "", fpkm_PCG$EnsID)

EarlyTargets_EGs_GO <- enrichGO(gene          = unique(EarlyTargets$EnsID_merge.y),
                             universe      = unique(fpkm_PCG$EnsID_merge.y),
                             keyType       = "ENSEMBL",
                             OrgDb         = org.Hs.eg.db,
                             ont           = "all",    
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.1,
                             qvalueCutoff  = 0.1,
                             readable      = TRUE)

EarlyTargets_df <- as.data.frame(EarlyTargets_EGs_GO)
#as above, no additional info

#list of names for STRING: 

write.csv(unique(EarlyTargets$EnsName.y), "EarlyTargets_STRING.csv", row.names = F)
#CXCL/CCL/IL6 in core, points out to HOX cluster as well as HMGA2/DNMT1/TET3/MT2A quartet, other direction to the MTs
#FOXL1 unconnected to these genes in STRING
