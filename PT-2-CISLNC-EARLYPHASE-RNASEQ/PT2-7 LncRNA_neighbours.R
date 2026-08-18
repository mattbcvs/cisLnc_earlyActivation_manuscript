library(dplyr)

#### import basic tables ####

fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026filt.csv", header = T)
length(unique(fpkm_allG$EnsID))#12740
table(unique(fpkm_allG[,c(2,55)])$V55)#597 lncs

#CAGE annotation info (to pinpoint TSS for novel lncs)
Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime_2026.csv", header = T)
Enhancer_lociII_DEsig_Enh <- Enhancer_lociII
length(unique(Enhancer_lociII_DEsig_Enh$EnsID))#597
length(unique(Enhancer_lociII_DEsig_Enh$MSTRG_Tx_ID))#1575 as expected
length(unique(filter(Enhancer_lociII_DEsig_Enh, !is.na(DiffExprs))$EnsID))#221

#get co-ords based on FANTOM TSS
Enhancer_lociII_DEsig_Enh$Enhancer_Coords <- paste(Enhancer_lociII_DEsig_Enh$chr, 
                                                   Enhancer_lociII_DEsig_Enh$Enh_Start, 
                                                   Enhancer_lociII_DEsig_Enh$Enh_Stop, sep = ",")

#### correct TSS for lncRNAs with a matched CAGE ####

#Get TSS from FANTOM and TSS from GENCODE in same column
Enhancer_lociII_DEsig_Enh$TSS_FANTOM_GENCODE <- Enhancer_lociII_DEsig_Enh$BestStart

#obtain TSS co-ords using these CAGE sites or just 5' limit from GENCODE/Stringtie transcripts for others
trial <- fpkm_allG
trial$Tx_start <- as.numeric(sapply(strsplit(sapply(strsplit(trial$Tx_Locus, ":"), "[[", 2), "-"), "[[", 1))
trial$Tx_stop <- as.numeric(gsub(" [+-]", "", sapply(strsplit(sapply(strsplit(trial$Tx_Locus, ":"), "[[", 2), "-"), "[[", 2)))

trial <- unique(trial[,c(2,5,59:60,8,47:48)])
#alternate/better TSS CAGE from FANTOM for these lncRNAs:
triali <- unique(filter(Enhancer_lociII_DEsig_Enh, CAGEvalidity == "Valid CAGE")[,c(2,42)])

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
length(unique(allGB$MSTRG_Tx_ID))#42511 TSS total (multiple TSS per gene now - this was a previous issue)


#### assign neighbours to lncRNAs ####

#isolate the 1575 lncRNA tx
allGB_LNCS <- filter(allGB, EnsID %in% filter(fpkm_allG, grepl("Lnc|fide", GeneClassUpdate))$EnsID)
length(unique(allGB_LNCS$EnsID))

allGB_LNCS$TSS <- allGB_LNCS$Tx_start
allGB_LNCS$TSS[allGB_LNCS$str == "-"] <- allGB_LNCS$Tx_stop[allGB_LNCS$str == "-"]

trial <- allGB_LNCS
trial <- split(trial, trial$EnsID)

#isolate PCGs
allGB_PCGs <- filter(allGB, EnsID %in% filter(fpkm_allG, EnsType == "protein_coding" & grepl("TF|CC|coding", V55)
                                              )$EnsID)

#for each lnc, subset the PCG table for close PCGs
triali <- lapply(trial, function(z){
  filter(allGB_PCGs[allGB_PCGs$str == "+",], 
         #for all + strand PCGs, find lncs on same chr
         (chr == unique(z$chr) & 
            #with a TSS (Tx_start for +ve) within given dist of a lncRNA TSS
            Tx_start > min(z$TSS)-1000000 & 
            Tx_start < max(z$TSS)+1000000)
         )
})

#equivalent for -ve strand
trialii <- lapply(trial, function(z){
  filter(allGB_PCGs[allGB_PCGs$str == "-",], 
         (chr == unique(z$chr) & 
            Tx_stop>min(z$TSS)-1000000 & 
            Tx_stop<max(z$TSS)+1000000)
         )
})

triali <- unique(bind_rows(triali, .id = "lnc_id"))
triali <- unique(triali)

trialii <- unique(bind_rows(trialii, .id = "lnc_id"))
trialii <- unique(trialii)

AllLNC_AllPCG_1 <- rbind(triali, trialii)

trial <- unique(merge(AllLNC_AllPCG_1[,c(1,3)], fpkm_allG[,c(2,3,58)], by.x = "lnc_id", by.y = "EnsID"))
trial <- unique(merge(trial, fpkm_allG[,c(2,3,58)], by = "EnsID"))

#2D pairs with TSS within 1Mbp table:
AllLNC_AllPCG_1 <- trial#8438 pairs

AllLNC_AllPCG_1$pairs <- paste(AllLNC_AllPCG_1$lnc_id, AllLNC_AllPCG_1$EnsID, sep="-")

#get format:
colnames(AllLNC_AllPCG_1)[1:2] <- c("EnsID.y", "EnsID")

AllLNC_AllPCG_1 <- AllLNC_AllPCG_1[,c(2:4,1,5:7)]

AllLNC_AllPCG_1$pairs <- paste(AllLNC_AllPCG_1$EnsID, AllLNC_AllPCG_1$EnsID.y, sep="-")

#save after adding some annotation info, distance, pair orientation etc
#write.csv(AllLNC_AllPCG_1, "\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_1_Q3.csv", row.names =  F)


#### extra annotation for 2d pairs ####

#plot distance between and spearman
#add in TSS and distance between, allTSS:
allGB$TSS_FANTOM_GENCODE <- allGB$Tx_start
allGB$TSS_FANTOM_GENCODE[allGB$str == "-"] <- allGB$Tx_stop[allGB$str == "-"]

#n.b. multiple TSS per gene
triali <- unique(merge(AllLNC_AllPCG_1, allGB[,c(2,8)], by.x = "EnsID", by.y = "EnsID"))
triali <- unique(merge(triali, allGB[,c(2,8)], by.x = "EnsID.y", by.y = "EnsID"))
triali$AbsDistLnc_PCG <- abs(triali$TSS_FANTOM_GENCODE.x - triali$TSS_FANTOM_GENCODE.y)/1000
triali$DistLnc_PCG <- (triali$TSS_FANTOM_GENCODE.x - triali$TSS_FANTOM_GENCODE.y)/1000

#shortest distance per pair:
trialii <- split(triali, triali$pairs)

trialii[[27]]
trialiii <- lapply(trialii, function(x){
  unique(x[x$AbsDistLnc_PCG == min((x$AbsDistLnc_PCG)),-c(8:9)])
})

trial <- bind_rows(trialiii)

#lnc-PCG relationship, need strand
colnames(fpkm_allG)
trial <- unique(merge(trial, fpkm_allG[,c(2,8)], by = "EnsID"))
trial <- unique(merge(trial, fpkm_allG[,c(2,8)], by.x = "EnsID.y", by.y = "EnsID"))

#need to import overlaps 
#(table contains genes where the maximum overlap between a pair of major isoforms from diff genes is 10% or above)
FPKM_CQV_OVERLAP_fpkm <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SameoverlapsG_lncs.csv")
FPKM_CQV_OVERLAP_fpkm <- merge(FPKM_CQV_OVERLAP_fpkm, fpkm_allG[,1:2], by.x = "MSTRG_ID.y", by.y = "MSTRG_ID")
FPKM_CQV_OVERLAP_fpkm$pairs <- paste(FPKM_CQV_OVERLAP_fpkm$EnsID.x, FPKM_CQV_OVERLAP_fpkm$EnsID, sep = "-")

trial$LNC_PCG_Type <- NA
#same strand overlaps - none!
trial$LNC_PCG_Type[trial$pairs %in% FPKM_CQV_OVERLAP_fpkm$pairs & trial$str.x == trial$str.y] <- "Sense Overlap"
#opposite strand overlaps
trial$LNC_PCG_Type[trial$pairs %in% FPKM_CQV_OVERLAP_fpkm$pairs & 
                     !trial$str.x == trial$str.y] <- "Antisense Overlap"
#intergenic - no overlaps
trial$LNC_PCG_Type[!trial$pairs %in% FPKM_CQV_OVERLAP_fpkm$pairs] <- "Intergenic"
#where the lncRNA is enhancer annotated (leave as sep column)
#trial$LNC_PCG_Type[!trial$pairs %in% triali$pairs & trial$EnsID %in% filter(fpkm_allGDE,  grepl("ELnc", GeneClassUpdate))$EnsID] <- "Intergenic Enhancer"

#transcription appears divergent from same promoter region:
trial$LNC_PCG_Type[!trial$str.x == trial$str.y & trial$AbsDistLnc_PCG <3] <- "Bidirectional"
table(trial$LNC_PCG_Type)#64 AS, 274 bidir, 8097 intergenic, 3 sense overlap

colnames(trial)[12] <- "LncRNA-PCG Relationship"

AllLNC_AllPCG_info <- trial

length(unique(AllLNC_AllPCG_info$EnsID)) #579 lncs near
length(unique(AllLNC_AllPCG_info$EnsID.y)) #5351 pcgs

#write.csv(AllLNC_AllPCG_info, "AllLNC_AllPCG_info_2026.csv", row.names = F)


#### old code - wrong/imprecise gene co-ords ####
#table of all lncRNAs + CAGE sites if available + TSS based on CAGE if available:
Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime.csv", header = T)
Enhancer_lociII_DEsig_Enh <- Enhancer_lociII
length(unique(Enhancer_lociII_DEsig_Enh$EnsID))#selecting all lncs = 558 genes, if just enhancer = 77 (7/2021)
#get co-ords based on FANTOM TSS
Enhancer_lociII_DEsig_Enh$Enhancer_Coords <- paste(Enhancer_lociII_DEsig_Enh$chr, 
                                                   Enhancer_lociII_DEsig_Enh$Enh_Start, 
                                                   Enhancer_lociII_DEsig_Enh$Enh_Stop, sep = ",")

#some stats (pre-AS artefact filtering):
#558 lncs expressed
length(unique(Enhancer_lociII$EnsID))
#234 DE lncs over 24 hours:
length(unique(filter(Enhancer_lociII, !is.na(DiffExprs))$EnsID))
#77 enhancer lncs expressed
length(unique(filter(Enhancer_lociII, !is.na(EnhancerVerdict))$EnsID))
#41 enhancer lncs expressed and DE
length(unique(filter(Enhancer_lociII, !is.na(DiffExprs), !is.na(EnhancerVerdict))$EnsID))

fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGClassUpdate.csv", header = T)
length(unique(fpkm_allG$EnsID))#11815 genes


fpkm_PCG_hSVSMC <- filter(fpkm_allG, `V57` == "Protein coding" | grepl("TF", V57) | grepl("CC", V57) )
fpkm_PCG_hSVSMC <- unique(fpkm_PCG_hSVSMC[,c(2,5,7:10,27,29,31,33,35,60)])
length(unique(fpkm_PCG_hSVSMC$EnsID))#10213 PC genes
length(unique(fpkm_allG$EnsID))#11815 genes
fpkm_PLARG <- filter(fpkm_allG, `V57` == "Bona fide lncRNA" | grepl("ELnc", V57) | grepl("VLnc", V57))
length(unique(fpkm_PLARG$EnsID))#558 lnc genes

trial <- filter(Enhancer_lociII_DEsig_Enh, EnsID %in% fpkm_PLARG$EnsID)
trial <- unique(trial[,c(19,17,18,43)])
trial <- split(trial, trial$Enhancer_Coords)


triali <- lapply(trial, function(z){
  filter(fpkm_PCG_hSVSMC[fpkm_PCG_hSVSMC$str == "+",], 
         (chr == z$chr & start>z$Enh_Start-250000 & start<z$Enh_Stop+250000))
})

trialii <- lapply(trial, function(z){
  filter(fpkm_PCG_hSVSMC[fpkm_PCG_hSVSMC$str == "-",], 
         (chr == z$chr & stop>z$Enh_Start-250000 & stop<z$Enh_Stop+250000))
})

triali <- unique(bind_rows(triali, .id = "Enhancer_Coords"))
triali <- unique(triali)

trialii <- unique(bind_rows(trialii, .id = "Enhancer_Coords"))
trialii <- unique(trialii)

ProximalDEBO <- rbind(triali, trialii)

trial <- merge(Enhancer_lociII_DEsig_Enh, ProximalDEBO, by = "Enhancer_Coords", all.x = T)

AllLNC_AllPCG_1 <- trial
colnames(AllLNC_AllPCG_1)
#lncRNA name missing for some reason:
AllLNC_AllPCG_1 <- unique(merge(fpkm_allG[,c(2,5)], AllLNC_AllPCG_1, by.x = "EnsID", by.y = "EnsID.x", all.y = T))
length(unique(AllLNC_AllPCG_1$EnsID))
length(unique(filter(AllLNC_AllPCG_1, !is.na(EnsID.y))$EnsID))#497 in range of a PCG
length(unique(AllLNC_AllPCG_1$EnsID.y))#1962 partners

#simplify:
#genes with any elnc tx should all be called elnc:
trial <- AllLNC_AllPCG_1
elncs <- unique(filter(trial, EnhancerVerdict == "Enhancer")$EnsID)
trial$EnhancerVerdict <- "Non-enhancer"
trial$EnhancerVerdict[trial$EnsID %in% elncs] <- "Enhancer"
#excess columns off
triali <- unique(trial[,c(1,2,43,45,46,56)])

#2D pairs with TSS within 250kbp table:
AllLNC_AllPCG_1 <- triali

#write.csv(AllLNC_AllPCG_1, "\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_1.csv", row.names =  F)

#n.b. some final edits needed after importing the above object:
#double check no AS artefacts, no PLAR only "PCGs", must have ENCODE annotation
AllLNC_AllPCG_1$pairs <- paste(AllLNC_AllPCG_1$EnsID, AllLNC_AllPCG_1$EnsID.y, sep="-")

AllLNC_AllPCG_1 <- filter(AllLNC_AllPCG_1, EnsID %in% fpkm_allG$EnsID, EnsID.y %in% fpkm_allG$EnsID) #all good - none removed
AllLNC_AllPCG_1 <- filter(AllLNC_AllPCG_1, EnsID.y %in% filter(fpkm_allG, EnsType == "protein_coding")$EnsID) #removes about 5
#this last step is bit more stringent than previous PCG analysis where just had to be called PCG by PLAR (not PLAR and Ens)
#1684 pairings of lnc-PCG expressed near to each other
