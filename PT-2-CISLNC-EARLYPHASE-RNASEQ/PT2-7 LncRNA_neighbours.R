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

#isolate the lncRNA tx
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

