#Quality control - unstranded data

#### key info, FPKM table ####
samplenames <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/rsem_4timepoints(NovelTx) calc-ci/", pattern = "*genes.results", full.names = TRUE)[1:16]
#samplenames <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/rsem_4timepoints(NovelTx)/", pattern = "*genes.results", full.names = TRUE)[1:16]

#sample info from Amira
actualnames <- c("1) Pt151 - 0h 2) Pt151 - 4h 3) Pt151 - 8h 4) Pt151 - 24h 
                 5) Pt157 - 0h 6) Pt157 - 4h 7) Pt157 - 8h 8) Pt157 - 24h 
                 9) Pt134 - 0h 10) Pt134 - 4h 11) Pt134 - 8h 12) Pt134 - 24h 
                 13) Pt2279 - 0h14) Pt2279 - 4h 15) Pt2279 - 8h 16) Pt2279 - 24h")
trial <- unlist(strsplit(actualnames, "[0-9]) "))[2:17]
actualnames <- strsplit(trial, "h")
actualnames <- sapply(actualnames, "[[", 1)
actualnames <- actualnames[c(seq(1,16,4), seq(1,16,4)+1, seq(1,16,4)+2, seq(1,16,4)+3)]

samplenames_order <- samplenames[c(1,12,16,5,
                                   9,13,2,6,
                                   10,14,3,7,
                                   11,15,4,8)]
filenames <- samplenames_order
rm(samplenames_order)
#generate counts table for downstream use
genes<-read.table(filenames[1],header=TRUE,sep="\t",stringsAsFactors = FALSE)[,1]
fpkm <-do.call(cbind,lapply(filenames,function(fn)read.table(fn,header=TRUE,sep="\t",stringsAsFactors = FALSE)[,7]))
fpkm <- data.frame(genes,fpkm,stringsAsFactors = FALSE)
colnames(fpkm)<-c("ENSEMBL",actualnames)

fpkm_list_ctrl <- as.list(as.data.frame(t(fpkm[,2:5])))
fpkm_list_pd <- as.list(as.data.frame(t(fpkm[,6:9])))
fpkm_list_il <- as.list(as.data.frame(t(fpkm[,10:13])))
fpkm_list_bo <- as.list(as.data.frame(t(fpkm[,14:17])))

fpkm_mean_ctrl <- sapply(fpkm_list_ctrl, mean)
fpkm_mean_pd <- sapply(fpkm_list_pd, mean)
fpkm_mean_il <- sapply(fpkm_list_il, mean)
fpkm_mean_bo <- sapply(fpkm_list_bo, mean)

std <- function(x) sd(x)/sqrt(length(x))

fpkm_se_basal <- sapply(fpkm_list_ctrl, std)
fpkm_se_chol <- sapply(fpkm_list_pd, std)
fpkm_se_mig <- sapply(fpkm_list_il, std)
fpkm_se_pdgf <- sapply(fpkm_list_bo, std)

fpkm_mean_treatment <- data.frame("Hour0_meanFPKM" = fpkm_mean_ctrl, "Hour0_seFPKM" = fpkm_se_basal,
                                  "Hour4_meanFPKM" = fpkm_mean_pd, "Hour4_seFPKM" = fpkm_se_chol,
                                  "Hour8_meanFPKM" = fpkm_mean_il, "Hour8_seFPKM" = fpkm_se_mig,
                                  "Hour24_meanFPKM" = fpkm_mean_bo, "Hour24_seFPKM" = fpkm_se_pdgf)

trial <- as.list(as.data.frame(t(fpkm_mean_treatment)))
fpkm_max_treatment <- as.numeric(sapply(trial, max))

fpkm <- cbind(fpkm, fpkm_mean_treatment, fpkm_max_treatment)

#### FPKM CQV ####

samplenames <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/rsem_4timepoints(NovelTx) calc-ci/", pattern = "*genes.results", full.names = TRUE)[1:16]

#sample info from Amira
actualnames <- c("1) Pt151 - 0h 2) Pt151 - 4h 3) Pt151 - 8h 4) Pt151 - 24h 
                 5) Pt157 - 0h 6) Pt157 - 4h 7) Pt157 - 8h 8) Pt157 - 24h 
                 9) Pt134 - 0h 10) Pt134 - 4h 11) Pt134 - 8h 12) Pt134 - 24h 
                 13) Pt2279 - 0h14) Pt2279 - 4h 15) Pt2279 - 8h 16) Pt2279 - 24h")
trial <- unlist(strsplit(actualnames, "[0-9]) "))[2:17]
actualnames <- strsplit(trial, "h")
actualnames <- sapply(actualnames, "[[", 1)
actualnames <- actualnames[c(seq(1,16,4), seq(1,16,4)+1, seq(1,16,4)+2, seq(1,16,4)+3)]

samplenames_order <- samplenames[c(1,12,16,5,
                                   9,13,2,6,
                                   10,14,3,7,
                                   11,15,4,8)]
filenames <- samplenames_order
rm(samplenames_order)
#generate counts table for downstream use
genes<-read.table(filenames[1],header=TRUE,sep="\t", stringsAsFactors = FALSE)[,1]

colnames(read.table(filenames[1],header=TRUE,sep="\t", nrows = 10, stringsAsFactors = FALSE))

fpkmCQV <-do.call(cbind,lapply(filenames,function(fn)read.table(fn,header=TRUE,sep="\t",stringsAsFactors = FALSE)[,17]))
fpkmCQV <- data.frame(genes,fpkmCQV,stringsAsFactors = FALSE)
colnames(fpkmCQV)<-c("ENSEMBL",actualnames)

fpkmCQV_list_ctrl <- as.list(as.data.frame(t(fpkmCQV[,2:5])))
fpkmCQV_list_pd <- as.list(as.data.frame(t(fpkmCQV[,6:9])))
fpkmCQV_list_il <- as.list(as.data.frame(t(fpkmCQV[,10:13])))
fpkmCQV_list_bo <- as.list(as.data.frame(t(fpkmCQV[,14:17])))

fpkmCQV_mean_ctrl <- sapply(fpkmCQV_list_ctrl, mean)
fpkmCQV_mean_pd <- sapply(fpkmCQV_list_pd, mean)
fpkmCQV_mean_il <- sapply(fpkmCQV_list_il, mean)
fpkmCQV_mean_bo <- sapply(fpkmCQV_list_bo, mean)

std <- function(x) sd(x)/sqrt(length(x))

fpkmCQV_se_basal <- sapply(fpkmCQV_list_ctrl, std)
fpkmCQV_se_chol <- sapply(fpkmCQV_list_pd, std)
fpkmCQV_se_mig <- sapply(fpkmCQV_list_il, std)
fpkmCQV_se_pdgf <- sapply(fpkmCQV_list_bo, std)

fpkmCQV_mean_treatment <- data.frame("Hour0_meanfpkmCQV" = fpkmCQV_mean_ctrl, "Hour0_sefpkmCQV" = fpkmCQV_se_basal,
                                     "Hour4_meanfpkmCQV" = fpkmCQV_mean_pd, "Hour4_sefpkmCQV" = fpkmCQV_se_chol,
                                     "Hour8_meanfpkmCQV" = fpkmCQV_mean_il, "Hour8_sefpkmCQV" = fpkmCQV_se_mig,
                                     "Hour24_meanfpkmCQV" = fpkmCQV_mean_bo, "Hour24_sefpkmCQV" = fpkmCQV_se_pdgf)

trial <- as.list(as.data.frame(t(fpkmCQV_mean_treatment)))
fpkmCQV_max_treatment <- as.numeric(sapply(trial, max))

fpkmCQV <- cbind(fpkmCQV, fpkmCQV_mean_treatment, fpkmCQV_max_treatment)


#### Overlapping extent ####

#overlaps:
#consider only major isoforms >10% of genes output:
samplenames <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/rsem_4timepoints(NovelTx)_NonMergedEnsID/", 
                          pattern = "*isoforms.results", full.names = TRUE)[1:16]

samplenames_order <- samplenames[c(1,12,16,5,
                                   9,13,2,6,
                                   10,14,3,7,
                                   11,15,4,8)]
iso_filenames <- samplenames_order

actualnames <- c("1) Pt151 - 0h 2) Pt151 - 4h 3) Pt151 - 8h 4) Pt151 - 24h 
                 5) Pt157 - 0h 6) Pt157 - 4h 7) Pt157 - 8h 8) Pt157 - 24h 
                 9) Pt134 - 0h 10) Pt134 - 4h 11) Pt134 - 8h 12) Pt134 - 24h 
                 13) Pt2279 - 0h14) Pt2279 - 4h 15) Pt2279 - 8h 16) Pt2279 - 24h")
trial <- unlist(strsplit(actualnames, "[0-9]) "))[2:17]
actualnames <- strsplit(trial, "h")
actualnames <- sapply(actualnames, "[[", 1)

sampleInfo_Timecourse <- data.frame("SampleNo" = c(1,10:16,2:9),
                                    "IsoFileName" = samplenames,
                                    "ActualNames" = actualnames[c(1,10:16,2:9)])

#sample names re-ordered into the same
sampleInfo_Timecourse <- sampleInfo_Timecourse[c(1,12,16,5,
                                                 9,13,2,6,
                                                 10,14,3,7,
                                                 11,15,4,8),]

#available info per transcript
colnames(read.delim(sampleInfo_Timecourse$IsoFileName[1], nrow = 10))

#iso table with per sample percentage:
isoforms<-read.table(sampleInfo_Timecourse$IsoFileName[1],
                     header=TRUE,
                     sep="\t",
                     stringsAsFactors = FALSE,
                     colClasses = c("character", "character", rep("NULL", 6)))

#iso_pct
trial <- lapply(sampleInfo_Timecourse$IsoFileName, function(x){
  read.table(x, header = T, stringsAsFactors = F, 
             colClasses = c(rep("NULL", 7), "numeric")) %>% as_tibble()
})

trialii <- do.call(cbind, trial)

iso_pct <- data.frame(isoforms, trialii, stringsAsFactors = FALSE)

#mean percentage across all samples (ignore if 0)
trial <- as.matrix(iso_pct[,3:18])
trial[trial ==0] <- NA
iso_pct$meanTxPerc <-  rowMeans(trial, na.rm = T)
iso_pct$meanTxPerc[rowSums(iso_pct[,3:18]) == 0] <- 0

#sort through for major isoforms, rough approximation here:
trial <- split(iso_pct[c(1,2,19)], iso_pct$gene_id)
iso_dominant <- sapply(trial, function(x){
  if(sum(x$meanTxPerc >(1/3*100)) >0) { 
    x$transcript_id[which(x$meanTxPerc >(1/3*100))] 
  } #if any are over 33 then take them
  else if(sum(x$meanTxPerc >(1/4*100)) >0) { 
    x$transcript_id[which(x$meanTxPerc >(1/4*100))] 
  } #if any are over 25 then take them
  else if(sum(x$meanTxPerc >(1/10*100)) >0) { 
    x$transcript_id[which(x$meanTxPerc >(1/10*100))] 
  } #if any are over 10 then take them
  else x$transcript_id[which(x$meanTxPerc >(1/100*100))]
})
iso_pct$iso_dominant <- NA
iso_pct$iso_dominant[iso_pct$transcript_id %in% unlist(iso_dominant)] <- "Dominant_iso"

iso_pct_10 <- filter(iso_pct, iso_dominant == "Dominant_iso")
#some weird cases like NR2F2-AS1 to review in future, not consistent across reps per condition (maybe remove large ICV tx?)

#shouldn't be losing any genes:
length(unique(iso_pct$gene_id))
length(unique(iso_pct_10$gene_id))
#limitation: a v. highly expressed PCG may have a low expressed tx that overlaps a low expressed lnc

#tx/gene tables built in Timecourse4Timepointsmapping - now using the above iso_pct_10 to filter
#fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGClassUpdate.csv", header = T)
#fpkm_allTx <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allTx.csv", header = T)
length(unique(fpkm_allTx$EnsID))
#fpkm_allTx <- filter(fpkm_allTx, MSTRG_Tx_ID %in% iso_pct_10$transcript_id, EnsID %in% fpkm_allG$EnsID)
length(unique(fpkm_allTx$EnsID)) #11771 lose a few genes somewhere, some fpkm <1 genes in fpkm_allTx probs

Tx_coordsChr <- sapply(sapply(fpkm_allTx$Tx_Locus, strsplit, "\\:"), "[[", 1)
Tx_coords <- sapply(sapply(sapply(fpkm_allTx$Tx_Locus, strsplit, "\\:"), "[[", 2), strsplit, "[[", 2)
Tx_coordsStart <- sapply(sapply(Tx_coords, strsplit, "\\-"), "[[", 1)
Tx_coordsStop <- sapply(sapply(Tx_coords, strsplit, "\\-"), "[[", 2)
Tx_coordsStop <- gsub(" ", "", Tx_coordsStop)
Tx_coordsStop <- gsub("\\+", "", Tx_coordsStop)

fpkm_allTx$Txchr <- Tx_coordsChr
fpkm_allTx$TxStart <- Tx_coordsStart
fpkm_allTx$TxStop <- Tx_coordsStop

Txcoords <- makeGRangesFromDataFrame(filter(fpkm_allTx)[,c(47,55:57)],
                                     start.field = "TxStart", 
                                     end.field = "TxStop", 
                                     keep.extra.columns = T, 
                                     seqnames.field = "Txchr", 
                                     ignore.strand = T)

#find overlaps to other tx (then filter to just those from other genes)
SameOverlapsindex <- findOverlaps(Txcoords, Txcoords)
Sameoverlaps <- data.frame("Tx1" = Txcoords$MSTRG_Tx_ID[queryHits(SameOverlapsindex)],
                           "Tx2" = Txcoords$MSTRG_Tx_ID[subjectHits(SameOverlapsindex)])

#bp of overlap per lncRNA transcript
overlaps <- pintersect(Txcoords[queryHits(SameOverlapsindex)], Txcoords[subjectHits(SameOverlapsindex)])
Sameoverlaps$percentOverlap<- width(overlaps)/width(Txcoords[queryHits(SameOverlapsindex)])*100
#additional info
Sameoverlaps <- merge(Sameoverlaps, fpkm_allG[,c(1,2,3,33,47)], by.x = "Tx1", by.y = "MSTRG_Tx_ID")
Sameoverlaps <- merge(Sameoverlaps, fpkm_allG[,c(1,2,3,33,47)], by.x = "Tx2", by.y = "MSTRG_Tx_ID")
Sameoverlaps$pairs <- paste(Sameoverlaps$EnsID.x, Sameoverlaps$EnsID.y, sep="-")

#remove same gene overlaps:
Sameoverlaps <- filter(Sameoverlaps, !EnsID.x == EnsID.y)
length(unique(Sameoverlaps$Tx1))#3035 tx that overlap another gene tx to some degree (10% approach)
length(unique(Sameoverlaps$EnsID.x))#2124 genes expressed where major iso overlaps another gene's major iso (10% approach)

#useful ballpark for how much one gene outweighs another:
Sameoverlaps$AbundanceMismatch <- Sameoverlaps$fpkm_max_treatment.x/Sameoverlaps$fpkm_max_treatment.y

#summary per pair:
#max percent overlap:
trial <- split(Sameoverlaps, Sameoverlaps$pairs)

triali <- lapply(trial, function(x){
  x$MaxPercentOverlap <- max(x$percentOverlap)
  return(x)
})

triali <- bind_rows(triali)

#by gene:
SameoverlapsG <- unique(triali[,c(5,4,6,7,8,10,11,13,14)])

SameoverlapsG_lncs <- unique(filter(SameoverlapsG, EnsID.x %in% unique(filter(fpkm_allG, V55 == "Bona fide lncRNA")$EnsID)))

#add in any IGV fails from previous checking (some changes but manual checking so fine - check a sample of pass/fails):
FPKM_CQV_OVERLAP_fpkm <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SameoverlapsG_lncs.csv")

SameoverlapsG_lncs <- merge(SameoverlapsG_lncs, FPKM_CQV_OVERLAP_fpkm[,c(1,10)], by.x = "EnsID.x", by.y = "EnsID.x", all.x = T)

#manual checks on all:
#write.csv(SameoverlapsG_lncs, "SameoverlapsG_lncs_151225.csv", row.names = F)


missingFromPrior <- filter(SameoverlapsG_lncs, is.na(IGV))

length(unique(SameoverlapsG_lncs$EnsID.x))
length(unique(missingFromPrior$EnsID.x))

summary(FPKM_CQV_OVERLAP_fpkm$MaxPercentOverlap)
summary(SameoverlapsG_lncs$MaxPercentOverlap)

#334 missing... some maybe don't need checking ...
#one or other gene not expressed?
fpkm_allG_old <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGClassUpdate.csv", header = T)

#203 because one or other not expressed in previous
missingDueToNonExpression <- filter(missingFromPrior, !EnsID.x %in% fpkm_allG_old$EnsID | !EnsID.y %in% fpkm_allG_old$EnsID)


#all lncRNAs checked manually and passed/failed via visual inspeection of IGV locus
# require clear signal over exons that can only be explained by lncRNA annotaiton
# borderline called as fails

# expected to be about 4-5hrs of checking to get through... did ~40 in 30mins

#fail rate justifies, out of 40, 13 fails found, unacceptably high

# additional lncRNAs from
# a) lowering threshold to FPKM >0.8
# b) splitting up ens genes
# c) max pct overlap was previously set at 10%

#### Combine into table of data for assessing overlap artefacts ####

FPKM_CQV_OVERLAP <- data.frame("gene" = fpkm$ENSEMBL,
                               "fpkm_max" = fpkm$fpkm_max_treatment#,
                               #"fpkmCQV_max" = fpkmCQV$fpkmCQV_max_treatment
                               )
trial <- unique(merge(FPKM_CQV_OVERLAP, SameoverlapsG, by.x = "gene", by.y = "MSTRG_ID.x"))
FPKM_CQV_OVERLAP <- unique(merge(trial, fpkm_allG[,c(1,2)], by.x = "gene", by.y = "MSTRG_ID"), all.x = T)
#can now assess all expressed genes in the timecourse for their expression in hour 0
length(unique(FPKM_CQV_OVERLAP$EnsID))#11815 genes
length(unique(fpkm_allG$EnsID))

write.csv(FPKM_CQV_OVERLAP, "FPKM_CQV_OVERLAP_26.csv")

#take to manually check lncRNAs:
FPKM_CQV_OVERLAP_lncs <- unique(filter(FPKM_CQV_OVERLAP, EnsID %in% unique(filter(fpkm_allG, V55 == "Bona fide lncRNA")$EnsID)))



