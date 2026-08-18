library(dplyr)
#library(clusterProfiler)
#library(org.Hs.eg.db)

##### old code ####
##### proof that it works ####
library(ReMapEnrich)

#bit dodgy and not maintained... but fairly basic overall
#https://remap-cisreg.github.io/ReMapEnrich/vignettes/basic_use.html
#https://remap-cisreg.github.io/ReMapEnrich/vignettes/advanced_use.html

#tutorial shows how to find TFs that overlap sox2 binding site on chr22

#query bed file, here all SOX2 peaks on chr22
query <- bedToGranges(system.file("extdata",
                                  "ReMap_nrPeaks_public_chr22_SOX2.bed",
                                  package = "ReMapEnrich"))
#catalog, all ReMap peaks on chr22
catalog <- bedToGranges(system.file("extdata",
                                    "ReMap_nrPeaks_public_chr22.bed",
                                    package = "ReMapEnrich"))
query
catalog

#this part is a black box, but some sort of shuffling (of the catalog?) to get a background rate
enrichment.df <- enrichment(query, catalog, byChrom = TRUE)

View(enrichment.df)

enrichmentBarPlot(enrichment.df, sigDisplayQuantile = 0.5, top = 20, aRisk = 0.00001)

#this one doesn't work anymore, could reconstruct from table easily enoug tho
#enrichmentVolcanoPlot(na.omit(enrichment.df), sigDisplayQuantile = 0.9, aRisk = 0.00001)

enrichmentDotPlot(enrichment.df)


#you can use the whole catalog, but they only have 2018 as the latest:
# Create a local directory for the tutorial
demo.dir <- "~/ReMapEnrich_demo"
dir.create(demo.dir, showWarnings = FALSE, recursive = TRUE)

# Use the function DowloadRemapCatalog
remapCatalog2018hg38 <- downloadRemapCatalog(demo.dir)

#fix on the github doesn't seem to improve:
downloadRemapCatalog <- function(targetDir,
                                 fileName = "", 
                                 version = "2018",
                                 assembly = "hg38",
                                 force = FALSE,
                                 store = TRUE) {
  
  if (version != "2018" && version != "2015" && version != "2020" && version != "2022") {
    message("Invalid version of catalog, choose between 2015, 2018, 2020 and 2022.")
    stop()
  }
  if (assembly != "hg38" && assembly != "hg19" && assembly != "dm6" && assembly != "tair10" && assembly != "mm10") {
    message("Invalid assembly, choose between hg19 and hg38.")
    stop()
  }
  size <- "Multi"
  
  if (version == "2015") {
    size <- "0.5"
  }
  
  
  url <- "https://remap.univ-amu.fr/storage/"
  
  if (version == "2018" && assembly == "hg38") {
    url <- paste(url, "remap2018/hg38/MACS/remap2018_nr_macs2_hg38_v1_2.bed.gz", sep = "")
  }
  else if (version == "2018" && assembly == "hg19") {
    url <- paste(url, "remap2018/hg19/MACS/remap2018_nr_macs2_hg19_v1_2.bed.gz", sep = "")
  } 
  else if (version == "2015" && assembly == "hg38") {
    url <- paste(url, "remap2015/hg38/MACS/remap2015_nr_macs2_hg38_v1.bed.gz", sep = "")
  }
  else if (version == "2015" && assembly == "hg19") {
    url <- paste(url, "remap2015/hg19/MACS/All/remap2015_nr_macs2_hg19_v1.bed.gz", sep = "")
  }
  
  else if (version == '2020' && assembly == 'tair10'){
    url <- paste(url, "remap2020/tair10/tf/MACS2/remap2020_nr_macs2_TAIR10_v1_0.bed.gz", sep = "")
  }
  
  else if (version == '2020' && assembly == 'hg38'){
    url <- paste(url, "remap2020/hg38/MACS2/remap2020_nr_macs2_hg38_v1_0.bed.gz", sep = "")
  }
  
  else if (version == "2022" && assembly == "dm6") {
    url <- paste(url, "remap2022/dm6/MACS2/remap2022_nr_macs2_dm6_v1_0.bed.gz", sep = "")
  }
  else if (version == '2022' && assembly == "mm10"){
    url <- paste(url, "remap2022/mm10/MACS2/remap2022_nr_macs2_mm10_v1_0.bed.gz", sep = "")
  }
  
  else if (version == '2022' && assembly == 'hg38'){
    url <- paste(url, "remap2022/hg38/MACS2/remap2022_nr_macs2_hg38_v1_0.bed.gz" , sep = "")
  }
  
  else if (version == '2022' && assembly == 'hg19'){
    url <- paste(url, 'remap2022/hg19/MACS2/remap2022_nr_macs2_hg19_v1_0.bed.gz', sep = '')
    
  }
  
  if (fileName == "") {
    splits <- strsplit(url, "/")
    fileName <- gsub(".gz", "", tail(unlist(splits), n = 1))
  }
  filePath <-file.path(targetDir,fileName)
  fileExists <- file.exists(filePath)
  input <- "Y"
  if (!force && !fileExists) {
    input <- readline(prompt=paste("A ", size, " GB file will be downloaded. 
                          Do you want to continue Y/N : "))
    while (input != "Y" && input != "N") {
      input <- readline(prompt="Please type Y or N and press Enter : ")
    }
  }
  if (input == "Y") {
    if (fileExists && !force) {
      message("The file ", fileName, " already exists. You may want to use
                    'force = TRUE' to overwrite this file.")
      if (store) {
        return(filePath) 
      } else {
        return(bedToGranges(filePath))
      }
    } else {
      tempZipFile <- paste(tempfile(),".bed.gz", sep = "")
      utils::download.file(url, tempZipFile, method="curl")
      R.utils::gunzip(tempZipFile, filePath, overwrite = force)
      unlink(tempZipFile)
      message("A file has been created at ", filePath)
      if (store) {
        return(filePath)
      } else {
        remapCatalog <- bedToGranges(filePath)
        unlink(filePath)
        return(remapCatalog)
      }
    }
  }
}

#however, a bed file is available online for ReMap 2022:
#https://remap.univ-amu.fr/download_page
#also available (new in 2022) are "non-redundant" peaks and CRMs

##### downloading and filtering all of ReMap2022 ####
##### import key tables ####

#filter down (e.g. to just promoter regions or proximal regions of expressed genes) then re-save

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

#### correct TSS for lncRNAs ####

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

#
##### overlap with remap ####
#remap:
#non-redundant peaks bed file, read from local (5GB) is impossible on laptop, took from 14:05-14:37 on PC:
#process a bit at a time (5 mill peaks for now)
remapCatalog2022hg38 <- read.delim("remap2022_nr_macs2_hg38_v1_0.bed", nrows = 5000000,
                                   header = F)
dim(remapCatalog2022hg38) #5M as expected

trial <- makeGRangesFromDataFrame(remapCatalog2022hg38, 
                                  start.field = "V2", 
                                  end.field = "V3", seqnames.field = "V1", keep.extra.columns = T)
rm(remapCatalog2022hg38)

#remove non standards
seqlevels(trial)
table(as.character(seqnames(trial)))
trial <- keepStandardChromosomes(trial, pruning.mode = "coarse")
seqlevels(trial)
table(as.character(seqnames(trial)))

#query all genes granges with remap:
remapindex <- findOverlaps(query = trial, subject = allGB_GR_promoters_reduce)

#batch1, 737,042 NR promoter peaks:
batch1_SVSMC_match <- unique(trial[queryHits(remapindex)])
write.table(batch1_SVSMC_match, "batch1_SVSMC_matchi.bed", row.names = F, col.names = F, quote = F)
dim(batch1_SVSMC_match)


#next 5mil-20mil
remapCatalog2022hg38_2 <- read.delim("remap2022_nr_macs2_hg38_v1_0.bed", skip = 5000000, nrows = 15000000,
                                   header = F)
dim(remapCatalog2022hg38_2)

trial <- makeGRangesFromDataFrame(remapCatalog2022hg38_2, 
                                  start.field = "V2", 
                                  end.field = "V3", seqnames.field = "V1", keep.extra.columns = T)
rm(remapCatalog2022hg38_2)

#remove non standards
seqlevels(trial)
table(as.character(seqnames(trial)))
trial <- keepStandardChromosomes(trial, pruning.mode = "coarse")
seqlevels(trial)
table(as.character(seqnames(trial)))

#query all genes granges with remap:
remapindex <- findOverlaps(query = trial, subject = allGB_GR_promoters_reduce)

#batch2, 1.941091M
batch2_SVSMC_match <- unique(trial[queryHits(remapindex)])
write.table(batch2_SVSMC_match, "batch2_SVSMC_matchi.bed", row.names = F, col.names = F, quote = F)

#what is the re-import like?
batch2_SVSMC_match_reimp <- read.table("batch2_SVSMC_matchi.bed", nrows = 10)


#if using +/- 500kbp then seems like 95% of ReMap peaks are being picked up...
4.7/5
13.86/15
#promoters is far more discerning, 12-14% so far
0.74/5
1.94/15


#next 20mil-55mil
remapCatalog2022hg38_3 <- read.delim("remap2022_nr_macs2_hg38_v1_0.bed", 
                                     skip  = 20000000, 
                                     nrows = 35000000,
                                     header = F)
dim(remapCatalog2022hg38_3)

trial <- makeGRangesFromDataFrame(remapCatalog2022hg38_3, 
                                  start.field = "V2", 
                                  end.field = "V3", seqnames.field = "V1", keep.extra.columns = T)
rm(remapCatalog2022hg38_3)

#remove non standards
seqlevels(trial)
table(as.character(seqnames(trial)))
trial <- keepStandardChromosomes(trial, pruning.mode = "coarse")
seqlevels(trial)
table(as.character(seqnames(trial)))

#query all genes granges with remap:
remapindex <- findOverlaps(query = trial, subject = allGB_GR_promoters_reduce)

#batch3
batch3_SVSMC_match <- unique(trial[queryHits(remapindex)])
write.table(batch3_SVSMC_match, "batch3_SVSMC_matchi.bed", row.names = F, col.names = F, quote = F)

4.734/35


#final 18.2mil
remapCatalog2022hg38_4 <- read.delim("remap2022_nr_macs2_hg38_v1_0.bed", 
                                     skip  = 55000000,
                                     header = F)
dim(remapCatalog2022hg38_4)

trial <- makeGRangesFromDataFrame(remapCatalog2022hg38_4, 
                                  start.field = "V2", 
                                  end.field = "V3", seqnames.field = "V1", keep.extra.columns = T)
rm(remapCatalog2022hg38_4)

#remove non standards
seqlevels(trial)
table(as.character(seqnames(trial)))
trial <- keepStandardChromosomes(trial, pruning.mode = "coarse")
seqlevels(trial)
table(as.character(seqnames(trial)))

#query all genes granges with remap:
remapindex <- findOverlaps(query = trial, subject = allGB_GR_promoters_reduce)

#batch4
batch4_SVSMC_match <- unique(trial[queryHits(remapindex)])
write.table(batch4_SVSMC_match, "batch4_SVSMC_matchi.bed", row.names = F, col.names = F, quote = F)
1561657/13655741

#doesn't quite add up, should be 68.2... not 68.6mil

#nvm, presumably could either be new peaks but not updated no. on website?
#or only counting on std chromosomes maybe

remap_SVSMC_prom <- list(batch1_SVSMC_match, batch2_SVSMC_match, batch3_SVSMC_match, batch4_SVSMC_match)
class(remap_SVSMC_prom)

remap_SVSMC_prom <- do.call(c, as(remap_SVSMC_prom, "GRangesList"))

remap_SVSMC_prom <- unique(remap_SVSMC_prom)

#8,973,873 peaks, matches:
737042+1941091 + 4734083 + 1561657


##### reimport + enrichment testing ####

#reimport - not doable on laptop
#batch1_SVSMC_match <- read.table("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/batch1_SVSMC_match.bed")
#batch2_SVSMC_match <- read.table("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/batch2_SVSMC_match.bed")
#batch3_SVSMC_match <- read.table("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/batch3_SVSMC_match.bed")
#batch4_SVSMC_match <- read.table("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/batch4_SVSMC_match.bed")

#remap_SVSMC_prom <- unique(rbind(batch1_SVSMC_match, batch2_SVSMC_match, batch3_SVSMC_match, batch4_SVSMC_match))

#remap_SVSMC_prom_GR <- makeGRangesFromDataFrame(remap_SVSMC_prom, start.field = "V2", end.field = "V3", 
#                                                seqnames = "V1", keep.extra.columns = T)
#rm(remap_SVSMC_prom)
#rm(batch1_SVSMC_match)
#rm(batch2_SVSMC_match)
#rm(batch3_SVSMC_match)
#rm(batch4_SVSMC_match)
#gc()

#saveRDS(remap_SVSMC_prom_GR, "remap_SVSMC_prom_GR.rds")

remap_SVSMC_prom_GR <- readRDS("remap_SVSMC_prom_GR.rds")

#identify peaks for certain promoters

#early lncRNAs:
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

table(fpkm_allGDE$RegulationStart)
table(fpkm_allGDE$GeneClassUpdate)

fpkm_allGDE_Upwithin_4_lncs <- filter(fpkm_allGDE_Upwithin_4, grepl("Lnc|fide", GeneClassUpdate))

#subset the allGB object for prom coords for early induced lncs
ERlnc_GR_promoters <- allGB_GR_promoters[allGB_GR_promoters$EnsID %in% fpkm_allGDE_Upwithin_4_lncs$EnsID]


##### find all peaks within the DF lncRNA promoters ####

remapindex <- findOverlaps(query = remap_SVSMC_prom_GR, subject = ERlnc_GR_promoters)

#col to match on
remap_SVSMC_prom_GR$peakID <- paste(as.character(seqnames(remap_SVSMC_prom_GR)),
                                 start(remap_SVSMC_prom_GR), end(remap_SVSMC_prom_GR), remap_SVSMC_prom_GR$V4, sep = "_")

#return matched peaks:
ERlnc_remap2022 <- data.frame("PeakType" = remap_SVSMC_prom_GR$peakID[queryHits(remapindex)],
                              "PeakTF" = remap_SVSMC_prom_GR$id[queryHits(remapindex)],
                                        "LncRNA" = ERlnc_GR_promoters$MSTRG_Tx_ID[subjectHits(remapindex)])

#if needed matched peaks as a GR:
ERlnc_remap2022_GR <- remap_SVSMC_prom_GR[remap_SVSMC_prom_GR$peakID %in% ERlnc_remap2022$PeakType]


#column format, just the TF, called "id"
ERlnc_remap2022_GR$id <- sapply(strsplit(ERlnc_remap2022_GR$V6, ":"), "[[", 1)
remap_SVSMC_prom_GR$id <- sapply(strsplit(remap_SVSMC_prom_GR$V6, ":"), "[[", 1)


###
##### attempt 1 remap peaks vs. remap peaks in lnc promoters (incorrect but kept in for explanation) ####
#the test example is overlapping the peaks of SOX2 with peaks of other TFs
#this will overlap peaks found in promoters with peaks found in other promoters - not the desired format
enrichment.df <- enrichment(ERlnc_remap2022_GR, remap_SVSMC_prom)

enrichment.df.test <- enrichment(query, catalog, byChrom = TRUE)
#>1000s of TFs with p val 0
#presumably because it is v. likely for these peaks to overlap


#set up catalog correctly?
catalogFile2 <- system.file("extdata",
                           "ReMap_nrPeaks_public_chr22.bed",
                           package = "ReMapEnrich")
catalog2 <- bedToGranges(catalogFile2)
#can't see any obvious diff

#score column influences?
catalog3 <- catalog
catalog3$score <- NA
enrichment.df.test <- enrichment(query, catalog3, byChrom = TRUE)
#no

#not seeing any other obvious factors to be aware of

##### attempt 2, remap peaks vs. lnc promoters #####

#so maybe providing just as promoter regions (rather than peaks within promoters) will make more sense
enrichment.df.prom <- enrichment(ERlnc_GR_promoters, remap_SVSMC_prom)

#still 100s of TFs but p values lower now

#notable are SMARCA4, MED1, BRD4, EP300, EZH2, KMT2A (MLL1)

#this is what would be expected if lots of chrom remodellers are binding the lncs
#a first indication that it is working (but need some comparisons to be sure)

#check numbers manually to help understand:
#SMARCA4 201 overlaps in table:
#so between ER lnc promoters and remap peaks
#but in remap2022 there are 246...
#for 92 promoters
ERlnc_remap2022_SMARCA4 <- filter(ERlnc_remap2022, grepl("SMARCA4", PeakTF))


##### attempt 3, remap peaks vs. lnc promoters - any overlap extent #####

#it looks like defaults when enrichment function is overlapping 
#are for peaks to have >10% coverage with the promoter (and vice versa)
#less overlap stringency may get the match right
enrichment.df.prom <- enrichment(ERlnc_GR_promoters, remap_SVSMC_prom_GR, fractionQuery = 0.0001, fractionCatalog = 0.0001)


#now there are 246 overlaps, the correct number of SMARCA4 peak-promoter pairs
#makes sense to use this, allow binding sites at the edge of the region
ERlnc_remap2022_SMARCA4 <- filter(ERlnc_remap2022, grepl("SMARCA4", PeakTF))

#in total there are 29846 SMARCA4 peaks overlapping a expressed promoter in ReMap
lengths(split(remap_SVSMC_prom_GR@elementMetadata$id, 
              remap_SVSMC_prom_GR@elementMetadata$id))["SMARCA4"]

#there are 33593 peaks overlapping a early lnc promoter
sum(ERlnc_remap2022_GR$id == "SMARCA4")#154 of the 33593 peaks in an early lnc promoter region are SMARCA4
154/33593

#there are 29846 SMARCA4 peaks over a promoter generally
sum(remap_SVSMC_prom$id == "SMARCA4")
29846/8973873

#so somewhat enriched

#doesn't match... the 0.0082 number...
154/29846
#but this does
246/29846

#think the query peaks need to be non-overlapping otherwise the enrichment thinks there are 246 catalog peak matches
#when in reality only 154


##### attempt 4, remap peaks vs. non-redundant lnc promoters - any overlap extent ####

#reduced regions probs makes most sense otherwise getting some duplicates:
ERlnc_GR_promoters_reduce <- reduce(ERlnc_GR_promoters)
#115 TSS to 66 (61 lncRNAs)

enrichment.df.prom.reduce <- enrichment(ERlnc_GR_promoters_reduce, remap_SVSMC_prom_GR,
                                        fractionQuery = 0.0001, fractionCatalog = 0.0001)
#this is effectively, TSS per gene, rather than TSS per transcript
#only difference is any close together TSS will be merged

#in total there are 29846 SMARCA4 peaks overlapping an expressed promoter in ReMap
lengths(split(remap_SVSMC_prom_GR@elementMetadata$id, 
              remap_SVSMC_prom_GR@elementMetadata$id))["SMARCA4"]

#there are 33593 peaks overlapping a early lnc promoter
sum(ERlnc_remap2022_GR$id == "SMARCA4")#154 of the 33593 peaks in an early lnc promoter region are SMARCA4
154/33593# 154 is reported in this version of the table

#there are 29846 SMARCA4 peaks over a promoter generally
sum(remap_SVSMC_prom_GR$id == "SMARCA4")
29846/8973873

154/29846 #mapped ratio - the % of SMARCA4 peaks in promoter in target set

#possible enrichment via Fisher's? simpler than the poiss approach?

#hits, SMARCA4, selection ER lnc proms, background all DEG proms
a <- 154
b <- 33593
c <- 29846
d <- 8973873

fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                       "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")))

#BRD4, MED1, SMARCA4 are the top 3

#there are still 100s of significant peaks from the q value though

#the stat test is ppois: https://www.geeksforgeeks.org/a-guide-to-dpois-ppois-qpois-and-rpois-in-r/

#This function finds the probability that a certain number of successes or less occur based on an 
#average rate of success, In other words, we can say as this function returns the value of the inverse 
#Poisson cumulative density function.

#https://www.scribbr.co.uk/stats/poisson-distribution-meaning/

#see notes on Poisson in Timecourse folder

#for SMARCA4 this is the source of the p value:
ppois(154, 2.666667, lower.tail = F)

#the 2.66666 number "random average" comes from
#a) taking the query region widths
#b) shuffling them amongst the catalog

#https://github.com/remap-cisreg/ReMapEnrich/blob/master/R/compute_enrichment.R

#this is the code to make background:
shuffles <- replicate(shufflesNumber, shuffle(query, chromSizes, universe, included, byChrom))
#in this function: https://github.com/remap-cisreg/ReMapEnrich/blob/master/R/shuffle.R
#if the universe is NULL, then all chromosome regions are used
#the query region widths (early lnc promoters) are taken randomly across the genome and returned
#the default is to do it 6x
#you can change the universe so that it is only done within set chromosomal regions


# The theorical means are calculated from the shuffles overlaps.
shuffleCatCount <- sapply(shuffles, intersect, catalog = catalog, fractionQuery = fractionQuery,
                          fractionCatalog = fractionCatalog, categories = categories)


##### attempt 5, as above, restricted just to universe of all DEG promoters ####

DEGB_GR_promoters <- allGB_GR_promoters[allGB_GR_promoters$EnsID %in% fpkm_allGDE$EnsID]
DEGB_GR_promoters_reduce <- reduce(DEGB_GR_promoters)

#here the shuffling is attempted purely within DE promoter regions (previously genome wide?)
enrichment.df.prom.reduce.DEGuni <- enrichment(ERlnc_GR_promoters_reduce, remap_SVSMC_prom_GR,
                                        universe = DEGB_GR_promoters_reduce,
                                        fractionQuery = 0.0001, fractionCatalog = 0.0001)
#much more restrained enrichment

#SMARCA4:
#154 vs. random average of 148


##### attempt 6, as above, restricted just to universe of all promoters ####


#here the shuffling is attempted purely within all promoter regions
enrichment.df.prom.reduce.allPromuni <- enrichment(ERlnc_GR_promoters_reduce, remap_SVSMC_prom_GR,
                                            universe = allGB_GR_promoters_reduce,
                                            fractionQuery = 0.0001, fractionCatalog = 0.0001)
#much more restrained enrichment

#SMARCA4:
#154 vs. random average of 140

#this seems unfair, reminder - 154 of the 33593 peaks found in early induced lnc promoters are SMARCA4
ERlnc_remap2022_GR
sum(ERlnc_remap2022_GR$id == "SMARCA4")
154/33593# as is still reported in this version of the table

#there are 29846 SMARCA4 peaks from the 8973873 peaks overlapping a promoter generally
remap_SVSMC_prom_GR
sum(remap_SVSMC_prom_GR$id == "SMARCA4")
29846/8973873

#0.46 vs. 0.33 seems a clear enrichment, the Fisher test done above is significant

#so how is the background rate so similar??


##### closer look at shuffling function #####

#shuffling function:
trial <- shuffle(ERlnc_GR_promoters_reduce, universe = allGB_GR_promoters_reduce)

#66 shuffled ranges to match that of the query sequence

#precisely matching regions taken
table(width(trial))
table(width(ERlnc_GR_promoters_reduce))
width(trial) == width(ERlnc_GR_promoters_reduce)

#presumably these correspond to other promoters in the allGB table? kind of...
allGB_GR_promoters_reduce_df <- as.data.frame(allGB_GR_promoters_reduce)
#this example (n.b. doesn't work if doing a re-shuffle)
trial[66]
#starts midway through a promoter and ends earlier
filter(allGB_GR_promoters_reduce_df, start >= 108000934, end <= 108006111, seqnames == "chr7")
#presumably to match on the distance

#this could be a key reason to use their poisson method then, so the promoter distances can be matched
#and the background shuffling is fairer way
#(is this permutation?)

#so is this trial region filled with SMARCA4 peaks?
remapindex <- findOverlaps(query = remap_SVSMC_prom_GR, subject = trial)

trial$locID <- paste(trial@seqnames@values, trial@ranges@start, trial@ranges@width, sep = "_")

#return matched peaks:
shuffled_remap2022 <- data.frame("PeakType" = remap_SVSMC_prom_GR$peakID[queryHits(remapindex)],
                              "region" = trial$locID[subjectHits(remapindex)])

#if needed matched peaks as a GR:
shuffled_remap2022_GR <- remap_SVSMC_prom_GR[remap_SVSMC_prom_GR$peakID %in% shuffled_remap2022$PeakType]
length(shuffled_remap2022_GR)
sum(shuffled_remap2022_GR$id == "SMARCA4")
146/49281
151/49183
157/47195
#etc

#so overall the lower % in background IS maintained
#the overall peak count in the shuffled regions is much higher tho (49k vs. 33k)
#unknown why, potentially due to ChIPseq peaks accumulating more at PCGs (more ubiquitous, higher expressed etc)
#seems this issue isn't considered in their approach
#so the conclusion is probs true (more SMARCA4 peaks at promoters) 
#but with a massive caveat (not in terms of overall proportion)
#DO NOT USE FOR THIS PURPOSE (WITHOUT ADAPTION)


##### n.b. issue with some of the Fisher's approaches above and next (by peaks or by gene question) ####

#the following code uses an approach of 
#"enrichment of peaks for TF-x" in "target promoters" vs. "all promoters"
#so % of TF-x peaks vs. all peaks
#but the "all peaks" is problematic, ReMap2022 is obvs incomplete/not all TFs/in all cell types/states

#better approach is
#"enrichment of TF-bound promoters" amongst "targets" vs. "all promoters"
#in that way, each TF is considered on it's own merit
#obvs still influenced by ReMap2022 "completeness", weaker profiled TFs will have less info
#but means the hits are more reliable/reflective of biology
#(less false positives)

##### Fishers attempt - SMARCA4 ####

#33593 peaks in ER lncs, of which 154 are SMARCA4
ERlnc_remap2022_GR
sum(ERlnc_remap2022_GR$id == "SMARCA4")
154/33593# as is still reported in this version of the table


#background of all peaks overlapping a promoter:
remap_SVSMC_prom_GR
sum(remap_SVSMC_prom_GR$id == "SMARCA4")
29846/8973873


#background of all peaks overlapping a DE promoter:
#subset the allGB object for prom coords for DEGs 
DEGB_GR_promoters <- allGB_GR_promoters[allGB_GR_promoters$EnsID %in% fpkm_allGDE$EnsID]

remapindex <- findOverlaps(query = remap_SVSMC_prom_GR, subject = DEGB_GR_promoters)

#col to match on
remap_SVSMC_prom_GR$peakID <- paste(as.character(seqnames(remap_SVSMC_prom_GR)),
                                    start(remap_SVSMC_prom_GR), end(remap_SVSMC_prom_GR), remap_SVSMC_prom_GR$V4, sep = "_")

#return matched peaks:
DEG_remap2022 <- data.frame("PeakType" = remap_SVSMC_prom_GR$peakID[queryHits(remapindex)],
                              "DEG" = DEGB_GR_promoters$MSTRG_Tx_ID[subjectHits(remapindex)])

#if needed matched peaks as a GR:
DEG_remap2022_GR <- remap_SVSMC_prom_GR[remap_SVSMC_prom_GR$peakID %in% DEG_remap2022$PeakType]

#column format, just the TF, called "id"
DEG_remap2022_GR$id <- sapply(strsplit(DEG_remap2022_GR$V6, ":"), "[[", 1)


#3690750 peaks
DEG_remap2022_GR
sum(DEG_remap2022_GR$id == "SMARCA4")
12895/3690750

a <- sum(ERlnc_remap2022_GR$id == "SMARCA4")
b <- length(ERlnc_remap2022_GR)
c <- sum(DEG_remap2022_GR$id == "SMARCA4")
d <- length(DEG_remap2022_GR)

fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                       "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")))

#SMARCA4 peaks enriched in ER lnc promoters relative to all promoters, DE promoters (OR 1.3)

###
##### re-run for early lncs + all TFs ####

TFs_inTargetProms <- unique(ERlnc_remap2022_GR$id)

TFs_inTargetProms_Fish <- list()

for(i in 1:length(TFs_inTargetProms)){
  a <- sum(ERlnc_remap2022_GR$id == TFs_inTargetProms[i])
  b <- length(ERlnc_remap2022_GR)
  c <- sum(DEG_remap2022_GR$id == TFs_inTargetProms[i])
  d <- length(DEG_remap2022_GR)
  
  TFs_inTargetProms_Fish[[i]] <- c(
    a,b,c,d,
    fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                           "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")))$p,
    fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                           "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")))$est
  )
}

names(TFs_inTargetProms_Fish) <- TFs_inTargetProms

TFs_inTargetProms_Fish_df <- as.data.frame(t(bind_rows(TFs_inTargetProms_Fish)))

TFs_inTargetProms_Fish_df$id <- TFs_inTargetProms

#looks quite similar to attempt 5

#wonder if this means that the lnc promoters are not too disimilar from the other promoters

#and if I'm missing something on why the poiss version makes more sense? shuffles make more robust?

#and if the cisLnc TFs are actually just quite generic transcriptional activators

#corrected p vals for context:
#remove TFBS with low coverage:
summary(TFs_inTargetProms_Fish_df$V1)
TFs_inTargetProms_Fish_df <- filter(TFs_inTargetProms_Fish_df, V1 > 10)
TFs_inTargetProms_Fish_df$BH_pVal <- p.adjust(TFs_inTargetProms_Fish_df$V5, method = "BH")

#SMARCA4, MED1, BRD4, BICRA all enriched in lnc promoters vs. other DEG promoters

##### re-run for early lncs - all promoter background ####

TFs_inTargetProms <- unique(ERlnc_remap2022_GR$id)

TFs_inTargetProms_all_Fish <- list()

for(i in 1:length(TFs_inTargetProms)){
  a <- sum(ERlnc_remap2022_GR$id == TFs_inTargetProms[i])
  b <- length(ERlnc_remap2022_GR)
  c <- sum(remap_SVSMC_prom_GR$id == TFs_inTargetProms[i])
  d <- length(remap_SVSMC_prom_GR)
  
  TFs_inTargetProms_all_Fish[[i]] <- c(
    a,b,c,d,
    fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                           "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")))$p,
    fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                           "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")))$est
  )
}

names(TFs_inTargetProms_all_Fish) <- TFs_inTargetProms

TFs_inTargetProms_all_Fish_df <- as.data.frame(t(bind_rows(TFs_inTargetProms_all_Fish)))

TFs_inTargetProms_all_Fish_df$id <- TFs_inTargetProms

#corrected p vals for context:
#remove TFBS with low coverage:
summary(TFs_inTargetProms_all_Fish_df$V1)
TFs_inTargetProms_all_Fish_df <- filter(TFs_inTargetProms_all_Fish_df, V1 > 10)
TFs_inTargetProms_all_Fish_df$BH_pVal <- p.adjust(TFs_inTargetProms_all_Fish_df$V5, method = "BH")

#looks similar to other background? very much so, little diff
filter(TFs_inTargetProms_Fish_df, BH_pVal < 0.05)$id
filter(TFs_inTargetProms_all_Fish_df, BH_pVal < 0.05)$id

filter(TFs_inTargetProms_Fish_df, BH_pVal < 0.05)$id %in%
filter(TFs_inTargetProms_all_Fish_df, BH_pVal < 0.05)$id


##### re-run for early vs other DEGs ####

#background of all peaks overlapping a DE promoter:
#subset the allGB object for prom coords for DEGs 
DE4up_GB_GR_promoters <- allGB_GR_promoters[allGB_GR_promoters$EnsID %in% fpkm_allGDE_Upwithin_4$EnsID]

remapindex <- findOverlaps(query = remap_SVSMC_prom_GR, subject = DE4up_GB_GR_promoters)

#return matched peaks:
DEG_4up_remap2022 <- data.frame("PeakType" = remap_SVSMC_prom_GR$peakID[queryHits(remapindex)],
                            "DEG" = DE4up_GB_GR_promoters$MSTRG_Tx_ID[subjectHits(remapindex)])

#if needed matched peaks as a GR:
DEG_4up_remap2022_GR <- remap_SVSMC_prom_GR[remap_SVSMC_prom_GR$peakID %in% DEG_4up_remap2022$PeakType]

#column format, just the TF, called "id"
DEG_4up_remap2022_GR$id <- sapply(strsplit(DEG_4up_remap2022_GR$V6, ":"), "[[", 1)


TFs_inTargetProms <- unique(DEG_4up_remap2022_GR$id)

TFs_in_4hrDEGProms_Fish <- list()

for(i in 1:length(TFs_inTargetProms)){
  a <- sum(DEG_4up_remap2022_GR$id == TFs_inTargetProms[i])
  b <- length(DEG_4up_remap2022_GR)
  c <- sum(DEG_remap2022_GR$id == TFs_inTargetProms[i])
  d <- length(DEG_remap2022_GR)
  
  TFs_in_4hrDEGProms_Fish[[i]] <- c(
    a,b,c,d,
    fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                           "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")))$p,
    fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                           "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")))$est
  )
}

names(TFs_in_4hrDEGProms_Fish) <- TFs_inTargetProms

TFs_in_4hrDEGProms_Fish_df <- as.data.frame(t(bind_rows(TFs_in_4hrDEGProms_Fish)))

TFs_in_4hrDEGProms_Fish_df$id <- TFs_inTargetProms

#corrected p vals for context:
#remove TFBS with low coverage:
summary(TFs_in_4hrDEGProms_Fish_df$V1)
TFs_in_4hrDEGProms_Fish_df <- filter(TFs_in_4hrDEGProms_Fish_df, V1 > 10)
TFs_in_4hrDEGProms_Fish_df$BH_pVal <- p.adjust(TFs_in_4hrDEGProms_Fish_df$V5, method = "BH")

#looks similar to other background? very much so, little diff
filter(TFs_in_4hrDEGProms_Fish_df, BH_pVal < 0.05)$id
filter(TFs_inTargetProms_Fish_df, BH_pVal < 0.05)$id

#somewhat, though a set of specific lncRNA TFs seems clear
filter(TFs_inTargetProms_Fish_df, BH_pVal < 0.05)$id[!filter(TFs_inTargetProms_Fish_df, BH_pVal < 0.05)$id %in%
                                                       filter(TFs_in_4hrDEGProms_Fish_df, BH_pVal < 0.05)$id]
#includes SMARCA2, BRD3, BRD2, STAT5, ONECUT2


##### re-run for early vs other DEGs - all promoter background ####

TFs_in_4hrDEGProms_all_Fish <- list()

for(i in 1:length(TFs_inTargetProms)){
  a <- sum(DEG_4up_remap2022_GR$id == TFs_inTargetProms[i])
  b <- length(DEG_4up_remap2022_GR)
  c <- sum(remap_SVSMC_prom_GR$id == TFs_inTargetProms[i])
  d <- length(remap_SVSMC_prom_GR)
  
  TFs_in_4hrDEGProms_all_Fish[[i]] <- c(
    a,b,c,d,
    fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                           "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")))$p,
    fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                           "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")))$est
  )
}

names(TFs_in_4hrDEGProms_all_Fish) <- TFs_inTargetProms

TFs_in_4hrDEGProms_all_Fish_df <- as.data.frame(t(bind_rows(TFs_in_4hrDEGProms_all_Fish)))

TFs_in_4hrDEGProms_all_Fish$id <- TFs_inTargetProms


##### conclusion Fisher's vs. poisson ####

#poisson approach from R package has been interesting but has a flaw in this approach 
#in that a higher number of peaks are being found in the shuffled backgrounds
#unclear why, possibly because contains PCGs with more obvious TF binding
#either way the approach is flawed for (promoter analysis at least) as doesn't account for more/less peaks in target vs. background
#Fisher's will work better to address this issue, e.g. a SMARCA4 increase at early lnc is not captured otherwise

#a set of TFs have enriched binding capacity amongst early induced DELs and DEGs relative to other DEGs

#the DEGs have genes known to function downstream of cisLncs:  
#SMARCA4
#BRD4
#MED1, MED12, MED25

#other associated chromatin remodellers in the top include
#MAML1
#SUPT16H (FACT complex)
#BICRA

#also IEG TFs:
#MYC, FOSL1, JUNB, JUN

#and immune TFs:
#RELA, NFKB2, STAT3

#but some are only found enriched at lnc sites:
#SMARCA2 (related to SWI/SNF)
#BRD2, BRD3 (also BET like BRD4 though diff function)

#does it help fig2? kind of... try a set of figures


##### all TFs, all clusters, Fisher's by peak - enrich of peaks in target proms vs. all used proms #####

#peaks for promoters in each cluster
#subset the allGB object for prom coords for DEGs
all_DEGs_list <- list(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Downwithin_4,
                      fpkm_allGDE_Upwithin_8, fpkm_allGDE_Downwithin_8,
                      fpkm_allGDE_Upwithin_24, fpkm_allGDE_Downwithin_24)

all_DEGs_Fish_list <- list()

#for each TF (filter these? e.g. exp level, e.g. lambert list):
#enrichment of TF-x peaks within target promoters compared to TF-x peaks within all promoters
#this may be skewed by a large number of peaks within a small subset of the promoters?
#e.g. in more ubiquitously found promoters or less specific promoters
#alt: enrichment of TF-x promoters within target promoters compared to all promoters
#this treats all promoters equally and gives a binary yes/no label for TF-x

#first option (by peaks) looks like this
for (j in 1:length(all_DEGs_list)){

  #promoter regions for given cluster of interest 
  DE4up_GB_GR_promoters <- allGB_GR_promoters[allGB_GR_promoters$EnsID %in% all_DEGs_list[[j]]$EnsID]
  
  remapindex <- findOverlaps(query = remap_SVSMC_prom_GR, subject = DE4up_GB_GR_promoters)
  
  #return matched peaks from ReMap2022:
  DEG_4up_remap2022 <- data.frame("PeakType" = remap_SVSMC_prom_GR$peakID[queryHits(remapindex)],
                                  "DEG" = DE4up_GB_GR_promoters$MSTRG_Tx_ID[subjectHits(remapindex)])
  
  #matched peaks in target genes as a GR:
  DEG_4up_remap2022_GR <- remap_SVSMC_prom_GR[remap_SVSMC_prom_GR$peakID %in% DEG_4up_remap2022$PeakType]
  
  TFs_inTargetProms <- unique(DEG_4up_remap2022_GR$id)
  TFs_in_4hrDEGProms_all_Fish <- list()
  
  
  for(i in 1:length(TFs_inTargetProms)){
    
    a <- sum(DEG_4up_remap2022_GR$id == TFs_inTargetProms[i])
    b <- length(DEG_4up_remap2022_GR)
    c <- sum(remap_SVSMC_prom_GR$id == TFs_inTargetProms[i])
    d <- length(remap_SVSMC_prom_GR)
    
    TFs_in_4hrDEGProms_all_Fish[[i]] <- c(
      a,b,c,d,
      fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                             "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p,
      fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                             "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est
    )
  }
  
  names(TFs_in_4hrDEGProms_all_Fish) <- TFs_inTargetProms
  
  TFs_in_4hrDEGProms_all_Fish_df <- as.data.frame(t(bind_rows(TFs_in_4hrDEGProms_all_Fish)))
  
  TFs_in_4hrDEGProms_all_Fish_df$id <- TFs_inTargetProms
  
  all_DEGs_Fish_list[[j]] <- TFs_in_4hrDEGProms_all_Fish_df
}

names(all_DEGs_Fish_list) <- c("4hrUp", "4hrDown", "8hrUp", "8hrDown", "24hrUp", "24hrDown")

all_DEGs_Fish_df <- as.data.frame((bind_rows(all_DEGs_Fish_list, .id = "cluster")))

#some TFs presumably have v. low number of peaks:
sapply(all_DEGs_Fish_list, function(x){summary(x$V1)})

#correct p val
all_DEGs_Fish_list_correct <- list()

#remove TFs which are not expressed or classed as TF (Lambert):
TFsExpressed <- filter(fpkm_allG, grepl("TF", GeneClassUpdate))
TFsExpressed <- unique(c(TFsExpressed$EnsName, 
                         gsub(" ", "", unlist(strsplit(filter(TFsExpressed, !is.na(AllNames))$AllNames, ",")))))
#could also remove anything not in the Lambert list too..

for(j in 1:length(all_DEGs_Fish_list)){
  #require at least 20 TF peaks to be in the selected group (otherwise may over-adjust p due to low interest TFs)
  trial <- filter(all_DEGs_Fish_list[[j]], V1 >20, id %in% TFsExpressed)
  trial$BH_p <- p.adjust(trial$V5, method = "BH")
  all_DEGs_Fish_list_correct[[j]] <- trial
}

names(all_DEGs_Fish_list_correct) <- c("4hrUp", "4hrDown", "8hrUp", "8hrDown", "24hrUp", "24hrDown")

all_DEGs_Fish_df_correct <- as.data.frame((bind_rows(all_DEGs_Fish_list_correct, .id = "cluster")))

#enrichments over time:
all_DEGs_Fish_df_correct_sig <- filter(all_DEGs_Fish_df_correct, BH_p < 0.05)

#write.csv(all_DEGs_Fish_df_correct, "TFenriched_DEGs_byPeak.csv", row.names = F)


##### all TFs, all clusters, Fisher's by gene - enrich of TF-targeted genes in target genes vs. all used genes #####

all_DEGs_list <- list(fpkm_allGDE_Upwithin_4, fpkm_allGDE_Downwithin_4,
                      fpkm_allGDE_Upwithin_8, fpkm_allGDE_Downwithin_8,
                      fpkm_allGDE_Upwithin_24, fpkm_allGDE_Downwithin_24)

all_DEGs_Fish_list <- list()

#for each TF (filter these? e.g. exp level, e.g. lambert list):
#enrichment of TF-x peaks within target promoters compared to TF-x peaks within all promoters
#this may be skewed by a large number of peaks within a small subset of the promoters?
#e.g. in more ubiquitously found promoters or less specific promoters
#alt: enrichment of TF-x promoters within target promoters compared to all promoters
#this treats all promoters equally and gives a binary yes/no label for TF-x

#second option (by gene) looks like this

remap_SVSMC_prom_GR$peakID <- paste(as.character(seqnames(remap_SVSMC_prom_GR)),
                                    start(remap_SVSMC_prom_GR), end(remap_SVSMC_prom_GR), remap_SVSMC_prom_GR$V4, sep = "_")

remap_SVSMC_prom_GR$id <- sapply(strsplit(remap_SVSMC_prom_GR$V6, ":"), "[[", 1)


#prom ID for promoter co-ordinates:
allGB_GR_promoters_reduce$promID <- paste(as.character(seqnames(allGB_GR_promoters_reduce)),
                                      start(allGB_GR_promoters_reduce), end(allGB_GR_promoters_reduce), 
                                      allGB_GR_promoters_reduce$MSTRG_Tx_ID, sep = "_")

#add in gene ID per reduced prom:
remapindex <- findOverlaps(query = allGB_GR_promoters, subject = allGB_GR_promoters_reduce)

#return matched peaks from ReMap2022:
allPromRed_EnsID <- unique(data.frame("EnsID" = allGB_GR_promoters$EnsID[queryHits(remapindex)],
                                "promID" = allGB_GR_promoters_reduce$promID[subjectHits(remapindex)]))

#now make match of reduced proms to remap peaks
remapindex <- findOverlaps(query = remap_SVSMC_prom_GR, subject = allGB_GR_promoters_reduce)

#return matched peaks from ReMap2022:
allProm_remap2022 <- data.frame("peakID" = remap_SVSMC_prom_GR$peakID[queryHits(remapindex)],
                                "peakTF" = remap_SVSMC_prom_GR$id[queryHits(remapindex)],
                                "peakTF_cell" = remap_SVSMC_prom_GR$V6[queryHits(remapindex)],
                                "promID" = allGB_GR_promoters_reduce$promID[subjectHits(remapindex)])

#now can combine all above to get for each gene ID, add any bound TFs in ReMap 
trial <- unique(merge(unique(allProm_remap2022[,c(2,4)]), allPromRed_EnsID, by = "promID", all.y = T))

length(unique(trial$promID))#14693 promoters
length(unique(allProm_remap2022$promID))#14442
length(unique(allGB_GR_promoters_reduce$promID))#14693 (probs inc some without a remap overlap)

length(unique(trial$EnsID))#10761 genes
length(unique(allGB_GR_promoters$EnsID))#10761 genes
length(unique(fpkm_allG$EnsID))#10761 genes
length(unique(allPromRed_EnsID$EnsID))#10761 genes

allG_remap2022_TFs <- trial
rm(trial)

#use this to find ReMap regulators which are expressed / with high potential to be active in each window:
fpkm_thresh <- 10
ReMaps_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >fpkm_thresh | Hour4_meanFPKM >fpkm_thresh))
ReMaps_04 <- unique(c(ReMaps_04$EnsName, 
                      gsub(" ", "", unlist(strsplit(filter(ReMaps_04, !is.na(AllNames))$AllNames, ",")))))

ReMaps_48 <- unique(filter(fpkm_allG, Hour4_meanFPKM >fpkm_thresh | Hour8_meanFPKM >fpkm_thresh))
ReMaps_48 <- unique(c(ReMaps_48$EnsName, 
                      gsub(" ", "", unlist(strsplit(filter(ReMaps_48, !is.na(AllNames))$AllNames, ",")))))

ReMaps_824 <- unique(filter(fpkm_allG, Hour8_meanFPKM >fpkm_thresh | Hour24_meanFPKM >fpkm_thresh))
ReMaps_824 <- unique(c(ReMaps_824$EnsName, 
                       gsub(" ", "", unlist(strsplit(filter(ReMaps_824, !is.na(AllNames))$AllNames, ",")))))

#why higher FPKM e.g. 10/20? high and arbitrary...
#TFs lower expressed less likely to bind across a wide group of genes across diff loci, e.g. ~60 up regulated lncRNAs
#the method is using mixed/non-matched data, so increased stringency best to remove uninteresting TFs that would dilute p after MHT

TF_perWindow <- list(ReMaps_04, ReMaps_04, ReMaps_48, ReMaps_48, ReMaps_824, ReMaps_824)


#background, promoters per window:
expr_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >1 | Hour4_meanFPKM >1)$EnsID)
expr_48 <- unique(filter(fpkm_allG, Hour4_meanFPKM >1 | Hour8_meanFPKM >1)$EnsID)
expr_824 <- unique(filter(fpkm_allG, Hour8_meanFPKM >1 | Hour24_meanFPKM >1)$EnsID)

exprs_perWindow <- list(expr_04, expr_04, expr_48, expr_48, expr_824, expr_824)


for (j in 1:length(all_DEGs_list)){
  
  #TFs expressed in window of interest:
  allG_remap2022_TFs_window <- filter(allG_remap2022_TFs, EnsID %in% exprs_perWindow[[j]],
                                      peakTF %in% TF_perWindow[[j]])
  
  #promoter regions for lncRNAs in given cluster of interest 
  TargetGenes <- allG_remap2022_TFs_window[allG_remap2022_TFs_window$EnsID %in% all_DEGs_list[[j]]$EnsID,]
  
  TFs_inTargetProms <- unique(TargetGenes$peakTF)
  TFs_inTargetProms <- TFs_inTargetProms[!is.na(TFs_inTargetProms)]
  TFs_in_4hrDEGProms_all_Fish <- list()
  
  for(i in 1:length(TFs_inTargetProms)){
    
    a <- length(unique(filter(TargetGenes, peakTF == TFs_inTargetProms[i])$EnsID))
    b <- length(unique(TargetGenes$EnsID))
    c <- length(unique(filter(allG_remap2022_TFs_window, peakTF == TFs_inTargetProms[i])$EnsID))
    d <- length(unique(allG_remap2022_TFs_window$EnsID))
    
    TFs_in_4hrDEGProms_all_Fish[[i]] <- c(
      a,b,c,d,
      fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                             "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p,
      fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                             "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est
    )
  }
  
  names(TFs_in_4hrDEGProms_all_Fish) <- TFs_inTargetProms
  
  TFs_in_4hrDEGProms_all_Fish_df <- as.data.frame(t(bind_rows(TFs_in_4hrDEGProms_all_Fish, .id = "TFs")))
  
  TFs_in_4hrDEGProms_all_Fish_df$id <- TFs_inTargetProms
  
  all_DEGs_Fish_list[[j]] <- TFs_in_4hrDEGProms_all_Fish_df
}

names(all_DEGs_Fish_list) <- c("04hrUp", "04hrDown", "48hrUp", "48hrDown", "824hrUp", "824hrDown")

all_DEGs_Fish_df <- as.data.frame((bind_rows(all_DEGs_Fish_list, .id = "cluster")))


#correct p val
all_DEGs_Fish_list_correct <- list()

#remove TFs which are not expressed or classed as TF (Lambert):
#TFsExpressed <- filter(fpkm_allG, grepl("TF", GeneClassUpdate))
#TFsExpressed <- unique(c(TFsExpressed$EnsName, 
#                         gsub(" ", "", unlist(strsplit(filter(TFsExpressed, !is.na(AllNames))$AllNames, ",")))))
#could also remove anything not in the Lambert list too..

for(j in 1:length(all_DEGs_Fish_list)){
  
  #remove any TFs which do not target at least 25 genes (rough area where high sig starts SRSF9)
  
  trial <- filter(all_DEGs_Fish_list[[j]], V3 >25)
  trial$BH_p <- p.adjust(trial$V5, method = "BH")
  all_DEGs_Fish_list_correct[[j]] <- trial
}

names(all_DEGs_Fish_list_correct) <- c("04hrUp", "04hrDown", "48hrUp", "48hrDown", "824hrUp", "824hrDown")

all_DEGs_Fish_df_correct <- as.data.frame((bind_rows(all_DEGs_Fish_list_correct, .id = "cluster")))

#enrichments over time:
all_DEGs_Fish_df_correct_sig <- filter(all_DEGs_Fish_df_correct, BH_p < 0.05)

table(all_DEGs_Fish_df_correct_sig$cluster)

#note - works really well for 48hr up - larger and really specific set of TFs bound there
#not many found for others unexpectedly
#the promoters have a mix of up and down TFBS
#the peaks are aggregated across broad cell types
#probs works well for 48hrUp as cell cycle machinery is more universal?

#write.csv(all_DEGs_Fish_df_correct, "all_DEGs_Fish_df_byGene.csv", row.names = F)


#### plot p vals/ORs against each other - expecting broad agreement though higher power with peak test? ####



#### all TFs, all clusters, Fisher's by gene - focus on just lncRNAs ####

#requires pre-amble of above by gene test

lnc_DEGs_Fish_list <- list()

#option for greater tf stringency:
#fpkm_allG_10 <- filter(fpkm_allG, fpkm_max_treatment >10)
#ReMap2022_inSVSMC <- unique(c(fpkm_allG_10$EnsName, 
#                              gsub(" ", "", unlist(strsplit(filter(fpkm_allG_10, !is.na(AllNames))$AllNames, ",")))))
#allG_remap2022_TFs_10FPKM <- filter(allG_remap2022_TFs, peakTF %in% ReMap2022_inSVSMC)
#
#fpkm_allG_20 <- filter(fpkm_allG, fpkm_max_treatment >20)
#ReMap2022_inSVSMC <- unique(c(fpkm_allG_20$EnsName, 
#                              gsub(" ", "", unlist(strsplit(filter(fpkm_allG_20, !is.na(AllNames))$AllNames, ",")))))
#allG_remap2022_TFs_20FPKM <- filter(allG_remap2022_TFs, peakTF %in% ReMap2022_inSVSMC)##
#
#fpkm_allG_25 <- filter(fpkm_allG, fpkm_max_treatment >25)
#ReMap2022_inSVSMC <- unique(c(fpkm_allG_25$EnsName, 
#                              gsub(" ", "", unlist(strsplit(filter(fpkm_allG_25, !is.na(AllNames))$AllNames, ",")))))
#allG_remap2022_TFs_25FPKM <- filter(allG_remap2022_TFs, peakTF %in% ReMap2022_inSVSMC)

#use this to find ReMap binders which are expressed / with high potential to be active in each window:
fpkm_thresh <- 10
ReMaps_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >fpkm_thresh | Hour4_meanFPKM >fpkm_thresh))
ReMaps_04 <- unique(c(ReMaps_04$EnsName, 
                              gsub(" ", "", unlist(strsplit(filter(ReMaps_04, !is.na(AllNames))$AllNames, ",")))))
                    
ReMaps_48 <- unique(filter(fpkm_allG, Hour4_meanFPKM >fpkm_thresh | Hour8_meanFPKM >fpkm_thresh))
ReMaps_48 <- unique(c(ReMaps_48$EnsName, 
                      gsub(" ", "", unlist(strsplit(filter(ReMaps_48, !is.na(AllNames))$AllNames, ",")))))

ReMaps_824 <- unique(filter(fpkm_allG, Hour8_meanFPKM >fpkm_thresh | Hour24_meanFPKM >fpkm_thresh))
ReMaps_824 <- unique(c(ReMaps_824$EnsName, 
                      gsub(" ", "", unlist(strsplit(filter(ReMaps_824, !is.na(AllNames))$AllNames, ",")))))

#why higher FPKM e.g. 10/20? high and arbitrary...
#TFs lower expressed less likely to bind across a wide group of genes across diff loci, e.g. ~60 up regulated lncRNAs
#the method is using mixed/non-matched data, so increased stringency best to remove uninteresting TFs that would dilute p after MHT

TF_perWindow <- list(ReMaps_04, ReMaps_04, ReMaps_48, ReMaps_48, ReMaps_824, ReMaps_824)


#background, promoters per window:
expr_04 <- unique(filter(fpkm_allG, Hour0_meanFPKM >1 | Hour4_meanFPKM >1)$EnsID)
expr_48 <- unique(filter(fpkm_allG, Hour4_meanFPKM >1 | Hour8_meanFPKM >1)$EnsID)
expr_824 <- unique(filter(fpkm_allG, Hour8_meanFPKM >1 | Hour24_meanFPKM >1)$EnsID)

exprs_perWindow <- list(expr_04, expr_04, expr_48, expr_48, expr_824, expr_824)


for (j in 1:length(all_DEGs_list)){
  
  #TFs expressed in window of interest:
  allG_remap2022_TFs_window <- filter(allG_remap2022_TFs, EnsID %in% exprs_perWindow[[j]],
                                      peakTF %in% TF_perWindow[[j]])
  
  #promoter regions for lncRNAs in given cluster of interest 
  TargetGenes <- allG_remap2022_TFs_window[allG_remap2022_TFs_window$EnsID %in% all_DEGs_list[[j]]$EnsID &
                                             allG_remap2022_TFs_window$EnsID %in% filter(fpkm_allG, grepl("Lnc|fide", GeneClassUpdate))$EnsID,]
  
  TFs_inTargetProms <- unique(TargetGenes$peakTF)
  TFs_inTargetProms <- TFs_inTargetProms[!is.na(TFs_inTargetProms)]
  TFs_in_4hrDEGProms_all_Fish <- list()
  
  for(i in 1:length(TFs_inTargetProms)){
    
    a <- length(unique(filter(TargetGenes, peakTF == TFs_inTargetProms[i])$EnsID))
    b <- length(unique(TargetGenes$EnsID))
    c <- length(unique(filter(allG_remap2022_TFs_window, peakTF == TFs_inTargetProms[i])$EnsID))
    d <- length(unique(allG_remap2022_TFs_window$EnsID))
    
    TFs_in_4hrDEGProms_all_Fish[[i]] <- c(
      a,b,c,d,
      fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                             "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$p,
      fisher.test(data.frame("ER_lnc_proms" = c(a,b-a),
                             "Other_DEG_proms" = c(c-a,d-c-(b-a)), row.names = c("Cluster", "other")), alternative = "greater")$est
    )
  }
  
  names(TFs_in_4hrDEGProms_all_Fish) <- TFs_inTargetProms
  
  TFs_in_4hrDEGProms_all_Fish_df <- as.data.frame(t(bind_rows(TFs_in_4hrDEGProms_all_Fish, .id = "TFs")))
  
  TFs_in_4hrDEGProms_all_Fish_df$id <- TFs_inTargetProms
  
  lnc_DEGs_Fish_list[[j]] <- TFs_in_4hrDEGProms_all_Fish_df
}

names(lnc_DEGs_Fish_list) <- c("04hrUp", "04hrDown", "48hrUp", "48hrDown", "824hrUp", "824hrDown")

lnc_DEGs_Fish_list_df <- as.data.frame((bind_rows(lnc_DEGs_Fish_list, .id = "cluster")))


#correct p val
lnc_DEGs_Fish_list_correct <- list()

#remove TFs which are not expressed or classed as TF (Lambert):
#TFsExpressed <- filter(fpkm_allG, grepl("TF", GeneClassUpdate))
#TFsExpressed <- unique(c(TFsExpressed$EnsName, 
#                         gsub(" ", "", unlist(strsplit(filter(TFsExpressed, !is.na(AllNames))$AllNames, ",")))))
#could also remove anything not in the Lambert list too..

for(j in 1:length(lnc_DEGs_Fish_list)){
  
  #remove any TFs which do not target at least 25 genes (based on all G)
  #180 is the point where more than 2 hits start appearing
  
  trial <- filter(lnc_DEGs_Fish_list[[j]], V3 >25)
  trial$BH_p <- p.adjust(trial$V5, method = "BH")
  lnc_DEGs_Fish_list_correct[[j]] <- trial
}

names(lnc_DEGs_Fish_list_correct) <- c("04hrUp", "04hrDown", "48hrUp", "48hrDown", "824hrUp", "824hrDown")

lnc_DEGs_Fish_list_correct <- as.data.frame((bind_rows(lnc_DEGs_Fish_list_correct, .id = "cluster")))

#HMGB2 does not survive MHT
#could subset to just regulators of certain function? hmmm
#or use the results from all DEGs then apply to lncRNAs to support hypothesis

#enrichments over time:
lnc_DEGs_Fish_list_correct_sig <- filter(lnc_DEGs_Fish_list_correct, BH_p < 0.1)

table(lnc_DEGs_Fish_list_correct_sig$cluster)

#no strong enrichments after multiple hypothesis correcting - FPKM >20 TFs only means can get HMGB2 to 0.08-0.09

#write.csv(lnc_DEGs_Fish_list_correct, "lnc_DEGs_Fish_df_byGene.csv", row.names = F)


#### HMGB2 bound lncRNAs ####

#interesting signal - role in DNA flexibility

EarlylncsTFs <- allG_remap2022_TFs_25FPKM[allG_remap2022_TFs_25FPKM$EnsID %in% all_DEGs_list[[1]]$EnsID &
                                           allG_remap2022_TFs_25FPKM$EnsID %in% filter(fpkm_allG, grepl("Lnc|fide", GeneClassUpdate))$EnsID,]

EarlylncsTFs_HMGB2 <- filter(EarlylncsTFs, peakTF == "HMGB2")

#any cclncs?

#### Influential partners of cis-acting lncs ####
library(dplyr)

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

CoRegPairs_04_48_24_extended <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CoRegPairs_04_48_24_extendedIII.csv")

length(unique(CoRegPairs_04_48_24_extended$EnsID))

CoRegPairs_04_48_24_extended_T2 <- filter(CoRegPairs_04_48_24_extended, corSig == "Yes" | 
                                            eQTL_supported == "Yes" | 
                                            !loopMethod == "Neither")

length(unique(CoRegPairs_04_48_24_extended$EnsID))#115
length(unique(CoRegPairs_04_48_24_extended$EnsID.y))#236

length(unique(CoRegPairs_04_48_24_extended_T2$EnsID))#82
length(unique(CoRegPairs_04_48_24_extended_T2$EnsID.y))#125


CoRegPairs_04_48_24_extended_TF <- filter(CoRegPairs_04_48_24_extended, GeneClassUpdate.y %in% c("TF", "TF + CC"))
length(unique(CoRegPairs_04_48_24_extended_TF$EnsName.y))#32 TFs

CoRegPairs_04_48_24_extended_T2_TF <- filter(CoRegPairs_04_48_24_extended_T2, GeneClassUpdate.y %in% c("TF", "TF + CC"))
length(unique(CoRegPairs_04_48_24_extended_T2_TF$EnsName.y))#22 TFs


# add influential info:
HubGenesMM <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/HubGenesMMall.csv")
ISMARA_influential <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/ISMARA_influential.csv")
LISA_influential <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/LISA_influential.csv")


# add hub info for lncs:
trial <- merge(CoRegPairs_04_48_24_extended, HubGenesMM[,c(1,6:8,11)], by.x = "EnsID", by.y = "EnsID", all.x = T)
colnames(trial)[34:37] <- paste("Lnc_", colnames(trial)[34:37], sep = "")
# add hub info for PCGs:
trial <- merge(trial, HubGenesMM[,c(1,6:8,11)], by.x = "EnsID.y", by.y = "EnsID", all.x = T)
colnames(trial)[38:41] <- paste("PCG_", colnames(trial)[38:41], sep = "")
# add ISMARA info:
trial <- merge(trial, ISMARA_influential, by.x = "EnsName.y", by.y = "Motif", all.x = T)
#add LISA info:
trial <- merge(trial, LISA_influential[,c(1,8:14)], by.x = "EnsName.y", by.y = "TF", all.x = T)
CoRegPairs_04_48_24_extended_inf <- trial

CoRegPairs_04_48_24_extended_inf_T2 <- filter(CoRegPairs_04_48_24_extended_inf, corSig == "Yes" | 
                                                eQTL_supported == "Yes" | 
                                                !loopMethod == "Neither")

length(unique(CoRegPairs_04_48_24_extended_inf$EnsID))#115
length(unique(CoRegPairs_04_48_24_extended_inf$EnsID.y))#236

length(unique(CoRegPairs_04_48_24_extended_inf_T2$EnsID))#82
length(unique(CoRegPairs_04_48_24_extended_inf_T2$EnsID.y))#125


#### HMGB2 amongst cclncs ####

HMGB2_CClncPairs <- filter(CoRegPairs_04_48_24_extended_inf, EnsID %in% EarlylncsTFs_HMGB2$EnsID)

length(unique(HMGB2_CClncPairs$EnsID))

9/14 #most are cclncs

#strong cclncs?
HMGB2_CClncPairs_T2 <- filter(CoRegPairs_04_48_24_extended_inf_T2, EnsID %in% EarlylncsTFs_HMGB2$EnsID)

length(unique(HMGB2_CClncPairs_T2$EnsID))

7/14
7/9
#78% of these are tier 2 (HiC/eQTL confirmed)

#targeting influential PCGs
HMGB2_CClncPairs_inf <- filter(CoRegPairs_04_48_24_extended_inf, PCG_moduleOfInterest == "Yes" | 
                                 Lnc_moduleOfInterest == "Yes" | 
                                 !is.na(Zscore) | !is.na(bestRank),
                               EnsID %in% EarlylncsTFs_HMGB2$EnsID)

length(unique(HMGB2_CClncPairs_inf$EnsID))
5/14
#36% are connected to predicted influential genes

#includes MSTRG.12915 (FOXL1) as well as AC002480.4 (adjacent to IPL-IL6) and GMEB1 (not taken to lab in the end)


#how many targets are also HMGB2-bound?

DEGs_HMGB2 <- allG_remap2022_TFs_25FPKM[allG_remap2022_TFs_25FPKM$peakTF == "HMGB2",]

HMGB2_CClncPairs_doubleHMGB2 <- filter(HMGB2_CClncPairs, EnsID.y %in% DEGs_HMGB2$EnsID)

#only FOXL1 and HAS2

#significance? evidence that both can be bound and activated by HMGB2?



#### Quick GPS2 check ####

EarlylncsTFs_GPS2 <- filter(EarlylncsTFs, peakTF == "GPS2")

GPS2_CClncPairs <- filter(CoRegPairs_04_48_24_extended_inf, EnsID %in% EarlylncsTFs_GPS2$EnsID)

length(unique(GPS2_CClncPairs$EnsID))

7/11 #most are cclncs

#strong cclncs?
GPS2_CClncPairs_T2 <- filter(CoRegPairs_04_48_24_extended_inf_T2, EnsID %in% EarlylncsTFs_GPS2$EnsID)

length(unique(GPS2_CClncPairs_T2$EnsID))

7/11
7/7
#100% of these are tier 2 (HiC/eQTL confirmed)

#targeting influential PCGs
GPS2_CClncPairs_inf <- filter(CoRegPairs_04_48_24_extended_inf, PCG_moduleOfInterest == "Yes" | 
                                 Lnc_moduleOfInterest == "Yes" | 
                                 !is.na(Zscore) | !is.na(bestRank),
                               EnsID %in% EarlylncsTFs_GPS2$EnsID)

length(unique(GPS2_CClncPairs_inf$EnsID))
2/14
#14%, low percentage are connected to predicted influential genes

#includes CXCL8-INKILN


#how many targets are also GPS2-bound?

DEGs_GPS2 <- allG_remap2022_TFs_25FPKM[allG_remap2022_TFs_25FPKM$peakTF == "GPS2",]

GPS2_CClncPairs_doubleGPS2 <- filter(GPS2_CClncPairs, EnsID.y %in% DEGs_GPS2$EnsID)

#3x genes CXCL1, HIVEP2, KLF7

#significance? evidence that both can be bound and activated by HMGB2?

#both GSP2 and HMGB2?
HMGB2_GPS2_CClncPairs <- filter(CoRegPairs_04_48_24_extended_inf, EnsID %in% c(EarlylncsTFs_HMGB2$EnsID, EarlylncsTFs_GPS2$EnsID))

length(unique(HMGB2_GPS2_CClncPairs$EnsID)) 

#total of 13

length(unique(HMGB2_CClncPairs$EnsID))#9 
length(unique(GPS2_CClncPairs$EnsID)) #7


#so they co-localise at 3 loci, mostly distinct


##### characterise/visualisation ####


#number per cluster:
table(all_DEGs_Fish_df_correct_sig$cluster)

#GO terms per TF set:

clusters2check <- c("4hrUp", "4hrDown", "8hrUp", "8hrDown", "24hrUp", "24hrDown")

all_DEGs_Fish_df_correct_sig_GO2 <- list()

for(i in 1:length(clusters2check)){
  all_DEGs_Fish_df_correct_sig_GO2[[i]] <- 
    enrichGO(gene         = filter(all_DEGs_Fish_df_correct_sig, cluster %in% clusters2check[i])$id,
             universe      = unique(remap_SVSMC_prom_GR$id),
             keyType       = "SYMBOL",
             OrgDb         = org.Hs.eg.db,
             ont           = "all",
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.05,
             readable      = TRUE)
}

#saveRDS(all_DEGs_Fish_df_correct_sig_GO2, "TFenriched_DEGs_byPeak_GO.rds")

View(as.data.frame(all_DEGs_Fish_df_correct_sig_GO2[[1]])) #miRNA metabolism, inflamm, defence response, apopt
View(as.data.frame(all_DEGs_Fish_df_correct_sig_GO2[[2]])) #empty
View(as.data.frame(all_DEGs_Fish_df_correct_sig_GO2[[3]])) #empty
View(as.data.frame(all_DEGs_Fish_df_correct_sig_GO2[[4]])) #empty
View(as.data.frame(all_DEGs_Fish_df_correct_sig_GO2[[5]])) #2 terms
View(as.data.frame(all_DEGs_Fish_df_correct_sig_GO2[[6]])) #miRNA metabolism, 

#strong concerted bio signal amongst 4hr up, not amongst others
#chromatin remodelling and organisation is not necessarily specific signature to any timeframe cluster
#(though the flavour of genes would maybe change)

dotplot(simplify(all_DEGs_Fish_df_correct_sig_GO[[1]]), showCategory = 30) + theme(axis.text.y = element_text(size = 10))
dotplot(all_DEGs_Fish_df_correct_sig_GO[[1]], showCategory = 30) + theme(axis.text.y = element_text(size = 10))


#top20 p per cluster:
trial <- split(all_DEGs_Fish_df_correct_sig, all_DEGs_Fish_df_correct_sig$cluster)

triali <- lapply(trial, function(x){
  x <- x[order(x$BH_p, decreasing = F),][1:20,8]
})

#they're all super significant, heatmap of p

trialii <- filter(all_DEGs_Fish_df_correct_sig, id %in% unlist(triali))
trialii$BH_p <- -log10(trialii$BH_p)

trialiii <- trialii[,c(1,8:9)] %>% 
  pivot_wider(names_from = id, values_from = BH_p, values_fill = 0)

trialiii <- as.data.frame(trialiii)
rownames(trialiii) <- trialiii$cluster
trialiii <- trialiii[,-1]

trialiv <- as.data.frame(t(trialiii))
trialiv <- trialiv[,c(1,3,5,2,4,6)]

#cap at -log10 p of 12 (3rd Q max is 8):
trialv <- as.matrix(trialiv)
max(trialv)
summary(trialv)
trialv[trialv >12] <- 12

myColor <- colorRampPalette(c("white", "red"))(50)
myBreaks <- seq(min(as.matrix(trialv)), max(as.matrix(trialv)), length.out=50)

pheatmap::pheatmap(trialv, #scale = "row", 
                   fontsize_row = 8, cluster_cols = F,
                   color = myColor, breaks = myBreaks)


#4hr up includes SMARCA4 and some other chrom remodellers like BICRA, MAML, SUPT16H

#Also MED1, CTCF

#no BRD4 or KMT2A


#just chromatin remodelling TFs:
#can track from GO terms
chromRemod <- unique(unlist(lapply(all_DEGs_Fish_df_correct_sig_GO2, function(x){
  unlist(strsplit(filter(as.data.frame(x), grepl("chromatin remodeling", Description))$geneID, "/"))
}
)))


#just chromatin organisation TFs:
#can track from GO terms
chromOrg <- unique(unlist(lapply(all_DEGs_Fish_df_correct_sig_GO2, function(x){
  unlist(strsplit(filter(as.data.frame(x), grepl("chromatin organization", Description))$geneID, "/"))
}
)))

#chromOrg is a slightly more encompassing term (includes BRD4, MED)
sum(chromOrg %in% chromRemod)

trialii <- filter(all_DEGs_Fish_df_correct_sig, id %in% chromOrg)
trialii$BH_p <- -log10(trialii$BH_p)

trialiii <- trialii[,c(1,8:9)] %>% 
  pivot_wider(names_from = id, values_from = BH_p, values_fill = 0)

trialiii <- as.data.frame(trialiii)
rownames(trialiii) <- trialiii$cluster
trialiii <- trialiii[,-1]

trialiv <- as.data.frame(t(trialiii))
trialiv <- trialiv[,c(1,3,5,2,4,6)]

#cap at -log10 p of 12 (3rd Q max is 8):
trialv <- as.matrix(trialiv)
max(trialv)
summary(trialv)
trialv[trialv >12] <- 12

myColor <- colorRampPalette(c("white", "red"))(50)
myBreaks <- seq(min(as.matrix(trialv)), max(as.matrix(trialv)), length.out=50)

pheatmap::pheatmap(trialv, #scale = "row", 
                   fontsize_row = 8, cluster_cols = F,
                   color = myColor, breaks = myBreaks)


#just histone binding TFs:
#can track from GO terms
histBind <- unique(unlist(lapply(all_DEGs_Fish_df_correct_sig_GO2, function(x){
  unlist(strsplit(filter(as.data.frame(x), grepl("histone binding", Description))$geneID, "/"))
}
)))

#chromOrg is a slightly more encompassing term (includes BRD4, MED)
sum(chromOrg %in% histBind)

#32/35 are in the chromOrg list

trialii <- filter(all_DEGs_Fish_df_correct_sig, id %in% histBind)
trialii$BH_p <- -log10(trialii$BH_p)

trialiii <- trialii[,c(1,8:9)] %>% 
  pivot_wider(names_from = id, values_from = BH_p, values_fill = 0)

trialiii <- as.data.frame(trialiii)
rownames(trialiii) <- trialiii$cluster
trialiii <- trialiii[,-1]

trialiv <- as.data.frame(t(trialiii))
trialiv <- trialiv[,c(1,3,5,2,4,6)]

#cap at -log10 p of 12 (3rd Q max is 8):
trialv <- as.matrix(trialiv)
max(trialv)
summary(trialv)
trialv[trialv >12] <- 12

myColor <- colorRampPalette(c("white", "red"))(50)
myBreaks <- seq(min(as.matrix(trialv)), max(as.matrix(trialv)), length.out=50)

pheatmap::pheatmap(trialv, #scale = "row", 
                   fontsize_row = 8, cluster_cols = F,
                   color = myColor, breaks = myBreaks)


#just co-activators
coactivator <- unique(unlist(lapply(all_DEGs_Fish_df_correct_sig_GO2, function(x){
  unlist(strsplit(filter(as.data.frame(x), grepl("coactivator activity", Description))$geneID, "/"))
}
)))

#chromOrg is a slightly more encompassing term (includes BRD4, MED)
sum(chromOrg %in% coactivator)
sum(coactivator %in% histBind)

19/27 #are in the chromOrg list
8/27 #in the histBind list

trialii <- filter(all_DEGs_Fish_df_correct_sig, id %in% coactivator)
trialii$BH_p <- -log10(trialii$BH_p)

trialiii <- trialii[,c(1,8:9)] %>% 
  pivot_wider(names_from = id, values_from = BH_p, values_fill = 0)

trialiii <- as.data.frame(trialiii)
rownames(trialiii) <- trialiii$cluster
trialiii <- trialiii[,-1]

trialiii <- rbind(trialiii, rep(0, 15), rep (0,15))
rownames(trialiii)[5:6] <- c("8hrUp", "24hrUp")

trialiv <- as.data.frame(t(trialiii))
trialiv <- trialiv[,c(1,5,6,2,3,4)]

#cap at -log10 p of 12 (3rd Q max is 8):
trialv <- as.matrix(trialiv)
max(trialv)
summary(trialv)
trialv[trialv >12] <- 12

myColor <- colorRampPalette(c("white", "red"))(50)
myBreaks <- seq(min(as.matrix(trialv)), max(as.matrix(trialv)), length.out=50)

pheatmap::pheatmap(trialv, #scale = "row", 
                   fontsize_row = 8, cluster_cols = F,
                   color = myColor, breaks = myBreaks)


#co-repressor
corepressor <- unique(unlist(lapply(all_DEGs_Fish_df_correct_sig_GO2, function(x){
  unlist(strsplit(filter(as.data.frame(x), grepl("corepressor activity", Description))$geneID, "/"))
}
)))

#chromOrg is a slightly more encompassing term (includes BRD4, MED)
sum(chromOrg %in% corepressor)
sum(corepressor %in% histBind)
sum(corepressor %in% coactivator)

11/16# are in the chromOrg list
5/16 #in the histBind list

trialii <- filter(all_DEGs_Fish_df_correct_sig, id %in% corepressor)
trialii$BH_p <- -log10(trialii$BH_p)

trialiii <- trialii[,c(1,8:9)] %>% 
  pivot_wider(names_from = id, values_from = BH_p, values_fill = 0)

trialiii <- as.data.frame(trialiii)
rownames(trialiii) <- trialiii$cluster
trialiii <- trialiii[,-1]

trialiii <- rbind(trialiii, rep(0, 16))#, rep (0,15))
rownames(trialiii)[6] <- c("24hrUp")

trialiv <- as.data.frame(t(trialiii))
trialiv <- trialiv[,c(1,3,6,2,4,5)]

#cap at -log10 p of 12 (3rd Q max is 8):
trialv <- as.matrix(trialiv)
max(trialv)
summary(trialv)
trialv[trialv >12] <- 12

myColor <- colorRampPalette(c("white", "red"))(50)
myBreaks <- seq(min(as.matrix(trialv)), max(as.matrix(trialv)), length.out=50)

pheatmap::pheatmap(trialv, #scale = "row", 
                   fontsize_row = 8, cluster_cols = F,
                   color = myColor, breaks = myBreaks)

  
##### to do #####

#change colour scale
#consider other terms like "chrom organisation" "SWI/SNF" alongside
#ensure consistent with lit
