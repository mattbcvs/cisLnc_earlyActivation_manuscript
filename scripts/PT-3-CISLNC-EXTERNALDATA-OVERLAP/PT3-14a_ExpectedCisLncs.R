#A test set of lncRNAs with described or strong evidence suggesting cis-acting lncRNA 

#a) TransCistor - finds lncRNAs which when knocked down have a bias to hit surrounding genes
#b) TransCistor literature search - lncRNAs they found in lit with previously described with cis (not a good, detailed or reliable table so needs double-checking)
#c) Fanucchi et al.(UMLILO) - lncRNAs which when knocked down have a role in changing surrounding chromatin near immune genes
#d) ASO iPSC + HDF study - lncRNAs with a consistent molecular phenotype, do their lncRNA-controlled genes match SVSMC co-reg neighbours?
#e) ASO iPSC + HDF study - lncRNAs with a consistent mol. pheno. that are identified as cis-acting through overlap with HiC (18x iPSC, 1x HDF)
#f) low throughput lit: e.g. AMANZI(not annotated in GEN - or de novo)
#g) other sources as they become known

#be aware that most of this will be based on ASO in iPSC or HDF i.e. unstimulated

library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(rcompanion)
library(ggbeeswarm)


#From TransCistor:
Trancistor <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/TransCistor.csv", skip = 2)
Trancistor <- filter(Trancistor, !File.Name == "", Species == "Human")

#identify their cis lncs, FDR of 0.25:
Trancistor_cis <- filter(Trancistor, Biobase::rowMin(as.matrix(Trancistor[,12:15])) <= 0.25)
unique(Trancistor_cis$Symbol)
unique(Trancistor_cis$ENSEMBL.ID)


#from HDF only one: RP11-398K22.12
"ENSG00000229852"


#From Yip, 18x lncs
ipsc_cis <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/IpSC_Cis_compare_S4.csv")
length(unique(ipsc_cis$target.geneID))
length(unique(ipsc_cis$target.geneID))

ipsc_cis_concordant <- filter(ipsc_cis, molecular.phenotype.concordant == "yes")

length(unique(ipsc_cis_concordant$target.geneID))
length(unique(ipsc_cis_concordant$target.geneID))#32 (35 reported in paper)

ipsc_cis_concordant_DEGs <- filter(ipsc_cis, abs(EdgeR.Zscore) >1.645, EdgeR.FDR <0.1, molecular.phenotype.concordant == "yes")
length(unique(ipsc_cis_concordant_DEGs$target.geneID))#18 (18 reported in paper)


#Agrawal et al. 2024, Nuclear annotation of lncRNAs
Agrawal_HiC_ASO <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Agrawal_HiC_ASO_lncs.csv")
length(unique(Agrawal_HiC_ASO$Nuclear.lncRNA))#81 lncRNAs...not 83
Agrawal_HiC_ASO_targets <- filter(Agrawal_HiC_ASO, Is.candidate.target.gene == "Yes")
length(unique(Agrawal_HiC_ASO_targets$Nuclear.lncRNA))#67 have one DEG and HiC link
#number of ASO hit HiC targets:
table(sapply(lapply(split(Agrawal_HiC_ASO_targets, Agrawal_HiC_ASO_targets$Nuclear.lncRNA), dim), "[[", 1))
#27 are less than 3, 40x have 3x DEG/HiC targets, not 33...

#table contains matches TO THE SAME LNCRNA GENE...needs some cautious filtering:
Agrawal_HiC_ASO_targetsDistinct <- filter(Agrawal_HiC_ASO_targets, !Nuclear.lncRNA == Differentially.Expressed.gene..geneID.)
#unique lnc to PCG pairs (exclude diff ASO pairs)
Agrawal_HiC_ASO_targetsDistinct <- Agrawal_HiC_ASO_targetsDistinct[,c(1,4,5)]
Agrawal_HiC_ASO_targetsDistinct <- unique(Agrawal_HiC_ASO_targetsDistinct)
length(unique(Agrawal_HiC_ASO_targetsDistinct$Nuclear.lncRNA))#47 have one DEG and HiC link to an alternative gene
table(sapply(lapply(split(Agrawal_HiC_ASO_targetsDistinct, Agrawal_HiC_ASO_targetsDistinct$Nuclear.lncRNA), dim), "[[", 1))
#now managed to get to the 33 lncRNAs...

#should they be used though? many ignored in Ramilowski as no second ASO to back up
#however, they do have HiC to connect them which is some validation
#it includes HMGA2-RP11 so is a nice bit of overlap with my work too
trial <- sapply(lapply(split(Agrawal_HiC_ASO_targetsDistinct, Agrawal_HiC_ASO_targetsDistinct$Nuclear.lncRNA), dim), "[[", 1)
Agrawal_cisLncs <- names(trial)[trial >2]


#Fanucchi, proven effect on locus and WDR5,histone etc (UMLILO already covered)
"ENSG00000232949" #IPL-IL6 aka AC002480.3
"ENSG00000232759" #AC002480.4 which is upstream similar exons, possibly another isoform
"ENSG00000228420" #IPL-CSF1
"ENSG00000225331" #IPL-ICOSLG, LINC01768
"ENSG00000228277" #UMLILO


#Other literature
"HOTAIR"


#all loci:
ControlCisLncs <- list("TransCistor" = unique(Trancistor_cis$ENSEMBL.ID),
                       "HDF_paper" = "ENSG00000229852",
                       "Agrawal" = Agrawal_cisLncs,
                       "Yip" = unique(ipsc_cis_concordant_DEGs$target.geneID),
                       "Fanucchi" = c("ENSG00000232949",
                                      "ENSG00000232759",
                                      "ENSG00000228420",
                                      "ENSG00000225331",
                                      "ENSG00000228277"),
                       #these mostly from Igor's review after a quick check of their papers
                       #must have knockdown and specific regulatory effect on neighbour via a cis-mechanism (recruitment/looping etc) established
                       "OtherLit" = c("ENSG00000229807",#XIST
                                      "ENSG00000228630",#HOTAIR
                                      "ENSG00000243766",#HOTTIP
                                      "ENSG00000281358"#RASSF1-AS1
                       ) 
)

ControlCisLncs <- bind_rows(lapply(ControlCisLncs, as.data.frame), .id = "source")
colnames(ControlCisLncs) <- c("source", "Ens_ID")
#write.csv(ControlCisLncs, "ControlCisLncs.csv")


#ControlCisLncs <- unique(c(Trancistor_cis$ENSEMBL.ID,
#                           ipsc_cis_concordant_DEGs$target.geneID,
#                           "ENSG00000232949", "ENSG00000232759", "ENSG00000228420", "ENSG00000225331"))
#write.csv(ControlCisLncs, "ControlCisLncs.csv")


##### SVSMC table import ####

fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGClassUpdate.csv", header = T)

FPKM_CQV_OVERLAP_fpkm <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/FPKM_CQV_OVERLAP_fpkm.csv")
table(FPKM_CQV_OVERLAP_fpkm$IGV)#413 pass, 168 fail

#remove artefacts from here already
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


#from neighbour associations, all lncRNA-PCGs within a 250kbp window
AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_1.csv", header = T)

# no AS artefacts, no PLAR only "PCGs", must have ENCODE annotation
AllLNC_AllPCG_1$pairs <- paste(AllLNC_AllPCG_1$EnsID, AllLNC_AllPCG_1$EnsID.y, sep="-")

AllLNC_AllPCG_1 <- filter(AllLNC_AllPCG_1, EnsID %in% fpkm_allG$EnsID, EnsID.y %in% fpkm_allG$EnsID) #removes about 800
AllLNC_AllPCG_1 <- filter(AllLNC_AllPCG_1, EnsID.y %in% filter(fpkm_allG, EnsType == "protein_coding")$EnsID) #removes about 5
#1653 pairings of lnc-PCG expressed near to each other


##### TranCistor import ####

#TransCistor:
Trancistor <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/TransCistor.csv", skip = 2)
Trancistor <- filter(Trancistor, !File.Name == "", Species == "Human")

#identify their cis lncs, FDR of 0.25:
Trancistor_cis <- filter(Trancistor, Biobase::rowMin(as.matrix(Trancistor[,12:15])) <= 0.25)

#TransCistor also has a "curated" group of cis-acting lncRNAs in supp table...
#not convinced tho as MYOSLID reported - source is a review which links to a paper (Baker group) which states it is not cis-acting
#therefore the curation is highly suspect
#some are reported as cis-acting via "another perturbation" - but not clear what this means
#in principle it makes sense though - look for cis activity amongst the found genes

#identify any in ours, by symbol (TransCistor has used ENSEMBL thankfully):
unique(filter(fpkm_allGDE, EnsName %in% Trancistor$Symbol)[,c(1,2)]) #32x DELs are assessed by TransCistor
unique(filter(fpkm_allGDE, EnsName %in% Trancistor_cis$Symbol)[,c(1,2,50)]) #4x DELs are detected as cis-acting by TransCistor

#ID might also be a good way to merge?
fpkm_allGDE$Ens_ID_merge <- gsub("\\.[0-9]*", "", fpkm_allGDE$EnsID)
unique(filter(fpkm_allGDE, EnsName %in% Trancistor$Symbol | Ens_ID_merge %in% Trancistor$ENSEMBL.ID)[,c(1,2)]) #39x DELs are assessed by TransCistor
unique(filter(fpkm_allGDE, EnsName %in% Trancistor_cis$Symbol | Ens_ID_merge %in% Trancistor_cis$ENSEMBL.ID)[,c(1,2,50)]) #5x DELs are called cis-acting by TransCistor

#by ID, 5x lncRNAs with predicted cis-acting activity are seen to be DEL in SVSMC


##### FANTOM ASO import ####

#HDF ASO screen, significant DEGs: https://fantom.gsc.riken.jp/6/datafiles/Core_FANTOM6/RELEASE_latest/analysis/DEGs/
HDF_ASO <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Ram_lncRNA_ASO_results.tsv", )
dim(HDF_ASO)#returns all 340 ASOs
colnames(HDF_ASO)[1:20]#against a set of genes
length(unique(HDF_ASO$geneID))
HDF_ASO[1:5,1:5]
#e.g. for Gapmer1, looks like 464 repressed, 385 induced
table(HDF_ASO[,3])
#search for them here to compare to more info:
#https://fantom.gsc.riken.jp/zenbu/reports/#FANTOM6_DESeq_ASO

#total of 340 high efficiency ASOs for 154 lncRNAs, from table S2 (sub-cellular fractionation also available in another S2 sub-table):
HDF_ASO_lnc_ID <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Ram_lncRNA_ASO_match.csv")
#154 ASO-targeted lncRNAs assayed for a molecular fingerprint identified thru CAGEseq in HDFs:
HDF_lncSymbols <- unique(filter(HDF_ASO_lnc_ID, perturb_id %in% colnames(HDF_ASO))$KD.geneSymbol)
table((filter(HDF_ASO_lnc_ID, perturb_id %in% colnames(HDF_ASO))$KD.geneSymbol))

#their approach is to find >5 common genes between ASO pairs - and make sure the DEG similarity is above random pairs
#this leaves only 13x lncRNAs
13/119 #10.9% of those assayed had a reproducible molecular phenotype

#treating others as untrustworthy

#will take some code to figure out DEGs per lnc but do-able
#ASO pairs are listed in a supp table so don't need to figure out again
#Pairs are in supp table 5:
HDF_ASO_pairs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Ram_lncRNA_ASO_pairs.csv")
length(unique(HDF_ASO_pairs$KD.geneSymbol))#119 lncRNAs with at least 2x working ASOs
#I added this column based on the paper to get to the 13x lncRNAs:
HDF_ASO_pairs_rep <- filter(HDF_ASO_pairs, Reproducible == "yes")
HDF_lncSymbols_rep <- unique(HDF_ASO_pairs_rep$KD.geneSymbol)


#iPSC ASO screen (123 lncRNAs tested, excluded ~95 lncRNAs without reproducible knockdown): Table S3:
Yip_ASO <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Yip_lncRNA_ASO_results.csv")
iPSC_lncID <- unique(Yip_ASO$target.geneID)#28 lncRNAs which have a consistent transcriptomic effect from 2x ASOs
length(unique(Yip_ASO$DEG.HGNC))#707 confirmed affected targets (hit by 2 ASOs)

#less stringent iPSC table: https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2022/data/processed_raw_data/


#confirm if the 13 HDF and 28 iPSC lncs are in SVSMC
#SVSMC lncs already matched-up to FANTOM ID
Enhancer_lociII <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Enhancer_lociIItime.csv", header = T)
length(unique(Enhancer_lociII$EnsID))#selecting all lncs = 558 genes, if just enhancer = 77 (7/2021)
#remove AS overlap artefacts:
Enhancer_lociII <- filter(Enhancer_lociII, EnsID %in% fpkm_allG$EnsID)

#need a FANTOM symbol column to compare to SVSMC:
allLncs_BestCAGE_NameFinder <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/FANTOMNameFinder.csv", stringsAsFactors = F,
                                          header = F, sep = "\t")

#contains an ID for all the HDF/iPSC lncs?
sum(HDF_lncSymbols_rep %in% gsub(";", "", allLncs_BestCAGE_NameFinder$V2))# yes
sum(iPSC_lncID %in% gsub("\\.[0-9]*;", "", allLncs_BestCAGE_NameFinder$V1))# yes

#need a table to convert between the mess of IDs - get a stable ID between HDF, iPSC and the SVSMC data:
allLncs_BestCAGE_NameFinder$symbol <- gsub(";", "", allLncs_BestCAGE_NameFinder$V2)
allLncs_BestCAGE_NameFinder$fantIDmerge <- gsub("\\.[0-9]*;", "", allLncs_BestCAGE_NameFinder$V1)
allLncs_BestCAGE_NameFinder$fantID <- gsub(";", "", allLncs_BestCAGE_NameFinder$V1)

#39x lncRNAs with repeatable ASO profile in one or other screen

SVSMC_lncRNAs_wKD <- filter(Enhancer_lociII, FANTOM_ID %in% ASO_lncRNAs$fantID)
length(unique(SVSMC_lncRNAs_wKD$EnsID)) #17 of them expressed in the SVSMC timecourse

#10x are DELs
SVSMC_DElncRNAs_wKD <- unique(filter(SVSMC_lncRNAs_wKD, EnsID %in% fpkm_allGDE$EnsID))
#all ID formats for the 10 genes:
colnames(SVSMC_DElncRNAs_wKD)
SVSMC_DElncRNAs_wKD <- unique(merge(SVSMC_DElncRNAs_wKD[,c(1,14)], allLncs_BestCAGE_NameFinder[,c(3:5)], by.x = "FANTOM_ID", by.y = "fantID"))
SVSMC_DElncRNAs_wKD <- unique(merge(SVSMC_DElncRNAs_wKD, fpkm_allG[,c(2,5)], by = "EnsID"))

#add stable ID column (FANTOM symbol) to HDF table of repeatable ASO pairs:
#5 of the 13 genes are DELs in SV data 
HDF_ASO_pairs_repSV <- merge(SVSMC_DElncRNAs_wKD, HDF_ASO_pairs_rep, by.y = "KD.geneSymbol", by.x = "symbol")
unique(HDF_ASO_pairs_repSV$symbol)

#add stable ID column to Yip table:
Yip_ASO_repSV <- merge(SVSMC_DElncRNAs_wKD, Yip_ASO, by.y = "target.geneID", by.x = "fantIDmerge")
unique(Yip_ASO_repSV$symbol)

#10x genes to check out, do they have a consistently-ASO affected gene that is seen in the SVSMC data?
unique(c(HDF_ASO_pairs_repSV$symbol, Yip_ASO_repSV$symbol))


#HDF lncRNAs-affected genes first
#identify columns for each pair in the HDF table:
unique(HDF_ASO_pairs_repSV$pair_id)

#their naming schema changes in the ASO table, need to match it back:
HDF_ASO_lnc_ID_targets <- filter(HDF_ASO_lnc_ID, Description == "Targeted lncRNA")
HDF_ASO_lnc_ID_targets <- filter(HDF_ASO_lnc_ID_targets,  KD.geneSymbol %in% HDF_ASO_pairs_repSV$symbol)
HDF_ASO_lnc_ID_targets$GapNo <- sapply(strsplit(HDF_ASO_lnc_ID_targets$perturb_id, "_"), "[[", 3)
#need slightly diff section for the "AD" ASOs
HDF_ASO_lnc_ID_targets$GapNo[grepl("_AD_", HDF_ASO_lnc_ID_targets$perturb_id)] <- 
  paste(sapply(strsplit(HDF_ASO_lnc_ID_targets$perturb_id[grepl("_AD_", HDF_ASO_lnc_ID_targets$perturb_id)], "_"), "[[", 3),
        sapply(strsplit(HDF_ASO_lnc_ID_targets$perturb_id[grepl("_AD_", HDF_ASO_lnc_ID_targets$perturb_id)], "_"), "[[", 4), sep = "_")

HDF_ASO_lnc_ID_targets$GapName <- paste(HDF_ASO_lnc_ID_targets$KD.geneSymbol, HDF_ASO_lnc_ID_targets$GapNo, sep = "-ASO_")
HDF_ASO_lnc_ID_targets$GapName

#take the second ASO, add the gene ID to get a proper usable ID to match up pairs
#use pair_id column NOT the other one where the gene symbol changes...
HDF_ASO_pairs_repSV$pair_idmerge <- sapply(split(HDF_ASO_pairs_repSV, HDF_ASO_pairs_repSV$pair_id), function(x){
  gsub("\\|", paste("\\|", x$symbol, "-", sep = ""), x$pair_id)
})
HDF_ASO_pairs_repSV$partnerA <- sapply(strsplit(HDF_ASO_pairs_repSV$pair_idmerge, "\\|"), "[[", 1)
HDF_ASO_pairs_repSV$partnerB <- sapply(strsplit(HDF_ASO_pairs_repSV$pair_idmerge, "\\|"), "[[", 2)

#add the ASO ID, what a faff this is:
HDF_ASO_pairs_repSV <- merge(HDF_ASO_pairs_repSV, HDF_ASO_lnc_ID_targets[,c(5,10)], by.x = "partnerA", by.y = "GapName", all.x = T)
HDF_ASO_pairs_repSV <- merge(HDF_ASO_pairs_repSV, HDF_ASO_lnc_ID_targets[,c(5,10)], by.x = "partnerB", by.y = "GapName", all.x = T)

#now know which columns in the big table refer to which ASO
#for each pair, find common up or down genes, e.g. first 2:
HDF_ASO[,1][HDF_ASO[,3] == 1 & HDF_ASO[,4] == 1] #matches first pair in table 5 for common up
HDF_ASO[,1][HDF_ASO[,3] == -1 & HDF_ASO[,4] == -1] #matches for common down

#for our pairs of interest:
ASO_1 <- HDF_ASO_pairs_repSV$perturb_id.x[1]
table(HDF_ASO[,ASO_1])
ASO_2 <- HDF_ASO_pairs_repSV$perturb_id.y[1]

HDF_ASO[,1][HDF_ASO[,ASO_1] == 1 & HDF_ASO[,ASO_2] == 1] #matching up genes for pair1
HDF_ASO[,1][HDF_ASO[,ASO_1] == -1 & HDF_ASO[,ASO_2] == -1] #matching down genes for pair1

HDF_DEGs_pair <- list()

for (i in 1:dim(HDF_ASO_pairs_repSV)[1]){
  ASO_1 <- HDF_ASO_pairs_repSV$perturb_id.x[i]
  ASO_2 <- HDF_ASO_pairs_repSV$perturb_id.y[i]
  
  HDF_DEGs_pair[[i]] <- list(HDF_ASO[,1][HDF_ASO[,ASO_1] == 1 & HDF_ASO[,ASO_2] == 1], #matching up genes for pair i
                             HDF_ASO[,1][HDF_ASO[,ASO_1] == -1 & HDF_ASO[,ASO_2] == -1]) #matching down genes for pair i
}

names(HDF_DEGs_pair) <- HDF_ASO_pairs_repSV$EnsName

HDF_DEGs_pair[2]

sapply(HDF_DEGs_pair[[1]], length)#match to pairs table DEG numbers
sapply(HDF_DEGs_pair[[2]], length)#match to pairs table DEG numbers
sapply(HDF_DEGs_pair[[3]], length)#match to pairs table DEG numbers
sapply(HDF_DEGs_pair[[4]], length)#match to pairs table DEG numbers
sapply(HDF_DEGs_pair[[5]], length)#match to pairs table DEG numbers

#all 5 of these lnc-regulated genes in HDF show up in the SVSMC:
HDF_ASO_pairs_repSV_CoReg <- filter(HDF_ASO_pairs_repSV, EnsID %in% AllLNC_AllPCG_1$EnsID)
HDF_ASO_pairs_repSV_CoReg$EnsName

names(HDF_DEGs_pair)

AllLNC_AllPCG_1i <- AllLNC_AllPCG_1
AllLNC_AllPCG_1i$EnsID.y.merge <- gsub("\\.[0-9]*", "", AllLNC_AllPCG_1i$EnsID.y)

#AC005592.2 check:
names(HDF_DEGs_pair)
HDF_DEGs_pair[1]
HDF_ASO_pairs_repSV_CoReg[,3:7]

filter(AllLNC_AllPCG_1i, EnsID == HDF_ASO_pairs_repSV_CoReg$EnsID[1],
       EnsID.y.merge %in% HDF_DEGs_pair[[1]][[1]])
filter(AllLNC_AllPCG_1i, EnsID == HDF_ASO_pairs_repSV_CoReg$EnsID[1],
       EnsID.y.merge %in% HDF_DEGs_pair[[1]][[2]])
#no genes in range are confirmed via ASO

#MSTRG.17405 check:
names(HDF_DEGs_pair)
HDF_DEGs_pair[2]
HDF_ASO_pairs_repSV_CoReg[,3:7]

filter(AllLNC_AllPCG_1i, EnsID == HDF_ASO_pairs_repSV_CoReg$EnsID[2],
       EnsID.y.merge %in% HDF_DEGs_pair[[2]][[1]])
filter(AllLNC_AllPCG_1i, EnsID == HDF_ASO_pairs_repSV_CoReg$EnsID[2],
       EnsID.y.merge %in% HDF_DEGs_pair[[2]][[2]])
#no genes in range are confirmed via ASO

#FGD5-AS1 check:
names(HDF_DEGs_pair)
HDF_DEGs_pair[3]
HDF_ASO_pairs_repSV_CoReg[,3:7]

filter(AllLNC_AllPCG_1i, EnsID == HDF_ASO_pairs_repSV_CoReg$EnsID[3],
       EnsID.y.merge %in% HDF_DEGs_pair[[3]][[1]])
filter(AllLNC_AllPCG_1i, EnsID == HDF_ASO_pairs_repSV_CoReg$EnsID[3],
       EnsID.y.merge %in% HDF_DEGs_pair[[3]][[2]])
#no genes in range are confirmed via ASO

#RP11-417E7.1 check:
names(HDF_DEGs_pair)
HDF_DEGs_pair[4]
HDF_ASO_pairs_repSV_CoReg[,3:7]

filter(AllLNC_AllPCG_1i, EnsID == HDF_ASO_pairs_repSV_CoReg$EnsID[4],
       EnsID.y.merge %in% HDF_DEGs_pair[[4]][[1]])
filter(AllLNC_AllPCG_1i, EnsID == HDF_ASO_pairs_repSV_CoReg$EnsID[4],
       EnsID.y.merge %in% HDF_DEGs_pair[[4]][[2]])
#no genes in range are confirmed via ASO

#RP11-422J8.1 check:
names(HDF_DEGs_pair)
HDF_DEGs_pair[5]
HDF_ASO_pairs_repSV_CoReg[,3:7]

filter(AllLNC_AllPCG_1i, EnsID == HDF_ASO_pairs_repSV_CoReg$EnsID[5]
       ,EnsID.y.merge %in% HDF_DEGs_pair[[5]][[1]]
       )
filter(AllLNC_AllPCG_1i, EnsID == HDF_ASO_pairs_repSV_CoReg$EnsID[5],
       EnsID.y.merge %in% HDF_DEGs_pair[[5]][[2]])
#no genes in range are confirmed via ASO

#no ID of any cis-acting potential for 5x SVSMC lncRNAs with data from HDF ASO screen


#iPSC
length(unique(Yip_ASO_repSV$EnsID))
Yip_ASO_repSV_CoReg <- filter(Yip_ASO_repSV, EnsID %in% CoRegPairs_04_48_24_extended$EnsID)
length(unique(Yip_ASO_repSV_CoReg$EnsID))#4x genes to check

filter(CoRegPairs_04_48_24_extended, EnsID %in% Yip_ASO_repSV_CoReg$EnsID, EnsName.y %in% Yip_ASO_repSV_CoReg$DEG.HGNC)


#7x CoRegs that have been assigned a consistent molecular profile by ASO screening in either HDF or iPSC
#the predicted target in SVSMC is not regulated in any of these screens 


##### FANTOM ASO papers- identified/predicted cis lncs ####

#from HDF only one: RP11-398K22.12(putative lncRNA)
filter(fpkm_allG, grepl("ENSG00000229852", EnsID))

#from iPSC table of 18x cis-suggestive lncRNAs
#consistent ASO (36x lncs total this time rather than 28, some subtle diff) and a DEG that is HiC contact
ipsc_cis <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/IpSC_Cis_compare_S4.csv")

length(unique(ipsc_cis$target.geneID))

ipsc_cis_concordant <- filter(ipsc_cis, molecular.phenotype.concordant == "yes")

length(unique(ipsc_cis_concordant$target.geneID))
length(unique(ipsc_cis_concordant$target.geneID))#32 (35 reported in paper)

ipsc_cis_concordant_DEGs <- filter(ipsc_cis, abs(EdgeR.Zscore) >1.645, EdgeR.FDR <0.1, molecular.phenotype.concordant == "yes")
length(unique(ipsc_cis_concordant_DEGs$target.geneID))#18 (18 reported in paper)

#total pairs
dim(unique(ipsc_cis_concordant_DEGs[,c(2,10)]))#69 (69 reported in paper)

#they have ~350k pairs in background... insane

unique(ipsc_cis_concordant_DEGs$target.geneID)

#presence in the SVSMC:
Enhancer_lociII_DEsig_Enh$FANTOM_ID_merge <- gsub("\\.[0-9]*", "", Enhancer_lociII_DEsig_Enh$FANTOM_ID)

#lncRNAs that have a consistent molecular phenotype and cis contact a DEG during kd
SVSMC_lncs_iPSC_cis <- unique(filter(Enhancer_lociII_DEsig_Enh, FANTOM_ID_merge %in% unique(ipsc_cis_concordant_DEGs$target.geneID))[,c(1,16,44)])

unique(filter(fpkm_allG, EnsID %in% SVSMC_lncs_iPSC_cis$EnsID)[,c(2,4,5,60)])
unique(filter(fpkm_allGDE, EnsID %in% SVSMC_lncs_iPSC_cis$EnsID)[,c(2,49,50)])
filter(AllLNC_AllPCG_1, EnsID %in% SVSMC_lncs_iPSC_cis$EnsID)

#in addition: chr6:28448457-28489790 is the region for CATG00000087927.1 and some peaks there in IGV but no annotated gene


#### Agrawal et al. 2024, Nuclear annotation of lncRNAs ####

#big range of HiC connected lncRNAs across multiple studies
#16x cell types and lines, candidate target-lncRNA pairs can be joined by "2 degrees of seperation"
#i.e 2x loops in a chain
#done for nuclear lncRNAs only
#HiC targets which are DE after ASO in Ramilowski HDF set identified for 33 out of 83 genes
#not as stringent as Ramilowski, but may give some useful lncRNA-target links (linked by both ASO and HiC)
#so basically, HiC connected nuclear lncRNAs across multiple cell lines, tested PERMISSIVELY (one ASO only) for an effect in HDF

#Supp table 4a
Agrawal_HiC_ASO <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Agrawal_HiC_ASO_lncs.csv")
length(unique(Agrawal_HiC_ASO$Nuclear.lncRNA))#81 lncRNAs...not 83
Agrawal_HiC_ASO_targets <- filter(Agrawal_HiC_ASO, Is.candidate.target.gene == "Yes")
length(unique(Agrawal_HiC_ASO_targets$Nuclear.lncRNA))#67 have one DEG and HiC link
#number of ASO hit HiC targets:
table(sapply(lapply(split(Agrawal_HiC_ASO_targets, Agrawal_HiC_ASO_targets$Nuclear.lncRNA), dim), "[[", 1))
#27 are less than 3, 40x have 3x DEG/HiC targets, not 33...

#table contains matches TO THE SAME LNCRNA GENE...needs some cautious filtering:
Agrawal_HiC_ASO_targetsDistinct <- filter(Agrawal_HiC_ASO_targets, !Nuclear.lncRNA == Differentially.Expressed.gene..geneID.)
#unique lnc to PCG pairs (exclude diff ASO pairs)
Agrawal_HiC_ASO_targetsDistinct <- Agrawal_HiC_ASO_targetsDistinct[,c(1,4,5)]
Agrawal_HiC_ASO_targetsDistinct <- unique(Agrawal_HiC_ASO_targetsDistinct)
length(unique(Agrawal_HiC_ASO_targetsDistinct$Nuclear.lncRNA))#47 have one DEG and HiC link to an alternative gene
table(sapply(lapply(split(Agrawal_HiC_ASO_targetsDistinct, Agrawal_HiC_ASO_targetsDistinct$Nuclear.lncRNA), dim), "[[", 1))
#now managed to get to the 33 lncRNAs...

#should they be used though? many ignored in Ramilowski as no second ASO to back up
#however, they do have HiC to connect them which is some validation
#it includes HMGA2-RP11 so is a nice bit of overlap with my work too
trial <- sapply(lapply(split(Agrawal_HiC_ASO_targetsDistinct, Agrawal_HiC_ASO_targetsDistinct$Nuclear.lncRNA), dim), "[[", 1)
Agrawal_cisLncs <- names(trial)[trial >2]

#FANTOM ID, no decimal
Enhancer_lociII_DEsig_Enh$FANTOM_ID_merge <- gsub("\\.[0-9]*", "", Enhancer_lociII_DEsig_Enh$FANTOM_ID)
unique(filter(fpkm_allG, EnsID %in% filter(Enhancer_lociII_DEsig_Enh, FANTOM_ID_merge %in% Agrawal_cisLncs)$EnsID)[,c(2,4:5,60)])
#16x lncRNAs confirmed
unique(filter(fpkm_allG, EnsID %in% filter(Enhancer_lociII_DEsig_Enh, FANTOM_ID_merge %in% Agrawal_cisLncs)$EnsID)$EnsID)
unique(filter(fpkm_allG, EnsID %in% filter(Enhancer_lociII_DEsig_Enh, FANTOM_ID_merge %in% Agrawal_cisLncs)$EnsID)$EnsName)

#6 of which are DE
unique(filter(fpkm_allGDE, EnsID %in% filter(Enhancer_lociII_DEsig_Enh, FANTOM_ID_merge %in% Agrawal_cisLncs)$EnsID)$EnsID)


#### Additional Cis-acting lncRNAs from literature found in CoRegs ####

#TransCistor curated-list overlaps Table S3:
#HOTAIR (potential cis-acting but controversial: recent reasonable looking paper https://www.sciencedirect.com/science/article/pii/S2589004220301929)

#TransCistor "Cis based on a different perturbation experiment" - no idea what this other experiment is tho....
#BOLA3-AS1 (ASO, HDF, not in the list with repeatable ASO pairs)
#FGD5-AS1

#4x lncRNAs checked out by Fannuchi UMLILO, LINC01678 and AC063976.3(ENSG00000224015) are missing from ELs
filter(fpkm_allG, grepl("LINC01678", AllNames)| grepl("LINC01678", EnsName))
filter(fpkm_allG, EnsID == "ENSG00000225331")
filter(fpkm_allG, grepl("AC063976.3", AllNames)| grepl("AC063976.3", EnsName))

#but:
#AC002480.4 (transcript match to IL6-IPL lncRNA gene, the ASO targets the second exon and depletes IL6 + changes epi locus too)
#positive control? and way to de-risk project


#### Total control cis-acting lncRNAs to check up on ####

#From TransCistor:
Trancistor <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/TransCistor.csv", skip = 2)
Trancistor <- filter(Trancistor, !File.Name == "", Species == "Human")

#identify their cis lncs, FDR of 0.25:
Trancistor_cis <- filter(Trancistor, Biobase::rowMin(as.matrix(Trancistor[,12:15])) <= 0.25)
unique(Trancistor_cis$Symbol)
unique(Trancistor_cis$ENSEMBL.ID)


#from HDF only one: RP11-398K22.12
"ENSG00000229852"


#From Yip, 18x lncs
ipsc_cis <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/IpSC_Cis_compare_S4.csv")
length(unique(ipsc_cis$target.geneID))
length(unique(ipsc_cis$target.geneID))

ipsc_cis_concordant <- filter(ipsc_cis, molecular.phenotype.concordant == "yes")

length(unique(ipsc_cis_concordant$target.geneID))
length(unique(ipsc_cis_concordant$target.geneID))#32 (35 reported in paper)

ipsc_cis_concordant_DEGs <- filter(ipsc_cis, abs(EdgeR.Zscore) >1.645, EdgeR.FDR <0.1, molecular.phenotype.concordant == "yes")
length(unique(ipsc_cis_concordant_DEGs$target.geneID))#18 (18 reported in paper)


#Agrawal
Agrawal_HiC_ASO <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/Agrawal_HiC_ASO_lncs.csv")
length(unique(Agrawal_HiC_ASO$Nuclear.lncRNA))#81 lncRNAs...not 83
Agrawal_HiC_ASO_targets <- filter(Agrawal_HiC_ASO, Is.candidate.target.gene == "Yes")
length(unique(Agrawal_HiC_ASO_targets$Nuclear.lncRNA))#67 have one DEG and HiC link
#number of ASO hit HiC targets:
table(sapply(lapply(split(Agrawal_HiC_ASO_targets, Agrawal_HiC_ASO_targets$Nuclear.lncRNA), dim), "[[", 1))
#27 are less than 3, 40x have 3x DEG/HiC targets, not 33...

#table contains matches TO THE SAME LNCRNA GENE...needs some cautious filtering:
Agrawal_HiC_ASO_targetsDistinct <- filter(Agrawal_HiC_ASO_targets, !Nuclear.lncRNA == Differentially.Expressed.gene..geneID.)
#unique lnc to PCG pairs (exclude diff ASO pairs)
Agrawal_HiC_ASO_targetsDistinct <- Agrawal_HiC_ASO_targetsDistinct[,c(1,4,5)]
Agrawal_HiC_ASO_targetsDistinct <- unique(Agrawal_HiC_ASO_targetsDistinct)
length(unique(Agrawal_HiC_ASO_targetsDistinct$Nuclear.lncRNA))#47 have one DEG and HiC link to an alternative gene
table(sapply(lapply(split(Agrawal_HiC_ASO_targetsDistinct, Agrawal_HiC_ASO_targetsDistinct$Nuclear.lncRNA), dim), "[[", 1))
#now managed to get to the 33 lncRNAs...

#should they be used though? many ignored in Ramilowski as no second ASO to back up
#however, they do have HiC to connect them which is some validation
#it includes HMGA2-RP11 so is a nice bit of overlap with my work too
trial <- sapply(lapply(split(Agrawal_HiC_ASO_targetsDistinct, Agrawal_HiC_ASO_targetsDistinct$Nuclear.lncRNA), dim), "[[", 1)
Agrawal_cisLncs <- names(trial)[trial >2]


#Fanucchi, proven effect on locus and WDR5,histone etc (UMLILO already covered)
"ENSG00000232949" #IPL-IL6 aka AC002480.3
"ENSG00000232759" #AC002480.4 which is upstream similar exons, possibly another isoform
"ENSG00000228420" #IPL-CSF1
"ENSG00000225331" #IPL-ICOSLG, LINC01768
"ENSG00000228277" #UMLILO


#Other literature
"HOTAIR"


#all loci:
ControlCisLncs <- list("TransCistor" = unique(Trancistor_cis$ENSEMBL.ID),
                       "HDF_paper" = "ENSG00000229852",
                       "Agrawal" = Agrawal_cisLncs,
                       "Yip" = unique(ipsc_cis_concordant_DEGs$target.geneID),
                       "Fanucchi" = c("ENSG00000232949",
                                      "ENSG00000232759",
                                      "ENSG00000228420",
                                      "ENSG00000225331",
                                      "ENSG00000228277"),
                       #these mostly from Igor's review after a quick check of their papers
                       #must have knockdown and specific regulatory effect on neighbour via a cis-mechanism (recruitment/looping etc) established
                       "OtherLit" = c("ENSG00000229807",#XIST
                                      "ENSG00000228630",#HOTAIR
                                      "ENSG00000243766",#HOTTIP
                                      "ENSG00000281358"#RASSF1-AS1
                                      ) 
                       )

ControlCisLncs <- bind_rows(lapply(ControlCisLncs, as.data.frame), .id = "source")
colnames(ControlCisLncs) <- c("source", "Ens_ID")
#write.csv(ControlCisLncs, "ControlCisLncs.csv")


ControlCisLncs <- unique(c(Trancistor_cis$ENSEMBL.ID,
                 ipsc_cis_concordant_DEGs$target.geneID,
                 "ENSG00000232949", "ENSG00000232759", "ENSG00000228420", "ENSG00000225331"))
#write.csv(ControlCisLncs, "ControlCisLncs.csv")

### extra code ####
#for hSVSMC
ControlCisLnc <- 
#TransCistor
c(unique(filter(fpkm_allGDE, EnsName %in% Trancistor_cis$Symbol | Ens_ID_merge %in% Trancistor_cis$ENSEMBL.ID)[,c(1,2,50)])$EnsName, #5x DELs are called cis-acting by TransCistor
#Yip iPSC lncRNAs with strong cis potential
unique(filter(fpkm_allG, EnsID %in% SVSMC_lncs_iPSC_cis$EnsID)[,c(2,4,5,60)])$EnsName,
#Agrawal HiC/ASO confirmed lncRNAs (more permissive)
#unique(filter(fpkm_allG, EnsID %in% filter(Enhancer_lociII_DEsig_Enh, FANTOM_ID_merge %in% Agrawal_cisLncs)$EnsID)$EnsID),
#+ additional from lit:
"AC002480.4", "HOTAIR")

#Amongst DELs, 9x control cis-acting lncs
filter(fpkm_allGDE, EnsName %in% ControlCisLnc)

#Amongst CoRegs, 3x TransCis, 1x from Fanucchi as well as HOTAIR

#so 3x predicted cis-acting(but either only one ASO or multiple ASOs that disagree for each), 1x established (but only one ASO) and 1x controversial

#notably only one of these is found by Spearman's - justifying a better way to narrow down the coRegs
