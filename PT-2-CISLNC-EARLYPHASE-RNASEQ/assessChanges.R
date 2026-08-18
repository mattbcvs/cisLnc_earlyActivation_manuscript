#big discrepancy in pairing distances:

#possible culprits are:
#splitting up the gencode genes, now ~1000 PCGs more in the pile that were previously lumped together by stringtie
#also means more lncRNAs, gencode etc, potentially of a particular type (more bidirectionals?)
#lowering the FPKM thresh to 0.8, more lncRNAs in the mix now
#different transcript filtering

AllLNC_AllPCG_1i <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_info_2026.csv", header = T)
AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_info_aug25_redone.csv", header = T)

#mistake in this one
#AllLNC_AllPCG_1ii <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_info_aug25.csv", header = T)

missingInAug <- filter(AllLNC_AllPCG_1, !pairs %in% AllLNC_AllPCG_1ii$pairs)

length(unique(AllLNC_AllPCG_1i$pairs)) #8438 pairs now
length(unique(AllLNC_AllPCG_1$pairs)) #5441 pairs previously

#still all found in a 1Mbp range, similar neighbours
summary(AllLNC_AllPCG_1i$DisLnc_PCG)
summary(AllLNC_AllPCG_1$DisLnc_PCG)

#can see a v. slight tendency to find more neighbours now
summary(as.numeric(table(AllLNC_AllPCG_1i$EnsID)))
summary(as.numeric(table(AllLNC_AllPCG_1$EnsID)))

#there are 3614 pairs in the new but not old
newPairs <- filter(AllLNC_AllPCG_1i, !pairs %in% AllLNC_AllPCG_1$pairs)

#so presumably these are from lncRNAs - or PCGs - that were not present in last analysis:
fpkm_allG_old <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGClassUpdate.csv", header = T)
FPKM_CQV_OVERLAP_fpkm_old <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/FPKM_CQV_OVERLAP_fpkm.csv")
table(FPKM_CQV_OVERLAP_fpkm_old$IGV)#413 pass, 168 fail

#remove artefacts (push back to step 7?)
fpkm_allG_old <- filter(fpkm_allG_old, grepl("chr", chr), !grepl("artefacts", GeneClassUpdate), !is.na(GeneClassUpdate))
#remove manual or Thresh4 fails
fpkm_allG_old <- filter(fpkm_allG_old, 
                                !EnsID %in% filter(FPKM_CQV_OVERLAP_fpkm_old, IGV == "fail")$EnsID, #remove manual fails
                        )

length(unique(newPairs$EnsID))#546 lncRNAs
length(unique(filter(newPairs, EnsID %in% fpkm_allG_old$EnsID)$EnsID))#273 found in previous
length(unique(filter(newPairs, EnsID %in% fpkm_allG_old$EnsID,
                     EnsID.y %in% fpkm_allG_old$EnsID)$EnsID))#only 23 of these with the neighbouring PCG labelled

#so why not picked up on previously?
newPairs_shouldBeInOld <- filter(newPairs, EnsID %in% fpkm_allG_old$EnsID,
                                 EnsID.y %in% fpkm_allG_old$EnsID)

#transcript diffs? probably
newPairs_shouldBeInOld[1,]
filter(fpkm_allG_old, EnsID %in% "ENSG00000219665.8")#1 tx
filter(fpkm_allG, EnsID %in% "ENSG00000219665.8")#5 tx

filter(fpkm_allG_old, EnsName %in% "NFIX")#3 tx
filter(fpkm_allG, EnsName %in% "NFIX")#5 tx


#and a small number of lost pairs interestingly
lostPairs <- filter(AllLNC_AllPCG_1, !pairs %in% AllLNC_AllPCG_1i$pairs)

#presumably igv fails?
length(unique(lostPairs$EnsID))#75 lncRNAs
length(unique(filter(lostPairs, EnsID %in% fpkm_allG$EnsID)$EnsID))#35 found in the update
length(unique(filter(lostPairs, EnsID %in% fpkm_allG$EnsID,
                     EnsID.y %in% fpkm_allG$EnsID)$EnsID))#only 11 of these with the neighbouring PCG labelled

#remainder? annotation changes
oldPairs_shouldBeInNew <- filter(lostPairs, EnsID %in% fpkm_allG$EnsID,
                                 EnsID.y %in% fpkm_allG$EnsID)

oldPairs_shouldBeInNew[1,]
filter(fpkm_allG_old, EnsID %in% "ENSG00000230615.6")#1 tx
filter(fpkm_allG, EnsID %in% "ENSG00000230615.6")#more tx, some low express are PCG

filter(fpkm_allG_old, EnsName %in% "KDM4A")#2 tx
filter(fpkm_allG, EnsName %in% "KDM4A")#3 tx
