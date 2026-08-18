#is the split of PCG lnc bit simplistic?

#e.g. where are pseudogenes etc...

#for fpkm_allG:
PCG <- filter(fpkm_allG, grepl("coding|TF|CC", GeneClassUpdate))

PCG <- unique(PCG[,c(2,6,60)])

table(PCG$EnsType)#contains some pseudogenes
table(PCG$GeneClassUpdate)

#lncs
lncs <- filter(fpkm_allG, grepl("fide|Lnc", GeneClassUpdate))

lncs <- unique(lncs[,c(2,6,60)])

table(lncs$EnsType)#contains some pseudogenes
table(lncs$GeneClassUpdate)

#remaining
other <- filter(fpkm_allG, !EnsID %in% c(PCG$EnsID, lncs$EnsID))

other <- unique(other[,c(2,6,60)])

table(other$EnsType)#contains a lot of pseudogenes
table(other$GeneClassUpdate)#just the putative lncs (>50% are pseudogenes)
