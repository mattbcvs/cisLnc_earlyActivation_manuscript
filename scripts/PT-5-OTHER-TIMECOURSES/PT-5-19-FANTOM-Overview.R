#### FANTOM timecourses ####
library(dplyr)

#identify appropriate timecourse data for running similar stats
#similar timepoints
#quality data

#Supp Table 18 has the library IDs across the timecourses used in FANTOM CAT paper:
Supp18 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/nature21374-s2/SuppTab18.csv")

unique(Supp18$series_name)#25x timecourses
Supp18_comps <- as.data.frame(table(Supp18$series_name))#indicates number of comparisons

#final timepoint:
Supp18list <- split(Supp18, Supp18$series_name)
trial <- sapply(Supp18list, function(x){
  x$qry_replicate_ID[length(x[,1])]
})

#first timepoint:
triali <- sapply(Supp18list, function(x){
  x$qry_replicate_ID[1]
})

Supp18_comps <- data.frame(Supp18_comps, "firstTime" = triali, "lastTime" = trial)

#filter to those of interest:
Supp18_compsII <- filter(Supp18_comps, Freq > 1, #not timecourse
                         grepl("hr", firstTime), #first timepoint should be hrs not days or a non-time variable
                         !grepl("Aortic|MCF7|VEGFC", Var1) #finish too soon
                         )

#full timepoints available for remaining:
Supp18II <- filter(Supp18, series_name %in% Supp18_compsII$Var1)
Supp18list <- split(Supp18II, Supp18II$series_name)
trial <- sapply(Supp18list, function(x){
  paste(x$qry_replicate_ID, collapse = "-")
})
Supp18_compsII$allTimes <- trial

#all are appropriate so far -4x contain some very early
#6x are immune based, 4 of these are infection, 2 are direct stims for immune pathways
#4x extend to 48hrs and 2x to beyond a week

#n.b. 4x others may be of use for other purposes, confirming co-regulation in early timeframes:
remainingComps <- filter(Supp18_comps, !Var1 %in% Supp18_compsII$Var1)
remainingComps <- filter(remainingComps, grepl("hr", firstTime), !grepl("VEGFC", Var1))
#the VEGFC has v. few DEGs (prelim) so is not counted:


#extract counts for all from big FANT table:
#Human Early response data, including some additional timecourses in similar timeframe to SVSMC but classed as "activation":
Supp18_ER <- filter(Supp18, series_name %in% Supp18_compsII$Var1 | series_name %in% remainingComps$Var1)

#number of comparisons per ER
table(Supp18_ER$series_tag)

#all libraries across ERs:
Supp18_ER_libs <- unique(c(unlist(strsplit(Supp18_ER$ref_sample_ID, ",")), unlist(strsplit(Supp18_ER$qry_sample_ID, ","))))

Libraries_bigTable <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Libraries", header = F)
Libraries_bigTable <- as.character(Libraries_bigTable)

#are all libraries in this big table? yes (421):
sum(Supp18_ER_libs %in% Libraries_bigTable)

#subsetted table from FANTOM:
#subsetted fantom atlas, now diverging from previous by getting counts aiming to re-run QC, timepoint selection + DE analysis:
genes <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/genes_robust", 
                  header = T)
ER_counts <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/FANTOM_CAT.expression_ERtime.genes.lv3_robust.count.csv", 
                      header = T)
ER_counts <- cbind(genes, ER_counts)
sum(Supp18_ER_libs %in% colnames(ER_counts))


#now assess each of 0-24hr datasets in terms of:
#QC for 0hr -24 hr timepoints
#PCA shift over key timepoints (beyond 48hr too)
#DEG number, DEL number

#then decide which are useful, and how:
#confirming SVSMC effects are wider
#confirming type of stimuli that induce
#confirming what are the drivers of the phenomenon

#same for shorter but these are less priority rn