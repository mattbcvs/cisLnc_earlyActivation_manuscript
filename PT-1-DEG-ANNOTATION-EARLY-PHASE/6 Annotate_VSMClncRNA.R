#identify VSMC enriched lncRNAs in the timecourse dataset
allLncs_BestCAGE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/allLncs_FANTOMCAGE_3PLARtimecourse.csv",
                             header = T, stringsAsFactors = F)

FANTOM_OntEnrich <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/nature21374-s2/table11.csv", header = T, stringsAsFactors = F)

allLncs_BestCAGE_FANTOM <- merge(allLncs_BestCAGE, FANTOM_OntEnrich[,c(1,3,5,8)], by.x = "FANTOM_ID", by.y = "CAT_geneID")

#Extract FANTOM enrichments:
trial <- sapply(filter(allLncs_BestCAGE_FANTOM, grepl("CL_0000359", associated_sample_ontology))$associated_sample_ontology, strsplit, ",")
triali <- sapply(trial, function(x){x[grepl("CL_0000359", x)]})
trialii <- sapply(triali, strsplit, "\\[")
trial <- sapply(trialii, "[[", 2)
triali <- sapply(trial, strsplit, "\\]")
trialii <- sapply(triali, "[[", 1)

allLncs_BestCAGE_FANTOM$vSMC_FANTenrich[grepl("CL_0000359", allLncs_BestCAGE_FANTOM$associated_sample_ontology)] <- as.numeric(unlist(trialii))
allLncs_BestCAGE_FANTOM$vSMC_FANTenrich <- as.numeric(allLncs_BestCAGE_FANTOM$vSMC_FANTenrich)

trial <- sapply(filter(allLncs_BestCAGE_FANTOM, grepl("CL_0002593", associated_sample_ontology))$associated_sample_ontology, strsplit, ",")
triali <- sapply(trial, function(x){x[grepl("CL_0002593", x)]})
trialii <- sapply(triali, strsplit, "\\[")
trial <- sapply(trialii, "[[", 2)
triali <- sapply(trial, strsplit, "\\]")
trialii <- sapply(triali, "[[", 1)

allLncs_BestCAGE_FANTOM$ITA_SMC_FANTenrich[grepl("CL_0002593", allLncs_BestCAGE_FANTOM$associated_sample_ontology)] <- as.numeric(unlist(trialii))
allLncs_BestCAGE_FANTOM$ITA_SMC_FANTenrich <- as.numeric(allLncs_BestCAGE_FANTOM$ITA_SMC_FANTenrich)

trial <- sapply(filter(allLncs_BestCAGE_FANTOM, grepl("CL_0002596", associated_sample_ontology))$associated_sample_ontology, strsplit, ",")
triali <- sapply(trial, function(x){x[grepl("CL_0002596", x)]})
trialii <- sapply(triali, strsplit, "\\[")
trial <- sapply(trialii, "[[", 2)
triali <- sapply(trial, strsplit, "\\]")
trialii <- sapply(triali, "[[", 1)

allLncs_BestCAGE_FANTOM$Carotid_SMC_FANTenrich[grepl("CL_0002596", allLncs_BestCAGE_FANTOM$associated_sample_ontology)] <- as.numeric(unlist(trialii))
allLncs_BestCAGE_FANTOM$Carotid_SMC_FANTenrich <- as.numeric(allLncs_BestCAGE_FANTOM$Carotid_SMC_FANTenrich)

trial <- sapply(filter(allLncs_BestCAGE_FANTOM, grepl("CL_0002592", associated_sample_ontology))$associated_sample_ontology, strsplit, ",")
triali <- sapply(trial, function(x){x[grepl("CL_0002592", x)]})
trialii <- sapply(triali, strsplit, "\\[")
trial <- sapply(trialii, "[[", 2)
triali <- sapply(trial, strsplit, "\\]")
trialii <- sapply(triali, "[[", 1)

allLncs_BestCAGE_FANTOM$ca_SMC_FANTenrich[grepl("CL_0002592", allLncs_BestCAGE_FANTOM$associated_sample_ontology)] <- as.numeric(unlist(trialii))
allLncs_BestCAGE_FANTOM$ca_SMC_FANTenrich <- as.numeric(allLncs_BestCAGE_FANTOM$ca_SMC_FANTenrich)

trial <- sapply(filter(allLncs_BestCAGE_FANTOM, grepl("CL_0002595", associated_sample_ontology))$associated_sample_ontology, strsplit, ",")
triali <- sapply(trial, function(x){x[grepl("CL_0002595", x)]})
trialii <- sapply(triali, strsplit, "\\[")
trial <- sapply(trialii, "[[", 2)
triali <- sapply(trial, strsplit, "\\]")
trialii <- sapply(triali, "[[", 1)

allLncs_BestCAGE_FANTOM$Sub_SMC_FANTenrich[grepl("CL_0002595", allLncs_BestCAGE_FANTOM$associated_sample_ontology)] <- as.numeric(unlist(trialii))
allLncs_BestCAGE_FANTOM$Sub_SMC_FANTenrich <- as.numeric(allLncs_BestCAGE_FANTOM$Sub_SMC_FANTenrich)

trial <- sapply(filter(allLncs_BestCAGE_FANTOM, grepl("CL_0002539", associated_sample_ontology))$associated_sample_ontology, strsplit, ",")
triali <- sapply(trial, function(x){x[grepl("CL_0002539", x)]})
trialii <- sapply(triali, strsplit, "\\[")
trial <- sapply(trialii, "[[", 2)
triali <- sapply(trial, strsplit, "\\]")
trialii <- sapply(triali, "[[", 1)

allLncs_BestCAGE_FANTOM$ao_SMC_FANTenrich[grepl("CL_0002539", allLncs_BestCAGE_FANTOM$associated_sample_ontology)] <- as.numeric(unlist(trialii))
allLncs_BestCAGE_FANTOM$ao_SMC_FANTenrich <- as.numeric(allLncs_BestCAGE_FANTOM$ao_SMC_FANTenrich)

trial <- sapply(filter(allLncs_BestCAGE_FANTOM, grepl("CL_0002588", associated_sample_ontology))$associated_sample_ontology, strsplit, ",")
triali <- sapply(trial, function(x){x[grepl("CL_0002588", x)]})
trialii <- sapply(triali, strsplit, "\\[")
trial <- sapply(trialii, "[[", 2)
triali <- sapply(trial, strsplit, "\\]")
trialii <- sapply(triali, "[[", 1)

allLncs_BestCAGE_FANTOM$uv_SMC_FANTenrich[grepl("CL_0002588", allLncs_BestCAGE_FANTOM$associated_sample_ontology)] <- as.numeric(unlist(trialii))
allLncs_BestCAGE_FANTOM$uv_SMC_FANTenrich <- as.numeric(allLncs_BestCAGE_FANTOM$uv_SMC_FANTenrich)

trial <- sapply(filter(allLncs_BestCAGE_FANTOM, grepl("CL_0002594", associated_sample_ontology))$associated_sample_ontology, strsplit, ",")
triali <- sapply(trial, function(x){x[grepl("CL_0002594", x)]})
trialii <- sapply(triali, strsplit, "\\[")
trial <- sapply(trialii, "[[", 2)
triali <- sapply(trial, strsplit, "\\]")
trialii <- sapply(triali, "[[", 1)

allLncs_BestCAGE_FANTOM$ua_SMC_FANTenrich[grepl("CL_0002594", allLncs_BestCAGE_FANTOM$associated_sample_ontology)] <- as.numeric(unlist(trialii))
allLncs_BestCAGE_FANTOM$ua_SMC_FANTenrich <- as.numeric(allLncs_BestCAGE_FANTOM$ua_SMC_FANTenrich)

trial <- sapply(filter(allLncs_BestCAGE_FANTOM, grepl("CL_0002590", associated_sample_ontology))$associated_sample_ontology, strsplit, ",")
triali <- sapply(trial, function(x){x[grepl("CL_0002590", x)]})
trialii <- sapply(triali, strsplit, "\\[")
trial <- sapply(trialii, "[[", 2)
triali <- sapply(trial, strsplit, "\\]")
trialii <- sapply(triali, "[[", 1)

allLncs_BestCAGE_FANTOM$bv_SMC_FANTenrich[grepl("CL_0002590", allLncs_BestCAGE_FANTOM$associated_sample_ontology)] <- as.numeric(unlist(trialii))
allLncs_BestCAGE_FANTOM$bv_SMC_FANTenrich <- as.numeric(allLncs_BestCAGE_FANTOM$bv_SMC_FANTenrich)

trial <- sapply(filter(allLncs_BestCAGE_FANTOM, grepl("CL_0002589", associated_sample_ontology))$associated_sample_ontology, strsplit, ",")
triali <- sapply(trial, function(x){x[grepl("CL_0002589", x)]})
trialii <- sapply(triali, strsplit, "\\[")
trial <- sapply(trialii, "[[", 2)
triali <- sapply(trial, strsplit, "\\]")
trialii <- sapply(triali, "[[", 1)

allLncs_BestCAGE_FANTOM$bc_SMC_FANTenrich[grepl("CL_0002589", allLncs_BestCAGE_FANTOM$associated_sample_ontology)] <- as.numeric(unlist(trialii))
allLncs_BestCAGE_FANTOM$bc_SMC_FANTenrich <- as.numeric(allLncs_BestCAGE_FANTOM$bc_SMC_FANTenrich)

trial <- as.list(data.frame(t(allLncs_BestCAGE_FANTOM[,19:28])))
triali <- sapply(trial, function(x){length(x[!is.na(x)])})

allLncs_BestCAGE_FANTOM$vSMCFacetNo <- unlist(triali)

table(allLncs_BestCAGE_FANTOM$vSMCFacetNo)

#add max. enrichment
trial <- sapply(data.frame(t(allLncs_BestCAGE_FANTOM[allLncs_BestCAGE_FANTOM$vSMCFacetNo >0,19:28])), max, na.rm = T)
allLncs_BestCAGE_FANTOM$maxEnrich[allLncs_BestCAGE_FANTOM$vSMCFacetNo >0] <- unlist(trial)

trial <- sapply(data.frame(t(allLncs_BestCAGE_FANTOM[allLncs_BestCAGE_FANTOM$vSMCFacetNo >0,19:28])), function(x){
  colnames(allLncs_BestCAGE_FANTOM[,19:28])[which(x == max(x, na.rm = T))]})
allLncs_BestCAGE_FANTOM$maxEnrichtype[allLncs_BestCAGE_FANTOM$vSMCFacetNo >0] <- unlist(trial)

#write.csv(allLncs_BestCAGE_FANTOM, "allLncs_BestCAGE_FANTOMtime.csv", row.names = F)
