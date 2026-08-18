#aim is to replace genes inappropriately merged by stringtie with original genv26 annotation

library(dplyr)

#GTF created through stringtie --merge on the 3x .gtfs in the 2021 IJMS publication
#GENCODEv26 primary assembly merged with all 3 stranded .gtfs as in 3PLARVSMC_explore in the "R scripts" folder of that manuscript
#There is a filtered version of this .gtf that was put into RSEM round II, but here we import the full merged .gtf (either would work?):
stringtie_gtf <- read.delim("3PLAR_allgenv26_Timecourse_ff.gtf", 
                            header= F, stringsAsFactors=F)
head(stringtie_gtf)


#
#### start with tx only ####
#transcripts only
stringtie_transcripts <- filter(stringtie_gtf, V3 == "transcript")

#seperate gene id (from stringtie)
stringtie_transcripts$gene_id <- gsub("gene_id ", "", sapply(stringtie_transcripts$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^gene_id", y)]
    })
  }))

#seperate ref gene id (from ensembl if present)
stringtie_transcripts$ref_gene_id <- gsub(";", "", gsub("ref_gene_id ", "", sapply(stringtie_transcripts$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("ref_gene_id", y)]
  })
})))

#isolate rows for multi gene clusters here:
stringtie_transcripts$ref_gene_id[stringtie_transcripts$ref_gene_id == "character(0)"] <- NA

stringtie_transcripts$transcript_id <- gsub(";", "", gsub("transcript_id ", "", sapply(stringtie_transcripts$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^transcript_id", y)]
  })
})))

trial <- split(stringtie_transcripts, stringtie_transcripts$gene_id)
triali <- sapply(trial, function(x){
  #should have 2 ENSG IDs in the ref col to be considered for splitting up
  length(unique(filter(x, grepl("ENSG", ref_gene_id))$ref_gene_id)) >1
})
multiGeneClusters <- trial[triali]

multiGeneClusters[[1]]
multiGeneClusters[[2]]
multiGeneClusters[[3]]
multiGeneClusters[[10]]

#these rows will be removed, and replaced with GENv26 entries
multiGeneClusters <- bind_rows(multiGeneClusters)

length(unique(multiGeneClusters$transcript_id))

#12789 tx in total to remove
#the novels will be lost, the enst will be switched with unaltered .gtf
multiGeneClusters_ensTx <- filter(multiGeneClusters, grepl("ENST", transcript_id))

#10021 ens tx to replace
length(unique(multiGeneClusters_ensTx$transcript_id))

#so after the process, should lose some tx:
12789-10021

#2768 transcripts (all MSTRG) will go

#genv26 transcripts
genv26_gtf <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker_group-bioinformatics/annotation/human/GENCODE/gencode.v26.primary_assembly.annotation.gtf", 
                         header= F, stringsAsFactors=F, skip = 5)
head(genv26_gtf, 10)
#transcripts only
genv26_transcripts <- filter(genv26_gtf, V3 == "transcript")
rm(genv26_gtf)
gc()

genv26_transcripts$transcript_id <- gsub("transcript_id ", "", sapply(genv26_transcripts$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^transcript_id", y)]
  })
}))

head(genv26_transcripts$transcript_id)

#transcripts to replace
genv26_transcripts_to_append <- filter(genv26_transcripts, transcript_id %in% multiGeneClusters_ensTx$transcript_id)

#n.b. if wanting to re-use the annotation on new samples (e.g. RNAseq on lncRNA knockdowns) then also need the other gencodev26 genes not found in initial samples:
genv26_transcripts_to_appendi <- filter(genv26_transcripts, !transcript_id %in% c(stringtie_transcripts$transcript_id, multiGeneClusters_ensTx$transcript_id))

genv26_transcripts_to_append <- rbind(genv26_transcripts_to_append, genv26_transcripts_to_appendi)

#check for format diffs, e.g. "ref_gene_id"
head(genv26_transcripts_to_append[,1:9], 10)

head(multiGeneClusters_ensTx[,1:9], 10)

#so need to make a new "V9" column
#isolate gene_id, transcript_id, gene_name and put ref_gene_id as the same one
genv26_transcripts_to_append$gene_id <- gsub("gene_id ", "", sapply(genv26_transcripts_to_append$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^gene_id", y)]
  })
}))

genv26_transcripts_to_append$gene_name <- gsub("gene_name ", "", sapply(genv26_transcripts_to_append$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^gene_name", y)]
  })
}))

genv26_transcripts_to_append$newV9 <- paste0("gene_id ", genv26_transcripts_to_append$gene_id, "; transcript_id ", genv26_transcripts_to_append$transcript_id,
                                             "; gene_name ", genv26_transcripts_to_append$gene_name, "; ref_gene_id ", genv26_transcripts_to_append$gene_id, ";")

head(multiGeneClusters_ensTx[,1:9], 10)
head(genv26_transcripts_to_append[,c(1:8,13)], 10)

#replace the 12789 merged gene entries
#with the 10021 tx (or more if making for a new RNAseq) for these genes from genv26:

colnames(genv26_transcripts_to_append)[9] <- "oldV9"
colnames(genv26_transcripts_to_append)[13] <- "V9"

stringtie_transcripts_ <- rbind(filter(stringtie_transcripts, !V9 %in% multiGeneClusters$V9)[,1:9],
                                genv26_transcripts_to_append[,c(1:8,13)])
#for RNAseq on same samples:
#based on number of novel tx to lose, this object should be 
85369 - 2768
#which it is

#### now make a full .gtf ####

#need to have a full .gtf, exons etc for RSEM:
stringtie_gtf <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/PLAR_4timepoints/3PLAR_allgenv26_Timecourse_ff.gtf", 
                            header= F, stringsAsFactors=F)
head(stringtie_gtf)

genv26_gtf <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker_group-bioinformatics/annotation/human/GENCODE/gencode.v26.primary_assembly.annotation.gtf", 
                         header= F, stringsAsFactors=F, skip = 5)
head(genv26_gtf, 10)

#pull out rows for offending transcripts:
stringtie_gtf$transcript_id <- gsub(";", "", gsub("transcript_id ", "", sapply(stringtie_gtf$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^transcript_id", y)]
  })
})))

#filtered version of gtf:
stringtie_gtf_ <- filter(stringtie_gtf, !transcript_id %in% multiGeneClusters$transcript_id)

head(stringtie_gtf_)

#obtain relevant rows from GENv26
genv26_gtf$transcript_id <- gsub("transcript_id ", "", sapply(genv26_gtf$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^transcript_id", y)]
  })
}))

genv26_gtf_ <- filter(genv26_gtf, transcript_id %in% genv26_transcripts_to_append$transcript_id)

#note, the stringtie .gtf does not have gene or other entries, just tx and exon
genv26_gtf_ <- filter(genv26_gtf_, V3 %in% c("transcript", "exon"))

#format for transcript rows:
head(filter(stringtie_gtf_, V3 == "transcript"))
#format for exon rows:
head(filter(stringtie_gtf_, V3 == "exon"))

#requires gene id, name as well as tx id:
genv26_gtf_$gene_id <- gsub("gene_id ", "", sapply(genv26_gtf_$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^gene_id", y)]
  })
}))

genv26_gtf_$gene_name <- gsub("gene_name ", "", sapply(genv26_gtf_$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^gene_name", y)]
  })
}))

genv26_gtf_$exon_number <- gsub("exon_number ", "", sapply(genv26_gtf_$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^exon_number", y)]
  })
}))

#may be useful to have a "0" entry for ordering later:
genv26_gtf_$exon_number[genv26_gtf_$exon_number == "character(0)"] <- 0

genv26_gtf_tx <- filter(genv26_gtf_, V3 == "transcript")
genv26_gtf_tx$newV9 <- paste0("gene_id ", genv26_gtf_tx$gene_id, "; transcript_id ", genv26_gtf_tx$transcript_id,
                              "; gene_name ", genv26_gtf_tx$gene_name, "; ref_gene_id ", genv26_gtf_tx$gene_id, ";")
head(genv26_gtf_tx[,c(1:8,14)])

genv26_gtf_ex <- filter(genv26_gtf_, V3 == "exon")
genv26_gtf_ex$newV9 <- paste0("gene_id ", genv26_gtf_ex$gene_id, "; transcript_id ", genv26_gtf_ex$transcript_id,
                              "; exon_number ", genv26_gtf_ex$exon_number,
                              "; gene_name ", genv26_gtf_ex$gene_name, "; ref_gene_id ", genv26_gtf_ex$gene_id, ";")

head(genv26_gtf_ex[,c(1:8,14)])

genv26_gtf_f <- rbind(genv26_gtf_tx, genv26_gtf_ex)

#append
colnames(genv26_gtf_f)
colnames(genv26_gtf_f)[9] <- "oldV9"
colnames(genv26_gtf_f)[14] <- "V9"

trial <- rbind(stringtie_gtf_[,c(1:9)], genv26_gtf_f[,c(1:8,14)])


#### sanity checking ####

#confirm correct outcome
#the ref_gene_id column (not the mstrg id col), should have same no. of Ens genes
trial$ref_gene_id <- gsub("ref_gene_id ", "", sapply(trial$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^ref_gene_id", y)]
  })
}))
trial$ref_gene_id[trial$ref_gene_id == "character(0)"] <- NA
length(unique(trial$ref_gene_id))#17600 genes in the processed gtf (58279 but one NA for novel gene entries)

#if adding in all gencode then should match the gtf:
genv26_transcripts$gene_id <- gsub("gene_id ", "", sapply(genv26_transcripts$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^gene_id", y)]
  })
}))
length(unique(genv26_transcripts$gene_id))#58278 genes in the original (matches)

#otherwise should be same as stringtie gtf
stringtie_gtf$ref_gene_id <- gsub("ref_gene_id ", "", sapply(stringtie_gtf$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^ref_gene_id", y)]
  })
}))
stringtie_gtf$ref_gene_id[stringtie_gtf$ref_gene_id == "character(0)"] <- NA
length(unique(stringtie_gtf$ref_gene_id))#17600 genes here too

#n.b. one of the above is an NA for the novel genes
length(unique(filter(trial, grepl("ENSG", ref_gene_id))$ref_gene_id))
length(unique(filter(stringtie_gtf, grepl("ENSG", ref_gene_id))$ref_gene_id))

#if adding in all of gencode then should have the tx in this gtf, plus novels from stringtie
#otherwise should have 82601 transcripts as calculated above (full .gtf minus the novel tx as merged ens genes):
trial$transcript_id <- gsub(";", "", gsub("transcript_id ", "", sapply(trial$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^transcript_id", y)]
  })
})))
length(unique(trial$transcript_id))#82601 as expected (or 214,363 otherwise for all genv26 version)

length(unique(stringtie_gtf$transcript_id))#85369 tx in the gtf used for the timecourse
length(unique(multiGeneClusters$transcript_id))#12789 removed
length(unique(multiGeneClusters_ensTx$transcript_id))#10021 re-entered
genv26_transcripts_to_append_check2 <- filter(genv26_transcripts, !transcript_id %in% c(stringtie_transcripts$transcript_id))
length(unique(genv26_transcripts_to_append_check2$transcript_id))#131762 added
85369-12789+10021+131762

#there should be a larger number of genes now that the merged Ens have been split:
trial$gene_id <- gsub("gene_id ", "", sapply(trial$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^gene_id", y)]
  })
}))
trial$gene_id[trial$gene_id == "character(0)"] <- NA
length(unique(trial$gene_id))#18710 genes in the processed gtf

length(unique(stringtie_gtf$gene_id))#17327 genes in the original gtf



#### sort order #### 
trial$transcript_id <- gsub(";", "", gsub("transcript_id ", "", sapply(trial$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^transcript_id", y)]
  })
})))

trial$gene_id <- gsub("gene_id ", "", sapply(trial$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^gene_id", y)]
  })
}))

trial$exon_number <- gsub(";", "", gsub("exon_number ", "", sapply(trial$V9, function(x){
  sapply(strsplit(x, "; "), function(y){
    y[grepl("^exon_number", y)]
  })
})))
trial$exon_number[trial$exon_number == "character(0)"] <- 0
trial$exon_number <- as.numeric(trial$exon_number)

trial <- trial[order(trial$V1, trial$gene_id, trial$transcript_id, trial$exon_number),]


#### save, with quotes ####

triali <- sapply(trial$V9, strsplit, " ")
QuotedChars <- sapply(triali, gsub, pattern = ";", replacement = "")
trialii <- sapply(QuotedChars, function(x){
  x[seq(2,length(x), 2)] <- paste('"', x[seq(2,length(x), 2)], '"', ';', sep = "")
  paste(x, collapse = " ")
})
trialiii <- trial
trialiii$V9 <- unlist(trialii)

write.table(trialiii[,1:9], sep = "\t", row.names= F, quote = F, col.names = F,
            "3PLAR_allgenv26_Timecourse_ff_nonMerge_fullgenv26.gtf")
