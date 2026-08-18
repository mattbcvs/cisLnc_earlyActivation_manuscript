#Add annotation and other fine detail to table, from GENv26 and PLAR

#GTF created through stringtie --merge on the 3x .gtfs in the 2021 IJMS publication
#then merged ensembl genes sorted in PT-0 script
#output gtf here:
stringtie_gtf <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/3PLAR_allgenv26_Timecourse_ff_nonMerge.gtf", 
                            header= F, stringsAsFactors=F)
test_transcript<-stringtie_gtf[,"V3"]=="transcript"
transcript_table<-stringtie_gtf[test_transcript,]
head(transcript_table)

transcript_table_locus1<-paste(transcript_table[,1],transcript_table[,4],sep=":")
transcript_table_locus2<-paste(transcript_table_locus1,transcript_table[,5],sep="-")
transcript_table_locus<-paste(transcript_table_locus2,transcript_table[,7],sep=" ")
head(transcript_table_locus)

# use column 9 to extract name
column_name<-sapply(transcript_table[,9],strsplit,split=";")
gene_id_entry<- sapply(column_name, "[[", 1)
gene_id_split<-sapply(gene_id_entry,strsplit,split=" ")
gene_id<-sapply(gene_id_split, "[[", 2)
transcript_id_entry<- sapply(column_name, "[[", 2)
transcript_id_split<-sapply(transcript_id_entry,strsplit,split=" ")
transcript_id<-sapply(transcript_id_split, "[[", 3)

#subset those rows that have an ensembl id
test_ref<-grep("ref_gene",column_name)
#create list of locus or gene name for each row if novel or ensembl respectively
ref_id_list<-c()
nb<-nrow(transcript_table)
for (i in 1:nb) {
  if (i %in% test_ref) {
    ref_id_line<- column_name[[i]]
    ref_id_entry<-ref_id_line[4] #edited to extract ensembl id instead
    ref_id_split<-strsplit(ref_id_entry,split=" ")[[1]]
    ref_id<-ref_id_split[3]
    ref_id_list<-c(ref_id_list,ref_id)
  } else {
    ref_id_list<-c(ref_id_list,NA)
  }
}

ref_name_list<-c()
for (i in 1:nb) {
  if (i %in% test_ref) {
    ref_name_line<- column_name[[i]]
    ref_name_entry<-ref_name_line[3] #and the original loop for the name too
    ref_name_split<-strsplit(ref_name_entry,split=" ")[[1]]
    ref_name <-ref_name_split[3]
    ref_name_list<-c(ref_name_list,ref_name)
  } else {
    ref_name_list<-c(ref_name_list,NA)
  }
}

#edited to include ensembl id and name too
convert_file<-unique(data.frame(transcript_table[,c(1,4,5,7)],
                                gene_id=gene_id,
                                ref_id=ref_id_list,
                                ref_name=ref_name_list,
                                transcript_id=transcript_id,
                                transcript_table_locus=transcript_table_locus,
                                stringsAsFactors = FALSE))

# for each gene id - several locus coordinates
# I want to know if each Id can correspond to ENSEMBL Id 
# if there is a reference Id, it should be use by default


#new transcripts for known genes, assign correct EnsID and EnsName
trial <- split(convert_file, as.factor(convert_file$gene_id))

#add in EnsID if it is there
triali <- lapply(trial, function(x){
  #this line previously probs an issue: some genes (~9%) contained multiple EnsIDs
  #assigning the first was bit misleading, implying only one gene is here
  #have now split any ensIDs up manually and re-ran
  
  #works cos sorted alphabetically
  x$ref_id <- x$ref_id[1]
  
  x})

#add in EnsName if it is there (n.b. same issue as above)
trialii <- lapply(triali, function(x){
  x$ref_name <- x$ref_name[1]
  x})

#if no EnsID then add the stringtie-assigned gene id (which would be MSTRG...)
trialiii <- lapply(trialii, function(x){
  x$ref_id[!grepl("ENS", x$ref_id)] <- x$gene_id[!grepl("ENS", x$ref_id)]
  x})

#if no EnsID then remove the gene name...?
trialiv <- lapply(trialiii, function(x){
  x$ref_name[grepl("MSTRG", x$ref_id)] <- NA
  x})
trialv <- bind_rows(trialiv, .id = "column_label")

#taks 100 years to import
gencode_v26_gft_table <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker_group-bioinformatics/annotation/human/GENCODE/gencode.v26.primary_assembly.annotation.gtf", 
                                    header=FALSE, stringsAsFactors=FALSE, sep = "\t", skip = 5)
gencode_v26_gft_table_G <- filter(gencode_v26_gft_table, V3 == "gene")
gencode_v26_v9 <- strsplit(gencode_v26_gft_table_G$V9, ";")

annotation <- data.frame("ENSEMBL" = gsub(".*\\ ", "", sapply(gencode_v26_v9, "[[", 1)),
                         "type" = gsub(".*\\ ", "", sapply(gencode_v26_v9, "[[", 2)),
                         "chr" = gencode_v26_gft_table_G$V1,
                         "start" = gencode_v26_gft_table_G$V4,
                         "end" = gencode_v26_gft_table_G$V5,
                         "str" = gencode_v26_gft_table_G$V7)

#annotation<-gencode_v26_gft_table[,c(4,6,1:3,5)]
colnames(annotation)<-c("ENSEMBL","type","chr","start","end","str")

#add the genv26 annotation
file3 <- merge(trialv,annotation,by.x=7,by.y=1,all.x=TRUE)
#adjust the gene locus (to fit genv26)
file3[grepl("ENS", file3$ref_id),3:5] <- file3[grepl("ENS", file3$ref_id),12:14]
#then can cut off the additional gene locus cols (and "column label" col)
file3 <- file3[,-c(2,12:15)]

length(unique(file3$gene_id))#18710 genes measured by rsem
length(unique(file3$ref_id))#18691 in the created ref_id col

#some stringtie genes are only part of ensembl genes:
idChecks <- unique(file3[,c(1,6)])
idChecks_refDups <- unique(idChecks[which(duplicated(idChecks$ref_id)),1])#19 ensIDs are split into 2 
idChecks_refDups <- filter(idChecks, ref_id %in% idChecks_refDups)

idChecks_refDups$label <- c("a","b")
idChecks_refDups$ref_id2 <- paste0(idChecks_refDups$ref_id, idChecks_refDups$label)

file3_split <- filter(file3, ref_id %in% idChecks_refDups$ref_id)
file3_split <- merge(file3_split, idChecks_refDups, by = "gene_id")
file3_split$ref_id.x <- file3_split$ref_id2
colnames(file3_split)[2] <- "ref_id"

file3_ <- filter(file3, !ref_id %in% idChecks_refDups$ref_id)
file3_ <- rbind(file3_, file3_split[,c(2:6,1,7:10)])
file3 <- file3_

length(unique(file3$gene_id))#18710 genes measured by rsem
length(unique(file3$ref_id))#18710 in the created ref_id col

#think the cols are now split a better than previously
table(table(file3$ref_id))
table(table(file3$gene_id))

noEns <- filter(file3, is.na(ref_name))
length(unique(noEns$ref_id))#1092 novel genes measured by rsem, no entry in genv26

file3 <- file3[,c(6,1,7,10,2:5,8:9)]
colnames(file3) <- c("MSTRG_ID", "EnsID", "EnsName", "EnsType", "chr", "start", "stop",
                     "str", "MSTRG_Tx_ID", "Tx_Locus")

#will do for now, bear in mind the gene locus columns don't really make sense:
#taken from genv26

#integrate outputs from PLAR
#first batch of annotations, genes in first 4 timepoints
Timecourse.annot.bed.assemblyStats <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/PLAR_4timepoints/Timecourse3PLAR.annot.bed.assemblyStats.txt", header=FALSE, stringsAsFactors=FALSE)
Timecourse.lincs.f1 <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/PLAR_4timepoints/Timecourse3PLAR.lincs.f1.bed", header=FALSE, stringsAsFactors=FALSE)
Timecourse.lincs.f1.clean <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/PLAR_4timepoints/Timecourse3PLAR.lincs.f1.clean.bed", header=FALSE, stringsAsFactors=FALSE)
Timecourse.lincs.f1.pure <- read.delim("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/PLAR_4timepoints/Timecourse3PLAR.lincs.f1.pure.bed", header=FALSE, stringsAsFactors=FALSE)

#append tx_id
trial <- sapply(as.list(Timecourse.annot.bed.assemblyStats$V1), strsplit, "\\|")
triali <- unlist(sapply(trial, function(x){
  x[length(x)]
}))
trial <- sapply(as.list(triali), strsplit, "\\:")
Timecourse.annot.bed.assemblyStats$MSTRG_Tx_ID <- unlist(sapply(trial, "[[", 1))

PLAR_Prefilter <- merge(file3, Timecourse.annot.bed.assemblyStats[,c(9,4,5,7,8)], by = "MSTRG_Tx_ID", all.x = T)

trial <- sapply(as.list(Timecourse.lincs.f1$V4), strsplit, "\\|")
triali <- unlist(sapply(trial, function(x){
  x[length(x)]
}))
trial <- sapply(as.list(triali), strsplit, "\\:")
Timecourse.lincs.f1$MSTRG_Tx_ID <- unlist(sapply(trial, "[[", 1))

trial <- sapply(as.list(Timecourse.lincs.f1.clean$V4), strsplit, "\\|")
triali <- unlist(sapply(trial, function(x){
  x[grepl("\\:", x)]
}))
trial <- sapply(as.list(triali), strsplit, "\\:")
Timecourse.lincs.f1.clean$MSTRG_Tx_ID <- unlist(sapply(trial, "[[", 1))

trial <- sapply(as.list(Timecourse.lincs.f1.pure$V4), strsplit, "\\|")
triali <- unlist(sapply(trial, function(x){
  x[grepl("\\:", x)]
}))
trial <- sapply(as.list(triali), strsplit, "\\:")
Timecourse.lincs.f1.pure$MSTRG_Tx_ID <- unlist(sapply(trial, "[[", 1))

PLAR_Prefilter$linc_pred_level[PLAR_Prefilter$MSTRG_Tx_ID %in% Timecourse.lincs.f1$MSTRG_Tx_ID] <- "potential_linc"
PLAR_Prefilter$linc_pred_level[PLAR_Prefilter$MSTRG_Tx_ID %in% Timecourse.lincs.f1.clean$MSTRG_Tx_ID] <- "clean_linc"
PLAR_Prefilter$linc_pred_level[PLAR_Prefilter$MSTRG_Tx_ID %in% Timecourse.lincs.f1.pure$MSTRG_Tx_ID] <- "pure_linc"

colnames(PLAR_Prefilter)[11:14] <- c("Exon_count", "Spliced_length", "PLAR_Prefilter", "Tx_FPKM")
trial <- PLAR_Prefilter[,c(2:9,1,10:12,14,13,15)]

#write.csv(trial, "PLAR_Timecourse4points_2025.csv", row.names = F)

#### add in DE info ####

trial <- read.csv("PLAR_Timecourse4points_2025.csv")

PLAR_Timecourse_4Timepoints_DEnonDE <- merge(trial, Timecourse_4Timepoints_DEnonDE, by.x = "MSTRG_ID", by.y = "ENSEMBL")
PLAR_Timecourse_4Timepoints_DEnonDE <- PLAR_Timecourse_4Timepoints_DEnonDE[,c(1:8,16:53,9:15)] 

#add average isoform fpkm:
isoforms<-read.table(iso_filenames[1],header=TRUE,sep="\t",stringsAsFactors = FALSE)[,1]

#extract isoform FPKM info
iso_fpkm <-do.call(cbind,lapply(iso_filenames,function(fn)read.table(fn,header=TRUE,sep="\t",stringsAsFactors = FALSE)[,7]))
iso_fpkm <- data.frame(isoforms, iso_fpkm, stringsAsFactors = FALSE)
colnames(iso_fpkm)<-c("ENSEMBL",actualnames)

fpkm_list_ctrl <- as.list(as.data.frame(t(iso_fpkm[,2:5])))
fpkm_list_pd <- as.list(as.data.frame(t(iso_fpkm[,6:9])))
fpkm_list_il <- as.list(as.data.frame(t(iso_fpkm[,10:13])))
fpkm_list_bo <- as.list(as.data.frame(t(iso_fpkm[,14:17])))

fpkm_mean_ctrl <- sapply(fpkm_list_ctrl, mean)
fpkm_mean_pd <- sapply(fpkm_list_pd, mean)
fpkm_mean_il <- sapply(fpkm_list_il, mean)
fpkm_mean_bo <- sapply(fpkm_list_bo, mean)

fpkm_mean_treatment <- data.frame(fpkm_mean_ctrl, fpkm_mean_pd, fpkm_mean_il, fpkm_mean_bo)

trial <- as.list(as.data.frame(t(fpkm_mean_treatment)))
fpkm_max_treatment <- as.numeric(sapply(trial, max))

iso_fpkm <- cbind(iso_fpkm,fpkm_mean_treatment, fpkm_max_treatment)
colnames(iso_fpkm)[22] <- "isoform_fpkm_max_treatment"

PLAR_Timecourse_4Timepoints_DEnonDE <- merge(PLAR_Timecourse_4Timepoints_DEnonDE, iso_fpkm[,c(1,22)], by.x = "MSTRG_Tx_ID", by.y = "ENSEMBL")
PLAR_Timecourse_4Timepoints_DEnonDE <- PLAR_Timecourse_4Timepoints_DEnonDE[,c(2:47,1,48:51,54,52:53)]
colnames(PLAR_Timecourse_4Timepoints_DEnonDE)[c(51:52)] <- c("Tx_Max", "Tx_Max_Average")

#quick check on lncs of interest:
filter(PLAR_Timecourse_4Timepoints_DEnonDE, EnsName == "SMILR")
filter(PLAR_Timecourse_4Timepoints_DEnonDE, EnsName == "AC002480.3")


#write.csv(PLAR_Timecourse_4Timepoints_DEnonDE, "PLAR_Timecourse_4Timepoints_DEnonDE_2026.csv", row.names = F)

