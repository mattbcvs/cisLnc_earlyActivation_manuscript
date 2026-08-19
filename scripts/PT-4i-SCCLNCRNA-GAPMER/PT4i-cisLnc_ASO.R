library(tximport)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(reshape2)
library(clusterProfiler)
library(org.Hs.eg.db)

#
##### Normalise counts using FPKM, generate table: ####

#FPKM is a standard way to normalise RNAseq reads, explainer here: http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

#not for identifying differentially expressed genes but necessary for later reference + clustering

#get file names
filenames <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/Matt cisLnc project/RNAseq_IPLIL6_MSTRG12913/RSEM_cisLnc", pattern = "*genes.results", full.names = TRUE)

samplenames <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/Matt cisLnc project/RNAseq_IPLIL6_MSTRG12913/RSEM_cisLnc", pattern = "*genes.results")
samplenames <- sapply(strsplit(samplenames, "_"), "[[", 3)
samplenames <- sapply(strsplit(samplenames, "\\."), "[[", 1)

#read in an object containing gene names
genes<-read.table(filenames[1],header=TRUE,sep="\t",stringsAsFactors = FALSE)[,1]
#for each sample, read in columns containing FPKM
fpkm <-do.call(cbind,lapply(filenames,function(fn)read.table(fn,header=TRUE,sep="\t",stringsAsFactors = FALSE)[,7]))
#build table out of the above
fpkm <- data.frame(genes,fpkm,stringsAsFactors = FALSE)
colnames(fpkm)<-c("ENSEMBL",samplenames)

#add in mean 
fpkm$averagefpkm_GapNeg <- rowMeans(fpkm[,grep("GapNeg", colnames(fpkm))])
fpkm$averagefpkm_12913 <- rowMeans(fpkm[,grep("12913", colnames(fpkm))])
fpkm$averagefpkm_IPLIL6 <- rowMeans(fpkm[,grep("IPLIL6", colnames(fpkm))])
fpkm$fpkm_max_treatment <- rowMax(as.matrix(fpkm[,14:17]))

colnames(fpkm)

#provide no. genes expressed >1 FPKM
length(unique(filter(fpkm, fpkm_max_treatment >1)$ENSEMBL)) #13822 (prev. 16,589) genes expressed >1 FPKM


#
#### sample info table ####

cisLnc_meta <- data.frame("SampleNames" = samplenames, "Patient" = c(rep("P176", 3),
                                                                     rep("P182", 3),
                                                                     rep("P248", 3),
                                                                     rep("P323", 2),
                                                                     rep("P330", 2)),
                          "Treatment" = c(rep(c("GapNeg", "GapIPLIL6", "Gap12913"), 3),
                                          rep(c("GapNeg", "Gap12913"), 2))
                          )
rownames(cisLnc_meta) <- cisLnc_meta$SampleNames


#
#### deseq2 for 12913 knockdown ####

#get filenames for isoform results
fileenamesIso <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/Matt cisLnc project/RNAseq_IPLIL6_MSTRG12913/RSEM_cisLnc/", pattern = "*isoforms.results", full.names = TRUE)

#read in list of tx
tx2gene <- read.csv(fileenamesIso[1], sep = "\t", stringsAsFactors = F) #237079 transcripts

#create object to feed into DESeq2, note specific handling of rsem results
txi <- tximport(fileenamesIso[c(1,3,4,6,7,9:13)],
                type = "rsem", tx2gene = tx2gene[,1:2])

#design option from prelim analysis
dds <- DESeqDataSetFromTximport(txi, cisLnc_meta[c(1,3,4,6,7,9:13),]
                                , design= ~Patient + Treatment)
dds <- dds[rowSums(counts(dds))>10, ]
dds <- DESeq(dds)

vst <- rlog(dds, blind=FALSE)
vst_assay <- assay(vst)

data <- plotPCA(vst,intgroup=c("Treatment"),ntop=500,returnData=TRUE)
names = rownames(colData(vst))
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  geom_text(aes(label=Patient),hjust=0.25, vjust=-0.5, show.legend = F)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #geom_text(aes(label=names),hjust=0.25, vjust=-0.5, show.legend = F)+
  #scale_x_continuous(limits = c(-10, 11)) +
  theme_bw()

ggplot(data, aes(PC1, PC2, color=as.factor(Patient))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #geom_text(aes(label=names),hjust=0.25, vjust=-0.5, show.legend = F)+
  #scale_x_continuous(limits = c(-10, 11)) +
  theme_bw()
#grouped by patient primarily, but consistent effects from 12913

#test this as a batch effect, limma to remove
vst_limma <- vst
mm <- model.matrix(~Treatment, colData(vst_limma))
assay(vst_limma) <- limma::removeBatchEffect(assay(vst_limma), vst_limma$Patient, design = mm)
data$group <- as.character(data$group)
data$group[data$group == "GapNeg"] <- "Gap-ve"
data$group[data$group == "Gap12913"] <- "GapFOXERLNC1"

ggplot(data, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #geom_text(aes(label=names),hjust=0.25, vjust=-0.5, show.legend = F)+
  #scale_x_continuous(limits = c(-10, 11)) +
  theme_bw()

data <- plotPCA(vst_limma,intgroup=c("Treatment"),ntop=500,returnData=TRUE)
names = rownames(colData(vst_limma))
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(`Gap-ve` = "#00BA38", `GapFOXERLNC1` = "#619CFF")) +
  #geom_text(aes(label=names),hjust=0.25, vjust=-0.5, show.legend = F)+
  #scale_x_continuous(limits = c(-10, 11)) +
  theme_bw() +
  theme(text = element_text(size=24))

ggplot(data, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=Patient),hjust=0.25, vjust=-0.5, show.legend = F)+
  #scale_x_continuous(limits = c(-10, 11)) +
  theme_bw()

results(dds)
resultsNames(dds)

ds_12913 <- as.data.frame(results(dds, contrast = c("Treatment", "Gap12913", "GapNeg")))
ds_12913$EnsID <- rownames(ds_12913)

#gtf used for mapping genv26 + novels, merged EnsIDs split up
trial <- read.delim("gencodes.v26.primary_assembly.annotation.gtf", header=FALSE, stringsAsFactors=FALSE, sep = "\t", skip = 5)
trial <- filter(trial, V3 == "gene")
trial$EnsID <- sapply(strsplit(trial$V9, ";"), "[[", 1)
trial$EnsID <- gsub("gene_id ", "", trial$EnsID)

trial$EnsType <- sapply(strsplit(trial$V9, ";"), "[[", 2)
trial$EnsType <- gsub("gene_type ", "", trial$EnsType)

trial$EnsName <- sapply(strsplit(trial$V9, ";"), "[[", 3)
trial$EnsName <- gsub("gene_name ", "", trial$EnsName)

triali <- trial[,c(10:12,1,4,5,7)]
rm(trial)
length(unique(triali$EnsID))
#format into a useful .gtf

SMCPLAR_genv26_gtf_table <- read.csv("PLAR_3VSMC_fullGenv26_2026.csv", header = T)

trialii <- merge(triali, fpkm, by.x = "EnsID", by.y = "ENSEMBL")
trialii <- merge(trialii, ds_12913[c(2,6,7)], by.x = "EnsID", by.y = "EnsID", all.x = T)
ccLnc_12913_SMCpheno <- trialii
colnames(ccLnc_12913_SMCpheno)

#write.csv(ccLnc_12913_SMCpheno, "ccLnc_12913_SMCpheno_gen26.csv", row.names = F)


#
#### (start from here to make figures) expressed and DE genes ####

ccLnc_12913_SMCpheno <- read.csv("ccLnc_12913_SMCpheno_gen26.csv")

ccLnc_12913_fpkm_allG <- unique(filter(ccLnc_12913_SMCpheno, fpkm_max_treatment >1))

ccLnc_12913_fpkm_allG$EnsName <- gsub(" ", "", ccLnc_12913_fpkm_allG$EnsName)

#163x induced genes:
ccLnc_12913_fpkm_allG_up <- filter(ccLnc_12913_fpkm_allG, padj <0.05, log2FoldChange > log2(1.25))
dim(ccLnc_12913_fpkm_allG_up)
#150x repressed genes:
ccLnc_12913_fpkm_allG_down <- filter(ccLnc_12913_fpkm_allG, padj <0.05, log2FoldChange < -log2(1.25))
dim(ccLnc_12913_fpkm_allG_down)
#relaxing slightly as HMGB2 is of high interest as preferentially binding to SCClnc proms
ccLnc_12913_fpkm_allG_DE <- filter(ccLnc_12913_fpkm_allG, padj <0.05, abs(log2FoldChange) > log2(1.25))

#write.csv(ccLnc_12913_fpkm_allG_DE, "ccLnc_12913_fpkm_allG_DE.csv", row.names = F)


#
#### volcano + outliers ####

#highlight top FC/CPM
#select interesting genes in up
ccLnc_12913_fpkm_allG_up <- ccLnc_12913_fpkm_allG_up[order(ccLnc_12913_fpkm_allG_up$log2FoldChange, decreasing = T),]
ccLnc_12913_fpkm_allG_up$FCrank <- 1:length(ccLnc_12913_fpkm_allG_up$EnsName)
ccLnc_12913_fpkm_allG_up <- ccLnc_12913_fpkm_allG_up[order(ccLnc_12913_fpkm_allG_up$fpkm_max_treatment, decreasing = T),]
ccLnc_12913_fpkm_allG_up$CPMrank <- 1:length(ccLnc_12913_fpkm_allG_up$EnsName)

#average rank in both:
ccLnc_12913_fpkm_allG_up$DErank <- rowMeans(ccLnc_12913_fpkm_allG_up[,27:28])

#sum rank from both:
#ccLnc_12913_fpkm_allG_up$DErank <- ccLnc_12913_fpkm_allG_up[,27]*ccLnc_12913_fpkm_allG_up[,28]
ccLnc_12913_fpkm_allG_up <- ccLnc_12913_fpkm_allG_up[order(ccLnc_12913_fpkm_allG_up$DErank, decreasing = F),]
ccLnc_12913_fpkm_allG_up$DErank2 <- 1:length(ccLnc_12913_fpkm_allG_up$EnsName)
ccLnc_12913_fpkm_allG_up$upDown <- "UpRegWithGapmer(Upstream_M.12913_activation_represses)"

#now down
ccLnc_12913_fpkm_allG_down <- ccLnc_12913_fpkm_allG_down[order(ccLnc_12913_fpkm_allG_down$log2FoldChange, decreasing = F),]#note change this to F for repressed
ccLnc_12913_fpkm_allG_down$FCrank <- 1:length(ccLnc_12913_fpkm_allG_down$EnsName)
ccLnc_12913_fpkm_allG_down <- ccLnc_12913_fpkm_allG_down[order(ccLnc_12913_fpkm_allG_down$fpkm_max_treatment, decreasing = T),]
ccLnc_12913_fpkm_allG_down$CPMrank <- 1:length(ccLnc_12913_fpkm_allG_down$EnsName)

#average rank in both:
ccLnc_12913_fpkm_allG_down$DErank <- rowMeans(ccLnc_12913_fpkm_allG_down[,27:28])

#sum rank in both:
#ccLnc_12913_fpkm_allG_down$DErank <- ccLnc_12913_fpkm_allG_down[,24]*ccLnc_12913_fpkm_allG_down[,25]
ccLnc_12913_fpkm_allG_down <- ccLnc_12913_fpkm_allG_down[order(ccLnc_12913_fpkm_allG_down$DErank, decreasing = F),]
ccLnc_12913_fpkm_allG_down$DErank2 <- 1:length(ccLnc_12913_fpkm_allG_down$EnsName)
ccLnc_12913_fpkm_allG_down$upDown <- "DownRegWithGapmer(Upstream_M.12913_activation_induces)"

ccLnc_12913_fpkm_allG_DE <- rbind(ccLnc_12913_fpkm_allG_up, ccLnc_12913_fpkm_allG_down)

#highlight top FC/CPM
library(ggrepel)

ggplot(filter(ccLnc_12913_fpkm_allG_DE)) + aes(x = fpkm_max_treatment, y = log2FoldChange) +
  geom_point(alpha = 0.5) +
  geom_label_repel(data = filter(ccLnc_12913_fpkm_allG_DE, DErank2 <=15 #| fpkm_max_treatment >200
                                 ), 
                   aes(x = fpkm_max_treatment, y = log2FoldChange, label = EnsName), size = 3, force = 20, nudge_x = 0.4, nudge_y = 0.2) +
  theme_bw() +
  theme(text = element_text(size = 22)) +
  scale_x_log10() +
  #coord_cartesian(xlim = c(4,12)) +
  #geom_label_repel(data = filter(kuppe_postMI_edgeR_CMECs_DE, !is.na(label)),
  #                 aes(x = logCPM, y = log2FoldChange, label = label), force = 150, nudge_x = 13, nudge_y = 4) +
  xlab("fpkm_max_treatment") +
  ylab("LFC Gap12913 vs.\nGapNeg")


#highlight top FC/p
#select interesting genes in up
ccLnc_12913_fpkm_allG_up <- ccLnc_12913_fpkm_allG_up[order(ccLnc_12913_fpkm_allG_up$log2FoldChange, decreasing = T),]
ccLnc_12913_fpkm_allG_up$FCrank <- 1:length(ccLnc_12913_fpkm_allG_up$EnsName)
ccLnc_12913_fpkm_allG_up <- ccLnc_12913_fpkm_allG_up[order(ccLnc_12913_fpkm_allG_up$padj, decreasing = F),]
ccLnc_12913_fpkm_allG_up$Prank <- 1:length(ccLnc_12913_fpkm_allG_up$EnsName)

#average rank in both:
ccLnc_12913_fpkm_allG_up$DErank <- rowMeans(ccLnc_12913_fpkm_allG_up[,c(27,32)])

#sum rank from both:
#ccLnc_12913_fpkm_allG_up$DErank <- ccLnc_12913_fpkm_allG_up[,27]*ccLnc_12913_fpkm_allG_up[,28]
ccLnc_12913_fpkm_allG_up <- ccLnc_12913_fpkm_allG_up[order(ccLnc_12913_fpkm_allG_up$DErank, decreasing = F),]
ccLnc_12913_fpkm_allG_up$DErank2 <- 1:length(ccLnc_12913_fpkm_allG_up$EnsName)
ccLnc_12913_fpkm_allG_up$upDown <- "UpRegWithGapmer(Upstream_M.12913_activation_represses)"

#now down
ccLnc_12913_fpkm_allG_down <- ccLnc_12913_fpkm_allG_down[order(ccLnc_12913_fpkm_allG_down$log2FoldChange, decreasing = F),]#note change this to F for repressed
ccLnc_12913_fpkm_allG_down$FCrank <- 1:length(ccLnc_12913_fpkm_allG_down$EnsName)
ccLnc_12913_fpkm_allG_down <- ccLnc_12913_fpkm_allG_down[order(ccLnc_12913_fpkm_allG_down$padj, decreasing = F),]
ccLnc_12913_fpkm_allG_down$Prank <- 1:length(ccLnc_12913_fpkm_allG_down$EnsName)

#average rank in both:
ccLnc_12913_fpkm_allG_down$DErank <- rowMeans(ccLnc_12913_fpkm_allG_down[,c(27,32)])

#sum rank in both:
#ccLnc_12913_fpkm_allG_down$DErank <- ccLnc_12913_fpkm_allG_down[,24]*ccLnc_12913_fpkm_allG_down[,25]
ccLnc_12913_fpkm_allG_down <- ccLnc_12913_fpkm_allG_down[order(ccLnc_12913_fpkm_allG_down$DErank, decreasing = F),]
ccLnc_12913_fpkm_allG_down$DErank2 <- 1:length(ccLnc_12913_fpkm_allG_down$EnsName)
ccLnc_12913_fpkm_allG_down$upDown <- "DownRegWithGapmer(Upstream_M.12913_activation_induces)"

ccLnc_12913_fpkm_allG_DE <- rbind(ccLnc_12913_fpkm_allG_up, ccLnc_12913_fpkm_allG_down)

#highlight top FC/p
library(ggrepel)

ggplot(filter(ccLnc_12913_fpkm_allG_DE)) + aes(x = log2FoldChange, y = -log10(padj)) +
  geom_point(alpha = 0.5) +
  geom_label_repel(data = filter(ccLnc_12913_fpkm_allG_DE, DErank2 <=15 #| fpkm_max_treatment >200
  ), 
  aes(x = log2FoldChange, y = -log10(padj), label = EnsName), size = 3, force = 80, nudge_y = 6.2) +
  theme_bw() +
  theme(text = element_text(size = 22)) +
  #scale_x_log10() +
  #coord_cartesian(xlim = c(4,12)) +
  #geom_label_repel(data = filter(kuppe_postMI_edgeR_CMECs_DE, !is.na(label)),
  #                 aes(x = logCPM, y = log2FoldChange, label = label), force = 150, nudge_x = 13, nudge_y = 4) +
  xlab("fpkm_max_treatment") +
  ylab("LFC Gap12913 vs.\nGapNeg")


 #activation of IFN seems likely, are they on in 0-24hr anyway?
fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE_2026filt.csv", header = T)

#induced in 0-24hr
fpkm_allGDE_24hr_On <- filter(fpkm_allGDE, LogFC_0_24 > log2(1.5), preadj_0_24 < 0.05)
fpkm_allGDE_24hr_Off <- filter(fpkm_allGDE, LogFC_0_24 < -log2(1.5), preadj_0_24 < 0.05)
fpkm_allGDE_4hr_On <- filter(fpkm_allGDE, LogFC_0_4 > log2(1.5), preadj_0_4 < 0.05)
fpkm_allGDE_4hr_Off <- filter(fpkm_allGDE, LogFC_0_4 < -log2(1.5), preadj_0_4 < 0.05)

ccLnc_12913_fpkm_allG_DE$IPon <- "NoChange_Timecourse"
ccLnc_12913_fpkm_allG_DE$IPon[ccLnc_12913_fpkm_allG_DE$EnsID %in% fpkm_allGDE$EnsID] <- "OtherChange_Timecourse"
ccLnc_12913_fpkm_allG_DE$IPon[ccLnc_12913_fpkm_allG_DE$EnsID %in% fpkm_allGDE_24hr_On$EnsID] <- "Induction_0v24h"
ccLnc_12913_fpkm_allG_DE$IPon[ccLnc_12913_fpkm_allG_DE$EnsID %in% fpkm_allGDE_24hr_Off$EnsID] <- "Repression_0v24h"

table(ccLnc_12913_fpkm_allG_DE$IPon)
table(ccLnc_12913_fpkm_allG_DE$IPon, ccLnc_12913_fpkm_allG_DE$upDown)

ggplot(ccLnc_12913_fpkm_allG_DE) + aes(x = fpkm_max_treatment, y = log2FoldChange, color = IPon) +
  geom_point(alpha = 0.5, size= 3) +
  theme_bw() +
  theme(text = element_text(size = 22)) +
  scale_color_manual(values = c(NoChange_Timecourse = "grey70", Induction_0v24h = "firebrick4", Repression_0v24h = "steelblue4",
                                Induction_0v4h = "firebrick1", Repression_0v4h = "steelblue1")) +
  scale_x_log10() +
  xlab("fpkm_max_treatment") +
  ylab("LFC Gap12913 vs.\nGapNeg")

table(filter(fpkm_allGDE_24hr_On)$GeneClassUpdate)
table(filter(fpkm_allGDE_24hr_On, EnsID %in% ccLnc_12913_fpkm_allG_down$EnsID)$GeneClassUpdate)

ccLnc_12913_fpkm_allG_DE$IPon <- "NoChange_Timecourse"
ccLnc_12913_fpkm_allG_DE$IPon[ccLnc_12913_fpkm_allG_DE$EnsID %in% fpkm_allGDE_4hr_On$EnsID] <- "Induction_0v4h"
ccLnc_12913_fpkm_allG_DE$IPon[ccLnc_12913_fpkm_allG_DE$EnsID %in% fpkm_allGDE_4hr_Off$EnsID] <- "Repression_0v4h"

#takeaways so far
#lots of outlying genes suggest IFN response, but not IFNs are upregulated
#others suggest apoptosis: GAS5 and cell cycle repression
#for down reg outliers: lots of ox stress protection, ECM degradation

#ANXA2 and MGST1 are v notable - both in FOXL1 harmonizome/JASPAR set of 5.5k genes
#other downregs in JASPAR include STMN1, NQO1, GREM1
#but also IFIT1, 3, OAS1

#enrichment of JASPAR FOXL1, FOXC2 or FOXC1 motifs amongst up or down reg genes should be checked...

#of 0-4hr reg genes found here, bias for upreg to be "more up reg"

#less clear for 0-24hr genes found here

#in general most of the genes regged here are not found in the 0-24 reg set


#### Figures for the above (entrez match sorted) ####

#series of bar plots of OR + sig

#for all tables (2025 + earlier)

#for various DEGs: 0.05/1.25, 0.01/1.25, 0.05/1.25/5fpkm

FOXL1_JASPAR <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/Matt cisLnc project/RNAseq_IPLIL6_MSTRG12913/Harmonizome_JASPAR_FOXL1.txt",
                         header = F, sep = "\t")
FOXL1_JASPAR <- bitr(FOXL1_JASPAR$V3, fromType = "ENTREZID", toType = "ENSEMBL", OrgDb = "org.Hs.eg.db")

#2025 set appears a little more specific, only 2500:
FOXL1jasp25 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/Matt cisLnc project/RNAseq_IPLIL6_MSTRG12913/Harmonizome_JASPAR25_FOXL1.txt",
                         header = F, sep = "\t")[,-1]
colnames(FOXL1jasp25)[1] <- "V1"
FOXL1jasp25 <- bitr(FOXL1jasp25$V3, fromType = "ENTREZID", toType = "ENSEMBL", OrgDb = "org.Hs.eg.db")

FOXC1_JASPAR <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/Matt cisLnc project/RNAseq_IPLIL6_MSTRG12913/Harmonizome_JASPAR_FOXC1.txt",
                         header = F, sep = "\t")
FOXC1_JASPAR <- bitr(FOXC1_JASPAR$V3, fromType = "ENTREZID", toType = "ENSEMBL", OrgDb = "org.Hs.eg.db")

#2025 set appears way more specific:
FOXC1jasp25 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/Matt cisLnc project/RNAseq_IPLIL6_MSTRG12913/Harmonizome_JASPAR25_FOXC1.txt",
                         header = F, sep = "\t")[,-1]
colnames(FOXC1jasp25)[1] <- "V1"
FOXC1jasp25 <- bitr(FOXC1jasp25$V3, fromType = "ENTREZID", toType = "ENSEMBL", OrgDb = "org.Hs.eg.db")

FOXC2jasp25 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/Matt cisLnc project/RNAseq_IPLIL6_MSTRG12913/Harmonizome_JASPAR25_FOXC2.txt",
                         header = F, sep = "\t")[,-1]
colnames(FOXC2jasp25)[1] <- "V1"
FOXC2jasp25 <- bitr(FOXC2jasp25$V3, fromType = "ENTREZID", toType = "ENSEMBL", OrgDb = "org.Hs.eg.db")

test_array <- data.frame("TF_table" = c(rep("FOXL1", 3), rep("FOXL1jasp25", 3), 
                                        rep("FOXC1", 3), rep("FOXC1jasp25", 3),
                                        rep("FOXC2jasp25", 3)),
                         "DEG_type" = rep(c("0.05_1.25_1", "0.01_1.25_1", "0.05_1.25_5"),5)
                         )

test_TFs <- list("FOXL1" = FOXL1_JASPAR, "FOXL1jasp25" = FOXL1jasp25,
                 "FOXC1" = FOXC1_JASPAR, "FOXC1jasp25" = FOXC1jasp25,
                 "FOXC2jasp25" = FOXC2jasp25)

ccLnc_12913_fpkm_allG$EnsName <- gsub(" ", "", ccLnc_12913_fpkm_allG$EnsName)
ccLnc_12913_fpkm_allG_DE <- filter(ccLnc_12913_fpkm_allG, padj <0.05, abs(log2FoldChange) > log2(1.25))
ccLnc_12913_fpkm_allG_DE2 <- filter(ccLnc_12913_fpkm_allG, padj <0.01, abs(log2FoldChange) > log2(1.25))
ccLnc_12913_fpkm_allG_DE3 <- filter(ccLnc_12913_fpkm_allG, padj <0.05, abs(log2FoldChange) > log2(1.25), fpkm_max_treatment>5)

test_DEGs <- list("0.05_1.25_1" = ccLnc_12913_fpkm_allG_DE,
                  "0.01_1.25_1" = ccLnc_12913_fpkm_allG_DE2,
                  "0.05_1.25_5" = ccLnc_12913_fpkm_allG_DE3)

test_results <- list()

for (i in 1:length(test_array[,1])){
  
  x <- test_TFs[[test_array[i,1]]]
  y <- test_DEGs[[test_array[i,2]]]
  
  d <- length(unique(filter(ccLnc_12913_fpkm_allG)$EnsID))
  c <- length(unique(filter(ccLnc_12913_fpkm_allG, gsub("\\.[0-9]*", "", EnsID) %in% x$ENSEMBL)$EnsID))

  b <- length(unique(filter(y, log2FoldChange >0)$EnsID))
  a <- length(unique(filter(y, log2FoldChange >0, gsub("\\.[0-9]*", "", EnsID) %in% x$ENSEMBL)$EnsID))

  bi <- length(unique(filter(y, log2FoldChange <0)$EnsID))
  ai <- length(unique(filter(y, log2FoldChange <0, gsub("\\.[0-9]*", "", EnsID) %in% x$ENSEMBL)$EnsID))
  
  test_results[[i]] <- c(a,b,ai,bi,c,d,
                         a/b,ai/bi,c/d,
                         fisher.test(data.frame("DownRegGenes" = c(a,b-a),
                                                "AllGenes" = c(c-a,d-c-(b-a)), row.names = c("FOXL1-JASPAR", "other")), alternative = "greater")$est,
                         fisher.test(data.frame("DownRegGenes" = c(a,b-a),
                                                "AllGenes" = c(c-a,d-c-(b-a)), row.names = c("FOXL1-JASPAR", "other")), alternative = "greater")$p,
                         
                         fisher.test(data.frame("DownRegGenes" = c(ai,bi-ai),
                                                "AllGenes" = c(c-ai,d-c-(bi-ai)), row.names = c("FOXL1-JASPAR", "other")), alternative = "greater")$est,
                         fisher.test(data.frame("DownRegGenes" = c(ai,bi-ai),
                                                "AllGenes" = c(c-ai,d-c-(bi-ai)), row.names = c("FOXL1-JASPAR", "other")), alternative = "greater")$p
                         )
}

trial <- sapply(test_results, data.frame)
names(trial) <- paste(test_array[,1], test_array[,2], sep = "_")
FOX_TFs_12913 <- as.data.frame(t(bind_rows(trial)))
colnames(FOX_TFs_12913) <- c("a", "b", "ai", "bi", "c", "d","percUp", "percDown", "PercBackground", "OR_up", "p_up", "OR_down", "p_down")
  
FOX_TFs_12913$OR_up2 <- FOX_TFs_12913$OR_up
FOX_TFs_12913$OR_up2[FOX_TFs_12913$p_up > 0.05] <- NA

FOX_TFs_12913$OR_down2 <- FOX_TFs_12913$OR_down
FOX_TFs_12913$OR_down2[FOX_TFs_12913$p_down > 0.05] <- NA

FOX_TFs_12913$test <- rownames(FOX_TFs_12913)

plot_FOX_TFs_12913 <- reshape2::melt(FOX_TFs_12913[,c(14:16)])
plot_FOX_TFs_12913$p <- reshape2::melt(FOX_TFs_12913[,c(11,13,16)])[,3]

plot_FOX_TFs_12913$TF_table <- sapply(strsplit(plot_FOX_TFs_12913$test, "_"), "[[", 1)
plot_FOX_TFs_12913$DEGs <- paste(sapply(strsplit(plot_FOX_TFs_12913$test, "_"), "[[", 2), 
                                 sapply(strsplit(plot_FOX_TFs_12913$test, "_"), "[[", 3),
                                 sapply(strsplit(plot_FOX_TFs_12913$test, "_"), "[[", 4), sep = "_")

ggplot(plot_FOX_TFs_12913) + aes(x = value, y = -log10(p), color = variable) +
  geom_point() +
  facet_wrap(~TF_table)

ggplot(filter(plot_FOX_TFs_12913, TF_table == "FOXC1")) + aes(x = value, y = -log10(p), color = variable) +
  geom_point() +
  facet_wrap(~DEGs) + xlab("OR") +ylab("-log10(p)")

ggplot(filter(plot_FOX_TFs_12913, TF_table == "FOXL1")) + aes(x = value, y = -log10(p), color = variable) +
  geom_point() +
  facet_wrap(~DEGs) + xlab("OR") +ylab("-log10(p)")

ggplot(filter(plot_FOX_TFs_12913, TF_table == "FOXC1jasp25")) + aes(x = value, y = -log10(p), color = variable) +
  geom_point() +
  facet_wrap(~DEGs) + xlab("OR") +ylab("-log10(p)")

ggplot(filter(plot_FOX_TFs_12913, TF_table == "FOXL1jasp25")) + aes(x = value, y = -log10(p), color = variable) +
  geom_point() +
  facet_wrap(~DEGs) + xlab("OR") +ylab("-log10(p)")

ggplot(filter(plot_FOX_TFs_12913, TF_table == "FOXC2jasp25")) + aes(x = value, y = -log10(p), color = variable) +
  geom_point() +
  facet_wrap(~DEGs) + xlab("OR") +ylab("-log10(p)")

#dot plot of x axis UpReg Genes, + Hiconf, DownReg Genes, + HiConf, y axis FOXL1 motif, FOXC2 motif, FOXC1 motif (Jaspar '25)

plot_FOX_TFs_12913_dot <- filter(plot_FOX_TFs_12913, grepl("_1$", DEGs), grepl("jasp25", test))

plot_FOX_TFs_12913_dot$variable2[grepl("OR_up", plot_FOX_TFs_12913_dot$variable)] <- "Induced"
plot_FOX_TFs_12913_dot$variable2[grepl("OR_down", plot_FOX_TFs_12913_dot$variable)] <- "Repressed"
plot_FOX_TFs_12913_dot$variable2[grepl("OR_up", plot_FOX_TFs_12913_dot$variable) & grepl("0.01", plot_FOX_TFs_12913_dot$DEGs)] <- "Induced\n(p<0.01)"
plot_FOX_TFs_12913_dot$variable2[grepl("OR_down", plot_FOX_TFs_12913_dot$variable) & grepl("0.01", plot_FOX_TFs_12913_dot$DEGs)] <- "Repressed\n(p<0.01)"

plot_FOX_TFs_12913_dot$TF_table <- gsub("jasp25", " JASPAR motif", plot_FOX_TFs_12913_dot$TF_table)

plot_FOX_TFs_12913_dot$p_adjust <- p.adjust(plot_FOX_TFs_12913_dot$p, method = "BH")

plot_FOX_TFs_12913_dot$p[plot_FOX_TFs_12913_dot$p >0.05] <- NA

plot_FOX_TFs_12913_dot$p2 <- plot_FOX_TFs_12913_dot$p
plot_FOX_TFs_12913_dot$p2[plot_FOX_TFs_12913_dot$p < 0.05] <- "p<0.05"
plot_FOX_TFs_12913_dot$p2[plot_FOX_TFs_12913_dot$p < 0.01] <- "p<0.01"

ggplot(plot_FOX_TFs_12913_dot) + aes(x = variable2, y = TF_table, color = value, size = p2) +
  geom_point() + Seurat::RotatedAxis() +
  scale_size_discrete(range = c(10,15), limits = c("p<0.05", #"p<0.01","p<0.001", 
                                                  "p<0.01"#, "p<0.000001", "p<0.00000001"
  )) +
  xlab("GapFOXERLNC DEGs") +
  ylab("") +
  theme(text = element_text(size = 24))

#
#### GO/KEGG/REACTOME terms - 12913 induced ####

#for enrichr:
#write.csv(ccLnc_12913_fpkm_allG_DE, "ccLnc_12913_fpkm_allG_DE.csv", row.names = F)

library(clusterProfiler)
library(org.Hs.eg.db)
ccLnc_12913_fpkm_allG
ccLnc_12913_fpkm_allG$EnsID_merge <- gsub("\\.[0-9]*", "", ccLnc_12913_fpkm_allG$EnsID)
ccLnc_12913_fpkm_allG_up$EnsID_merge <- gsub("\\.[0-9]*", "", ccLnc_12913_fpkm_allG_up$EnsID)

fpkm_DEG_Gap12913up_GO <- enrichGO(gene          = ccLnc_12913_fpkm_allG_up$EnsID_merge,
                                   universe      = ccLnc_12913_fpkm_allG$EnsID_merge,
                                   keyType       = "ENSEMBL",
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   readable      = TRUE)
fpkm_DEG_Gap12913up_GO_df <- as.data.frame(fpkm_DEG_Gap12913up_GO)

dotplot(simplify(fpkm_DEG_Gap12913up_GO), showCategory = 10)

#strong set related to immune response bias, response to biotic stimulus interferon signalling

#KEGG
convertEnsEnt <- bitr(unique(ccLnc_12913_fpkm_allG$EnsName), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

fpkm_DEG_Gap12913up_K <- enrichKEGG(gene = filter(convertEnsEnt, SYMBOL %in% ccLnc_12913_fpkm_allG_up$EnsName)$ENTREZID,
                                     universe = filter(convertEnsEnt, SYMBOL %in% ccLnc_12913_fpkm_allG$EnsName)$ENTREZID,
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     qvalueCutoff  = 0.05)
fpkm_DEG_Gap12913up_K_df <- data.frame(fpkm_DEG_Gap12913up_K)
unique(fpkm_DEG_GapIPLIL6up_K_df$Description)
filter(convertEnsEnt, 
       ENTREZID %in% unlist(strsplit(filter(fpkm_DEG_GapIPLIL6up_K_df, Description == unique(fpkm_DEG_GapIPLIL6up_K_df$Description)[4])$geneID, 
                                     split = "/")))
#not IFN, mainly diseases

dotplot(fpkm_DEG_Gap12913up_K, showCategory = 20)

#REACTOME
library(ReactomePA)
fpkm_DEG_Gap12913up_R <- enrichPathway(gene          = unique(filter(convertEnsEnt,SYMBOL %in% ccLnc_12913_fpkm_allG_up$EnsName)$ENTREZID),
                                       universe      = unique(filter(convertEnsEnt,SYMBOL %in% ccLnc_12913_fpkm_allG$EnsName)$ENTREZID),
                                       organism = "human",
                                       pvalueCutoff = 0.05,
                                       qvalueCutoff  = 0.05,
                                       readable      = TRUE)
fpkm_DEG_Gap12913up_R_df <- as.data.frame(fpkm_DEG_Gap12913up_R)

dotplot(fpkm_DEG_Gap12913up_R, showCategory = 20)


#combined figure/table Go-S/K/R:
fpkm_DEG_Gap12913up_GO_S <- simplify(fpkm_DEG_Gap12913up_GO, cutoff = 0.6)#62 at 0.7 (2x virus right at top)
fpkm_DEG_Gap12913up_GO_S_df <- as.data.frame(fpkm_DEG_Gap12913up_GO_S)

#description, gene ratio, adj p val, genes
colnames(fpkm_DEG_Gap12913up_GO_S_df)
colnames(fpkm_DEG_Gap12913up_K_df)
colnames(fpkm_DEG_Gap12913up_R_df)

trial <- fpkm_DEG_Gap12913up_GO_S_df[,c(2,3,9,11)]
triali <- fpkm_DEG_Gap12913up_K_df[,c(4,5,11,13)]
trialii <- fpkm_DEG_Gap12913up_R_df[,c(2,3,9,11)]

trial$termType <- "GO"
triali$termType <- "KEGG"
trialii$termType <- "REACTOME"

combined_GOKR <- rbind(trial, triali, trialii)

#this can be saved for supplement

#plotting:
combined_GOKR$selectHits <- as.numeric(sapply(strsplit(combined_GOKR$GeneRatio, "\\/"), "[[", 1))
combined_GOKR$select <- as.numeric(sapply(strsplit(combined_GOKR$GeneRatio, "\\/"), "[[", 2))
combined_GOKR$GeneRatio <- combined_GOKR$selectHits/combined_GOKR$select*100

#top 10 GO, top 5 others
combined_GOKR_plot <- filter(combined_GOKR, (termType == "GO" & p.adjust < 1.883957e-12) | 
                               (termType == "KEGG" & p.adjust < 1e-5) | #tempting to cut all these "off topic" diseases 
                               (termType == "REACTOME" & p.adjust < 4.405474e-05))

combined_GOKR_plot$Description <- tolower(combined_GOKR_plot$Description)
combined_GOKR_plot$Description[grepl("modulation", combined_GOKR_plot$Description)] <- "host response modulation from ifn"
combined_GOKR_plot <- combined_GOKR_plot[order(combined_GOKR_plot$Description),]

combined_GOKR_plot$Description <- factor(combined_GOKR_plot$Description)
combined_GOKR_plot$Description <- factor(combined_GOKR_plot$Description,
                                         levels = levels(combined_GOKR_plot$Description)[order( 
                                           combined_GOKR_plot$termType, -combined_GOKR_plot$GeneRatio, decreasing = T)])

colnames(combined_GOKR_plot)[2] <- "% GapFOXERLNC1\ninduced Genes"

ggplot(combined_GOKR_plot) + aes(x =`% GapFOXERLNC1\ninduced Genes`, fill = -log10(p.adjust), y = Description) +
  geom_bar(stat = "identity", color = "grey40")+
  theme_bw() +
  scale_x_continuous(breaks = seq(0,30,10)) +
  scale_fill_continuous(limits = c(0, max(-log10(combined_GOKR_plot$p.adjust)))) +
  theme(text = element_text(size=24)) +
  ylab("")


#
#### GO/K/R 12913 repressed ####

ccLnc_12913_fpkm_allG_down$EnsID_merge <- gsub("\\.[0-9]*", "", ccLnc_12913_fpkm_allG_down$EnsID)

fpkm_DEG_Gap12913down_GO <- enrichGO(gene          = ccLnc_12913_fpkm_allG_down$EnsID_merge,
                                   universe      = ccLnc_12913_fpkm_allG$EnsID_merge,
                                   keyType       = "ENSEMBL",
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   readable      = TRUE)
fpkm_DEG_Gap12913down_GO_df <- as.data.frame(fpkm_DEG_Gap12913down_GO)

dotplot(simplify(fpkm_DEG_Gap12913down_GO), showCategory = 10)


#KEGG
convertEnsEnt <- bitr(unique(ccLnc_12913_fpkm_allG$EnsName), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

fpkm_DEG_Gap12913down_K <- enrichKEGG(gene = filter(convertEnsEnt, SYMBOL %in% ccLnc_12913_fpkm_allG_down$EnsName)$ENTREZID,
                                    universe = filter(convertEnsEnt, SYMBOL %in% ccLnc_12913_fpkm_allG$EnsName)$ENTREZID,
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.05,
                                    qvalueCutoff  = 0.05)
fpkm_DEG_Gap12913down_K_df <- data.frame(fpkm_DEG_Gap12913down_K)
unique(fpkm_DEG_Gap12913down_K_df$Description)
filter(convertEnsEnt, 
       ENTREZID %in% unlist(strsplit(filter(fpkm_DEG_Gap12913down_K_df, Description == unique(fpkm_DEG_Gap12913down_K_df$Description)[1])$geneID, 
                                     split = "/")))
#just cell cycle

dotplot(fpkm_DEG_Gap12913down_K, showCategory = 20)


#REACTOME
library(ReactomePA)
fpkm_DEG_Gap12913down_R <- enrichPathway(gene          = unique(filter(convertEnsEnt,SYMBOL %in% ccLnc_12913_fpkm_allG_down$EnsName)$ENTREZID),
                                       universe      = unique(filter(convertEnsEnt,SYMBOL %in% ccLnc_12913_fpkm_allG$EnsName)$ENTREZID),
                                       organism = "human",
                                       pvalueCutoff = 0.05,
                                       qvalueCutoff  = 0.05,
                                       readable      = TRUE)
fpkm_DEG_Gap12913down_R_df <- as.data.frame(fpkm_DEG_Gap12913down_R)

dotplot(fpkm_DEG_Gap12913down_R, showCategory = 20)

#cell cycle, some ECM too


#combined figure/table Go-S/K/R:
fpkm_DEG_Gap12913down_GO_S <- simplify(fpkm_DEG_Gap12913down_GO, cutoff = 0.6)#too many similar at default
fpkm_DEG_Gap12913down_GO_S_df <- as.data.frame(fpkm_DEG_Gap12913down_GO_S)

#description, gene ratio, adj p val, genes
colnames(fpkm_DEG_Gap12913down_GO_S_df)
colnames(fpkm_DEG_Gap12913down_K_df)
colnames(fpkm_DEG_Gap12913down_R_df)

trial <- fpkm_DEG_Gap12913down_GO_S_df[,c(2,3,9,11)]
triali <- fpkm_DEG_Gap12913down_K_df[,c(4,5,11,13)]
trialii <- fpkm_DEG_Gap12913down_R_df[,c(2,3,9,11)]

trial$termType <- "GO"
triali$termType <- "KEGG"
trialii$termType <- "REACTOME"

combined_GOKR <- rbind(trial, triali, trialii)

#this can be saved for supplement

#plotting:
combined_GOKR$selectHits <- as.numeric(sapply(strsplit(combined_GOKR$GeneRatio, "\\/"), "[[", 1))
combined_GOKR$select <- as.numeric(sapply(strsplit(combined_GOKR$GeneRatio, "\\/"), "[[", 2))
combined_GOKR$GeneRatio <- combined_GOKR$selectHits/combined_GOKR$select*100

#top 10 GO, top 5 others
combined_GOKR_plot <- filter(combined_GOKR, (termType == "GO" & p.adjust < 7.864304e-08) | 
                               (termType == "KEGG" & p.adjust < 1.243611e-03) | #tempting to cut all these "off topic" diseases 
                               (termType == "REACTOME" & p.adjust < 2.345271e-08))

combined_GOKR_plot$Description <- tolower(combined_GOKR_plot$Description)
combined_GOKR_plot <- combined_GOKR_plot[order(combined_GOKR_plot$Description),]

combined_GOKR_plot$Description <- factor(combined_GOKR_plot$Description)
combined_GOKR_plot$Description <- factor(combined_GOKR_plot$Description,
                                         levels = levels(combined_GOKR_plot$Description)[order( 
                                           combined_GOKR_plot$termType, -combined_GOKR_plot$GeneRatio, decreasing = T)])

colnames(combined_GOKR_plot)[2] <- "% GapFOXERLNC1\nrepressed Genes"

ggplot(combined_GOKR_plot) + aes(x =`% GapFOXERLNC1\nrepressed Genes`, fill = -log10(p.adjust), y = Description) +
  geom_bar(stat = "identity", color = "grey40")+
  theme_bw() +
  scale_x_continuous(breaks = seq(0,30,10)) +
  scale_fill_continuous(limits = c(0, max(-log10(combined_GOKR_plot$p.adjust)))) +
  theme(text = element_text(size=24)) +
  ylab("")


#now with FOXL1 and FOXC2 JASPAR genes:



#
#### final volcano with some highlighted genes ####

#stick to those found in GO/K/R + outliers on the plot

#CDH13 - massive outlier, anti-apoptotic, adiponectin receptor, loss indicates SMC damage - needs a ref tho as not in GO/K/R
Notable <- c("RSAD2", "IFIT3", "IFIT2", "IL23A", "GBP4",
             
             #genes indicating loss of redox protection
             "NQO1", "APOE", "TXNRD1",
             
             #key genes indicating loss of proliferation
             "LGMN", "MKI67", "RRM2", "CDK1", "TOP2A")

ccLnc_12913_fpkm_allG_DE$Notable <- NA
ccLnc_12913_fpkm_allG_DE$Notable[ccLnc_12913_fpkm_allG_DE$EnsName %in% c("RSAD2", "IFIT3", "IFIT2", "IL23A", "GBP4")] <- "IFN-related"
ccLnc_12913_fpkm_allG_DE$Notable[ccLnc_12913_fpkm_allG_DE$EnsName %in% c("LGMN", "MKI67", "RRM2", "CDK1", "TOP2A")] <- "Proliferation"

ggplot(filter(ccLnc_12913_fpkm_allG_DE)) + aes(x = log2FoldChange, y = -log10(padj), color = Notable) +
  geom_point(data = ccLnc_12913_fpkm_allG_DE[is.na(ccLnc_12913_fpkm_allG_DE$Notable),], aes(colour = Notable), alpha = 0.2) +
  geom_point(data = ccLnc_12913_fpkm_allG_DE[!is.na(ccLnc_12913_fpkm_allG_DE$Notable),], aes(colour = Notable), alpha = 1,size=2) +
  #geom_label_repel(data = filter(ccLnc_12913_fpkm_allG_DE, !is.na(Notable)), 
  #                 aes(x = log2FoldChange, y = -log10(padj), label = EnsName), size = 3.5, force = 10, nudge_y = 10.2) +
  theme_bw() +
  theme(text = element_text(size = 22)) +
  #scale_x_log10() +
  #coord_cartesian(xlim = c(4,12)) +
  #geom_label_repel(data = filter(kuppe_postMI_edgeR_CMECs_DE, !is.na(label)),
  #                 aes(x = logCPM, y = log2FoldChange, label = label), force = 150, nudge_x = 13, nudge_y = 4) +
  xlab("LFC GapFOXERLNC1 vs.\nGapNeg") +
  ylab("-log10(p)")


ccLnc_12913_fpkm_allG_DE$Motif <- NA
ccLnc_12913_fpkm_allG_DE$Motif[gsub("\\.[0-9]*", "", ccLnc_12913_fpkm_allG_DE$EnsID) %in% FOXL1jasp25$ENSEMBL] <- "FOXL1"
ccLnc_12913_fpkm_allG_DE$Motif[gsub("\\.[0-9]*", "", ccLnc_12913_fpkm_allG_DE$EnsID) %in% FOXC2jasp25$ENSEMBL] <- "FOXC2"
#ccLnc_12913_fpkm_allG_DE$Motif[gsub("\\.[0-9]*", "", ccLnc_12913_fpkm_allG_DE$EnsID) %in% FOXL1jasp25$ENSEMBL & 
#                                 gsub("\\.[0-9]*", "", ccLnc_12913_fpkm_allG_DE$EnsID) %in% FOXC2jasp25$ENSEMBL] <- "Both"

table(ccLnc_12913_fpkm_allG_DE$Motif)
table(filter(ccLnc_12913_fpkm_allG_DE, log2FoldChange <0)$Motif)

ggplot(filter(ccLnc_12913_fpkm_allG_DE)) + aes(x = log2FoldChange, y = -log10(padj), color = Motif) +
  geom_point(data = ccLnc_12913_fpkm_allG_DE[is.na(ccLnc_12913_fpkm_allG_DE$Motif),], aes(colour = Motif), alpha = 0.45) +
  geom_point(data = ccLnc_12913_fpkm_allG_DE[!is.na(ccLnc_12913_fpkm_allG_DE$Motif),], aes(colour = Motif), alpha = 0.8,size=2) +
  #geom_label_repel(data = filter(ccLnc_12913_fpkm_allG_DE, !is.na(Motif)), 
  #                 aes(x = log2FoldChange, y = -log10(padj), label = EnsName), size = 3.5, force = 100, nudge_y = 10.2) +
  theme_bw() +
  theme(text = element_text(size = 22)) +
  geom_hline(yintercept = -log10(c(0.05, 0.01)), linetype = "dashed", color = "grey40") +
  #scale_x_log10() +
  #coord_cartesian(xlim = c(4,12)) +
  #geom_label_repel(data = filter(kuppe_postMI_edgeR_CMECs_DE, !is.na(label)),
  #                 aes(x = logCPM, y = log2FoldChange, label = label), force = 150, nudge_x = 13, nudge_y = 4) +
  xlab("LFC GapFOXERLNC1 vs.\nGapNeg") +
  ylab("-log10(p)")

#
#### deseq2 for ipl-il6 knockdown ####

#get filenames for isoform results
fileenamesIso <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/Baker-lab/Matt cisLnc project/RNAseq_IPLIL6_MSTRGIPLIL6/RSEM_cisLnc/", pattern = "*isoforms.results", full.names = TRUE)

#read in list of tx
tx2gene <- read.csv(fileenamesIso[1], sep = "\t", stringsAsFactors = F) #237079 transcripts

#create object to feed into DESeq2, note specific handling of rsem results
txi <- tximport(fileenamesIso[c(1,2,4,5,7,8)],
                type = "rsem", tx2gene = tx2gene[,1:2])

#design option from prelim analysis
dds <- DESeqDataSetFromTximport(txi, cisLnc_meta[c(1,2,4,5,7,8),]
                                , design= ~Patient + Treatment)
dds <- dds[rowSums(counts(dds))>10, ]
dds <- DESeq(dds)

vst <- rlog(dds, blind=FALSE)
vst_assay <- assay(vst)

data <- plotPCA(vst,intgroup=c("Treatment"),ntop=500,returnData=TRUE)
names = rownames(colData(vst))
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=3) +
  geom_text(aes(label=Patient),hjust=0.25, vjust=-0.5, show.legend = F)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #geom_text(aes(label=names),hjust=0.25, vjust=-0.5, show.legend = F)+
  #scale_x_continuous(limits = c(-10, 11)) +
  theme_bw()

ggplot(data, aes(PC1, PC2, color=as.factor(Patient))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #geom_text(aes(label=names),hjust=0.25, vjust=-0.5, show.legend = F)+
  #scale_x_continuous(limits = c(-10, 11)) +
  theme_bw()
#grouped by patient primarily, but consistent effects from IPLIL6

#test this as a batch effect, limma to remove
vst_limma <- vst
mm <- model.matrix(~Treatment, colData(vst_limma))
assay(vst_limma) <- limma::removeBatchEffect(assay(vst_limma), vst_limma$Patient, design = mm)

data <- plotPCA(vst_limma,intgroup=c("Treatment"),ntop=500,returnData=TRUE)
names = rownames(colData(vst_limma))
percentVar <- round(100 * attr(data, "percentVar"))
data$group <- as.character(data$group)
data$group[data$group == "GapNeg"] <- "Gap-ve"
data$group[data$group == "GapIPLIL6"] <- "GapIPL-IL6"

ggplot(data, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #geom_text(aes(label=names),hjust=0.25, vjust=-0.5, show.legend = F)+
  #scale_x_continuous(limits = c(-10, 11)) +
  theme_bw()

ggplot(data, aes(PC1, PC2, color=group)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(`Gap-ve` = "#00BA38", `GapIPL-IL6` = "#619CFF")) +
  #geom_text(aes(label=Patient),hjust=0.25, vjust=-0.5, show.legend = F)+
  #scale_x_continuous(limits = c(-10, 11)) +
  theme_bw() +
  theme(text = element_text(size=24))

results(dds)
resultsNames(dds)

ds_IPLIL6 <- as.data.frame(results(dds, contrast = c("Treatment", "GapIPLIL6", "GapNeg")))
ds_IPLIL6$EnsID <- rownames(ds_IPLIL6)

#gtf used for mapping genv26 + novels, merged EnsIDs split up
trial <- read.delim("gencode.v26.primary_assembly.annotation.gtf", header=FALSE, stringsAsFactors=FALSE, sep = "\t", skip = 5)
trial <- filter(trial, V3 == "gene")
trial$EnsID <- sapply(strsplit(trial$V9, ";"), "[[", 1)
trial$EnsID <- gsub("gene_id ", "", trial$EnsID)

trial$EnsType <- sapply(strsplit(trial$V9, ";"), "[[", 2)
trial$EnsType <- gsub("gene_type ", "", trial$EnsType)

trial$EnsName <- sapply(strsplit(trial$V9, ";"), "[[", 3)
trial$EnsName <- gsub("gene_name ", "", trial$EnsName)

triali <- trial[,c(10:12,1,4,5,7)]
rm(trial)
length(unique(triali$EnsID))
#format into a useful .gtf

trialii <- merge(triali, fpkm, by.x = "EnsID", by.y = "ENSEMBL")
trialii <- merge(trialii, ds_IPLIL6[c(2,6,7)], by.x = "EnsID", by.y = "EnsID", all.x = T)
ccLnc_IPLIL6_SMCpheno <- trialii
colnames(ccLnc_IPLIL6_SMCpheno)

#write.csv(ccLnc_IPLIL6_SMCpheno, "ccLnc_IPLIL6_SMCpheno_gen26.csv", row.names = F)


#
#### (start here for figures) expressed and DE genes ####

ccLnc_IPLIL6_SMCpheno <- read.csv("ccLnc_IPLIL6_SMCpheno_gen26.csv")

ccLnc_IPLIL6_fpkm_allG <- unique(filter(ccLnc_IPLIL6_SMCpheno, fpkm_max_treatment >1))

ccLnc_IPLIL6_fpkm_allG$EnsName <- gsub(" ", "", ccLnc_IPLIL6_fpkm_allG$EnsName)

#140x induced genes:
ccLnc_IPLIL6_fpkm_allG_up <- filter(ccLnc_IPLIL6_fpkm_allG, padj <0.05, log2FoldChange > log2(1.25))
dim(ccLnc_IPLIL6_fpkm_allG_up)
#242x repressed genes:
ccLnc_IPLIL6_fpkm_allG_down <- filter(ccLnc_IPLIL6_fpkm_allG, padj <0.05, log2FoldChange < -log2(1.25))
dim(ccLnc_IPLIL6_fpkm_allG_down)

ccLnc_IPLIL6_fpkm_allG_DE <- filter(ccLnc_IPLIL6_fpkm_allG, padj <0.05, abs(log2FoldChange) > log2(1.25))

#write.csv(ccLnc_IPLIL6_fpkm_allG_DE, "ccLnc_IPLIL6_fpkm_allG_DE.csv", row.names = F)


#
#### volcano + outliers ####

#highlight top FC/CPM
#select interesting genes in up
ccLnc_IPLIL6_fpkm_allG_up <- ccLnc_IPLIL6_fpkm_allG_up[order(ccLnc_IPLIL6_fpkm_allG_up$log2FoldChange, decreasing = T),]
ccLnc_IPLIL6_fpkm_allG_up$FCrank <- 1:length(ccLnc_IPLIL6_fpkm_allG_up$EnsName)
ccLnc_IPLIL6_fpkm_allG_up <- ccLnc_IPLIL6_fpkm_allG_up[order(ccLnc_IPLIL6_fpkm_allG_up$fpkm_max_treatment, decreasing = T),]
ccLnc_IPLIL6_fpkm_allG_up$CPMrank <- 1:length(ccLnc_IPLIL6_fpkm_allG_up$EnsName)

#average rank in both:
ccLnc_IPLIL6_fpkm_allG_up$DErank <- rowMeans(ccLnc_IPLIL6_fpkm_allG_up[,27:28])

#sum rank from both:
#ccLnc_IPLIL6_fpkm_allG_up$DErank <- ccLnc_IPLIL6_fpkm_allG_up[,27]*ccLnc_IPLIL6_fpkm_allG_up[,28]
ccLnc_IPLIL6_fpkm_allG_up <- ccLnc_IPLIL6_fpkm_allG_up[order(ccLnc_IPLIL6_fpkm_allG_up$DErank, decreasing = F),]
ccLnc_IPLIL6_fpkm_allG_up$DErank2 <- 1:length(ccLnc_IPLIL6_fpkm_allG_up$EnsName)
ccLnc_IPLIL6_fpkm_allG_up$upDown <- "UpRegWithGapmer(Upstream_IPLIL6_activation_represses)"

#now down
ccLnc_IPLIL6_fpkm_allG_down <- ccLnc_IPLIL6_fpkm_allG_down[order(ccLnc_IPLIL6_fpkm_allG_down$log2FoldChange, decreasing = F),]#note change this to F for repressed
ccLnc_IPLIL6_fpkm_allG_down$FCrank <- 1:length(ccLnc_IPLIL6_fpkm_allG_down$EnsName)
ccLnc_IPLIL6_fpkm_allG_down <- ccLnc_IPLIL6_fpkm_allG_down[order(ccLnc_IPLIL6_fpkm_allG_down$fpkm_max_treatment, decreasing = T),]
ccLnc_IPLIL6_fpkm_allG_down$CPMrank <- 1:length(ccLnc_IPLIL6_fpkm_allG_down$EnsName)

#average rank in both:
ccLnc_IPLIL6_fpkm_allG_down$DErank <- rowMeans(ccLnc_IPLIL6_fpkm_allG_down[,27:28])

#sum rank in both:
#ccLnc_IPLIL6_fpkm_allG_down$DErank <- ccLnc_IPLIL6_fpkm_allG_down[,24]*ccLnc_IPLIL6_fpkm_allG_down[,25]
ccLnc_IPLIL6_fpkm_allG_down <- ccLnc_IPLIL6_fpkm_allG_down[order(ccLnc_IPLIL6_fpkm_allG_down$DErank, decreasing = F),]
ccLnc_IPLIL6_fpkm_allG_down$DErank2 <- 1:length(ccLnc_IPLIL6_fpkm_allG_down$EnsName)
ccLnc_IPLIL6_fpkm_allG_down$upDown <- "DownRegWithGapmer(Upstream_IPLIL6_activation_induces)"

ccLnc_IPLIL6_fpkm_allG_DE <- rbind(ccLnc_IPLIL6_fpkm_allG_up, ccLnc_IPLIL6_fpkm_allG_down)

#highlight top FC/CPM
library(ggrepel)

ggplot(filter(ccLnc_IPLIL6_fpkm_allG_DE)) + aes(x = fpkm_max_treatment, y = log2FoldChange) +
  geom_point(alpha = 0.5) +
  geom_label_repel(data = filter(ccLnc_IPLIL6_fpkm_allG_DE, DErank2 <=15 #| fpkm_max_treatment >200
  ), 
  aes(x = fpkm_max_treatment, y = log2FoldChange, label = EnsName), size = 3, force = 20, nudge_x = 0.4, nudge_y = 0.2) +
  theme_bw() +
  theme(text = element_text(size = 22)) +
  scale_x_log10() +
  #coord_cartesian(xlim = c(4,12)) +
  #geom_label_repel(data = filter(kuppe_postMI_edgeR_CMECs_DE, !is.na(label)),
  #                 aes(x = logCPM, y = log2FoldChange, label = label), force = 150, nudge_x = 13, nudge_y = 4) +
  xlab("fpkm_max_treatment") +
  ylab("LFC GapIPLIL6 vs.\nGapNeg")


#highlight top FC/p
#select interesting genes in up
ccLnc_IPLIL6_fpkm_allG_up <- ccLnc_IPLIL6_fpkm_allG_up[order(ccLnc_IPLIL6_fpkm_allG_up$log2FoldChange, decreasing = T),]
ccLnc_IPLIL6_fpkm_allG_up$FCrank <- 1:length(ccLnc_IPLIL6_fpkm_allG_up$EnsName)
ccLnc_IPLIL6_fpkm_allG_up <- ccLnc_IPLIL6_fpkm_allG_up[order(ccLnc_IPLIL6_fpkm_allG_up$padj, decreasing = F),]
ccLnc_IPLIL6_fpkm_allG_up$Prank <- 1:length(ccLnc_IPLIL6_fpkm_allG_up$EnsName)

#average rank in both:
ccLnc_IPLIL6_fpkm_allG_up$DErank <- rowMeans(ccLnc_IPLIL6_fpkm_allG_up[,c(27,32)])

#sum rank from both:
#ccLnc_IPLIL6_fpkm_allG_up$DErank <- ccLnc_IPLIL6_fpkm_allG_up[,27]*ccLnc_IPLIL6_fpkm_allG_up[,28]
ccLnc_IPLIL6_fpkm_allG_up <- ccLnc_IPLIL6_fpkm_allG_up[order(ccLnc_IPLIL6_fpkm_allG_up$DErank, decreasing = F),]
ccLnc_IPLIL6_fpkm_allG_up$DErank2 <- 1:length(ccLnc_IPLIL6_fpkm_allG_up$EnsName)
ccLnc_IPLIL6_fpkm_allG_up$upDown <- "UpRegWithGapmer(Upstream_M.IPLIL6_activation_represses)"

#now down
ccLnc_IPLIL6_fpkm_allG_down <- ccLnc_IPLIL6_fpkm_allG_down[order(ccLnc_IPLIL6_fpkm_allG_down$log2FoldChange, decreasing = F),]#note change this to F for repressed
ccLnc_IPLIL6_fpkm_allG_down$FCrank <- 1:length(ccLnc_IPLIL6_fpkm_allG_down$EnsName)
ccLnc_IPLIL6_fpkm_allG_down <- ccLnc_IPLIL6_fpkm_allG_down[order(ccLnc_IPLIL6_fpkm_allG_down$padj, decreasing = F),]
ccLnc_IPLIL6_fpkm_allG_down$Prank <- 1:length(ccLnc_IPLIL6_fpkm_allG_down$EnsName)

#average rank in both:
ccLnc_IPLIL6_fpkm_allG_down$DErank <- rowMeans(ccLnc_IPLIL6_fpkm_allG_down[,c(27,32)])

#sum rank in both:
#ccLnc_IPLIL6_fpkm_allG_down$DErank <- ccLnc_IPLIL6_fpkm_allG_down[,24]*ccLnc_IPLIL6_fpkm_allG_down[,25]
ccLnc_IPLIL6_fpkm_allG_down <- ccLnc_IPLIL6_fpkm_allG_down[order(ccLnc_IPLIL6_fpkm_allG_down$DErank, decreasing = F),]
ccLnc_IPLIL6_fpkm_allG_down$DErank2 <- 1:length(ccLnc_IPLIL6_fpkm_allG_down$EnsName)
ccLnc_IPLIL6_fpkm_allG_down$upDown <- "DownRegWithGapmer(Upstream_M.IPLIL6_activation_induces)"

ccLnc_IPLIL6_fpkm_allG_DE <- rbind(ccLnc_IPLIL6_fpkm_allG_up, ccLnc_IPLIL6_fpkm_allG_down)

#highlight top FC/p
library(ggrepel)

ggplot(filter(ccLnc_IPLIL6_fpkm_allG_DE)) + aes(x = log2FoldChange, y = -log10(padj)) +
  geom_point(alpha = 0.5) +
  geom_label_repel(data = filter(ccLnc_IPLIL6_fpkm_allG_DE, DErank2 <=5 #| fpkm_max_treatment >200
  ), 
  aes(x = log2FoldChange, y = -log10(padj), label = EnsName), size = 3, force = 80, nudge_y = 6.2) +
  theme_bw() +
  theme(text = element_text(size = 22)) +
  #scale_x_log10() +
  #coord_cartesian(xlim = c(4,12)) +
  #geom_label_repel(data = filter(kuppe_postMI_edgeR_CMECs_DE, !is.na(label)),
  #                 aes(x = logCPM, y = log2FoldChange, label = label), force = 150, nudge_x = 13, nudge_y = 4) +
  xlab("fpkm_max_treatment") +
  ylab("LFC GapIPLIL6 vs.\nGapNeg")


#activation of IFN seems likely, are they on in 0-24hr anyway?
fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE_2026filt.csv", header = T)

#induced in 0-24hr
fpkm_allGDE_24hr_On <- filter(fpkm_allGDE, LogFC_0_24 > log2(1.5), preadj_0_24 < 0.05)
fpkm_allGDE_24hr_Off <- filter(fpkm_allGDE, LogFC_0_24 < -log2(1.5), preadj_0_24 < 0.05)
fpkm_allGDE_4hr_On <- filter(fpkm_allGDE, LogFC_0_4 > log2(1.5), preadj_0_4 < 0.05)
fpkm_allGDE_4hr_Off <- filter(fpkm_allGDE, LogFC_0_4 < -log2(1.5), preadj_0_4 < 0.05)

ccLnc_IPLIL6_fpkm_allG_DE$IPon <- "NoChange_Timecourse"
ccLnc_IPLIL6_fpkm_allG_DE$IPon[ccLnc_IPLIL6_fpkm_allG_DE$EnsID %in% fpkm_allGDE$EnsID] <- "OtherChange_Timecourse"
ccLnc_IPLIL6_fpkm_allG_DE$IPon[ccLnc_IPLIL6_fpkm_allG_DE$EnsID %in% fpkm_allGDE_24hr_On$EnsID] <- "Induction_0v24h"
ccLnc_IPLIL6_fpkm_allG_DE$IPon[ccLnc_IPLIL6_fpkm_allG_DE$EnsID %in% fpkm_allGDE_24hr_Off$EnsID] <- "Repression_0v24h"

table(ccLnc_IPLIL6_fpkm_allG_DE$IPon)
table(ccLnc_IPLIL6_fpkm_allG_DE$IPon, ccLnc_IPLIL6_fpkm_allG_DE$upDown)

ggplot(ccLnc_IPLIL6_fpkm_allG_DE) + aes(x = fpkm_max_treatment, y = log2FoldChange, color = IPon) +
  geom_point(alpha = 0.5, size= 3) +
  theme_bw() +
  theme(text = element_text(size = 22)) +
  scale_color_manual(values = c(NoChange_Timecourse = "grey70", Induction_0v24h = "firebrick4", Repression_0v24h = "steelblue4",
                                Induction_0v4h = "firebrick1", Repression_0v4h = "steelblue1")) +
  scale_x_log10() +
  xlab("fpkm_max_treatment") +
  ylab("LFC GapIPLIL6 vs.\nGapNeg")

table(filter(fpkm_allGDE_24hr_On)$GeneClassUpdate)
table(filter(fpkm_allGDE_24hr_On, EnsID %in% ccLnc_IPLIL6_fpkm_allG_down$EnsID)$GeneClassUpdate)

ccLnc_IPLIL6_fpkm_allG_DE$IPon <- "NoChange_Timecourse"
ccLnc_IPLIL6_fpkm_allG_DE$IPon[ccLnc_IPLIL6_fpkm_allG_DE$EnsID %in% fpkm_allGDE_4hr_On$EnsID] <- "Induction_0v4h"
ccLnc_IPLIL6_fpkm_allG_DE$IPon[ccLnc_IPLIL6_fpkm_allG_DE$EnsID %in% fpkm_allGDE_4hr_Off$EnsID] <- "Repression_0v4h"

table(ccLnc_IPLIL6_fpkm_allG_DE$IPon)
table(ccLnc_IPLIL6_fpkm_allG_DE$IPon, ccLnc_IPLIL6_fpkm_allG_DE$upDown)

ggplot(ccLnc_IPLIL6_fpkm_allG_DE) + aes(x = fpkm_max_treatment, y = log2FoldChange, color = IPon) +
  geom_point(alpha = 0.5, size= 3) +
  theme_bw() +
  theme(text = element_text(size = 22)) +
  scale_color_manual(values = c(NoChange_Timecourse = "grey70", Induction_0v24h = "firebrick4", Repression_0v24h = "steelblue4",
                                Induction_0v4h = "firebrick1", Repression_0v4h = "steelblue1")) +
  scale_x_log10() +
  xlab("fpkm_max_treatment") +
  ylab("LFC GapIPLIL6 vs.\nGapNeg")


#takeaways so far
#


#
#### GO/KEGG/REACTOME terms - GapIPLIL6 induced ####

#for enrichr:
#write.csv(ccLnc_IPLIL6_fpkm_allG_DE, "ccLnc_IPLIL6_fpkm_allG_DE.csv", row.names = F)

library(clusterProfiler)
library(org.Hs.eg.db)
ccLnc_IPLIL6_fpkm_allG
ccLnc_IPLIL6_fpkm_allG$EnsID_merge <- gsub("\\.[0-9]*", "", ccLnc_IPLIL6_fpkm_allG$EnsID)
ccLnc_IPLIL6_fpkm_allG_up$EnsID_merge <- gsub("\\.[0-9]*", "", ccLnc_IPLIL6_fpkm_allG_up$EnsID)

fpkm_DEG_GapIPLIL6up_GO <- enrichGO(gene          = ccLnc_IPLIL6_fpkm_allG_up$EnsID_merge,
                                    universe      = ccLnc_IPLIL6_fpkm_allG$EnsID_merge,
                                    keyType       = "ENSEMBL",
                                    OrgDb         = org.Hs.eg.db,
                                    ont           = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.05,
                                    readable      = TRUE)
fpkm_DEG_GapIPLIL6up_GO_df <- as.data.frame(fpkm_DEG_GapIPLIL6up_GO)

dotplot(simplify(fpkm_DEG_GapIPLIL6up_GO), showCategory = 20)

#can add an FC in for DEGs:
trial <- ccLnc_IPLIL6_fpkm_allG_DE$log2FoldChange
names(trial) <- ccLnc_IPLIL6_fpkm_allG_DE$EnsName
convertEnsEnt2 <- trial

fpkm_DEG_GapIPLIL6up_GO_S <- simplify(fpkm_DEG_GapIPLIL6up_GO)
edox <- setReadable(fpkm_DEG_GapIPLIL6up_GO_S, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, categorySizeBy= ~-log10(edox@result$pvalue[1:10]),  showCategory = 10, foldChange=convertEnsEnt2, fc_threshold = 1)

fpkm_DEG_GapIPLIL6up_GO_S_df <- as.data.frame(fpkm_DEG_GapIPLIL6up_GO_S)


#strong set related to immune response bias, response to biotic stimulus interferon signalling

#KEGG
convertEnsEnt <- bitr(unique(ccLnc_IPLIL6_fpkm_allG$EnsName), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

fpkm_DEG_GapIPLIL6up_K <- enrichKEGG(gene = filter(convertEnsEnt, SYMBOL %in% ccLnc_IPLIL6_fpkm_allG_up$EnsName)$ENTREZID,
                          universe = filter(convertEnsEnt, SYMBOL %in% ccLnc_IPLIL6_fpkm_allG$EnsName)$ENTREZID,
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)
fpkm_DEG_GapIPLIL6up_K_df <- data.frame(fpkm_DEG_GapIPLIL6up_K)
unique(fpkm_DEG_GapIPLIL6up_K_df$Description)
filter(convertEnsEnt, 
       ENTREZID %in% unlist(strsplit(filter(fpkm_DEG_GapIPLIL6up_K_df, Description == unique(fpkm_DEG_GapIPLIL6up_K_df$Description)[4])$geneID, 
                                     split = "/")))
#ferroptosis and fatty acid metabolism put to the fore - no clear steer on biosynthesis or degradation or efflux
#PPAR repressive of IL6
#includes HOMX1, ACSL5 strong genes

dotplot(fpkm_DEG_GapIPLIL6up_K, showCategory = 20)

library(pathview)

#can add an FC in for DEGs:
convertEnsEnt2 <- merge(convertEnsEnt, ccLnc_IPLIL6_fpkm_allG_DE[,c(3,25)], by.x = "SYMBOL", by.y = "EnsName")
trial <- convertEnsEnt2$log2FoldChange
names(trial) <- convertEnsEnt2$ENTREZID
convertEnsEnt2 <- trial

#ferroptosis
hsa04216 <- pathview(gene.data  = convertEnsEnt2,
                     pathway.id = "hsa04216",
                     species    = "hsa")

#PPAR
hsa03320 <- pathview(gene.data  = convertEnsEnt2,
                    pathway.id = "hsa03320",
                    species    = "hsa")


#REACTOME
library(ReactomePA)
fpkm_DEG_GapIPLIL6up_R <- enrichPathway(gene          = unique(filter(convertEnsEnt,SYMBOL %in% ccLnc_IPLIL6_fpkm_allG_up$EnsName)$ENTREZID),
                             universe      = unique(filter(convertEnsEnt,SYMBOL %in% ccLnc_IPLIL6_fpkm_allG$EnsName)$ENTREZID),
                             organism = "human",
                             pvalueCutoff = 0.05,
                             qvalueCutoff  = 0.05,
                             readable      = TRUE)
fpkm_DEG_GapIPLIL6up_R_df <- as.data.frame(fpkm_DEG_GapIPLIL6up_R)

dotplot(fpkm_DEG_GapIPLIL6up_R, showCategory = 20)

viewPathway("PPARA activates gene expression", 
            readable = TRUE, 
            foldChange = convertEnsEnt2, )

#a clear PPARA signal emerges
#this is a lipid metabolising pathway also associated with reduction in VSMC proliferation which is quite handy
#proposed as atheroprotective via cholesterol homeostasis, fatty acid degradation reducing inflammation
#quite good references to support the above
#PPARA is a repressor of IL6, but unsure of vice-versa, doesn't seem well looked at in VSMC

#is it - and associated genes - downregulated normally?? is this avoided by IPLIL6?

#combined figure/table Go-S/K/R:
fpkm_DEG_GapIPLIL6up_GO_S <- simplify(fpkm_DEG_GapIPLIL6up_GO, cutoff = 0.6)#101 at 0.7 (2x ferroptosis tho)
fpkm_DEG_GapIPLIL6up_GO_S_df <- as.data.frame(fpkm_DEG_GapIPLIL6up_GO_S)

#description, gene ratio, adj p val, genes
colnames(fpkm_DEG_GapIPLIL6up_GO_S_df)
colnames(fpkm_DEG_GapIPLIL6up_K_df)[]
colnames(fpkm_DEG_GapIPLIL6up_R_df)

trial <- fpkm_DEG_GapIPLIL6up_GO_S_df[,c(2,3,9,11)]
triali <- fpkm_DEG_GapIPLIL6up_K_df[,c(4,5,11,13)]
trialii <- fpkm_DEG_GapIPLIL6up_R_df[,c(2,3,9,11)]

trial$termType <- "GO"
triali$termType <- "KEGG"
trialii$termType <- "REACTOME"

combined_GOKR <- rbind(trial, triali, trialii)

#this can be saved for supplement

#plotting:
combined_GOKR$selectHits <- as.numeric(sapply(strsplit(combined_GOKR$GeneRatio, "\\/"), "[[", 1))
combined_GOKR$select <- as.numeric(sapply(strsplit(combined_GOKR$GeneRatio, "\\/"), "[[", 2))
combined_GOKR$GeneRatio <- combined_GOKR$selectHits/combined_GOKR$select*100

#top 10 GO, top 5 others
combined_GOKR_plot <- filter(combined_GOKR, (termType == "GO" & p.adjust < 0.002) | 
                               (termType == "KEGG" & p.adjust < 0.03) |
                               (termType == "REACTOME"))

combined_GOKR_plot <- combined_GOKR_plot[order(combined_GOKR_plot$Description),]

#CoUp_DE2_df_MF$DescriptionII <- stringr::str_wrap(CoUp_DE2_df$Description, width = 20)

combined_GOKR_plot$Description <- tolower(combined_GOKR_plot$Description)
combined_GOKR_plot$Description <- factor(combined_GOKR_plot$Description)
combined_GOKR_plot$Description <- factor(combined_GOKR_plot$Description,
                                      levels = levels(combined_GOKR_plot$Description)[order( 
                                        combined_GOKR_plot$termType, -combined_GOKR_plot$GeneRatio, decreasing = T)])

colnames(combined_GOKR_plot)[2] <- "% GapIPL-IL6\ninduced Genes"

ggplot(combined_GOKR_plot) + aes(x =`% GapIPL-IL6\ninduced Genes`, fill = -log10(p.adjust), y = Description) +
  geom_bar(stat = "identity", color = "grey40")+
  theme_bw() +
  scale_x_continuous(breaks = seq(0,30,5)) +
  scale_fill_continuous(limits = c(0, max(-log10(combined_GOKR_plot$p.adjust)))) +
  theme(text = element_text(size=24)) +
  ylab("")


#
#### GO/KEGG/REACTOME terms - GapIPLIL6 repressed ####

ccLnc_IPLIL6_fpkm_allG_down$EnsID_merge <- gsub("\\.[0-9]*", "", ccLnc_IPLIL6_fpkm_allG_down$EnsID)

fpkm_DEG_GapIPLIL6down_GO <- enrichGO(gene          = ccLnc_IPLIL6_fpkm_allG_down$EnsID_merge,
                                      universe      = ccLnc_IPLIL6_fpkm_allG$EnsID_merge,
                                      keyType       = "ENSEMBL",
                                      OrgDb         = org.Hs.eg.db,
                                      ont           = "BP",
                                      pAdjustMethod = "BH",
                                      pvalueCutoff  = 0.05,
                                      readable      = TRUE)
fpkm_DEG_GapIPLIL6down_GO_df <- as.data.frame(fpkm_DEG_GapIPLIL6down_GO)

dotplot(simplify(fpkm_DEG_GapIPLIL6down_GO), showCategory = 20)

fpkm_DEG_GapIPLIL6down_GO_S <- simplify(fpkm_DEG_GapIPLIL6down_GO)
edox <- setReadable(fpkm_DEG_GapIPLIL6down_GO_S, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, categorySizeBy= ~-log10(edox@result$pvalue[1:10]),  showCategory = 10, foldChange=convertEnsEnt2, fc_threshold = 0.8)

#cell cycle, some ECM too

#KEGG
convertEnsEnt <- bitr(unique(ccLnc_IPLIL6_fpkm_allG$EnsName), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

fpkm_DEG_GapIPLIL6down_K <- enrichKEGG(gene = filter(convertEnsEnt, SYMBOL %in% ccLnc_IPLIL6_fpkm_allG_down$EnsName)$ENTREZID,
                                     universe = filter(convertEnsEnt, SYMBOL %in% ccLnc_IPLIL6_fpkm_allG$EnsName)$ENTREZID,
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     qvalueCutoff  = 0.05)
fpkm_DEG_GapIPLIL6down_K_df <- data.frame(fpkm_DEG_GapIPLIL6down_K)
unique(fpkm_DEG_GapIPLIL6down_K_df$Description)
filter(convertEnsEnt, 
       ENTREZID %in% unlist(strsplit(filter(fpkm_DEG_GapIPLIL6down_K_df, Description == unique(fpkm_DEG_GapIPLIL6down_K_df$Description)[2])$geneID, 
                                     split = "/")))
#motor proteins - but in the cell motility category
#also VSMC contraction and cytoskeletion

dotplot(fpkm_DEG_GapIPLIL6down_K, showCategory = 20)

#can add an FC in for DEGs:
convertEnsEnt2 <- merge(convertEnsEnt, ccLnc_IPLIL6_fpkm_allG_DE[,c(3,25)], by.x = "SYMBOL", by.y = "EnsName")
trial <- convertEnsEnt2$log2FoldChange
names(trial) <- convertEnsEnt2$ENTREZID
convertEnsEnt2 <- trial

#VSMC contraction
hsa04270 <- pathview(gene.data  = convertEnsEnt2,
                     pathway.id = "hsa04270",
                     species    = "hsa")


#REACTOME
library(ReactomePA)
fpkm_DEG_GapIPLIL6down_R <- enrichPathway(gene          = unique(filter(convertEnsEnt,SYMBOL %in% ccLnc_IPLIL6_fpkm_allG_down$EnsName)$ENTREZID),
                                        universe      = unique(filter(convertEnsEnt,SYMBOL %in% ccLnc_IPLIL6_fpkm_allG$EnsName)$ENTREZID),
                                        organism = "human",
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff  = 0.05,
                                        readable      = TRUE)
fpkm_DEG_GapIPLIL6down_R_df <- as.data.frame(fpkm_DEG_GapIPLIL6down_R)

dotplot(fpkm_DEG_GapIPLIL6down_R, showCategory = 20)


#combined figure/table Go-S/K/R:
fpkm_DEG_GapIPLIL6down_GO_S <- simplify(fpkm_DEG_GapIPLIL6down_GO, cutoff = 0.65)#28 at 0.7 (2x cytokinesis tho)
fpkm_DEG_GapIPLIL6down_GO_S_df <- as.data.frame(fpkm_DEG_GapIPLIL6down_GO_S)

#description, gene ratio, adj p val, genes
colnames(fpkm_DEG_GapIPLIL6down_GO_S_df)
colnames(fpkm_DEG_GapIPLIL6down_K_df)[]
colnames(fpkm_DEG_GapIPLIL6down_R_df)

trial <- fpkm_DEG_GapIPLIL6down_GO_S_df[,c(2,3,9,11)]
triali <- fpkm_DEG_GapIPLIL6down_K_df[,c(4,5,11,13)]
trialii <- fpkm_DEG_GapIPLIL6down_R_df[,c(2,3,9,11)]

trial$termType <- "GO"
triali$termType <- "KEGG"
trialii$termType <- "REACTOME"

combined_GOKR <- rbind(trial, triali, trialii)

#this can be saved for supplement

#plotting:
combined_GOKR$selectHits <- as.numeric(sapply(strsplit(combined_GOKR$GeneRatio, "\\/"), "[[", 1))
combined_GOKR$select <- as.numeric(sapply(strsplit(combined_GOKR$GeneRatio, "\\/"), "[[", 2))
combined_GOKR$GeneRatio <- combined_GOKR$selectHits/combined_GOKR$select*100

#top 10 GO, top 5 others
combined_GOKR_plot <- filter(combined_GOKR, (termType == "GO" & p.adjust < 0.021) | 
                               (termType == "KEGG" & p.adjust < 0.0083 & selectHits >=9) |
                               (termType == "REACTOME" & p.adjust < 0.03))

#CoUp_DE2_df_MF$DescriptionII <- stringr::str_wrap(CoUp_DE2_df$Description, width = 20)

combined_GOKR_plot$Description[grepl("Extracellular", combined_GOKR_plot$Description)] <- "Extracellular matrix organization (REACTOME)"
combined_GOKR_plot$Description[grepl("ECM proteoglycans", combined_GOKR_plot$Description)] <- "Extracellular matrix proteoglycans"

combined_GOKR_plot$Description[grepl("Regulation of MITF-M-dependent genes involved in extracellular matrix, focal adhesion and epithelial-to-mesenchymal transition", 
                                     combined_GOKR_plot$Description)] <- "MITF-M-dependent genes in extracellular matrix"
combined_GOKR_plot <- combined_GOKR_plot[order(combined_GOKR_plot$Description),]
combined_GOKR_plot$Description <- tolower(combined_GOKR_plot$Description)

combined_GOKR_plot$Description <- factor(combined_GOKR_plot$Description)
combined_GOKR_plot$Description <- factor(combined_GOKR_plot$Description,
                                         levels = levels(combined_GOKR_plot$Description)[order( 
                                           combined_GOKR_plot$termType, -combined_GOKR_plot$GeneRatio, decreasing = T)])

colnames(combined_GOKR_plot)[2] <- "% GapIPL-IL6\nrepressed Genes"

ggplot(combined_GOKR_plot) + aes(x =`% GapIPL-IL6\nrepressed Genes`, fill = -log10(p.adjust), y = Description) +
  geom_bar(stat = "identity", color = "grey40")+
  theme_bw() +
  scale_x_continuous(breaks = seq(0,30,5)) +
  scale_fill_continuous(limits = c(0, max(-log10(combined_GOKR_plot$p.adjust)))) +
  theme(text = element_text(size=24)) +
  ylab("")

#
#### final volcano with some highlighted genes ####

#stick to those found in GO/K/R + outliers on the plot

#CDH13 - massive outlier, anti-apoptotic, adiponectin receptor, loss indicates SMC damage - needs a ref tho as not in GO/K/R
Notable <- c("CDH13",

#genes for signs of lipid accumulation
"ABCA1", "AGT", "ANGPTL4", "G0S2", "AKR1C1",

#genes indicating increased oxidative damage/need for redox protection
"NQO1", "HMOX1", "FTL", "SLC7A11",

#key genes indicating loss of contractile program
"ACTA2", "TAGLN", "LMOD1",

#key genes indicating loss of ECM program
"FBLN1", "COL1A1", "COL1A3", "FLNB",

#key genes indicating loss of proliferation
"MKI67", "TOP2A", "CDK1")

ccLnc_IPLIL6_fpkm_allG_DE$Notable <- NA
ccLnc_IPLIL6_fpkm_allG_DE$Notable[ccLnc_IPLIL6_fpkm_allG_DE$EnsName == "CDH13"] <- "Anti-apoptotic"
ccLnc_IPLIL6_fpkm_allG_DE$Notable[ccLnc_IPLIL6_fpkm_allG_DE$EnsName %in% c("ABCA1", "AGT", "ANGPTL4", "G0S2", "AKR1C1")] <- "Response to lipid"
ccLnc_IPLIL6_fpkm_allG_DE$Notable[ccLnc_IPLIL6_fpkm_allG_DE$EnsName %in% c("NQO1", "HMOX1", "FTL", "SLC7A11")] <- "Oxidative Damage/Ferroptosis"
ccLnc_IPLIL6_fpkm_allG_DE$Notable[ccLnc_IPLIL6_fpkm_allG_DE$EnsName %in% c("ACTA2", "TAGLN", "LMOD1")] <- "Contractility"
ccLnc_IPLIL6_fpkm_allG_DE$Notable[ccLnc_IPLIL6_fpkm_allG_DE$EnsName %in% c("MKI67", "TOP2A", "CDK1")] <- "Proliferation"
ccLnc_IPLIL6_fpkm_allG_DE$Notable[ccLnc_IPLIL6_fpkm_allG_DE$EnsName %in% c("FBLN1", "COL1A1", "COL1A3", "FLNB")] <- "ECM organisation"

ggplot(filter(ccLnc_IPLIL6_fpkm_allG_DE)) + aes(x = log2FoldChange, y = -log10(padj), color = Notable) +
  geom_point(data = ccLnc_IPLIL6_fpkm_allG_DE[is.na(ccLnc_IPLIL6_fpkm_allG_DE$Notable),], aes(colour = Notable), alpha = 0.2) +
  geom_point(data = ccLnc_IPLIL6_fpkm_allG_DE[!is.na(ccLnc_IPLIL6_fpkm_allG_DE$Notable),], aes(colour = Notable), alpha = 1,size=2) +
  geom_label_repel(data = filter(ccLnc_IPLIL6_fpkm_allG_DE, !is.na(Notable)), 
                   aes(x = log2FoldChange, y = -log10(padj), label = EnsName), size = 3.5, force = 100, nudge_y = 10.2) +
  theme_bw() +
  theme(text = element_text(size = 22)) +
  #scale_x_log10() +
  #coord_cartesian(xlim = c(4,12)) +
  #geom_label_repel(data = filter(kuppe_postMI_edgeR_CMECs_DE, !is.na(label)),
  #                 aes(x = logCPM, y = log2FoldChange, label = label), force = 150, nudge_x = 13, nudge_y = 4) +
  xlab("LFC GapIPL-IL6 vs.\nGapNeg") +
  ylab("-log10(p)")


geom_point(data = TPM_DEG_hypV_DE[TPM_DEG_hypV_DE$vegf == "Other", ], aes(colour = vegf)) +
  geom_point(data = TPM_DEG_hypV_DE[TPM_DEG_hypV_DE$vegf != "Other", ], aes(colour = vegf), alpha = 0.5) +