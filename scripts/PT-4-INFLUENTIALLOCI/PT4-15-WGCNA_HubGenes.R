#### 17 WGCNA ####

#PartII-WGCNA
library(tximport)
library(DESeq2)
library(WGCNA)
library(dplyr)
library(ggplot2)

#chiefly for identification of hub genes, then see if cis-acting lncRNAs target them

#this is an easy tutorial for a first look: https://bioinformaticsworkbook.org/tutorials/wgcna.html
#since supplemented with this: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
#e.g. to suport batch correction beforehand

#also one example done in "Baker-lab" BioInfoGroupResources area

#approach can take all DE genes but some say invalidates idea of scale-free network...

#Check out Quertermous snATACseq paper for idea on better hub gene ID

#### import data to build network ####

samplenames <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/rsem_4timepoints(NovelTx)_NonMergedEnsID", pattern = "*genes.results", full.names = TRUE)[1:16]

#sample info copied/pasted from Amira (matches PCA characteristics later)
#she orderered by patient, then condition
actualnames <- c("1) Pt151 - 0h 2) Pt151 - 4h 3) Pt151 - 8h 4) Pt151 - 24h 
                 5) Pt157 - 0h 6) Pt157 - 4h 7) Pt157 - 8h 8) Pt157 - 24h 
                 9) Pt134 - 0h 10) Pt134 - 4h 11) Pt134 - 8h 12) Pt134 - 24h 
                 13) Pt2279 - 0h14) Pt2279 - 4h 15) Pt2279 - 8h 16) Pt2279 - 24h")
trial <- unlist(strsplit(actualnames, "[0-9]) "))[2:17]
actualnames <- strsplit(trial, "h")
actualnames <- sapply(actualnames, "[[", 1)

sampleInfo_Timecourse <- data.frame("SampleNo" = c(1,10:16,2:9),
                                    "FileName" = samplenames,
                                    "ActualNames" = actualnames[c(1,10:16,2:9)])

#sample names re-ordered into the same
sampleInfo_Timecourse <- sampleInfo_Timecourse[c(1,12,16,5,
                                                 9,13,2,6,
                                                 10,14,3,7,
                                                 11,15,4,8),]

#tximport, better use of isoform level info during DESeq2 norm see: https://support.bioconductor.org/p/94003/
samplenames <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/rsem_4timepoints(NovelTx)_NonMergedEnsID/", 
                          pattern = "*isoforms.results", full.names = TRUE)[1:16]

samplenames_order <- samplenames[c(1,12,16,5,
                                   9,13,2,6,
                                   10,14,3,7,
                                   11,15,4,8)]
iso_filenames <- samplenames_order

tx2gene <- read.csv(iso_filenames[1], sep = "\t", stringsAsFactors = F) #85369 transcripts

#build meta-data table
sampleInfo_Timecourse$condition <- c(rep("0hr",4),rep("4hr",4),rep("8hr",4),rep("24hr",4))
sampleInfo_Timecourse$patient <- rep(c("pt151","pt157","pt134","pt2279"),4)

sampleInfo_Timecourse$condition <- factor(sampleInfo_Timecourse$condition, levels(as.factor(sampleInfo_Timecourse$condition))[c(1,3,4,2)])
row.names(sampleInfo_Timecourse) <- sampleInfo_Timecourse$ActualNames

#import in format for DESeq2 input
txi <- tximport(iso_filenames, type = "rsem", tx2gene = tx2gene)

#LRT model
dds <- DESeqDataSetFromTximport(txi, sampleInfo_Timecourse, design=~patient+condition)
dds <-dds[rowSums(counts(dds))>10, ]
ddsLRT <- DESeq(dds, test = "LRT", reduced = ~patient)
resLRT <- results(ddsLRT)

#pca check (need to remove batch)
rld <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup = c("patient", "condition"), returnData = TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))

colnames(pcaData)[c(8,9)] <- c("Timepoint", "Patient")

pcaData[,c(1:3,8:9)]

ggplot(pcaData, aes(x = PC1, y = PC2, color = Timepoint, shape = Patient)) +
  geom_point(size =3) +
  xlab(paste0("\nPC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance\n")) +
  coord_fixed() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size =16.5),
        axis.title.x = element_text(size = 16.5),
        legend.text = element_text(size =15),
        legend.title = element_text(size = 16.5))

assay(rld) <- limma::removeBatchEffect(assay(rld), rld$patient)

pcaData <- plotPCA(rld, intgroup = c("patient", "condition"), returnData = TRUE, ntop = 1500)
percentVar <- round(100 * attr(pcaData, "percentVar"))

colnames(pcaData)[c(8,9)] <- c("Patient", "Timepoint")

ggplot(pcaData, aes(x = PC1, y = PC2, color = Timepoint, shape = Patient)) +
  geom_point(size =3) +
  xlab(paste0("\nPC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance\n")) +
  coord_fixed() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size =16.5),
        axis.title.x = element_text(size = 16.5),
        legend.text = element_text(size =15),
        legend.title = element_text(size = 16.5))

expr_normalized <- assay(rld)

#### prep for WGCNA ####
#WGCNA works on genes as cols:
#input_mat = t(expr_normalized)
#write.csv(input_mat, "input_mat_WGCNA_2026.csv", row.names = T)
input_mat <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/input_mat_WGCNA_2026.csv",  
                      header = T, row.names = 1)
fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026filt.csv")
length(unique(fpkm_allG$EnsID))#12740

#only FPKM 1 genes
input_mat <- input_mat[,colnames(input_mat) %in% fpkm_allG$MSTRG_ID]
dim(input_mat)

allowWGCNAThreads()          # allow multi-threading (optional)
#> Allowing multi-threading with up to 8 threads.

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.698, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power = 9
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue), shouldn't be issue but for safety

#
#### WGCNA first attempt ####

#below are settings obtained from the above prep script and various tutorials
#decided to run all 11k genes as one block
netwk <- blockwiseModules(input_mat,               # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed", #neg and pos rather than just any type
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,#supposedly the easiest way to increase clustering default is 2
                          pamRespectsDendro = F,
                          detectCutHeight = 0.75, #from a tutorial (prior attempt pushed this a bit)
                          minModuleSize = 30,
                          maxBlockSize = 13000,
                          #maxBlockSize = 4000, #original run was probs a bit too small on blocks, laptop can probs handle the 11k genes
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time, allows module reclassing later)
                          saveTOMs = T,
                          loadTOM = F,
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

#cor <- temp_cor
saveRDS(netwk, "oneBlock_netwk_2026.rds")
#netwk <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_netwk.rds")


#### WGCNA subsequent attempts (can load the TOM) ####

netwk_0775 <- blockwiseModules(input_mat,               # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed", #neg and pos rather than just any type
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,#supposedly the easiest way to increase clustering default is 2, but with this data/parameters detectCutHeight works better
                          pamRespectsDendro = F,
                          detectCutHeight = 0.775, #used 0.8 originally, 0.775 is first point where sustained immune appears
                          minModuleSize = 30,
                          maxBlockSize = 13000,
                          #maxBlockSize = 4000, #original run was probs a bit too small on blocks, laptop can probs handle the 11k genes
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time, allows module reclassing later)
                          saveTOMs = F,
                          loadTOM = T,
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

#cor <- temp_cor
saveRDS(netwk_0775, "oneBlock_netwk_0775_2026.rds")


#### WGCNA subsequent attempts (can load the TOM) ####

netwk_0785 <- blockwiseModules(input_mat,               # <= input here
                               
                               # == Adjacency Function ==
                               power = picked_power,                # <= power here
                               networkType = "signed", #neg and pos rather than just any type
                               
                               # == Tree and Block Options ==
                               deepSplit = 2,#supposedly the easiest way to increase clustering default is 2, but with this data/parameters detectCutHeight works better
                               pamRespectsDendro = F,
                               detectCutHeight = 0.785, #used 0.8 originally, 0.775 was first point where sustained immune appears
                               minModuleSize = 30,
                               maxBlockSize = 13000,
                               #maxBlockSize = 4000, #original run was probs a bit too small on blocks, laptop can probs handle the 11k genes
                               
                               # == Module Adjustments ==
                               reassignThreshold = 0,
                               mergeCutHeight = 0.25,
                               
                               # == TOM == Archive the run results in TOM file (saves time, allows module reclassing later)
                               saveTOMs = F,
                               loadTOM = T,
                               
                               # == Output Options
                               numericLabels = T,
                               verbose = 3)

#cor <- temp_cor
saveRDS(netwk_0785, "oneBlock_netwk_0785_2026.rds")


#### WGCNA subsequent attempts (can load the TOM) ####

netwk_07825 <- blockwiseModules(input_mat,               # <= input here
                               
                               # == Adjacency Function ==
                               power = picked_power,                # <= power here
                               networkType = "signed", #neg and pos rather than just any type
                               
                               # == Tree and Block Options ==
                               deepSplit = 2,#supposedly the easiest way to increase clustering default is 2, but with this data/parameters detectCutHeight works better
                               pamRespectsDendro = F,
                               detectCutHeight = 0.7825, #used 0.8 originally, 0.775 was first point where sustained immune appears
                               minModuleSize = 30,
                               maxBlockSize = 13000,
                               #maxBlockSize = 4000, #original run was probs a bit too small on blocks, laptop can probs handle the 11k genes
                               
                               # == Module Adjustments ==
                               reassignThreshold = 0,
                               mergeCutHeight = 0.25,
                               
                               # == TOM == Archive the run results in TOM file (saves time, allows module reclassing later)
                               saveTOMs = F,
                               loadTOM = T,
                               
                               # == Output Options
                               numericLabels = T,
                               verbose = 3)

#cor <- temp_cor
saveRDS(netwk_0785, "oneBlock_netwk_0785_2026.rds")

#0.780 does not provide a rapid,sustained induced cluster
#nor does 0.7825

#### (start here to analyse network) import netwks and key objects for analysis #####

#various saved netwks
#netwk <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_netwk_2026.rds")
#netwk <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_netwk_0775_2026.rds")

#this one identifies sustained early induced cluster
netwk <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_netwk_0785_2026.rds")

#all counts
input_mat <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/input_mat_WGCNA_2026.csv",  
                      header = T, row.names = 1)

fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allG_2026filt.csv")
length(unique(fpkm_allG$EnsID))


#
#### initial module correlation table + plots ####

#assign any netwk to "netwk"
mergedColors = labels2colors(netwk$colors) #e.g 5 = green cluster
# Plot the dendrogram and the module colors underneath, per block
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)

#clusters per gene:
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

#size of clusters
table(module_df$colors)

#build a heatmap, modules expression across samples as basic first look
#confirm only looking at assayed , shouldn't actually filter anything
input_mat <- input_mat[,colnames(input_mat) %in% fpkm_allG$MSTRG_ID]
dim(input_mat)

#calculate module eigengene - the general direction of expression dynamics through the samples per module
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
library(tidyr)
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME$treatment <- factor(mME$treatment)
mME$treatment <- factor(mME$treatment, levels = levels(mME$treatment)[c(1,5,9,13,
                                                                        3,7,11,15,
                                                                        4,8,12,16,
                                                                        2,6,10,14)])

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue3",
    high = "red3",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90),
        text = element_text(size=15)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="Correlation")


#### module bio relevance via GO/KEGG ####

#expected gene ontology?
library(clusterProfiler)
library(org.Hs.eg.db)

module_df <- unique(merge(module_df, fpkm_allG[,c(1,2)], by.x = "gene_id", by.y = "MSTRG_ID"))
module_df$Ens_ID_merge <- gsub("\\.[0-9]*", "", module_df$EnsID)
fpkm_allG$Ens_ID_merge <- gsub("\\.[0-9]*", "", fpkm_allG$EnsID)

modules <- split(module_df, module_df$colors)

#check an individual module:
individual <- enrichGO(gene          = unique(modules$greenyellow$Ens_ID_merge),
                       universe      = unique(fpkm_allG$Ens_ID_merge),
                       keyType       = "ENSEMBL",
                       OrgDb         = org.Hs.eg.db,
                       ont           = "all",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
View(data.frame(individual))

#check all modules
trial <- lapply(modules, function(x){
  enrichGO(gene          = unique(x$Ens_ID_merge),
           universe      = unique(fpkm_allG$Ens_ID_merge),
           keyType       = "ENSEMBL",
           OrgDb         = org.Hs.eg.db,
           ont           = "all",
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05,
           readable      = TRUE)
})

#saveRDS(trial, "oneBlock_ds2_0785_WGCNA_GOall_2026.rds")

#0.75 modules
#trial <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_WGCNA_GO.rds")
#0.8 modules
#trial <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_14mod_WGCNA_GO.rds")

#GO per module:
GO_per_Module <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_ds2_0785_WGCNA_GOall_2026.rds")

#check no. genes per modules
sapply(modules, function(x){dim(x)[1]})
#check no. GO terms per modules (how strong is the theme found for each set of genes)
sapply(GO_per_Module, function(x){dim(x)[1]})
#check how many genes have been assigned to a GO term (also indicates how strong theme is)
#remove empty if needed
sapply(GO_per_Module[1], function(x){length(unique(unlist(strsplit(x$geneID, "/"))))})
#as percentage
sapply(GO_per_Module[1], function(x){length(unique(unlist(strsplit(x$geneID, "/"))))})/sapply(modules[1], function(x){dim(x)[1]})


#plot GO themes found per module
names(GO_per_Module)
modNumber <- "black"

edox2 <- enrichplot::pairwise_termsim(GO_per_Module[[modNumber]])
enrichplot::treeplot(edox2, showCategory = 30, #fontSize =2, #extend = 0.1, #hilight =F #nWords = 3,
                     cluster.params = list(n = 6), 
                     #label_format_tiplab = function(x) stringr::str_wrap(x, width=40),
                     label_format = function(x) stringr::str_wrap(x, width=25))
View(data.frame(GO_per_Module[[modNumber]])) 

dotplot(GO_per_Module[[modNumber]], showCategory = 20)

#use to build a supplementary table
names(GO_per_Module)
#write.csv(as.data.frame(GO_per_Module[["black"]]), "black.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[["blue"]]), "blue.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[["brown"]]), "brown.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[["green"]]), "green.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[["greenyellow"]]), "greenyellow.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[["magenta"]]), "magenta.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[["pink"]]), "pink.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[["purple"]]), "purple.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[["red"]]), "red.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[["turquoise"]]), "turquoise.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[["yellow"]]), "yellow.csv", row.names = F)


#KEGG
convertEnsEnt <- bitr(unique(fpkm_allG$Ens_ID_merge), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

trialK <- lapply(modules, function(x){
  enrichKEGG(gene = unique(filter(convertEnsEnt, ENSEMBL %in% 
                                    x$Ens_ID_merge)$ENTREZID),
             universe = unique(convertEnsEnt$ENTREZID),
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.05)
})

#saveRDS(trialK, "oneBlock_ds2_0785_WGCNA_KEGGall_2026.rds")

KEGG_per_Module <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_ds2_0785_WGCNA_KEGGall_2026.rds")
names(KEGG_per_Module)

modNumber <- "yellow"
View(as.data.frame(KEGG_per_Module[[modNumber]]))

names(KEGG_per_Module)
#write.csv(as.data.frame(KEGG_per_Module[["black"]]), "blackKEGG.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[["blue"]]), "blueKEGG.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[["brown"]]), "brownKEGG.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[["green"]]), "greenKEGG.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[["greenyellow"]]), "greenyellowKEGG.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[["magenta"]]), "magentaKEGG.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[["pink"]]), "pinkKEGG.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[["purple"]]), "purpleKEGG.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[["red"]]), "redKEGG.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[["turquoise"]]), "turquoiseKEGG.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[["yellow"]]), "yellowKEGG.csv", row.names = F)


#
#### central genes within modules ####

#correlate each gene's expression with the eigengene per cluster

dim(MEs0)[2]
geneModuleMembership = as.data.frame(cor(input_mat, MEs0[,-dim(MEs0)[2]], use = "p"))
#each gene will correlate most with it's corresponding clusters eigengene?
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 16)) #nSamples here

#access modules top MM value genes
modInfo <- list()
modNames = substring(names(MEs0[,-dim(MEs0)[2]]), 3)

#geneModuleMembershipAll <- geneModuleMembership
geneModuleMembership <- geneModuleMembership[,paste("ME", modNames, sep="")]
#MMPvalueAll <- MMPvalue
MMPvalue <- MMPvalue[,paste("ME", modNames, sep="")]

Clusters <- modNames
#Clusters <- modNames[c(1:3)]#selected clusters of interest

for(i in 1:length(Clusters)){
  module <-  Clusters[[i]]
  moduleGenes = mergedColors==module
  moduleGenesMM <- data.frame("MM" = geneModuleMembership[,colnames(geneModuleMembership) == paste("ME", module, sep="")], 
                              "MM.p" = MMPvalue[,colnames(MMPvalue) == paste("ME", module, sep="")],
                              row.names =  rownames(geneModuleMembership))[moduleGenes,]
  modInfo[[i]] <- moduleGenesMM
}
names(modInfo) <- Clusters

moduleGenesMM <- bind_rows(modInfo, .id = "Module")
moduleGenesMM$MSTRG_ID <- rownames(moduleGenesMM)

#add better info, e.g. EnsID:
colnames(fpkm_allG)
moduleGenesMM <- unique(merge(fpkm_allG[,c(1,2,3,4,58)], moduleGenesMM, by = "MSTRG_ID"))

fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE_2026filt.csv", header = T)
length(unique(fpkm_allGDE$EnsID))

fpkm_allGDE_Upwithin_4 <- filter(fpkm_allGDE, RegulationStart == "Induced <4hrs")
fpkm_allGDE_Downwithin_4 <- filter(fpkm_allGDE, RegulationStart == "Repressed <4hrs")

fpkm_allGDE_Upwithin_8 <- filter(fpkm_allGDE, RegulationStart == "Induced 4-8hrs")
fpkm_allGDE_Downwithin_8 <- filter(fpkm_allGDE, RegulationStart == "Repressed 4-8hrs")

fpkm_allGDE_Upwithin_24 <- filter(fpkm_allGDE, RegulationStart == "Induced 8-24hrs")
fpkm_allGDE_Downwithin_24 <- filter(fpkm_allGDE, RegulationStart == "Repressed 8-24hrs")

moduleGenesMM$DEG_Cluster[!moduleGenesMM$EnsID %in% fpkm_allGDE$EnsID] <- "NonDEG"
moduleGenesMM$DEG_Cluster[moduleGenesMM$EnsID %in% fpkm_allGDE_Upwithin_4$EnsID] <- "Induced <4hrs"
moduleGenesMM$DEG_Cluster[moduleGenesMM$EnsID %in% fpkm_allGDE_Downwithin_4$EnsID] <- "Repressed <4hrs"
moduleGenesMM$DEG_Cluster[moduleGenesMM$EnsID %in% fpkm_allGDE_Upwithin_8$EnsID] <- "Induced 4-8hrs"
moduleGenesMM$DEG_Cluster[moduleGenesMM$EnsID %in% fpkm_allGDE_Downwithin_8$EnsID] <- "Repressed 4-8hrs"
moduleGenesMM$DEG_Cluster[moduleGenesMM$EnsID %in% fpkm_allGDE_Upwithin_24$EnsID] <- "Induced 8-24hrs"
moduleGenesMM$DEG_Cluster[moduleGenesMM$EnsID %in% fpkm_allGDE_Downwithin_24$EnsID] <- "Repressed 8-24hrs"

moduleGenesMM$DEG_Simple <- "NonDEG"
moduleGenesMM$DEG_Simple[moduleGenesMM$EnsID %in% fpkm_allGDE$EnsID] <- "DEG"


#
#### IP driven modules ####

#plot all modules + whether DEGs are enriched (plot the unassigned, grey seperately)
ggplot(filter(moduleGenesMM, !Module == "grey"
)) + aes(x = DEG_Simple, y = MM, color = DEG_Simple) +
  geom_violin() +
  #geom_jitter() +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  theme_minimal() +
  theme(strip.background = element_rect(),
        axis.text.x = element_blank(),
        text = element_text(size=24)) +
  #coord_cartesian(ylim = c(0.8, 1.1)) +
  #geom_hline(yintercept = min(filter(moduleGenesMM, MM.p <0.05)$MM), linetype = "dashed", color = "grey70") +
  facet_wrap(~Module, ncol = 4) +
  xlab("") +
  ylab("Module Membership")

ggplot(filter(moduleGenesMM, Module == "grey"
)) + aes(x = DEG_Simple, y = MM, color = DEG_Simple) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  theme_minimal() +
  theme(strip.background = element_rect(),
        axis.text.x = element_blank(),
        text = element_text(size=24)) +
  coord_cartesian(ylim = c(-1, 1.2)) +
  #geom_hline(yintercept = min(filter(moduleGenesMM, MM.p <0.05)$MM), linetype = "dashed", color = "grey70") +
  facet_wrap(~Module, ncol = 6) +
  xlab("") +
  ylab("Module Membership")

#stats, is MM higher for DEGs in a given cluster?
#pool sizes
#table(moduleGenesMM$Module)
#table(moduleGenesMM$DEG_Simple, moduleGenesMM$Module)

#runs stats
#is actually a t test based on high numbers but just a label here:
wilcox_p_MM_DEG <- list()
wilcox_means_MM_DEG <- list()

#normality/variance per cluster here (non-normal but equal var, generally way over 30 n per pool so t.test)
#ggplot() + 
#  stat_qq(aes(sample = filter(moduleGenesMM, Module == Clusters[2], DEG_Simple == "DEG")$MM), colour = "green") + 
#  stat_qq(aes(sample = filter(moduleGenesMM, Module == Clusters[2], !DEG_Simple == "DEG")$MM), colour = "red") +
#  geom_abline(aes(slope = 1, intercept = 0), linetype = 2)

wilcox_means_MM_DEG_val1 <- list()
wilcox_means_MM_DEG_val2 <- list()

for (i in 1:length(Clusters)){
  wilcox_p_MM_DEG[[i]] <- t.test(filter(moduleGenesMM, Module == Clusters[i], DEG_Simple == "DEG")$MM, 
                                 filter(moduleGenesMM, Module == Clusters[i], !DEG_Simple == "DEG")$MM, var.equal = T)$p.value
  wilcox_means_MM_DEG[[i]] <- c(median(filter(moduleGenesMM, Module == Clusters[i], DEG_Simple == "DEG")$MM)/
                                  median(filter(moduleGenesMM, Module == Clusters[i], !DEG_Simple == "DEG")$MM))  
  wilcox_means_MM_DEG_val1[[i]] <- median(filter(moduleGenesMM, Module == Clusters[i], DEG_Simple == "DEG")$MM)
  wilcox_means_MM_DEG_val2[[i]] <- median(filter(moduleGenesMM, Module == Clusters[i], !DEG_Simple == "DEG")$MM)
}

names(wilcox_p_MM_DEG) <- Clusters
#allows defining of IP-driven clusters, separate to patient effects/noise
p.adjust(unlist(wilcox_p_MM_DEG), method = "BH")
#diff in mean (possible should be mean diff)
names(wilcox_means_MM_DEG) <- Clusters
unlist(wilcox_means_MM_DEG)

names(wilcox_means_MM_DEG_val1) <- Clusters
unlist(wilcox_means_MM_DEG_val1)


#Themed/IP-driven clusters:
names(wilcox_p_MM_DEG)[p.adjust(unlist(wilcox_p_MM_DEG), method = "BH") < 0.05 & 
                         wilcox_means_MM_DEG >1.1 & !names(wilcox_p_MM_DEG) == "grey"]

#excluding grey despite higher IP as these are unassigned genes...
IPdriven <- names(wilcox_p_MM_DEG)[p.adjust(unlist(wilcox_p_MM_DEG), method = "BH") < 0.05 & 
                                     wilcox_means_MM_DEG >1.12 & !names(wilcox_p_MM_DEG) == "grey"]

dim(filter(moduleGenesMM, Module %in% IPdriven))
2985/12740 #23% if cutting turquoise/pink/black/brown @1.1 - dds 0.785

#other modules not of direct interest for SMC activation
moduleGenesMM$moduleOfInterest <- "No"
moduleGenesMM$moduleOfInterest[moduleGenesMM$Module %in% IPdriven] <- "Yes"

#replot modules according to DEG pattern:
ggplot(filter(moduleGenesMM, Module %in% IPdriven)) + aes(x = DEG_Simple, y = MM, color = DEG_Simple) +
  geom_violin() +
  #geom_jitter() +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  theme_minimal() +
  theme(strip.background = element_rect(),
        axis.text.x = element_blank(),
        text = element_text(size=24)) +
  #coord_cartesian(ylim = c(0.8, 1.1)) +
  #geom_hline(yintercept = min(filter(moduleGenesMM, MM.p <0.05)$MM), linetype = "dashed", color = "grey70") +
  facet_wrap(~Module, ncol = 6) +
  xlab("") +
  ylab("Module Membership")

#low responders
ggplot(filter(moduleGenesMM, Module %in% c("turquoise", "brown", "black", "pink"),
)) + aes(x = DEG_Simple, y = MM, color = DEG_Simple) +
  geom_violin() +
  #geom_jitter() +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  theme_minimal() +
  theme(strip.background = element_rect(),
        axis.text.x = element_blank(),
        text = element_text(size=24)) +
  #coord_cartesian(ylim = c(0.8, 1.1)) +
  #geom_hline(yintercept = min(filter(moduleGenesMM, MM.p <0.05)$MM), linetype = "dashed", color = "grey70") +
  facet_wrap(~Module, ncol = 6) +
  xlab("") +
  ylab("Module Membership")

#anti DEGs
ggplot(filter(moduleGenesMM, Module %in% c("blue", "greenyellow"))) + aes(x = DEG_Simple, y = MM, color = DEG_Simple) +
  geom_violin() +
  #geom_jitter() +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  theme_minimal() +
  theme(strip.background = element_rect(),
        axis.text.x = element_blank(),
        text = element_text(size=24)) +
  #coord_cartesian(ylim = c(0.8, 1.1)) +
  #geom_hline(yintercept = min(filter(moduleGenesMM, MM.p <0.05)$MM), linetype = "dashed", color = "grey70") +
  facet_wrap(~Module, ncol = 6) +
  xlab("") +
  ylab("Module Membership")


#annotate module label based on above:
moduleGenesMM$ModuleSummary[grepl("^yellow", moduleGenesMM$Module)] <-       "Immune/muscle proliferation"
moduleGenesMM$ModuleSummary[grepl("green$", moduleGenesMM$Module)] <-       "Unclear"
moduleGenesMM$ModuleSummary[grepl("red", moduleGenesMM$Module)] <-          "Muscle homeostasis"
moduleGenesMM$ModuleSummary[grepl("magenta", moduleGenesMM$Module)] <-          "Cell division"
moduleGenesMM$ModuleSummary[grepl("purple", moduleGenesMM$Module)] <-          "Immune(sustained)"


#
#### association of core cell cycle genes with hub-ness ####

moduleGenesMM_cellCycle <- filter(moduleGenesMM, Module %in% c("magenta"))
moduleGenesMM_cellCycle$Core_SG2M <- "Other"
moduleGenesMM_cellCycle$Core_SG2M[moduleGenesMM_cellCycle$EnsID %in% filter(fpkm_allG, grepl("CC", GeneClassUpdate))$EnsID] <- "Core SG2M"

#plot all modules + whether DEGs are enriched (plot the unassigned, grey seperately)
ggplot(filter(moduleGenesMM_cellCycle)) + aes(x = Core_SG2M, y = MM, color = Core_SG2M) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  theme_minimal() +
  theme(strip.background = element_rect(),
        axis.text.x = element_blank(),
        text = element_text(size=24)) +
  #coord_cartesian(ylim = c(0.4, 1.1)) +
  #geom_hline(yintercept = min(filter(moduleGenesMM, MM.p <0.05)$MM), linetype = "dashed", color = "grey70") +
  facet_wrap(~Module, ncol = 6) +
  xlab("") +
  ylab("Module Membership")

t.test(filter(moduleGenesMM_cellCycle, Module == "magenta", Core_SG2M == "Core SG2M")$MM, 
       filter(moduleGenesMM_cellCycle, Module == "magenta", !Core_SG2M == "Core SG2M")$MM, var.equal = T)


#
#### key canonical immune signals notable as hub genes of purple ####

moduleGenesMM_immune <- filter(moduleGenesMM, Module %in% c("purple"))
moduleGenesMM_immune$immuneMarkers <- NA
moduleGenesMM_immune$immuneMarkers[moduleGenesMM_immune$EnsName %in% c("CXCL8", "IL6", "IL1B", "CCL2", "MIR3142HG")] <- "ImmuneMarkers"

#plot all modules + whether DEGs are enriched (plot the unassigned, grey seperately)
ggplot(filter(moduleGenesMM_immune)) + aes(x = "", y = MM, color = immuneMarkers) +
  geom_jitter(width = 0.07) +
  theme_minimal() +
  ggrepel::geom_label_repel(data = filter(moduleGenesMM_immune, !is.na(immuneMarkers)),
                            aes(x = "", y = MM, color = immuneMarkers, label = EnsName), 
                            nudge_x = 1.5, nudge_y = -0.05, size =6, force = 10) +
  theme(strip.background = element_rect(),
        axis.text.x = element_blank(),
        legend.position = "none",
        text = element_text(size=24)) +
  facet_wrap(~Module, ncol = 1) +
  xlab("") +
  ylab("Module Membership")


#
#### finish table, define hub genes via a threshold, label modules ####

#annotate scclnc association
SCClncRNAs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SCClncRNAs.csv")

AllLNC_AllPCG_2d3d <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_2d3d_2026.csv")
CoRegPairs_04_48_24_extended <- filter(AllLNC_AllPCG_2d3d, AbsDistLnc_PCG <250,
                                       #AllLNC_AllPCG_1,
                                       (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                     fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                       fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                        fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                         (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                       fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
                                       )


moduleGenesMM$CCLnc_association <- "None"
moduleGenesMM$CCLnc_association[moduleGenesMM$EnsID %in% CoRegPairs_04_48_24_extended$EnsID] <- "CClncRNA"
moduleGenesMM$CCLnc_association[moduleGenesMM$EnsID %in% CoRegPairs_04_48_24_extended$EnsID.y] <- "CClncRNA target"
moduleGenesMM$CCLnc_association[moduleGenesMM$EnsID %in% SCClncRNAs$EnsID] <- "SCClncRNA"
moduleGenesMM$CCLnc_association[moduleGenesMM$EnsID %in% SCClncRNAs$EnsID.y] <- "SCClncRNA target"

Clusters <- unique(moduleGenesMM$Module)

#assign hub genes, give % cclncs and % PCG
hubs <- list()
for (i in 1:length(Clusters)){
  module <- filter(moduleGenesMM, Module == Clusters[i])
  module$Hubness <- "Other"
  module$Hubness[module$EnsID %in% module[order(module$MM, decreasing = T),][1:(length(module$EnsID)*0.3),]$EnsID] <- "Top30"
  module$Hubness[module$EnsID %in% module[order(module$MM, decreasing = T),][1:(length(module$EnsID)/5),]$EnsID] <- "Top20"
  module$Hubness[module$EnsID %in% module[order(module$MM, decreasing = T),][1:(length(module$EnsID)/10),]$EnsID] <- "Top10"
  hubs[[i]] <- module
}

moduleGenesMM <- bind_rows(hubs)

#write.csv(moduleGenesMM, "moduleGenesMM_Feb2026.csv", row.names = F)


#### plot hub genes as %, display alongside heatmap and module labels ####

#key info import
#moduleGenesMM <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/moduleGenesMM_Nov25.csv")
#moduleGenesMM$ModuleSummary <- gsub("\\)", "", moduleGenesMM$ModuleSummary)
#IPdriven <- c("red", "blue", "brown", "green", "black") #default 2026

#which modules have highest % of CClncs or CClnc targets in their hubs?
moduleGenesMM_ <- filter(moduleGenesMM, Module %in% IPdriven)

#using top20% here:
for (i in 1:length(IPdriven)){
  hubsCheck <- filter(moduleGenesMM_, Module %in% IPdriven[i], 
                      Hubness %in% c("Top10", "Top20", "Top30"))
  print(
    paste(round(dim(filter(hubsCheck, grepl("target$", CCLnc_association)))[1]/dim(hubsCheck)[1]*100, 2), 
          IPdriven[i], 
          unique(filter(moduleGenesMM_, Module %in% IPdriven[i])$Module),
          "Targets")
  )
}

for (i in 1:length(IPdriven)){
  hubsCheck <- filter(moduleGenesMM_, Module %in% IPdriven[i], 
                      Hubness %in% c("Top10", "Top20", "Top30"))
  print(
    paste(round(dim(filter(hubsCheck, grepl("CClncRNA$", CCLnc_association)))[1]/dim(hubsCheck)[1]*100, 2), IPdriven[i], "CClncRNAs")
    )
}

#as percentage of DEGs in hub genes:
HubGenesMM_top <- filter(moduleGenesMM_, Hubness %in% c("Top10", "Top20", "Top30"), DEG_Simple == "DEG")

ggplot(filter(HubGenesMM_top, !CCLnc_association %in% c("None"))) + 
  aes(y = Module, fill = CCLnc_association) +
  geom_bar(position = position_dodge(width = 0.8))

table(HubGenesMM_top$Module)
table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$Module)

trial <- data.frame("CClncRNA" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$Module)[1,]/table(HubGenesMM_top$Module)*100,
                    "SCClncRNA" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$Module)[4,]/table(HubGenesMM_top$Module)*100)[,-3]
trial <- reshape2::melt(trial)

triali <- data.frame("CClncRNAn" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$Module)[1,],
                     "SCClncRNAn" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$Module)[4,])
triali <- reshape2::melt(triali)

trial$valueN <- triali$value

trial$variable <- gsub("\\.", " ", trial$variable)
trial$variable <- gsub(" Freq", "", trial$variable)

modules_hubs_cclncs <- trial
modules_hubs_cclncs$CClncRNA.Var1 <- factor(modules_hubs_cclncs$CClncRNA.Var1)
#ordering here, options are:http://127.0.0.1:38213/graphics/plot_zoom_png?width=607&height=314
#a) by theme (makes most sense probs)
modules_hubs_cclncs$CClncRNA.Var1 <- factor(modules_hubs_cclncs$CClncRNA.Var1, levels(modules_hubs_cclncs$CClncRNA.Var1)[c(6,8,7,1,3,2,5,4)])
#b) by % in Hubs
#trial$CClncRNA.Var1 <- factor(trial$CClncRNA.Var1, levels(trial$CClncRNA.Var1)[c(6,3,2,8,7,1,4,5)])

triali <- aggregate(valueN ~ CClncRNA.Var1, data = modules_hubs_cclncs, sum)

#modules with lots of CClncs in their hub genes
ggplot(modules_hubs_cclncs#[6:10,]
       ) + aes(y = CClncRNA.Var1, fill = variable, x = value, label = valueN) +
  geom_bar(stat = "identity", position = position_stack()) +
  geom_label(data = triali, size = 8, color = "black", inherit.aes = F, label.size = 0.8,
             label.padding = unit(0.17, "lines"),
             aes(y = CClncRNA.Var1, x = 5, label = valueN)) +
  xlab("% of Hub DEGs") +
  ylab("") +
  scale_fill_manual(values = c("CClncRNA" = "olivedrab4", "SCClncRNA" = "olivedrab2")) +
  scale_x_continuous(limits = c(0,7), breaks = c(0,5)) +
  theme_minimal() +
  theme(text = element_text(size=24))

#targets:
trial <- data.frame("CClncRNA target" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$Module)[2,]/table(HubGenesMM_top$Module)*100,
                    "SCClncRNA target" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$Module)[5,]/table(HubGenesMM_top$Module)*100)[,-3]
trial <- reshape2::melt(trial)

triali <- data.frame("CClncRNA targetn" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$Module)[2,],
                     "SCClncRNA targetn" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$Module)[5,])
triali <- reshape2::melt(triali)

trial$valueN <- triali$value

trial$variable <- gsub("\\.", " ", trial$variable)
trial$variable <- gsub(" Freq", "", trial$variable)

modules_hubs_cclncTargets <- trial
modules_hubs_cclncTargets$CClncRNA.target.Var1 <- factor(modules_hubs_cclncTargets$CClncRNA.target.Var1)
#ordering here, options are:
#a) by theme (makes most sense probs)
modules_hubs_cclncTargets$CClncRNA.target.Var1 <- factor(modules_hubs_cclncTargets$CClncRNA.target.Var1, 
                                                   levels(modules_hubs_cclncTargets$CClncRNA.target.Var1)[c(6,8,7,1,3,2,5,4)])
#b) by % in Hubs
#trial$CClncRNA.Var1 <- factor(trial$CClncRNA.Var1, levels(trial$CClncRNA.Var1)[c(6,3,2,8,7,1,4,5)])

triali <- aggregate(valueN ~ CClncRNA.target.Var1, data = modules_hubs_cclncTargets, sum)

#modules with lots of CClncs in their hub genes
ggplot(modules_hubs_cclncTargets) + aes(y = CClncRNA.target.Var1, fill = variable, x = value, label = valueN) +
  geom_bar(stat = "identity", position = position_stack()) +
  geom_label(data = triali, size = 8, color = "black", inherit.aes = F, label.size = 0.8,
             label.padding = unit(0.17, "lines"),
             aes(y = CClncRNA.target.Var1, x = 16, label = valueN)) +
  xlab("% of Hub DEGs") +
  ylab("") +
  scale_fill_manual(values = c("CClncRNA target" = "mediumorchid4", "SCClncRNA target" = "mediumorchid1")) +
  scale_x_continuous(limits = c(0,18), breaks = c(0,10,20)) +
  theme_minimal() +
  theme(text = element_text(size=24))


#### do a % per all of module for comparison/enrichment ####

#as percentage of DEGs in total module:
moduleGenesMM_DEGs <- filter(moduleGenesMM_, DEG_Simple == "DEG")

ggplot(filter(moduleGenesMM_DEGs, !CCLnc_association %in% c("None"))) + 
  aes(y = Module, fill = CCLnc_association) +
  geom_bar(position = position_dodge(width = 0.8))

table(moduleGenesMM_DEGs$Module)
table(moduleGenesMM_DEGs$CCLnc_association, moduleGenesMM_DEGs$Module)

trial <- data.frame("CClncRNA" = table(moduleGenesMM_DEGs$CCLnc_association, moduleGenesMM_DEGs$Module)[1,]/table(moduleGenesMM_DEGs$Module)*100,
                    "SCClncRNA" = table(moduleGenesMM_DEGs$CCLnc_association, moduleGenesMM_DEGs$Module)[4,]/table(moduleGenesMM_DEGs$Module)*100)[,-3]
trial <- reshape2::melt(trial)

triali <- data.frame("CClncRNAn" = table(moduleGenesMM_DEGs$CCLnc_association, moduleGenesMM_DEGs$Module)[1,],
                     "SCClncRNAn" = table(moduleGenesMM_DEGs$CCLnc_association, moduleGenesMM_DEGs$Module)[4,])
triali <- reshape2::melt(triali)

trial$valueN <- triali$value

trial$variable <- gsub("\\.", " ", trial$variable)
trial$variable <- gsub(" Freq", "", trial$variable)

modules_cclncs <- trial
modules_cclncs$CClncRNA.Var1 <- factor(modules_cclncs$CClncRNA.Var1)
#ordering here, options are:
#a) by theme (makes most sense probs)
modules_cclncs$CClncRNA.Var1 <- factor(modules_cclncs$CClncRNA.Var1, levels(modules_cclncs$CClncRNA.Var1)[c(6,8,7,1,3,2,5,4)])
#b) by % in Hubs
#trial$CClncRNA.Var1 <- factor(trial$CClncRNA.Var1, levels(trial$CClncRNA.Var1)[c(6,3,2,8,7,1,4,5)])

triali <- aggregate(valueN ~ CClncRNA.Var1, data = modules_cclncs, sum)

#modules with lots of CClncs in their hub genes
ggplot(modules_cclncs) + aes(y = CClncRNA.Var1, fill = variable, x = value, label = valueN) +
  geom_bar(stat = "identity", position = position_stack()) +
  geom_label(data = triali, size = 8, color = "black", inherit.aes = F, label.size = 0.8,
             label.padding = unit(0.17, "lines"),
             aes(y = CClncRNA.Var1, x = 6, label = valueN)) +
  xlab("% of all module DEGs") +
  ylab("") +
  scale_fill_manual(values = c("CClncRNA" = "olivedrab4", "SCClncRNA" = "olivedrab2")) +
  scale_x_continuous(limits = c(0,7), breaks = c(0,5)) +
  theme_minimal() +
  theme(text = element_text(size=24))

#targets:
trial <- data.frame("CClncRNA target" = table(moduleGenesMM_DEGs$CCLnc_association, moduleGenesMM_DEGs$Module)[2,]/table(moduleGenesMM_DEGs$Module)*100,
                    "SCClncRNA target" = table(moduleGenesMM_DEGs$CCLnc_association, moduleGenesMM_DEGs$Module)[5,]/table(moduleGenesMM_DEGs$Module)*100)[,-3]
trial <- reshape2::melt(trial)

triali <- data.frame("CClncRNA targetn" = table(moduleGenesMM_DEGs$CCLnc_association, moduleGenesMM_DEGs$Module)[2,],
                     "SCClncRNA targetn" = table(moduleGenesMM_DEGs$CCLnc_association, moduleGenesMM_DEGs$Module)[5,])
triali <- reshape2::melt(triali)

trial$valueN <- triali$value

trial$variable <- gsub("\\.", " ", trial$variable)
trial$variable <- gsub(" Freq", "", trial$variable)

modules_cclncTargets <- trial
modules_cclncTargets$CClncRNA.target.Var1 <- factor(modules_cclncTargets$CClncRNA.target.Var1)
#ordering here, options are:
#a) by theme (makes most sense probs)
modules_cclncTargets$CClncRNA.target.Var1 <- factor(modules_cclncTargets$CClncRNA.target.Var1, 
                                                         levels(modules_cclncTargets$CClncRNA.target.Var1)[c(6,8,7,1,3,2,5,4)])
#b) by % in Hubs
#trial$CClncRNA.Var1 <- factor(trial$CClncRNA.Var1, levels(trial$CClncRNA.Var1)[c(6,3,2,8,7,1,4,5)])

triali <- aggregate(valueN ~ CClncRNA.target.Var1, data = modules_cclncTargets, sum)

#modules with lots of CClncs in their hub genes
ggplot(modules_cclncTargets) + aes(y = CClncRNA.target.Var1, fill = variable, x = value, label = valueN) +
  geom_bar(stat = "identity", position = position_stack()) +
  geom_label(data = triali, size = 8, color = "black", inherit.aes = F, label.size = 0.8,
             label.padding = unit(0.17, "lines"),
             aes(y = CClncRNA.target.Var1, x = 16, label = valueN)) +
  xlab("% of all module DEGs") +
  ylab("") +
  scale_fill_manual(values = c("CClncRNA target" = "mediumorchid4", "SCClncRNA target" = "mediumorchid1")) +
  scale_x_continuous(limits = c(0,18), breaks = c(0,10,20)) +
  theme_minimal() +
  theme(text = element_text(size=24))

#looks like some chance for enrichment in the top connected - check vs other DE PCGs
table(filter(moduleGenesMM_DEGs, grepl("TF|CC|coding", GeneClassUpdate), EnsType == "protein_coding")$Module)
table(filter(HubGenesMM_top, grepl("TF|CC|coding", GeneClassUpdate), EnsType == "protein_coding")$Module)
table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$Module)
table(moduleGenesMM_DEGs$CCLnc_association, moduleGenesMM_DEGs$Module)

#purple - cclnc targets
a <- 10
b <- 85
c <- 17
d <- 242
a/b
c/d
fisher.test(data.frame("DEL" = c(a,b-a),
                       "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")
#0.03/2.8

#yellow - scclnc targets
a <- 14
b <- 225
c <- 25
d <- 613
a/b
c/d
fisher.test(data.frame("DEL" = c(a,b-a),
                       "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")
#0.035/2.3 (better)

#all - scclnc targets
table(filter(moduleGenesMM_DEGs, grepl("TF|CC|coding", GeneClassUpdate), EnsType == "protein_coding")$CCLnc_association)
table(filter(HubGenesMM_top, grepl("TF|CC|coding", GeneClassUpdate), EnsType == "protein_coding")$CCLnc_association)
dim(filter(moduleGenesMM_DEGs, grepl("TF|CC|coding", GeneClassUpdate), EnsType == "protein_coding"))
dim(filter(HubGenesMM_top, grepl("TF|CC|coding", GeneClassUpdate), EnsType == "protein_coding"))

a <- 40
b <- 818
c <- 78
d <- 2269
a/b
c/d
fisher.test(data.frame("DEL" = c(a,b-a),
                       "EL" = c(c-a,d-c-(b-a)), row.names = c("CC-target", "Other")), alternative = "greater")
#works, great result - 0.0036/1.9

#good result at 30%, with scclncRNA targets, neatly joining together Fig3 with WGCNA

#targets of scclncRNAs are higher up in connectivity than others

trial <- data.frame(a/b*100, c/d*100)
colnames(trial) <- c("Top 30% most central", "All module")

trial <- reshape2::melt(trial)
trial$variable <- as.factor(trial$variable)
trial$variable <- factor(trial$variable, levels(trial$variable)[])

ggplot(trial) + aes(y = value, x = variable, fill = variable) +
  geom_bar(stat = "identity", color = "grey60") +
  xlab("All IP-driven\nmodules") +
  ylab("") +
  theme_minimal() +
  scale_fill_manual(values = c(`Top 30% most central` = "mediumorchid1", `All module` = "grey60")) +
  scale_y_continuous(breaks = seq(0,6,2)) +
  theme(text = element_text(size =20), axis.text.x = element_blank())


#
#### integrate above 2 figures - hub vs all module % for SCClncRNA targets ####


modules_cclncTargets$geneType <- "All module"
modules_hubs_cclncTargets$geneType <- "Top 30% most central"

plot_scclncTargets_hubEnrich <- rbind(modules_cclncTargets[6:10,],
                                      modules_hubs_cclncTargets[6:10,])

#order by strongest percentages:
#reset
plot_scclncTargets_hubEnrich$CClncRNA.target.Var1 <- as.factor(as.character(plot_scclncTargets_hubEnrich$CClncRNA.target.Var1))
plot_scclncTargets_hubEnrich$CClncRNA.target.Var1 <- factor(plot_scclncTargets_hubEnrich$CClncRNA.target.Var1, 
                                                            levels = levels(plot_scclncTargets_hubEnrich$CClncRNA.target.Var1)[
                                                              order(plot_scclncTargets_hubEnrich[6:10,]$value, decreasing = F)
                                                            ])

plot_scclncTargets_hubEnrich$geneType <- as.factor(as.character(plot_scclncTargets_hubEnrich$geneType))
plot_scclncTargets_hubEnrich$geneType <- factor(plot_scclncTargets_hubEnrich$geneType, 
                                                            levels = levels(plot_scclncTargets_hubEnrich$geneType)[2:1
                                                            ])

ggplot(plot_scclncTargets_hubEnrich) + aes(x = value, y = CClncRNA.target.Var1, fill = geneType) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  xlab("% SCClncRNA targets") +
  ylab("") +
  scale_fill_manual(values = c("All module" = "grey60", "Top 30% most central" = "mediumorchid1")) +
  theme_minimal()+
  theme(text = element_text(size=24)) +
  scale_x_continuous(limits = c(0,8), breaks = c(0,2,4,6,8))

#all hub targeted:
SCClncRNAs_hubs <- filter(SCClncRNAs, EnsID.y %in% filter(moduleGenesMM_, !Hubness == "Other", CCLnc_association == "SCClncRNA target")$EnsID)

#0-4hr induced lncs:
EarlyUp_SCClncRNAs_hubs <- filter(SCClncRNAs_hubs, Lnc_Cluster == "Induced <4hrs")

length(unique(SCClncRNAs_hubs$EnsID))
length(unique(EarlyUp_SCClncRNAs_hubs$EnsID))
13/31

length(unique(SCClncRNAs_hubs$EnsID.y))
length(unique(EarlyUp_SCClncRNAs_hubs$EnsID.y))
19/40

#0-4hr repressed lncs:
EarlyDown_SCClncRNAs_hubs <- filter(SCClncRNAs_hubs, Lnc_Cluster == "Repressed <4hrs")

length(unique(SCClncRNAs_hubs$EnsID))
length(unique(EarlyDown_SCClncRNAs_hubs$EnsID))
11/31

length(unique(SCClncRNAs_hubs$EnsID.y))
length(unique(EarlyDown_SCClncRNAs_hubs$EnsID.y))
14/40


#### focus on lncRNA hubs ####

table(HubGenesMM_top$Module)
table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$Module)

HubGenesMM_top$CCLnc_association2 <- HubGenesMM_top$CCLnc_association
HubGenesMM_top$CCLnc_association2[HubGenesMM_top$EnsID %in% filter(fpkm_allG, grepl("fide", GeneClassUpdate))$EnsID & 
                                    !HubGenesMM_top$CCLnc_association2 == "SCClncRNA"] <- "LncRNA"
table(HubGenesMM_top$CCLnc_association2)
table(HubGenesMM_top$Module)
table(HubGenesMM_top$CCLnc_association2, HubGenesMM_top$Module)

trial <- data.frame("lncRNA" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$Module)[2,]/table(HubGenesMM_top$Module)*100,
                    "SCClncRNA" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$Module)[4,]/table(HubGenesMM_top$Module)*100)[,-3]
trial <- reshape2::melt(trial)

triali <- data.frame("lncRNAn" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$Module)[2,],
                     "SCClncRNAn" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$Module)[4,])
triali <- reshape2::melt(triali)

trial$valueN <- triali$value

trial$variable <- gsub("\\.", " ", trial$variable)
trial$variable <- gsub(" Freq", "", trial$variable)

modules_hubs_cclncs <- trial

#add in module details
modules_hubs_cclncs <- merge(modules_hubs_cclncs, unique(HubGenesMM_top[,c(6,12)]), by.x = "lncRNA.Var1", by.y = "Module")

modules_hubs_cclncs$ModuleSummary <- factor(modules_hubs_cclncs$ModuleSummary)
modules_hubs_cclncs$ModuleSummary <- factor(modules_hubs_cclncs$ModuleSummary, levels(modules_hubs_cclncs$ModuleSummary)[c(5,4,3,1,2)])

modules_hubs_cclncs$variable <- factor(modules_hubs_cclncs$variable)
modules_hubs_cclncs$variable <- factor(modules_hubs_cclncs$variable, levels(modules_hubs_cclncs$variable)[c(2,1)])

triali <- aggregate(valueN ~ ModuleSummary, data = modules_hubs_cclncs, sum)

#modules with lots of CClncs in their hub genes
ggplot(modules_hubs_cclncs#[6:10,]
) + aes(y = ModuleSummary, fill = variable, x = value, label = valueN) +
  geom_bar(stat = "identity", position = position_stack()) +
  geom_label(data = triali, size = 5, color = "black", inherit.aes = F, label.size = 0.8,
             label.padding = unit(0.17, "lines"),
             aes(y = ModuleSummary, x = 6.6, label = valueN)) +
  xlab("% of Hub DEGs") +
  ylab("") +
  scale_fill_manual(values = c("lncRNA" = "grey60", "SCClncRNA" = "olivedrab2")) +
  scale_x_continuous(limits = c(0,7.5), breaks = c(0,5)) +
  theme_minimal() +
  theme(text = element_text(size=24))


#
#### closer look at SCClncRNA hub genes ####

colnames(HubGenesMM_top)

#annotate scclnc association
SCClncRNAs <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/SCClncRNAs.csv")

SCClncRNA_hubs <- merge(SCClncRNAs, HubGenesMM_top[,c(2,6,12,7,8,14)], by.x = "EnsID.y", by.y = "EnsID", all.x = T)
colnames(SCClncRNA_hubs)[26:30] <- paste("Target", colnames(SCClncRNA_hubs)[26:30], sep = "_")

SCClncRNA_hubs <- merge(SCClncRNA_hubs, HubGenesMM_top[,c(2,6,12,7,8,14)], by.x = "EnsID", by.y = "EnsID", all.x = T)
colnames(SCClncRNA_hubs)[31:35] <- paste("SCClncRNA", colnames(SCClncRNA_hubs)[31:35], sep = "_")

#useful later
SCClncRNA_hubs$Sources2 <- "Other"
SCClncRNA_hubs$Sources2[!is.na(SCClncRNAs$Sources)] <- "GapmeR evidence of cis-acting"

length(unique(SCClncRNA_hubs$EnsID))
SCClncRNA_hubs_top30 <- filter(SCClncRNA_hubs, !is.na(Target_Hubness) | 
                                 !is.na(SCClncRNA_Hubness))
length(unique(SCClncRNA_hubs_top30$EnsID))
40/87 #46% ScclncRNAs are hub-targeting

#number of hub genes targeted by lnc with prior described cis activity
filter(SCClncRNA_hubs, !is.na(Target_Hubness), !is.na(Sources))[,1:7]
unique(filter(SCClncRNA_hubs, !is.na(Target_Hubness), !is.na(Sources))$EnsName.x)
unique(filter(SCClncRNAs, !is.na(Sources))$EnsID)
9/14 #64% of SCClncRNAs with ASO-confirmed target effect
#cyrano, cd27-as1, IPL-IL6, etc

#early induced scclncRNAs:
SCClncRNA_hubs_Up4 <- filter(SCClncRNA_hubs, Lnc_Cluster == "Induced <4hrs")
length(unique(SCClncRNA_hubs_Up4$EnsID))
colnames(SCClncRNA_hubs)
SCClncRNA_hubs_Up4_top30 <- filter(SCClncRNA_hubs_Up4, !is.na(Target_Hubness) | !is.na(SCClncRNA_Hubness))
length(unique(SCClncRNA_hubs_Up4_top30$EnsID))
18/37 #48% of 0-4hr induced SCClncRNAs are hub-targeting (no greater chance amongst these)

#2x loci with ASO validation
filter(SCClncRNA_hubs_Up4_top30, !is.na(Target_Hubness), !is.na(Sources))[,1:7]

#down reg
SCClncRNA_hubs_Down4 <- filter(SCClncRNA_hubs, Lnc_Cluster == "Repressed <4hrs")
length(unique(SCClncRNA_hubs_Down4$EnsID))
colnames(SCClncRNA_hubs)
SCClncRNA_hubs_Down4_top30 <- filter(SCClncRNA_hubs_Down4, !is.na(Target_Hubness) | !is.na(SCClncRNA_Hubness))
length(unique(SCClncRNA_hubs_Down4_top30$EnsID))
13/21

#bar plot, % hub-targeting SCClncRNAs with a) lit. support b) early induction, early repression
SCClncRNA_detail <- data.frame("Trait" = c("Expected\ncis-acting", "0-4hr induced", "0-4hr repressed"),
                               "hubSCClncRNAs_perc" = c(9/40*100, 18/40*100, 13/40*100))

ggplot(SCClncRNA_detail) + aes(x = Trait, y = hubSCClncRNAs_perc) +
  geom_bar(stat = "identity") +
  ylab("Hub-targeting SCClncRNAs") +
  xlab("") +
  theme_minimal()

#color by ASo evidence
SCClncRNAs$Sources2 <- "Other"
SCClncRNAs$Sources2[!is.na(SCClncRNAs$Sources)] <- "GapmeR evidence of cis-acting"

table(unique(SCClncRNA_hubs[,c(1,9,36)])$Lnc_Cluster, unique(SCClncRNA_hubs[,c(1,9,36)])$Sources2)
trial <- table(unique(SCClncRNA_hubs_top30[,c(1,9,36)])$Lnc_Cluster, unique(SCClncRNA_hubs_top30[,c(1,9,36)])$Sources2)
trial <- as.data.frame(trial)

ggplot(trial) + aes(x = Var1, y = Freq, fill = Var2) +
  geom_bar(stat = "identity", color = "grey60") + 
  scale_fill_manual(values = c(`GapmeR evidence of cis-acting` = "olivedrab1", Other = "olivedrab")) +
  theme_minimal() +
  theme(text = element_text(size = 24),
        axis.title.y = element_text(size = 20, angle = 0, hjust = 1, vjust = 0.5)) +
  xlab("") +
  ylab("Module\ncentre-targeting\nSCClncRNAs") +
  Seurat::RotatedAxis()
