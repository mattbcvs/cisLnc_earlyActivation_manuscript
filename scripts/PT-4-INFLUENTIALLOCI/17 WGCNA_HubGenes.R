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

#tximport trial, better use of isoform level info see: https://support.bioconductor.org/p/94003/
samplenames <- list.files("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/rsem_4timepoints(NovelTx)/", pattern = "*isoforms.results", full.names = TRUE)[1:16]

samplenames_order <- samplenames[c(1,12,16,5,
                                   9,13,2,6,
                                   10,14,3,7,
                                   11,15,4,8)]
iso_filenames <- samplenames_order

tx2gene <- read.csv(iso_filenames[1], sep = "\t", stringsAsFactors = F) #85369 transcripts

condition<-c(rep("0hr",4),rep("4hr",4),rep("8hr",4),rep("24hr",4))
patient<-rep(c("pt151","pt157","pt134","pt2279"),4)
type<-c(rep("paired-end",16))
data<-data.frame(condition=condition,type=type,patient=patient, stringsAsFactors = T)
data$condition <- factor(data$condition, levels(as.factor(data$condition))[c(1,3,4,2)])

actualnames <- c("1) Pt151 - 0h 2) Pt151 - 4h 3) Pt151 - 8h 4) Pt151 - 24h 5) Pt157 - 0h 6) Pt157 - 4h 7) Pt157 - 8h 8) Pt157 - 24h 9) Pt134 - 0h 10) Pt134 - 4h 11) Pt134 - 8h 12) Pt134 - 24h 13) Pt2279 - 0h14) Pt2279 - 4h 15) Pt2279 - 8h 16) Pt2279 - 24h")
trial <- unlist(strsplit(actualnames, "[0-9]) "))[2:17]
actualnames <- strsplit(trial, "h")
actualnames <- sapply(actualnames, "[[", 1)
actualnames <- actualnames[c(seq(1,16,4), seq(1,16,4)+1, seq(1,16,4)+2, seq(1,16,4)+3)]

row.names(data) <- actualnames

txi <- tximport(iso_filenames, type = "rsem", tx2gene = tx2gene)

#LRT model
dds <- DESeqDataSetFromTximport(txi, data, design=~patient+condition)
dds <- dds[rowSums(counts(dds))>10, ]
ddsLRT <- DESeq(dds, test = "LRT", reduced = ~patient)
resLRT <- results(ddsLRT)

#pca check (need to remove batch)
rld <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup = c("patient", "condition"), returnData = TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))

colnames(pcaData)[c(4,5)] <- c("Patient", "Timepoint")

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

colnames(pcaData)[c(4,5)] <- c("Patient", "Timepoint")

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
#write.csv(input_mat, "\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/input_mat_WGCNA.csv", 
#          row.names = T)
input_mat <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/input_mat_WGCNA.csv",  
                      header = T, row.names = 1)
#fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGClassUpdate.csv",
#                      header = T)
length(unique(fpkm_allG$EnsID))#11815 genes pre-IGV checks
length(unique(fpkm_allG_filt_manual$EnsID))#10761 genes after artefact removal + post-IGV checks

#only FPKM 1 genes
input_mat <- input_mat[,colnames(input_mat) %in% fpkm_allG_filt_manual$MSTRG_ID]

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
                          detectCutHeight = 0.75, #from a tutorial
                          minModuleSize = 30,
                          maxBlockSize = 12000,
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
saveRDS(netwk, "oneBlock_netwk_nov_dsII.rds")
#netwk <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_netwk.rds")


#### WGCNA subsequent attempts ####

netwk <- blockwiseModules(input_mat,               # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed", #neg and pos rather than just any type
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,#supposedly the easiest way to increase clustering default is 2, but with this data/parameters detectCutHeight works better
                          pamRespectsDendro = F,
                          detectCutHeight = 0.775, #used 0.8 originally, 0.775 is first point where sustained immune appears
                          minModuleSize = 30,
                          maxBlockSize = 12000,
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
saveRDS(netwk, "oneBlock_netwk_nov_ds2_0775.rds")


#### (start here to analyse network) import netwks and key objects for analysis #####

#various saved netwks
#netwk <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_netwk_nov.rds")  #ds3
#netwk <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_netwk_nov_dsII.rds")
#netwk <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_netwk_nov_ds2_08.rds")
netwk <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_netwk_nov_ds2_0775.rds")

#all counts
input_mat <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/input_mat_WGCNA.csv",  
                      header = T, row.names = 1)
fpkm_allG <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGClassUpdate.csv",
                      header = T)
FPKM_CQV_OVERLAP_fpkm <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/FPKM_CQV_OVERLAP_fpkm.csv")
table(FPKM_CQV_OVERLAP_fpkm$IGV)#413 pass, 168 fail

#should have removed artefacts from here already
fpkm_allG_filt <- filter(fpkm_allG, grepl("chr", chr), !grepl("artefacts", GeneClassUpdate), !is.na(GeneClassUpdate))
#remove manual or Thresh4 fails
fpkm_allG_filt_manual <- filter(fpkm_allG_filt, 
                                !EnsID %in% filter(FPKM_CQV_OVERLAP_fpkm, IGV == "fail")$EnsID, #remove manual fails
)
fpkm_allG <- fpkm_allG_filt_manual

#with removed artefacts
input_mat <- input_mat[,colnames(input_mat) %in% fpkm_allG_filt_manual$MSTRG_ID]

fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE_clustered.csv", header = T)
fpkm_allGDE_filt <- filter(fpkm_allGDE, grepl("chr", chr), !grepl("artefacts", GeneClassUpdate), !is.na(GeneClassUpdate))
#remove manual/Thresh4 fails
fpkm_allGDE <- filter(fpkm_allGDE_filt, EnsID %in% fpkm_allG_filt_manual$EnsID)


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
#check cluster numbers match previous
#netwk dsIII 0.75:
#black      blue     brown     green      grey      pink       red turquoise    yellow 
#374      2248      1161       414      2930       314       410      2479       431 

#netwk dsII, 0.75:
#black      blue     brown     green      grey       red turquoise    yellow 
#372      2134      1164       515      2942       413      2475       746 

#netwk dsII, 0.76:
#black      blue     brown     green      grey      pink       red turquoise    yellow 
#325      2054      1309       597      2589       141       419      2419       908

#netwk dsII, 0.775: (SELECTED)
#black        blue       brown       green greenyellow(immune theme)grey     magenta        pink      purple         red   turquoise      yellow 
#517        1856         926         643         185        2074         299         448         243         604        2120         846 

#netwk dsII 0.8:
#black        blue       brown       green greenyellow(no theme)grey     magenta        pink      purple         red   turquoise      yellow 
#638        1917         836         744         163        1312         445         525         302         743        2386         750

#netwk dsII, 0.8 + defaults:
#black         blue        brown         cyan        green  greenyellow         grey       grey60    lightcyan      magenta midnightblue         pink 
#419         1294          711           98          471          338         3000           38           43          394           46          403 
#purple          red       salmon          tan    turquoise       yellow 
#373          443          213          302         1619          556 


#netwk dsII, 0.8 + defaults, PRD off:
#black          blue         brown          cyan     darkgreen       darkred darkturquoise         green   greenyellow          grey        grey60 
#523           771           765           291            49            51            37           741           392          1312           177 
#lightcyan    lightgreen   lightyellow       magenta  midnightblue          pink        purple           red     royalblue        salmon           tan 
#230           159           116           446           257           459           429           552           114           302           363 
#turquoise        yellow 
##1483           742 

#build a heatmap, modules expression across samples as basic first look
#confirm only looking at assayed , shouldn't actually filter anything
input_mat <- input_mat[,colnames(input_mat) %in% fpkm_allG$MSTRG_ID]

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
    low = "#67A9CF",
    high = "#D6604D",
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

#saveRDS(trial, "oneBlock_ds2_0775_WGCNA_GOall.rds")

#0.75 modules
#trial <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_WGCNA_GO.rds")
#0.8 modules
#trial <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_14mod_WGCNA_GO.rds")

#GO per module:
GO_per_Module <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_ds2_0775_WGCNA_GO.rds")

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
modNumber <- 3

edox2 <- enrichplot::pairwise_termsim(GO_per_Module[[modNumber]])
enrichplot::treeplot(edox2, showCategory = 30, #fontSize =2, #extend = 0.1, #hilight =F #nWords = 3,
                     cluster.params = list(n = 6), 
                     #label_format_tiplab = function(x) stringr::str_wrap(x, width=40),
                     label_format = function(x) stringr::str_wrap(x, width=25))
View(data.frame(GO_per_Module[[10]])) 

dotplot(GO_per_Module[[modNumber]], showCategory = 20)

#for supplementary (move to post module selection)
#write.csv(as.data.frame(GO_per_Module[[4]]), "green.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[[5]]), "greenyellow.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[[7]]), "magenta.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[[8]]), "pink.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[[9]]), "purple.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[[10]]), "red.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[[11]]), "turquoise.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[[12]]), "yellow.csv", row.names = F)

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

#saveRDS(trialK, "oneBlock_ds2_0775_WGCNA_KEGGall.rds")
KEGG_per_Module <- readRDS("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/oneBlock_ds2_0775_WGCNA_KEGGall.rds")
names(KEGG_per_Module)
View(as.data.frame(KEGG_per_Module[[10]]))

#make some plots:
ownPlots <- list()
trialKi <- KEGG_per_Module

sapply(trialKi, dim)

for (i in c(1,3,4,5,8,9:12)){
  trialKi[[i]] <- data.frame(trialKi[[i]])
  #View(trialKi[[i]])
  trialKi[[i]]$selectHits <- as.numeric(sapply(strsplit(trialKi[[i]]$GeneRatio, "\\/"), "[[", 1))
  trialKi[[i]]$select <- as.numeric(sapply(strsplit(trialKi[[i]]$GeneRatio, "\\/"), "[[", 2))
  trialKi[[i]]$geneRatio <- trialKi[[i]]$selectHits/trialKi[[i]]$select*100
  
  trialKi[[i]] <- trialKi[[i]][order(trialKi[[i]]$Description),]
  
  trialKi[[i]]$DescriptionII <- stringr::str_wrap(trialKi[[i]]$Description, width = 40)
  
  trialKi[[i]]$Description <- factor(trialKi[[i]]$Description, labels = trialKi[[i]]$DescriptionII)
  trialKi[[i]]$Description <- factor(trialKi[[i]]$Description,
                                     levels = levels(trialKi[[i]]$Description)[order(trialKi[[i]]$geneRatio, decreasing = F)])
  
  trialKi[[i]] <- trialKi[[i]][order(-trialKi[[i]]$geneRatio, trialKi[[i]]$p.adjust),]
  
  trialKi[[i]] <- trialKi[[i]][,c(4,11,17)]
  
  colnames(trialKi[[i]])[3] <- "% of module\n with KEGG term"
}

sizeTerms <- 10
sapply(trialKi, dim)

KEGGPlots <- list()

for (i in c(1,3,4,5,8,9:12)){
  KEGGPlots[[i]] <- ggplot(filter(trialKi[[i]][1:sizeTerms,], !is.na(Description))) + 
    aes(x = `% of module\n with KEGG term`, y = Description, color = -log10(p.adjust)) +
    geom_point(stat = "identity", size = 5)+
    theme_bw() +
    scale_x_continuous(limits = c(0,max(trialKi[[i]][1:sizeTerms,3])+0.2), breaks = seq(0,30,5)) +
    theme(text = element_text(size=24)) +
    #xlab("") +
    ylab("")
}

names(trialKi)
KEGGPlots[[1]]
View(as.data.frame(trialKi[[10]]))
#with no filter on MM, some kinda weird terms, but some good too:
#"muscle cytoskeleton, lipid and athero"


#### central genes within modules ####

#correlate each gene's expression with the eigengene per cluster
#closest correlations are most connected - i.e. hubs

#take out last col (13 in this case given the 12 modules)
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
moduleGenesMM <- unique(merge(fpkm_allG[,c(1,2,5,6,60)], moduleGenesMM, by = "MSTRG_ID"))

fpkm_allGDE <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/fpkm_allGDE_clustered.csv", header = T)
fpkm_allGDE_filt <- filter(fpkm_allGDE, grepl("chr", chr), !grepl("artefacts", GeneClassUpdate), !is.na(GeneClassUpdate))
#remove manual/Thresh4 fails
fpkm_allGDE <- filter(fpkm_allGDE_filt, EnsID %in% fpkm_allG_filt_manual$EnsID)

moduleGenesMM <- merge(moduleGenesMM, fpkm_allGDE[,c(1,50)], by = "EnsID", all.x = T)
moduleGenesMM$DEG_Cluster[is.na(moduleGenesMM$DEG_Cluster)] <- "NonDEG"
moduleGenesMM$DEG_Simple <- "NonDEG"
moduleGenesMM$DEG_Simple[!moduleGenesMM$DEG_Cluster == "NonDEG"] <- "DEG"


#### IP driven modules ####

#plot all modules + whether DEGs are enriched (plot the unassigned, grey seperately)
ggplot(filter(moduleGenesMM, !Module == "grey"
)) + aes(x = DEG_Simple, y = MM, color = DEG_Simple) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.2) +
  theme_minimal() +
  theme(strip.background = element_rect(),
        axis.text.x = element_blank(),
        text = element_text(size=24)) +
  coord_cartesian(ylim = c(0.4, 1.1)) +
  #geom_hline(yintercept = min(filter(moduleGenesMM, MM.p <0.05)$MM), linetype = "dashed", color = "grey70") +
  facet_wrap(~Module, ncol = 6) +
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
table(moduleGenesMM$Module)
table(moduleGenesMM$DEG_Simple, moduleGenesMM$Module)

#runs stats
#is actually a t test based on high numbers but just a label here:
wilcox_p_MM_DEG <- list()
wilcox_means_MM_DEG <- list()

#normality/variance per cluster here (non-normal but equal var, generally way over 30 n per pool so t.test)
ggplot() + 
  stat_qq(aes(sample = filter(moduleGenesMM, Module == Clusters[2], DEG_Simple == "DEG")$MM), colour = "green") + 
  stat_qq(aes(sample = filter(moduleGenesMM, Module == Clusters[2], !DEG_Simple == "DEG")$MM), colour = "red") +
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2)

for (i in 1:length(Clusters)){
  wilcox_p_MM_DEG[[i]] <- t.test(filter(moduleGenesMM, Module == Clusters[i], DEG_Simple == "DEG")$MM, 
                                 filter(moduleGenesMM, Module == Clusters[i], !DEG_Simple == "DEG")$MM, var.equal = T)$p.value
  wilcox_means_MM_DEG[[i]] <- c(mean(filter(moduleGenesMM, Module == Clusters[i], DEG_Simple == "DEG")$MM)/
                                  mean(filter(moduleGenesMM, Module == Clusters[i], !DEG_Simple == "DEG")$MM))  
}

names(wilcox_p_MM_DEG) <- Clusters
#allows defining of IP-driven clusters, separate to patient effects/noise
p.adjust(unlist(wilcox_p_MM_DEG), method = "BH")
#diff in mean (possible should be mean diff)
names(wilcox_means_MM_DEG) <- Clusters

#Themed/IP-driven clusters:
names(wilcox_p_MM_DEG)[p.adjust(unlist(wilcox_p_MM_DEG), method = "BH") < 0.05 & 
                         wilcox_means_MM_DEG >1 & !names(wilcox_p_MM_DEG) == "grey"]

#excluding grey despite higher IP as these are unassigned genes...
IPdriven <- names(wilcox_p_MM_DEG)[p.adjust(unlist(wilcox_p_MM_DEG), method = "BH") < 0.05 & 
                                     wilcox_means_MM_DEG >1 & !names(wilcox_p_MM_DEG) == "grey"]

dim(filter(moduleGenesMM, Module %in% IPdriven))
5388/10761 #50% in IP-driven module

#other modules not of direct interest for SMC activation
moduleGenesMM$moduleOfInterest <- "No"
moduleGenesMM$moduleOfInterest[moduleGenesMM$Module %in% IPdriven] <- "Yes"

#save GO and KEGG for these and put in supplementary
which(names(GO_per_Module) %in% IPdriven)
names(GO_per_Module)[which(names(GO_per_Module) %in% IPdriven)]

#write.csv(as.data.frame(GO_per_Module[[4]]), "green.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[[5]]), "greenyellow.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[[7]]), "magenta.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[[8]]), "pink.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[[9]]), "purple.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[[10]]), "red.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[[11]]), "turquoise.csv", row.names = F)
#write.csv(as.data.frame(GO_per_Module[[12]]), "yellow.csv", row.names = F)

which(names(KEGG_per_Module) %in% IPdriven)
names(KEGG_per_Module)[which(names(KEGG_per_Module) %in% IPdriven)]

#write.csv(as.data.frame(KEGG_per_Module[[4]]), "Kgreen.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[[5]]), "Kgreenyellow.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[[7]]), "Kmagenta.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[[8]]), "Kpink.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[[9]]), "Kpurple.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[[10]]), "Kred.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[[11]]), "Kturquoise.csv", row.names = F)
#write.csv(as.data.frame(KEGG_per_Module[[12]]), "Kyellow.csv", row.names = F)

#annotate module label based on above:
moduleGenesMM$ModuleSummary[grepl("turquoise", moduleGenesMM$Module)] <-    "Ribosome/respiration)"
moduleGenesMM$ModuleSummary[grepl("^yellow", moduleGenesMM$Module)] <-       "Innate immune-1"
moduleGenesMM$ModuleSummary[grepl("^green", moduleGenesMM$Module)] <-       "Wnt/muscle-1"
moduleGenesMM$ModuleSummary[grepl("red", moduleGenesMM$Module)] <-          "Cell division-1"
moduleGenesMM$ModuleSummary[grepl("^green", moduleGenesMM$Module)] <-       "Wnt/muscle-1"
moduleGenesMM$ModuleSummary[grepl("pink", moduleGenesMM$Module)] <-       "Wnt/muscle-2"
moduleGenesMM$ModuleSummary[grepl("^magenta", moduleGenesMM$Module)] <-       "Actin/cytoskeleton"
moduleGenesMM$ModuleSummary[grepl("purple", moduleGenesMM$Module)] <-       "Cell division-2"
moduleGenesMM$ModuleSummary[grepl("greenyellow", moduleGenesMM$Module)] <-  "Innate immune-2"


#### association of core cell cycle genes with hub-ness ####

moduleGenesMM_cellCycle <- filter(moduleGenesMM, Module %in% c("red", "purple"))
moduleGenesMM_cellCycle$Core_SG2M <- "Other"
moduleGenesMM_cellCycle$Core_SG2M[moduleGenesMM_cellCycle$EnsID %in% filter(fpkm_allG, grepl("CC", GeneClassUpdate))$EnsID] <- "Core_SG2M"

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
  facet_wrap(~ModuleSummary, ncol = 6) +
  xlab("") +
  ylab("Module Membership")

t.test(filter(moduleGenesMM_cellCycle, Module == "red", Core_SG2M == "Core_SG2M")$MM, 
       filter(moduleGenesMM_cellCycle, Module == "red", !Core_SG2M == "Core_SG2M")$MM, var.equal = T)$p.value
t.test(filter(moduleGenesMM_cellCycle, Module == "purple", Core_SG2M == "Core_SG2M")$MM, 
       filter(moduleGenesMM_cellCycle, Module == "purple", !Core_SG2M == "Core_SG2M")$MM, var.equal = T)$p.value

#### key canonical immune signals as hub genes ####

moduleGenesMM_immune <- filter(moduleGenesMM, Module %in% c("greenyellow"))
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
  facet_wrap(~ModuleSummary, ncol = 1) +
  xlab("") +
  ylab("Module Membership")


#### finish table, define hub genes via a threshold, label modules ####

#annotate cclnc association
CoRegPairs_04_48_24_extended <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CoRegPairs_04_48_24_extended_SCClnc_Nov25.csv")
StrongCoRegPairs_04_48_24_extended <- filter(CoRegPairs_04_48_24_extended, 
                                             (corSig == "Yes" | 
                                                !loopMethod == "Neither" | 
                                                eQTLvalidations >0 | 
                                                FANTOM_eQTL == "Yes"))

moduleGenesMM$CCLnc_association <- "None"
moduleGenesMM$CCLnc_association[moduleGenesMM$EnsID %in% CoRegPairs_04_48_24_extended$EnsID] <- "CClncRNA"
moduleGenesMM$CCLnc_association[moduleGenesMM$EnsID %in% CoRegPairs_04_48_24_extended$EnsID.y] <- "CClncRNA target"
moduleGenesMM$CCLnc_association[moduleGenesMM$EnsID %in% StrongCoRegPairs_04_48_24_extended$EnsID] <- "SCClncRNA"
moduleGenesMM$CCLnc_association[moduleGenesMM$EnsID %in% StrongCoRegPairs_04_48_24_extended$EnsID.y] <- "SCClncRNA target"

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

#write.csv(moduleGenesMM, "moduleGenesMM_Nov25.csv", row.names = F)


#### plot hub genes as %, display alongside heatmap and module labels ####

#key info import
moduleGenesMM <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/moduleGenesMM_Nov25.csv")
moduleGenesMM$ModuleSummary <- gsub("\\)", "", moduleGenesMM$ModuleSummary)
IPdriven <- c("greenyellow", "purple", "red", "turquoise", "yellow", "green", "magenta", "pink")

#which modules have highest % of CClncs or CClnc targets in their hubs?
moduleGenesMM_ <- filter(moduleGenesMM, Module %in% IPdriven)

#using top20% here:
for (i in 1:length(IPdriven)){
  hubsCheck <- filter(moduleGenesMM_, Module %in% IPdriven[i], 
                      Hubness %in% c("Top10", "Top20"))
  print(
    paste(round(dim(filter(hubsCheck, grepl("target$", CCLnc_association)))[1]/dim(hubsCheck)[1]*100, 2), 
          IPdriven[i], 
          unique(filter(moduleGenesMM_, Module %in% IPdriven[i])$ModuleSummary),
          "Targets")
  )
}

for (i in 1:length(IPdriven)){
  hubsCheck <- filter(moduleGenesMM_, Module %in% IPdriven[i], 
                      Hubness %in% c("Top10", "Top20"))
  print(
    paste(round(dim(filter(hubsCheck, grepl("CClncRNA$", CCLnc_association)))[1]/dim(hubsCheck)[1]*100, 2), IPdriven[i], "CClncRNAs")
    )
}

#as percentage of DEGs in hub genes:
HubGenesMM_top <- filter(moduleGenesMM_, Hubness %in% c("Top10", "Top20"), DEG_Simple == "DEG")

ggplot(filter(HubGenesMM_top, !is.na(CCLnc_association), !is.na(ModuleSummary), !CCLnc_association %in% c("OtherDE", "None"))) + 
  aes(y = ModuleSummary, fill = CCLnc_association) +
  geom_bar(position = position_dodge(width = 0.8))

table(HubGenesMM_top$ModuleSummary)
table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$ModuleSummary)

trial <- data.frame("CClncRNA" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$ModuleSummary)[1,]/table(HubGenesMM_top$ModuleSummary)*100,
                    "SCClncRNA" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$ModuleSummary)[4,]/table(HubGenesMM_top$ModuleSummary)*100)[,-3]
trial <- reshape2::melt(trial)

triali <- data.frame("CClncRNAn" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$ModuleSummary)[1,],
                     "CClncRNA targetn" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$ModuleSummary)[4,])
triali <- reshape2::melt(triali)

trial$valueN <- triali$value

trial$variable <- gsub("\\.", " ", trial$variable)
trial$variable <- gsub(" Freq", "", trial$variable)

modules_hubs_cclncs <- trial
modules_hubs_cclncs$CClncRNA.Var1 <- factor(modules_hubs_cclncs$CClncRNA.Var1)
#ordering here, options are:
#a) by theme (makes most sense probs)
modules_hubs_cclncs$CClncRNA.Var1 <- factor(modules_hubs_cclncs$CClncRNA.Var1, levels(modules_hubs_cclncs$CClncRNA.Var1)[c(6,8,7,1,3,2,5,4)])
#b) by % in Hubs
#trial$CClncRNA.Var1 <- factor(trial$CClncRNA.Var1, levels(trial$CClncRNA.Var1)[c(6,3,2,8,7,1,4,5)])

triali <- aggregate(valueN ~ CClncRNA.Var1, data = modules_hubs_cclncs, sum)

#modules with lots of CClncs in their hub genes
ggplot(modules_hubs_cclncs) + aes(y = CClncRNA.Var1, fill = variable, x = value, label = valueN) +
  geom_bar(stat = "identity", position = position_stack()) +
  geom_label(data = triali, size = 8, color = "black", inherit.aes = F, label.size = 0.8,
             label.padding = unit(0.17, "lines"),
             aes(y = CClncRNA.Var1, x = 7, label = valueN)) +
  xlab("% of DE Hub Genes") +
  ylab("") +
  scale_fill_manual(values = c("CClncRNA" = "olivedrab4", "SCClncRNA" = "olivedrab2")) +
  scale_x_continuous(limits = c(0,8), breaks = c(0,5)) +
  theme_minimal() +
  theme(text = element_text(size=24))

#targets:
trial <- data.frame("CClncRNA target" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$ModuleSummary)[2,]/table(HubGenesMM_top$ModuleSummary)*100,
                    "SCClncRNA target" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$ModuleSummary)[5,]/table(HubGenesMM_top$ModuleSummary)*100)[,-3]
trial <- reshape2::melt(trial)

triali <- data.frame("CClncRNA targetn" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$ModuleSummary)[2,],
                     "SCClncRNA targetn" = table(HubGenesMM_top$CCLnc_association, HubGenesMM_top$ModuleSummary)[5,])
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
             aes(y = CClncRNA.target.Var1, x = 25, label = valueN)) +
  xlab("% of DE Hub Genes") +
  ylab("") +
  scale_fill_manual(values = c("CClncRNA target" = "mediumorchid4", "SCClncRNA target" = "mediumorchid1")) +
  scale_x_continuous(limits = c(0,27), breaks = c(0,10,20)) +
  theme_minimal() +
  theme(text = element_text(size=24))


#### more detailed heatmap (run early code to remake) ####

#improved version of heatmap after annotating lncs/targets + GO/KEGG characterising:
colnames(mME)[1] <- "Patient"

mME$Patient <- as.character(mME$Patient)
mME$Patient <- gsub("Pt134", "Pt1", mME$Patient)
mME$Patient <- gsub("Pt151", "Pt2", mME$Patient)
mME$Patient <- gsub("Pt157", "Pt3", mME$Patient)
mME$Patient <- gsub("Pt2279", "Pt4", mME$Patient)

mME$Patient <- gsub(" - 0", "-0hr", mME$Patient)
mME$Patient <- gsub(" - 4", "-4hr", mME$Patient)
mME$Patient <- gsub(" - 8", "-8hr", mME$Patient)
mME$Patient <- gsub(" - 24", "-24hr", mME$Patient)

#module names again (re-align later to fit bar charts)

trial <- merge(mME, unique(moduleGenesMM_[,c(6,12)]), by.x = "name", by.y = "Module")

trial$Patient <- as.character(trial$Patient)
trial$Patient <- as.factor(trial$Patient)
trial$Patient <- factor(trial$Patient, levels = levels(trial$Patient)[c(1,5,9,13,
                                                                  3,7,11,15,
                                                                  4,8,12,16,
                                                                  2,6,10,14)])

trial <- filter(trial, name %in% IPdriven)

#can order by theme as here or order by no. cclncRNA targets (as in later figure)
trial$nameII <- as.character(trial$ModuleSummary)
trial$nameII <- as.factor(trial$nameII)
trial$nameII <- factor(trial$nameII, levels = levels(trial$nameII)[c(6,8,7,1,3,2,5,4)])

ggplot(trial, aes(x=Patient, y=nameII, fill=value)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_gradient2(
    low = "steelblue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(text = element_text(size=24),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=11)) +labs(y = "", x ="", fill="Correlation")



#### association of CClncRNA targets and hub-ness ####

#annotate cclnc targets
#just 2d as want to compare fairly to PCG
#not doing this anymore - there were scant few extra HiC
#AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_info_Aug25.csv", header = T)
#AllLNC_AllPCG_1 <- filter(AllLNC_AllPCG_1, DisLnc_PCG <400)

#clusters needed
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

CoRegPairs_04_48_24_extended <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CoRegPairs_04_48_24_extended_SCClnc_Nov25.csv")
CoRegPairs_04_48_24_extended$lncNameCol <- CoRegPairs_04_48_24_extended$EnsID
CoRegPairs_04_48_24_extended$lncNameCol[!is.na(CoRegPairs_04_48_24_extended$EnsName.x)] <- CoRegPairs_04_48_24_extended$EnsName.x[!is.na(CoRegPairs_04_48_24_extended$EnsName.x)]
#364 at 400kbp after HiC

#compare to genes with DE PCG nearby
AllPCG_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllPCG_AllPCG_1_Aug25.csv")

#2D pairs PCG
CoRegPairs_04_48_24_extendedPCG <- filter(AllPCG_AllPCG_1,
                                          #AllLNC_AllPCG_1,
                                          (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                        fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                            (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                          fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                           fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                            (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                          fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
)
#5288 pairs @250kbp
#@400kbp 7749
length(unique(CoRegPairs_04_48_24_extendedPCG$EnsID))#2743 PCGs with DE PCG neighbour

#annotate cclnc targets
moduleGenesMM$CoRegLnc <- "NotDE"
moduleGenesMM$CoRegLnc[moduleGenesMM$EnsID %in% c(fpkm_allGDE$EnsID)] <- "DE" #all DE for now, consider revisiting
table(moduleGenesMM$CoRegLnc, moduleGenesMM$Module)[,IPdriven]

moduleGenesMM$CoRegLnc[moduleGenesMM$EnsID %in% c(CoRegPairs_04_48_24_extendedPCG$EnsID)] <- "Target of CCPCG Only"
#moduleGenesMM$CoRegLnc[moduleGenesMM$EnsID %in% CoRegPairs_04_48_24_extended_naive$EnsID.y] <- "Target of CClncRNA"
moduleGenesMM$CoRegLnc[moduleGenesMM$EnsID %in% CoRegPairs_04_48_24_extended$EnsID.y] <- "Target of CClncRNA"

#reasonable number to test in all modules
table(moduleGenesMM$CoRegLnc) #340 targets of CClncs
table(moduleGenesMM$CoRegLnc, moduleGenesMM$Module)[,IPdriven]
sum(table(moduleGenesMM$CoRegLnc, moduleGenesMM$Module)[,IPdriven])

library(ggbeeswarm)

#main comparisons are between, no regulator, lncRNA regulator and PCG regulator
#does having a lnc as a regulator make it more likely to be hub than PCG or none?
#omit non-DE, don't need to do that test again... cos
#a)this was done earlier to establish whether a module was of interest or not
#b)all the lnc targeted are DE

#do co regulated PCGs have greater MM than normal DE PCG? only comparison interested in:
#MM per biotype per cluster:
moduleGenesMM$biotypeSimpl <- moduleGenesMM$GeneClassUpdate
moduleGenesMM$biotypeSimpl[grepl("Bona|VLnc|ELnc", moduleGenesMM$GeneClassUpdate)] <- "LncRNA"
moduleGenesMM$biotypeSimpl[grepl("coding|TF|CC", moduleGenesMM$GeneClassUpdate)] <- "PCG"

#just testing DE PCGs
test_CoRegMM <- filter(moduleGenesMM, biotypeSimpl == "PCG", !CoRegLnc %in% c("NotDE"))
table(test_CoRegMM$CoRegLnc, test_CoRegMM$Module)[,IPdriven]

#small n numbers, non-normal dist, equal variance?

ggplot(filter(test_CoRegMM, Module %in% IPdriven)) + aes(x = CoRegLnc, y = MM, color = CoRegLnc) +
  geom_quasirandom(alpha = 0.7) +
  #scale_y_log10() +
  coord_cartesian(ylim = c(0.6,1)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.5) +
  theme_minimal() +
  facet_wrap(~Module, ncol =4) +
  #  scale_y_log10() +
  xlab("") +
  theme(axis.text.x = element_blank())
#would probs work better doing individual scales per plot, but for now works
#in greenyellow and purple, being targeted by an early lncRNA gives more centrality than others
#sustained upregulated immune cluster, and purple cluster

#suggests a kruskal-wallis per module:
KW_CoReg_MM <- list()
for (i in c(1:length(IPdriven))){
  KW_CoReg_MM[[i]] <- kruskal.test(MM ~ CoRegLnc, filter(test_CoRegMM, Module == IPdriven[i]))$p.value
}
unlist(KW_CoReg_MM) <0.05
p.adjust(unlist(KW_CoReg_MM), method = "BH")
p.adjust(unlist(KW_CoReg_MM), method = "BH") <0.05
IPdriven
IPdriven[p.adjust(unlist(KW_CoReg_MM), method = "BH") <0.05]

#dunn test more appropriate than wilcox as a post-hoc apparently...
#https://www.theanalysisfactor.com/dunns-test-post-hoc-test-after-kruskal-wallis/
#https://www.reddit.com/r/statistics/comments/15fk7j1/q_kruskalwallis_multiple_testing/

#one sided dunn test, focus on comparison vs. cclncRNA only 
dt <- dunn.test::dunn.test(filter(test_CoRegMM, Module == IPdriven[1])$MM, 
                     filter(test_CoRegMM, Module == IPdriven[1])$CoRegLnc, method = "bonferroni")
dt$P[c(1,3)] *2

Dunn_TestOutPut <- list()
IPdriven_2test <- IPdriven[p.adjust(unlist(KW_CoReg_MM), method = "BH") <0.05]

for (i in c(1:length(IPdriven_2test))){
  dt <- dunn.test::dunn.test(filter(test_CoRegMM, Module == IPdriven_2test[i])$MM, 
                             filter(test_CoRegMM, Module == IPdriven_2test[i])$CoRegLnc, method = "bh")
  Dunn_TestOutPut[[i]] <- dt$P[c(1,3)] *2
}
names(Dunn_TestOutPut) <- IPdriven_2test
Dunn_TestOutPut
#significant vs. other but not PCG
#greenyellow closer
#but almost certainly driven by CXCL x3, not a thing across all lncs


#re-run with same timeframe targets only:
CoRegPairs_04_48_24_extendedSame <- filter(CoRegPairs_04_48_24_extended,
                                                 #AllLNC_AllPCG_1,
                                                 (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                               fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                                fpkm_allGDE_Downwithin_4$EnsID)) |
                                                   (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                 fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                                  fpkm_allGDE_Downwithin_8$EnsID)) |
                                                   (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                 fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
)#203 pairs @400kbp

CoRegPairs_04_48_24_samePCG <- filter(AllPCG_AllPCG_1,
                                      #AllLNC_AllPCG_1,
                                      (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                    fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                     fpkm_allGDE_Downwithin_4$EnsID)) |
                                        (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                      fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                       fpkm_allGDE_Downwithin_8$EnsID)) |
                                        (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                      fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                        fpkm_allGDE_Downwithin_24$EnsID))
)


moduleGenesMM$CoRegLncSame <- "NotDE"
moduleGenesMM$CoRegLncSame[moduleGenesMM$EnsID %in% c(fpkm_allGDE$EnsID)] <- "DE" #all DE for now, consider revisiting
table(moduleGenesMM$CoRegLncSame, moduleGenesMM$Module)[,IPdriven]

moduleGenesMM$CoRegLncSame[moduleGenesMM$EnsID %in% c(CoRegPairs_04_48_24_samePCG$EnsID)] <- "Target of CCPCG Only"
moduleGenesMM$CoRegLncSame[moduleGenesMM$EnsID %in% CoRegPairs_04_48_24_extendedSame$EnsID.y] <- "Target of CClncRNA"

#reasonable number to test in all modules
table(moduleGenesMM$CoRegLncSame) #195 targets of CClncs
table(moduleGenesMM$CoRegLncSame, moduleGenesMM$Module)[,IPdriven]

test_CoRegMM <- filter(moduleGenesMM, biotypeSimpl == "PCG", !CoRegLnc %in% c("NotDE"))
table(test_CoRegMM$CoRegLncSame, test_CoRegMM$Module)[,IPdriven]

ggplot(filter(test_CoRegMM, Module %in% IPdriven)) + aes(x = CoRegLncSame, y = MM, color = CoRegLncSame) +
  geom_quasirandom(alpha = 0.7) +
  #scale_y_log10() +
  coord_cartesian(ylim = c(0.6,1)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.5) +
  theme_minimal() +
  facet_wrap(~Module, ncol =4) +
  #  scale_y_log10() +
  xlab("") +
  theme(axis.text.x = element_blank())
#would probs work better doing individual scales per plot, but for now works
#in greenyellow and purple, being targeted by an early lncRNA gives more centrality than others
#sustained upregulated immune cluster, and purple cluster

#suggests a kruskal-wallis per module:
KW_CoReg_MMsame <- list()
for (i in c(1:length(IPdriven))){
  KW_CoReg_MMsame[[i]] <- kruskal.test(MM ~ CoRegLncSame, filter(test_CoRegMM, Module == IPdriven[i]))$p.value
}
unlist(KW_CoReg_MMsame) <0.05
p.adjust(unlist(KW_CoReg_MMsame), method = "BH")
p.adjust(unlist(KW_CoReg_MMsame), method = "BH") <0.05
IPdriven
IPdriven[p.adjust(unlist(KW_CoReg_MMsame), method = "BH") <0.05]
#another module comes in here

#dunn test more appropriate than wilcox as a post-hoc apparently...
#https://www.theanalysisfactor.com/dunns-test-post-hoc-test-after-kruskal-wallis/
#https://www.reddit.com/r/statistics/comments/15fk7j1/q_kruskalwallis_multiple_testing/

#one sided dunn test, focus on comparison vs. cclncRNA only 
Dunn_TestOutPutsame <- list()
IPdriven_2test <- IPdriven[p.adjust(unlist(KW_CoReg_MMsame), method = "BH") <0.05]

for (i in c(1:length(IPdriven_2test))){
  dt <- dunn.test::dunn.test(filter(test_CoRegMM, Module == IPdriven_2test[i])$MM, 
                             filter(test_CoRegMM, Module == IPdriven_2test[i])$CoRegLncSame, method = "bh")
  Dunn_TestOutPutsame[[i]] <- dt$P[c(1,3)] *2
}
names(Dunn_TestOutPutsame) <- IPdriven_2test
Dunn_TestOutPutsame
#significant vs. other and PCG for yellow
#greenyellow close at 0.08, but almost certainly driven by CXCL locus
#immune response targets of cclncRNAs have elevated hub-ness

#pushing the stats a bit but reasonable
#even if including the 3rd comparison still fairly sig
#no issue from HOX, CXCL etc

#could also compare just lncs/PCGs:
ggplot(filter(test_CoRegMM, !CoRegLncSame == "DE", Module %in% IPdriven)) + aes(x = CoRegLncSame, y = MM, color = CoRegLncSame) +
  geom_quasirandom(alpha = 0.7) +
  #scale_y_log10() +
  coord_cartesian(ylim = c(0.6,1)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.5) +
  theme_minimal() +
  facet_wrap(~Module, ncol =4) +
  #  scale_y_log10() +
  xlab("") +
  theme(axis.text.x = element_blank())

table(filter(test_CoRegMM, !CoRegLncSame == "DE", Module %in% IPdriven)$CoRegLncSame, 
      filter(test_CoRegMM, !CoRegLncSame == "DE", Module %in% IPdriven)$Module)

#simple wilcox:
wilcox_p_hubness <- list()
wilcox_means_hubness <- list()

for (i in 1:length(IPdriven)){
  wilcox_p_hubness[[i]] <- wilcox.test(filter(moduleGenesMM, Module %in% IPdriven[i], CoRegLncSame == "Target of CClncRNA")$MM, 
                                 filter(moduleGenesMM, Module %in% IPdriven[i], CoRegLncSame == "Target of CCPCG Only")$MM, var.equal = T)$p.value
  wilcox_means_hubness[[i]] <- c(mean(filter(moduleGenesMM, Module %in% IPdriven[i], CoRegLncSame == "Target of CClncRNA")$MM)/
                                  mean(filter(moduleGenesMM, Module %in% IPdriven[i], CoRegLncSame == "Target of CCPCG Only")$MM))  
}

names(wilcox_p_hubness) <- IPdriven
p.adjust(wilcox_p_hubness, method = "BH")
#actually less successful, kruskal step helps filter out weaker comparisons


#### association of SCClncRNA targets and hub-ness ####

#SCClncs possibly bit more sense to use here, better stats for GO/KEGG etc enrichment

#binned this off bcos too trivial to care about/explain
#annotate cclnc targets:
#just 2d as want to compare fairly to PCG
#AllLNC_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllLNC_AllPCG_info_Aug25.csv", header = T)
#AllLNC_AllPCG_1 <- filter(AllLNC_AllPCG_1, DisLnc_PCG <400)

#clusters needed for later
#regulated within 4 hours:
fpkm_allGDE_Upwithin_4 <- filter(fpkm_allGDE, (Hour0_meanFPKM >=1 | Hour4_meanFPKM >=1) &
                                   (LogFC_0_4 >= log2(1.5) & preadj_0_4 <0.05))
fpkm_allGDE_Downwithin_4 <- filter(fpkm_allGDE, (Hour0_meanFPKM >=1 | Hour4_meanFPKM >=1) &
                                     (LogFC_0_4 < -log2(1.5) & preadj_0_4 <0.05))
#of remaining, within 8 hours:
fpkm_allGDE_Upwithin_8 <- filter(fpkm_allGDE, !EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                                 ((Hour0_meanFPKM >=1 | Hour8_meanFPKM >=1) &
                                    (LogFC_0_8 >= log2(1.5) & preadj_0_8 <0.05))
                                 |
                                   ((Hour4_meanFPKM >=1 | Hour8_meanFPKM >=1) &
                                      (LogFC_4_8 >= log2(1.5) & preadj_4_8 <0.05)))
fpkm_allGDE_Downwithin_8 <- filter(fpkm_allGDE, !EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                                   ((Hour0_meanFPKM >=1 | Hour8_meanFPKM >=1) &
                                      (LogFC_0_8 < -log2(1.5) & preadj_0_8 <0.05))
                                   |
                                     ((Hour4_meanFPKM >=1 | Hour8_meanFPKM >=1) &
                                        (LogFC_4_8 < -log2(1.5) & preadj_4_8 <0.05)))
#of remaining, within 24 hours:
fpkm_allGDE_Upwithin_24 <- filter(fpkm_allGDE, !EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                                  !EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID),
                                  ((Hour0_meanFPKM >=1 | Hour24_meanFPKM >=1) &
                                     (LogFC_0_24 >= log2(1.5) & preadj_0_24 <0.05))
                                  |
                                    ((Hour4_meanFPKM >=1 | Hour24_meanFPKM >=1) &
                                       (LogFC_4_24 >= log2(1.5) & preadj_4_24 <0.05))
                                  |
                                    ((Hour8_meanFPKM >=1 | Hour24_meanFPKM >=1) &
                                       (LogFC_8_24 >= log2(1.5) & preadj_8_24 <0.05)))
fpkm_allGDE_Downwithin_24 <- filter(fpkm_allGDE, !EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, fpkm_allGDE_Downwithin_4$EnsID),
                                    !EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID),
                                    ((Hour0_meanFPKM >=1 | Hour24_meanFPKM >=1) &
                                       (LogFC_0_24 < -log2(1.5) & preadj_0_24 <0.05))
                                    |
                                      ((Hour4_meanFPKM >=1 | Hour24_meanFPKM >=1) &
                                         (LogFC_4_24 < -log2(1.5) & preadj_4_24 <0.05))
                                    |
                                      ((Hour8_meanFPKM >=1 | Hour24_meanFPKM >=1) &
                                         (LogFC_8_24 < -log2(1.5) & preadj_8_24 <0.05)))

CoRegPairs_04_48_24_extended <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/CoRegPairs_04_48_24_extended_SCClnc_Nov25.csv")
CoRegPairs_04_48_24_extended$lncNameCol <- CoRegPairs_04_48_24_extended$EnsID
CoRegPairs_04_48_24_extended$lncNameCol[!is.na(CoRegPairs_04_48_24_extended$EnsName.x)] <- CoRegPairs_04_48_24_extended$EnsName.x[!is.na(CoRegPairs_04_48_24_extended$EnsName.x)]

StrongCoRegPairs_04_48_24_extended <- filter(CoRegPairs_04_48_24_extended, 
                                             (corSig == "Yes" | 
                                                !loopMethod == "Neither" | 
                                                eQTLvalidations >0 | 
                                                FANTOM_eQTL == "Yes"))
#163 pairs at 400kbp after HiC

AllPCG_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllPCG_AllPCG_1_Aug25.csv")

#2D pairs PCG
CoRegPairs_04_48_24_extendedPCG <- filter(AllPCG_AllPCG_1,
                                          #AllLNC_AllPCG_1,
                                          (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                        fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                            (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                          fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                           fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                            (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                          fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID))
)
#5288 pairs @250kbp
#@400kbp 7749
length(unique(CoRegPairs_04_48_24_extendedPCG$EnsID))#2743 PCGs with DE PCG neighbour

#annotate cclnc targets - "Co-regulated neighbours include:"
moduleGenesMM$CoRegLnc <- "NotDE"
moduleGenesMM$CoRegLnc[moduleGenesMM$EnsID %in% c(fpkm_allGDE$EnsID)] <- "Neither" #all DE for now, consider revisiting
table(moduleGenesMM$CoRegLnc, moduleGenesMM$Module)[,IPdriven]

moduleGenesMM$CoRegLnc[moduleGenesMM$EnsID %in% c(CoRegPairs_04_48_24_extendedPCG$EnsID)] <- "PCGs only"
#moduleGenesMM$CoRegLnc[moduleGenesMM$EnsID %in% CoRegPairs_04_48_24_extended_naive$EnsID.y] <- "Target of CClncRNA"
moduleGenesMM$CoRegLnc[moduleGenesMM$EnsID %in% StrongCoRegPairs_04_48_24_extended$EnsID.y] <- "SCClncRNAs"

#reasonable number to test in all modules
table(moduleGenesMM$CoRegLnc) #340 targets of CClncs
table(moduleGenesMM$CoRegLnc, moduleGenesMM$Module)[,IPdriven]
sum(table(moduleGenesMM$CoRegLnc, moduleGenesMM$Module)[,IPdriven])

library(ggbeeswarm)

#main comparisons are between, no potential regulatory neighbour vs. a lncRNA one vs. a PCG regulator
#does having a lnc as a regulator make it more likely to be hub than PCG or none?
#omit non-DE, don't need to do that test again... cos
#a)this was done earlier to establish whether a module was of interest or not
#b)all the lnc targeted are DE

#do co regulated PCGs have greater MM than normal DE PCG? only comparison interested in:
#MM per biotype per cluster:
moduleGenesMM$biotypeSimpl <- moduleGenesMM$GeneClassUpdate
moduleGenesMM$biotypeSimpl[grepl("Bona|VLnc|ELnc", moduleGenesMM$GeneClassUpdate)] <- "LncRNA"
moduleGenesMM$biotypeSimpl[grepl("coding|TF|CC", moduleGenesMM$GeneClassUpdate)] <- "PCG"

test_CoRegMM <- filter(moduleGenesMM, biotypeSimpl == "PCG", !CoRegLnc %in% c("NotDE"))
table(test_CoRegMM$CoRegLnc, test_CoRegMM$Module)[,IPdriven]

test_CoRegMM$ModuleSummaryII <- paste("\n", test_CoRegMM$ModuleSummary, "\n", sep = "")
test_CoRegMM$ModuleSummaryII <- gsub("respiration", "resp.", test_CoRegMM$ModuleSummaryII)
test_CoRegMM$ModuleSummaryII <- gsub("cytoskeleton", "cytoskel.", test_CoRegMM$ModuleSummaryII)

test_CoRegMM$CoRegLnc <- as.factor(test_CoRegMM$CoRegLnc)
test_CoRegMM$CoRegLnc <- factor(test_CoRegMM$CoRegLnc, levels(test_CoRegMM$CoRegLnc)[c(3,2,1)])

#small n numbers, non-normal dist, equal variance?
ggplot(filter(test_CoRegMM, Module %in% IPdriven)) + aes(x = CoRegLnc, y = MM, fill = CoRegLnc, color = CoRegLnc) +
  #scale_y_log10() +
  #coord_cartesian(ylim = c(0.6,1)) +
  geom_quasirandom(alpha = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.35, alpha = 1, color = "black") +
  scale_color_manual(values = c(`Neither` = "grey70", `SCClncRNAs` = "olivedrab3", `PCGs only` = "mediumorchid1")) +
  scale_fill_manual(values = c(`Neither` = "grey70", `SCClncRNAs` = "olivedrab3", `PCGs only` = "mediumorchid1")) +
  labs(dictionary = c(CoRegLnc = "Co-regulated with")) +
  theme_minimal() +
  facet_wrap(~ModuleSummaryII, ncol =4, scales = "free") +
  #  scale_y_log10() +
  xlab("") +
  theme(text = element_text(size=24),
        strip.text = element_text(size=15),
        #strip.text = element_blank(),
        axis.text.y = element_text(size=15),
        axis.text.x = element_blank(),#legend.position = "none"
        )
#would probs work better doing individual scales per plot, but for now works

#suggests a kruskal-wallis per module:
KW_CoReg_MM <- list()
for (i in c(1:length(IPdriven))){
  KW_CoReg_MM[[i]] <- kruskal.test(MM ~ CoRegLnc, filter(test_CoRegMM, Module == IPdriven[i]))$p.value
}
unlist(KW_CoReg_MM) <0.05
p.adjust(unlist(KW_CoReg_MM), method = "BH")
p.adjust(unlist(KW_CoReg_MM), method = "BH") <0.05
IPdriven
IPdriven[p.adjust(unlist(KW_CoReg_MM), method = "BH") <0.05]

#dunn test more appropriate than wilcox as a post-hoc apparently...
#https://www.theanalysisfactor.com/dunns-test-post-hoc-test-after-kruskal-wallis/
#https://www.reddit.com/r/statistics/comments/15fk7j1/q_kruskalwallis_multiple_testing/

#one sided dunn test, focus on comparison vs. cclncRNA only 
Dunn_TestOutPut_T2 <- list()
IPdriven_2test <- IPdriven[p.adjust(unlist(KW_CoReg_MM), method = "BH") <0.05]

for (i in c(1:length(IPdriven_2test))){
  dt <- dunn.test::dunn.test(filter(test_CoRegMM, Module == IPdriven_2test[i])$MM, 
                             filter(test_CoRegMM, Module == IPdriven_2test[i])$CoRegLnc, method = "bh")
  Dunn_TestOutPut_T2[[i]] <- dt$P[c(2,3)] *2
}
names(Dunn_TestOutPut_T2) <- IPdriven_2test
Dunn_TestOutPut_T2
#significant vs. other but not PCG
#yellow close, seems like improved stats vs. others
Dunn_TestOutPut

#re-run with same timeframe targets only:
StrongCoRegPairs_04_48_24_extendedSame <- filter(StrongCoRegPairs_04_48_24_extended,
                                           #AllLNC_AllPCG_1,
                                           (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                         fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                          fpkm_allGDE_Downwithin_4$EnsID)) |
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                           fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                            fpkm_allGDE_Downwithin_8$EnsID)) |
                                             (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                           fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)))
#106 pairs @400kbp

CoRegPairs_04_48_24_samePCG <- filter(AllPCG_AllPCG_1,
                                      #AllLNC_AllPCG_1,
                                      (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                    fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                                                                     fpkm_allGDE_Downwithin_4$EnsID)) |
                                        (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                      fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                                                                       fpkm_allGDE_Downwithin_8$EnsID)) |
                                        (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                      fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                                                                        fpkm_allGDE_Downwithin_24$EnsID)))


moduleGenesMM$CoRegLncSame <- "NotDE"
moduleGenesMM$CoRegLncSame[moduleGenesMM$EnsID %in% c(fpkm_allGDE$EnsID)] <- "Neither" #all DE for now, consider revisiting
table(moduleGenesMM$CoRegLncSame, moduleGenesMM$Module)[,IPdriven]

moduleGenesMM$CoRegLncSame[moduleGenesMM$EnsID %in% c(CoRegPairs_04_48_24_samePCG$EnsID)] <- "PCGs Only"
moduleGenesMM$CoRegLncSame[moduleGenesMM$EnsID %in% StrongCoRegPairs_04_48_24_extendedSame$EnsID.y] <- "SCClncRNAs"

#reasonable number to test in all modules
table(moduleGenesMM$CoRegLncSame) #195 targets of CClncs
table(moduleGenesMM$CoRegLncSame, moduleGenesMM$Module)[,IPdriven]

test_CoRegMM <- filter(moduleGenesMM, biotypeSimpl == "PCG", !CoRegLnc %in% c("NotDE"))
table(test_CoRegMM$CoRegLncSame, test_CoRegMM$Module)[,IPdriven]

test_CoRegMM$ModuleSummaryII <- paste("\n", test_CoRegMM$ModuleSummary, "\n", sep = "")
test_CoRegMM$ModuleSummaryII <- gsub("respiration", "resp.", test_CoRegMM$ModuleSummaryII)
test_CoRegMM$ModuleSummaryII <- gsub("cytoskeleton", "cytoskel.", test_CoRegMM$ModuleSummaryII)

test_CoRegMM$CoRegLncSame <- as.factor(test_CoRegMM$CoRegLncSame)
test_CoRegMM$CoRegLncSame <- factor(test_CoRegMM$CoRegLncSame, levels(test_CoRegMM$CoRegLncSame)[c(3,2,1)])

#small n numbers, non-normal dist, equal variance?
ggplot(filter(test_CoRegMM, Module %in% IPdriven)) + aes(x = CoRegLncSame, y = MM, fill = CoRegLncSame, color = CoRegLncSame) +
  #scale_y_log10() +
  #coord_cartesian(ylim = c(0.6,1)) +
  geom_quasirandom(alpha = 1) +
  geom_boxplot(outlier.shape = NA, width = 0.35, alpha = 0.6, color = "black") +
  scale_color_manual(values = c(`Neither` = "grey70", `SCClncRNAs` = "olivedrab3", `PCGs Only` = "mediumorchid1")) +
  scale_fill_manual(values = c(`Neither` = "grey70", `SCClncRNAs` = "olivedrab3", `PCGs Only` = "mediumorchid1")) +
  labs(dictionary = c(CoRegLncSame = "Co-regulated with")) +
  theme_minimal() +
  facet_wrap(~ModuleSummaryII, ncol =4, scales = "free") +
  #  scale_y_log10() +
  xlab("") +
  theme(text = element_text(size=24),
        strip.text = element_text(size=15),
        #strip.text = element_blank(),
        axis.text.y = element_text(size=15),
        axis.text.x = element_blank(),legend.position = "bottom"
  )
#would probs work better doing individual scales per plot, but for now works

#suggests a kruskal-wallis per module:
KW_CoReg_MMsame <- list()
for (i in c(1:length(IPdriven))){
  KW_CoReg_MMsame[[i]] <- kruskal.test(MM ~ CoRegLncSame, filter(test_CoRegMM, Module == IPdriven[i]))$p.value
}
unlist(KW_CoReg_MMsame) <0.05
p.adjust(unlist(KW_CoReg_MMsame), method = "BH")
p.adjust(unlist(KW_CoReg_MMsame), method = "BH") <0.05
IPdriven
IPdriven[p.adjust(unlist(KW_CoReg_MMsame), method = "BH") <0.05]
#greenyellow gone

#dunn test more appropriate than wilcox as a post-hoc apparently...
#https://www.theanalysisfactor.com/dunns-test-post-hoc-test-after-kruskal-wallis/
#https://www.reddit.com/r/statistics/comments/15fk7j1/q_kruskalwallis_multiple_testing/

#one sided dunn test, focus on comparison vs. cclncRNA only 
Dunn_TestOutPutsame_T2 <- list()
IPdriven_2test <- IPdriven[p.adjust(unlist(KW_CoReg_MMsame), method = "BH") <0.05]

for (i in c(1:length(IPdriven_2test))){
  dt <- dunn.test::dunn.test(filter(test_CoRegMM, Module == IPdriven_2test[i])$MM, 
                             filter(test_CoRegMM, Module == IPdriven_2test[i])$CoRegLncSame, method = "bh")
  Dunn_TestOutPutsame_T2[[i]] <- dt$P[c(2,3)] *2
}
names(Dunn_TestOutPutsame_T2) <- IPdriven_2test
Dunn_TestOutPutsame_T2
#yellow module immune response targets of cclncRNAs have elevated hub-ness

#scclncRNAs, give a clear answer
#filtering for PCGs correlated with a lnc.. explains in part

#running without correlated PCGs gives greenyellow v. strongly significant (but only 4 genes)

#version for paper figure

#better labels:
test_CoRegMM$CoregulatedNeighbour <- test_CoRegMM$CoRegLncSame
test_CoRegMM$CoregulatedNeighbour[test_CoRegMM$CoRegLncSame_read == "Target of SCClncRNA"] <- "SCClncRNA"
test_CoRegMM$CoregulatedNeighbour[test_CoRegMM$CoRegLncSame_read == "Target of CCPCG Only"] <- "DE PCGs, no SCClncRNAs"
test_CoRegMM$CoregulatedNeighbour[test_CoRegMM$CoRegLncSame_read == "DE"] <- "Neither DE PCGs or SCClncRNAs"

ggplot(filter(test_CoRegMM, Module %in% IPdriven)) + aes(x = CoregulatedNeighbour, y = MM, color = CoregulatedNeighbour) +
  geom_quasirandom(alpha = 0.7) +
  #scale_y_log10() +
  coord_cartesian(ylim = c(0.5,1)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.5) +
  theme_minimal() +
  facet_wrap(~Module, ncol =4) +
  #  scale_y_log10() +
  xlab("") +
  ylab("Network Connectivity (MM)") +
  labs(color = "Co-regulated neighbours") +
  theme(axis.text.x = element_blank())

