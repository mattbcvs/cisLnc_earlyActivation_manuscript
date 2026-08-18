#### alternatives to WGCNA hub genes ####

#1 increased FC for targets vs PCGs

#### 1 compare 0-4hr FCs for SCClncRNA targets vs. co-reg targets with PCGs ####

#to keep fair use 250kbp distance
SCClncRNAs_250 <- filter(SCClncRNAs, AbsDistLnc_PCG <250)

#import PCG-PCG pairs within 250kbp
AllPCG_AllPCG_1 <- read.csv("\\\\cmvm.datastore.ed.ac.uk/cmvm/scs/groups/lncRNA_orthology/Timecourse/AllPCG_AllPCG_1_2026_250.csv")

#2D pairs PCG
CoRegPairs_04_48_24_extendedPCG <- filter(AllPCG_AllPCG_1,
                                          (EnsID %in% c(fpkm_allGDE_Upwithin_4$EnsID, 
                                                        fpkm_allGDE_Downwithin_4$EnsID) & EnsID.y %in% fpkm_allGDE$EnsID) |
                                            (EnsID %in% c(fpkm_allGDE_Upwithin_8$EnsID, 
                                                          fpkm_allGDE_Downwithin_8$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_8$EnsID, fpkm_allGDE_Downwithin_8$EnsID,
                                                                                                           fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)) |
                                            (EnsID %in% c(fpkm_allGDE_Upwithin_24$EnsID, 
                                                          fpkm_allGDE_Downwithin_24$EnsID) & EnsID.y %in% c(fpkm_allGDE_Upwithin_24$EnsID, fpkm_allGDE_Downwithin_24$EnsID)))
#7288 pairs @250kbp
length(unique(CoRegPairs_04_48_24_extendedPCG$EnsID))#2927 PCGs with DE PCG neighbour


#assign labels to the early induced PCGs
PCGDE_Upwithin_4 <- filter(fpkm_allGDE_Upwithin_4, EnsType == "protein_coding", grepl("coding|TF|CC", GeneClassUpdate))

#now only co-induced pairs
PCGDE_Upwithin_4$DEneighbourType <- "No"

CoInducedSCClncRNA_4hr <- filter(SCClncRNAs_250, 
                              EnsID %in% fpkm_allGDE_Upwithin_4$EnsID, 
                              EnsID.y %in% fpkm_allGDE_Upwithin_4$EnsID)

PCGDE_Upwithin_4$DEneighbourType[PCGDE_Upwithin_4$EnsID %in% CoInducedSCClncRNA_4hr$EnsID.y] <- "Yes"

table(PCGDE_Upwithin_4$DEneighbourType)#36 PCG coinduced with a SCCLnc in 4hrs, 315 are coinduced with a PCG

PCGDE_Upwithin_4$`Concordant lncRNA\nneighbour` <- PCGDE_Upwithin_4$DEneighbourType

ggplot(PCGDE_Upwithin_4) + aes(x = `Concordant lncRNA\nneighbour`, y = abs(LogFC_0_4), color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5, color = "grey50") +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10(breaks = c(0.3,1,3,10)) +
  scale_color_manual(values = c("Yes" = "#D6604D", "No" = "#67A9CF"))+
  xlab("Induced PCGs")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Log2FC (0-4hrs)")

#no increase in FC induction for SCClncRNA targets vs other amongst 0-4hr induced genes
t.test(PCGDE_Upwithin_4$LogFC_0_4 ~ PCGDE_Upwithin_4$`Concordant lncRNA\nneighbour`, PCGDE_Upwithin_4, var.equal = T)
t.test(PCGDE_Upwithin_4$LogFC_0_4 ~ PCGDE_Upwithin_4$`Concordant lncRNA\nneighbour`, PCGDE_Upwithin_4, var.equal = F)
wilcox.test(PCGDE_Upwithin_4$LogFC_0_4 ~ PCGDE_Upwithin_4$`Concordant lncRNA\nneighbour`, PCGDE_Upwithin_4)

#so no sense in adding the PCG view here

ggplot(PCGDE_Upwithin_4) + aes(x = `Concordant lncRNA\nneighbour`, y = abs(Hour4_meanFPKM), color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5, color = "grey50") +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10() +
  scale_color_manual(values = c("Yes" = "#D6604D", "No" = "#67A9CF"))+
  xlab("Induced PCGs")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Log2FC (0-4hrs)")

t.test(PCGDE_Upwithin_4$Hour4_meanFPKM ~ PCGDE_Upwithin_4$`Concordant lncRNA\nneighbour`, PCGDE_Upwithin_4, var.equal = F)

#driven by outlier?
PCGDE_Upwithin_4_lo <- filter(PCGDE_Upwithin_4, Hour4_meanFPKM <1000)
ggplot(PCGDE_Upwithin_4_lo) + aes(x = `Concordant lncRNA\nneighbour`, y = abs(Hour4_meanFPKM), color = `Concordant lncRNA\nneighbour`) +
  geom_quasirandom(alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.5, color = "grey50") +
  #geom_jitter(width = 0.3, alpha = 0.1) +
  #facet_wrap(~Cluster) +
  scale_y_log10() +
  scale_color_manual(values = c("Yes" = "#D6604D", "No" = "#67A9CF"))+
  xlab("Induced PCGs")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.background = element_rect()) +
  ylab("Log2FC (0-4hrs)")

#seems driven by the outliers
t.test(filter(PCGDE_Upwithin_4_lo, DEneighbourType == "No")$Hour4_meanFPKM, 
       filter(PCGDE_Upwithin_4_lo, DEneighbourType == "Yes")$Hour4_meanFPKM, var.equal = F)


#### 2 increased FC/FPKM ranking for targets vs PCGs ####

#but some FCs are also v high FPKMs:

ggplot(PCGDE_Upwithin_4) + aes(x = Hour4_meanFPKM, y = LogFC_0_4) +
  geom_point(alpha = 0.5) +
  scale_x_log10(labels = scales::label_comma())

#colored by co-induction

ggplot(PCGDE_Upwithin_4) + aes(x = Hour4_meanFPKM, y = LogFC_0_4, color = DEneighbourType) +
  geom_point(alpha = 0.5) +
  scale_x_log10(labels = scales::label_comma()) + 
  facet_wrap(~DEneighbourType) +
  theme_bw()

#could be that the centre of gravity shifts for scclncRNA targets? density plot:
ggplot(PCGDE_Upwithin_4) + aes(x = Hour4_meanFPKM, y = LogFC_0_4, color = DEneighbourType) +
  #geom_point(alpha = 0.25) +
  geom_density_2d(bins = 14) +
  scale_x_log10(labels = scales::label_comma()) + 
  theme_bw()

#suggests that is the case, nut not representative of the hour4 FPKM change (implies SCClncRNA targets generally lower)

#MANOVA, two variables at once:
dep_vars <- cbind(PCGDE_Upwithin_4$Hour4_meanFPKM, PCGDE_Upwithin_4$LogFC_0_4)
fit <- manova(dep_vars ~ DEneighbourType, data = PCGDE_Upwithin_4)
summary(fit)

library(effectsize)
effectsize::eta_squared(fit)
#would need to read up... suggests a very small effect but quite significant

#ELM suggests this is akin to "Hotellings T Test" when MANOVA is just two categories
hotelling.test(dep_vars ~ DEneighbourType, data = PCGDE_Upwithin_4)
#would need to investigate...

#driven by CXCL IL6 outliers? maybe a little but still seems strong
ggplot(filter(PCGDE_Upwithin_4, Hour4_meanFPKM < 500), aes(x = Hour4_meanFPKM, y = LogFC_0_4, color = DEneighbourType)) +
  geom_density_2d(bins = 12) +
  scale_x_log10(labels = scales::label_comma()) + 
  theme_bw()

#maybe a composite score is easiest to implement:
PCGDE_Upwithin_4$FC_FPKM_comp <- (PCGDE_Upwithin_4$Hour4_meanFPKM * PCGDE_Upwithin_4$LogFC_0_4)/2

ggplot(filter(PCGDE_Upwithin_4, Hour4_meanFPKM < 10000), aes(y = FC_FPKM_comp, x = DEneighbourType)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.3) +
  scale_y_log10(labels = scales::label_comma()) + 
  theme_bw()

t.test(FC_FPKM_comp ~ DEneighbourType, PCGDE_Upwithin_4)

#doesn't seem a fair rep of the distinctness on density though

#MANOVA with PCG co-regs as third group:
#now only co-induced pairs
PCGDE_Upwithin_4$DEneighbourType2 <- "No"

CoInducedSCClncRNA_4hr <- filter(SCClncRNAs_250, 
                                 EnsID %in% fpkm_allGDE_Upwithin_4$EnsID, 
                                 EnsID.y %in% fpkm_allGDE_Upwithin_4$EnsID)

CoInducedPCG_4hr <- filter(AllPCG_AllPCG_1, 
                                 EnsID %in% fpkm_allGDE_Upwithin_4$EnsID, 
                                 EnsID.y %in% fpkm_allGDE_Upwithin_4$EnsID)


PCGDE_Upwithin_4$DEneighbourType2[PCGDE_Upwithin_4$EnsID %in% CoInducedPCG_4hr$EnsID.y] <- "PCG"
PCGDE_Upwithin_4$DEneighbourType2[PCGDE_Upwithin_4$EnsID %in% CoInducedSCClncRNA_4hr$EnsID.y] <- "SCClncRNA"

ggplot(PCGDE_Upwithin_4) + aes(x = Hour4_meanFPKM, y = LogFC_0_4, color = DEneighbourType2) +
  geom_point() +
  #geom_density_2d(bins = 14) +
  scale_x_log10(labels = scales::label_comma()) + 
  theme_bw()

fit <- manova(dep_vars ~ DEneighbourType2, data = PCGDE_Upwithin_4)
summary(fit)

effectsize::eta_squared(fit)
library(MASS)
post_hoc <- lda(PCGDE_Upwithin_4$DEneighbourType2 ~ dep_vars, CV=F)
post_hoc

plot_lda <- data.frame(PCGDE_Upwithin_4[, "DEneighbourType2"], lda = predict(post_hoc)$x)
ggplot(plot_lda) + geom_point(aes(x = lda.LD1, y = lda.LD2, colour = PCGDE_Upwithin_4....DEneighbourType2..), size = 2)

