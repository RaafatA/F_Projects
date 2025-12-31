############################################
##               PACKAGES                 ##
############################################
library(WGCNA)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(Hmisc)
library(dplyr)
library(igraph)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

############################################
##           SAMPLE METADATA              ##
############################################

sample_ids <- colnames(final_xena_identifier_matched)

coldata <- data.frame(
  sample_id = sample_ids,
  Condition = ifelse(grepl("^TCGA", sample_ids), "Cancer", "Normal")
)
rownames(coldata) <- sample_ids
coldata$Condition <- factor(coldata$Condition, levels = c("Normal","Cancer"))

############################################
##      NORMALIZATION (DESeq2)            ##
############################################

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ Condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)

datExpr0 <- t(assay(vsd))   # samples × genes

############################################
##      FILTER LOW VARIANCE GENES         ##
############################################

geneVar <- apply(datExpr0, 2, var)
datExpr <- datExpr0[, geneVar > quantile(geneVar, 0.25)]

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

############################################
##       SAMPLE CLUSTERING (QC)           ##
############################################

sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main="Sample clustering (QC)", xlab="", sub="")

############################################
##           TRAIT MATRIX                ##
############################################

traitData <- data.frame(
  Condition = ifelse(coldata$Condition == "Cancer", 1, 0)
)
rownames(traitData) <- rownames(datExpr)

############################################
##     SOFT THRESHOLD SELECTION           ##
############################################

powers <- c(1:10, seq(12,20,2))
sft <- pickSoftThreshold(datExpr, powerVector=powers,
                         networkType="signed", verbose=5)

par(mfrow=c(1,2))
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Power", ylab="Scale Free R²", main="Scale-free topology")
abline(h=0.8, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Power", ylab="Mean connectivity")

softPower <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)

############################################
##         NETWORK CONSTRUCTION           ##
############################################

cor <- WGCNA::cor

net <- blockwiseModules(
  datExpr,
  power = softPower,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  verbose = 3
)

moduleColors <- labels2colors(net$colors)

plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors,
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03
)

############################################
##       MODULE EIGENGENES (ME)           ##
############################################

MEs <- orderMEs(net$MEs)

moduleTraitCor <- cor(MEs, traitData$Condition, use="p")
moduleTraitP   <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

moduleSummary <- data.frame(
  Module = rownames(moduleTraitCor),
  Correlation = moduleTraitCor[,1],
  Pvalue = moduleTraitP[,1]
)

sigModules <- subset(moduleSummary,
                     abs(Correlation) > 0.3 & Pvalue < 0.05)

############################################
##         GENE MEMBERSHIP (GM)           ##
############################################

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use="p"))
geneModuleMembershipP <- as.data.frame(
  corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr))
)

colnames(geneModuleMembership)  <- paste0("kME_", gsub("ME","",colnames(MEs)))
colnames(geneModuleMembershipP) <- paste0("p_kME_", gsub("ME","",colnames(MEs)))

geneInfo <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors
)

geneInfo$kME <- geneInfo$p_kME <- NA

for (i in 1:nrow(geneInfo)) {
  mod <- geneInfo$Module[i]
  geneInfo$kME[i]   <- geneModuleMembership[i, paste0("kME_",mod)]
  geneInfo$p_kME[i] <- geneModuleMembershipP[i, paste0("p_kME_",mod)]
}

############################################
##         GENE SIGNIFICANCE (GS)         ##
############################################

GS <- as.data.frame(cor(datExpr, traitData$Condition, use="p"))
GS_p <- as.data.frame(corPvalueStudent(as.matrix(GS), nrow(datExpr)))

colnames(GS) <- "GS.Condition"
colnames(GS_p) <- "p.GS.Condition"

geneInfo <- cbind(geneInfo, GS, GS_p)

############################################
##            HUB GENE TABLE              ##
############################################

hubGenes <- geneInfo %>%
  filter(abs(kME) > 0.8, abs(GS.Condition) > 0.3, p_kME < 0.05) %>%
  arrange(desc(abs(kME)))

############################################
##       mTOR / ANDROGEN ANALYSIS         ##
############################################

mtor_genes <- c("MTOR","RPTOR","RICTOR","MLST8","AKT1","AKT2","TSC1","TSC2","RHEB")
androgen_genes <- c("AR","SRD5A1","SRD5A2","HSD17B3","CYP17A1","STAR")

pathwayGenes <- unique(c(mtor_genes, androgen_genes))

pathwayHubGenes <- geneInfo %>%
  filter(Gene %in% pathwayGenes, abs(kME) > 0.7, abs(GS.Condition) > 0.3)

############################################
##       GS vs GM PLOT (MODULE)           ##
############################################

targetModule <- gsub("ME","", sigModules$Module[1])

df_plot <- geneInfo %>% filter(Module == targetModule)

ggplot(df_plot, aes(abs(kME), abs(GS.Condition))) +
  geom_point(alpha=0.7) +
  theme_classic(base_size=14) +
  labs(
    title=paste("GS vs GM — Module", targetModule),
    x="Module Membership |kME|",
    y="Gene Significance |GS|"
  )

############################################
##     MODULE EIGENGENE VS CONDITION      ##
############################################

eigengene <- MEs[, paste0("ME",targetModule)]

df_me <- data.frame(
  Eigengene = eigengene,
  Condition = coldata$Condition
)

ggplot(df_me, aes(Condition, Eigengene, fill=Condition)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.2) +
  theme_classic(base_size=14) +
  labs(title=paste("Module",targetModule,"Eigengene"))

############################################
##             SAVE OUTPUTS               ##
############################################

write.csv(geneInfo, "WGCNA_gene_membership_GS.csv", row.names=FALSE)
write.csv(hubGenes, "WGCNA_hub_genes.csv", row.names=FALSE)
write.csv(pathwayHubGenes, "Pathway_hub_genes.csv", row.names=FALSE)
write.csv(moduleSummary, "Module_trait_correlations.csv", row.names=FALSE)
