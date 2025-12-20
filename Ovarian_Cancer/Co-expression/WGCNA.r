# Load Packages
library(WGCNA)
library(pheatmap)
library(ggplot2)
# Correlation statistics
library(Hmisc)

# Important: enable WGCNA multi-threading
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

## Expression Data used is not raw count data, but normalized ##
# expr <- read.csv("normalized_expression.csv", row.names = 1, check.names = FALSE)

coldata <- data.frame(
  sample_id = colnames(final_xena_identifier_matched)
)

coldata$condition <- ifelse(
  grepl("^TCGA", coldata$sample_id),
  "tumor",
  ifelse(
    grepl("^GTEX|^GTEx|^gtex", coldata$sample_id),
    "normal",
    NA
  )
)

rownames(coldata) <- coldata$sample_id

coldata$condition <- factor(
  coldata$condition,
  levels = c("normal", "tumor")
)


# dds <- DESeqDataSetFromMatrix(
#  countData = counts,
#  colData = colData,
#  design = ~ Condition)
#
# dds <- dds[rowSums(counts(dds)) > 10, ]
## Run DESeq2 normalization
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
norm_expr <- assay(vsd)

geneVariance <- apply(norm_expr, 1, var) # Calculate variance for each gene
# Remove genes with very low variance (bottom 25%)
datExpr <- datExpr[, geneVariance > quantile(geneVariance, 0.25)]

# Transpose for WGCNA (samples as rows, genes as columns)
expr_matrix <- t(datExpr)

# Check for Missing Data and Clean Up
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# # Check sample clustering to identify outliers
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering", xlab = "", sub = "")


traitData <- as.data.frame(colData(dds))

# Convert Condition to numeric
traitData$Condition <- ifelse(
  traitData$Condition == "Cancer",
  1,
  0
)
# Set row names to match expression matrix
rownames(traitData) <- rownames(datExpr)
all(rownames(traitData) == rownames(datExpr))



# Extract sample names
sampleNames <- rownames(co_exp_df)

# Create binary disease status
traitData <- data.frame(
  Condition = ifelse(
    grepl("^TCGA", sampleNames),  
    1,                            
    0                             
  )
)

rownames(traitData) <- sampleNames

# Inspect
table(traitData$Condition)








################################################################################  WGCNA  ##################################################################
# Calculate network topology for various soft-thresholding powers
powers <- c(seq(1,10,1), seq(12,20,2))

# This function tests different power values
sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)

# Plot scale-free topology fit index as a function of power
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit R²",
     main="Scale-free topology fit")
abline(h = 0.8, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity",
     main="Mean connectivity")

# Choose the power: first value that reaches plateau in scale-free topology
# and gives reasonable mean connectivity
soft_power <- sft$powerEstimate
if (is.na(soft_power)) {
  soft_power <- 6  # Use 6 as default if automatic selection fails
  cat("Warning: Automatic power selection failed. Using default power =", soft_power, "\n")
} else {
  cat("Selected soft-thresholding power:", soft_power, "\n")
}



# CRITICAL: Fix for WGCNA cor() function bug
# WGCNA has its own cor() function that conflicts with stats::cor()
cor <- WGCNA::cor
 
# One-step network construction and module detection
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

# Plot the dendrogram with module colors
plotDendroAndColors(
  net$dendrograms[[1]], moduleColors,
  "Module colors",
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene Dendrogram and Module Colors"
)

mtor_genes <- c("MTOR","RPTOR","RICTOR","MLST8","AKT1","AKT2","TSC1","TSC2","RHEB","PIK3CA","PIK3CB")
testo_genes <- c("AR","SRD5A1","SRD5A2","HSD17B3","HSD17B2","CYP17A1","CYP11A1","STAR","NR3C1")




MEs <- orderMEs(net$MEs)

moduleTraitCor <- cor(MEs, traitData$Condition, use = "p")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Decide WHICH MODULES TO SHOW
moduleSummary <- data.frame(
  Module = rownames(moduleTraitCor),
  Correlation = moduleTraitCor[,1],
  Pvalue = moduleTraitP[,1]
)

# Significant Cancer-associated modules
sigModules <- moduleSummary |>
  subset(Pvalue < 0.05 & abs(Correlation) > 0.3)

sigModules


# Identify modules enriched for mTOR / androgen genes
geneModuleMap <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors
)

pathwayGenes <- unique(c(mtor_genes, androgen_genes))

pathwayModules <- geneModuleMap |>
  subset(Gene %in% pathwayGenes & Module %in% gsub("ME","", sigModules$Module))

unique(pathwayModules$Module)

######################################################## Visualization ########################################################################################



#-----------------------------------------------
# Pathway gene co-expression heatmap
#-----------------------------------------------
# Heatmap showing pairwise co-expression between mTOR and androgen signaling genes within a disease-associated WGCNA module. 
# Colors represent Pearson correlation coefficients across all samples. Genes cluster by pathway, indicating coordinated regulation.

library(pheatmap)

targetModule <- pathwayModules$Module[1]

genes_in_module <- geneModuleMap |>
  subset(Module == targetModule & Gene %in% pathwayGenes) |>
  pull(Gene)

expr_sub <- datExpr[, genes_in_module]

cor_mat <- cor(expr_sub, method = "pearson")

annotation <- data.frame(
  Pathway = ifelse(colnames(cor_mat) %in% mtor_genes, "mTOR", "Androgen")
)
rownames(annotation) <- colnames(cor_mat)

pheatmap(
  cor_mat,
  annotation_row = annotation,
  annotation_col = annotation,
  color = colorRampPalette(c("blue","white","red"))(100),
  border_color = NA,
  main = paste("Co-expression of mTOR and Androgen genes\nModule:", targetModule)
)


#-----------------------------------------------
# Module-restricted network
#-----------------------------------------------

# Network visualization of mTOR and androgen pathway genes within the selected WGCNA module. 
# Nodes represent genes and edges indicate strong co-expression relationships. 
# Node color reflects pathway membership.

library(igraph)

threshold <- 0.6
edges <- which(abs(cor_mat) > threshold & upper.tri(cor_mat), arr.ind = TRUE)

edge_list <- data.frame(
  from = rownames(cor_mat)[edges[,1]],
  to   = colnames(cor_mat)[edges[,2]],
  weight = cor_mat[edges]
)

g <- graph_from_data_frame(edge_list, directed = FALSE)

V(g)$pathway <- ifelse(V(g)$name %in% mtor_genes, "mTOR", "Androgen")

plot(
  g,
  vertex.color = ifelse(V(g)$pathway == "mTOR","firebrick","steelblue"),
  vertex.size = 22,
  vertex.label.cex = 0.8,
  edge.width = abs(E(g)$weight) * 3,
  main = "Module-restricted co-expression network"
)


#-----------------------------------------------
# Eigengene vs Condition
#-----------------------------------------------

# Comparison of module eigengene expression between normal (GTEx) and cancer (TCGA) samples. 
# The significant shift in eigengene values indicates disease-associated dysregulation of this co-expression module.


library(ggplot2)

eigengene <- MEs[, paste0("ME", targetModule)]

df <- data.frame(
  Eigengene = eigengene,
  Condition = factor(traitData$Condition, labels = c("Normal","Cancer"))
)

ggplot(df, aes(Condition, Eigengene, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5) +
  theme_classic(base_size = 14) +
  labs(
    title = paste("Module", targetModule, "eigengene vs condition"),
    y = "Module eigengene expression"
  )




















# Map Pathway Genes to Modules
geneModuleMap <- data.frame(Gene=colnames(datExpr), Module=moduleColors)
mtor_mods <- geneModuleMap %>% filter(Gene %in% mtor_genes)
testo_mods <- geneModuleMap %>% filter(Gene %in% testo_genes)


# Module Eigengene Correlation (Pathway-Level)
MEs <- orderMEs(net$MEs)
modCorMatrix <- cor(MEs, use="p")
modCorMatrix[
  paste0("ME", unique(mtor_mods$Module)),
  paste0("ME", unique(testo_mods$Module))
]
# Gene-Level Cross-Pathway Coexpression
expr_mtor <- datExpr[, colnames(datExpr) %in% mtor_genes]
expr_testo <- datExpr[, colnames(datExpr) %in% testo_genes]
cross_cor <- cor(expr_mtor, expr_testo, method="pearson")
pheatmap(cross_cor, main="mTOR–Testosterone Correlations")

#Statistical Significance of Cross-Correlations
cor_test <- rcorr(cbind(expr_mtor, expr_testo), type="pearson")
pval_matrix <- cor_test$P[colnames(expr_mtor), colnames(expr_testo)]
sig_pairs <- which(abs(cross_cor) > 0.6 & pval_matrix < 0.05, arr.ind=TRUE)


write.csv(cross_cor, "cross_correlation_matrix.csv")
write.csv(geneModuleMap, "gene_module_map.csv")

