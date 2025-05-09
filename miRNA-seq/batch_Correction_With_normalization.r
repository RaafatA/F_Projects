# Load required Bioconductor packages
library(DESeq2)    # differential expression and size-factor normalization :contentReference[oaicite:4]{index=4}
library(sva)       # ComBat for batch correction :contentReference[oaicite:5]{index=5}

# Read in your GTEx CPM matrix (rows = miRNAs, cols = samples)
gtex_cpm      <- as.matrix(read.csv("GTEx_ovary_CPM.csv", row.names=1))

# Read GTEx metadata (must contain a batch column e.g. 'Center')
meta_gtex     <- read.csv("GTEx_ovary_metadata.csv", row.names=1)

# Read TCGA raw counts (rows = miRNAs, cols = samples)
tcga_counts   <- as.matrix(read.csv("TCGA_ovary_counts.csv", row.names=1))

# Read TCGA metadata (must contain a batch column e.g. 'Plate' or 'Center')
meta_tcga     <- read.csv("TCGA_ovary_metadata.csv", row.names=1)



# Create DESeqDataSet for TCGA (no design needed for normalization)
dds_tcga <- DESeqDataSetFromMatrix(countData = tcga_counts,
                                   colData   = meta_tcga,
                                   design    = ~ 1)

# Estimate size factors (median-of-ratios normalization) :contentReference[oaicite:7]{index=7}
dds_tcga <- estimateSizeFactors(dds_tcga)

# Extract CPM‐style values ("fragments per million") :contentReference[oaicite:8]{index=8}
tcga_cpm <- fpm(dds_tcga, robust=TRUE)




# Extract batch labels
batch_gtex <- meta_gtex$Center     # e.g., sequencing center for GTEx :contentReference[oaicite:9]{index=9}  
batch_tcga <- meta_tcga$Plate      # e.g., plate ID for TCGA :contentReference[oaicite:10]{index=10}

# ComBat on GTEx CPM
gtex_combat <- ComBat(dat = gtex_cpm,
                      batch = batch_gtex,
                      mod   = NULL,         # no covariates
                      par.prior = TRUE)     # empirical Bayes adjustment :contentReference[oaicite:11]{index=11}

# ComBat on TCGA CPM
tcga_combat <- ComBat(dat = tcga_cpm,
                      batch = batch_tcga,
                      mod   = NULL,
                      par.prior = TRUE)     # preserves CPM scaling :contentReference[oaicite:12]{index=12}





# Merge corrected matrices
expr_all   <- cbind(gtex_combat, tcga_combat)
batch_all  <- factor(c(rep("GTEx",  ncol(gtex_combat)),
                       rep("TCGA", ncol(tcga_combat))))

# Final ComBat across projects
expr_final <- ComBat(dat       = expr_all,
                     batch     = batch_all,
                     mod       = NULL,
                     par.prior = TRUE)   # yields homogenized matrix :contentReference[oaicite:13]{index=13}

# Write out the batch‐corrected expression matrix
write.csv(expr_final,
          file = "Ovary_GTex_TCGA_batch_corrected_CPM.csv",
          quote = FALSE)
