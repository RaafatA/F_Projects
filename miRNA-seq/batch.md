Batch effects in RNA-Seq experiments can be addressed in two main ways within a DESeq2 workflow:

1. **Model-based correction**  
   Include **batch** as a covariate in the DESeq2 design formula (e.g. `~ batch + condition`).  
   This allows DESeq2 to account for batch during differential testing without altering the underlying counts.

2. **Data-driven correction**  
   Apply a batch-correction algorithm (e.g. **ComBat-seq** or `limma::removeBatchEffect`) **after** normalization (rlog/vst) but **before** downstream exploratory analyses like PCA and clustering.

---

### Ideal insertion points in your script

- **Option A (recommended)**  
  1. Add a `batch` column to **colData**.  
  2. Change:
     ```r
     design = ~ condition
     ```
     to
     ```r
     design = ~ batch + condition
     ```
     in `DESeqDataSetFromMatrix()`.

- **Option B (for visualization only)**  
  After the rlog transformation and before PCA/heatmap:
  ```r
  rld <- rlog(dds, blind = FALSE)
  library(limma)
  assay(rld) <- removeBatchEffect(
    assay(rld),
    batch = colData(rld)$batch
  )
  # Then proceed with plotPCA() or pheatmap()
```
