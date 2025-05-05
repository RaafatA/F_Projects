
Batch effects in RNA-Seq experiments can be addressed in two main ways within a DESeq2 workflow:

1. **Model-based correction**: Include `batch` as a covariate in the DESeq2 design formula (e.g. `~ batch + condition`).  This allows DESeq2 to **account for** batch during differential testing without altering the underlying counts (###Model-based correction).
2. **Data-driven correction**: Apply a batch-correction algorithm (e.g. **ComBat-seq** or `limma::removeBatchEffect`) **after** normalizing (rlog/vst) but **before** downstream exploratory analyses like PCA and clustering (###Data-driven correction).

In your script, the ideal insertion points are:

* **Option A (recommended)**: Add a `batch` column to **colData** and change

  ```r
  design = ~ condition
  ```

  to

  ```r
  design = ~ batch + condition
  ```

  in the call to `DESeqDataSetFromMatrix`.
* **Option B (for visualization only)**: After the rlog transformation (`rld <- rlog(dds, blind=FALSE)`) and **before** `plotPCA()` or `pheatmap()`, call

  ```r
  assay(rld) <- limma::removeBatchEffect(assay(rld), batch = colData(rld)$batch)
  ```
* **Option C (raw count correction)**: Run **ComBat-seq** on your raw `count_matrix` **before** creating the DESeq2 object, then proceed with DESeq2 on the batch-adjusted counts.

---

## 1. Model-based correction (DESeq2 design)

### When?

– **Always** for differential expression analysis when batch and condition are not completely confounded.
– You do **not** modify counts; instead the statistical model adjusts for batch during testing.

### How?

1. Add a `batch` column to your `colData` data.frame:

   ```r
   colData <- data.frame(
     row.names = sample_names,
     batch     = factor(c("B1","B1","B2","B2")),     # your actual batches
     condition = factor(c("cancer","cancer","normal","normal"))
   )
   ```
2. Create the DESeq2 dataset with batch in the design:

   ```r
   dds <- DESeqDataSetFromMatrix(
     countData = count_matrix,
     colData   = colData,
     design    = ~ batch + condition
   )
   ```
3. Run the standard DESeq2 pipeline:

   ```r
   dds <- estimateSizeFactors(dds)
   dds <- dds[rowSums(counts(dds)>0) >= 2, ]
   dds <- DESeq(dds)
   res <- results(dds, contrast=c("condition","normal","cancer"))
   ```

This approach is recommended in the DESeq2 vignette and by the community for properly accounting for known batch effects.

---

## 2. Data-driven correction (removeBatchEffect or ComBat)

### A. Using `limma::removeBatchEffect` on transformed data

* **Where**: Immediately **after**

  ```r
  rld <- rlog(dds, blind = FALSE)
  ```

  and **before** PCA or heatmap steps.
* **Why**: To produce batch-corrected expression values for **exploratory** plots (PCA/MDS/heatmap), since DESeq2 transformations do not remove batch.
* **Code**:

  ```r
  library(limma)
  # Correct the rlog matrix
  mat_rlog <- assay(rld)
  mat_corr <- removeBatchEffect(mat_rlog, batch=colData(rld)$batch)
  assay(rld) <- mat_corr
  # Then plot:
  pcaData <- plotPCA(rld, intgroup="condition", returnData=TRUE)
  ```

### B. Using **ComBat-seq** on raw counts

* **Where**: **Before** constructing the DESeq2 object
  i.e. right after you build `count_matrix`, but **before** `DESeqDataSetFromMatrix()`.
* **Why**: If you prefer to remove batch completely from the raw counts (e.g. for combined analyses across batches) and still maintain integer counts for DESeq2.
* **Code**:

  ```r
  library(sva)
  combat_counts <- ComBat_seq(counts = count_matrix,
                              batch  = colData$batch)
  dds <- DESeqDataSetFromMatrix(combat_counts, colData, design=~condition)
  ```

  Note: If you use ComBat-seq, you typically **omit** `batch` from the design, tested via `~ condition` alone.

---

## Recommended 

1. **After** you create `count_matrix` (end of step 8) and define **colData**, **insert or update**:

   ```r
   colData$batch <- factor(c(...))    # your batch assignments
   dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                  colData   = colData,
                                  design    = ~ batch + condition)
   ```
2. **Keep** the rest of your DESeq2 workflow unchanged (size factor estimation, filtering, `DESeq()`, `results()`).
3. For **visualizations** (PCA, heatmaps), apply `removeBatchEffect` on `rld` **before** `plotPCA()` or heatmap calls.

This ensures robust differential testing **and** accurate exploratory plots that reflect biological variation rather than technical batch effects.
