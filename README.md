# TensorPLS

An R package for exploring patterns in three-way omics data using N-PLS-DA.

---

## Why TensorPLS?

Multi-omics and longitudinal studies often generate data with a natural **multi-way structure**  
(e.g., *subjects × features × time*).  


---

# 1) Core questions addressed by TensorPLS

When applying N-PLS-DA with TensorPLS, you can answer:

1. **Is the model able to discriminate groups?**  

2. **Which variables/features mostly contribute to group discrimination?**  

3. **Is the model able to predict outside the sample (Q²)?**  

4. **How much variation is captured (R²)?**  

5. **Which time points or blocks contribute most?**  

---

# 2) Algorithms and data structures used

### Tensors
TensorPLS natively handles **3D arrays** (tensors) that encode multiple modes  
(e.g., *samples × variables × time/blocks*).  

- The **3D structure** is preserved during preprocessing and **imputation of missing values**, ensuring that the multi-way information is respected when filling gaps.  
- For the **modeling step (PLS-DA)**, the tensor is **flattened into a 2D matrix** (samples × unfolded features), since standard PLS-DA operates on matrix data.  

### Partial Least Squares (PLS) & N-PLS-DA
PLS finds latent components by **maximizing the covariance between `X` (predictors) and `Y` (response)**, not just variance in `X`.  
<p align="center">
  <img src="https://github.com/alejanner/TensorPLS/blob/main/man/figures/PLSDA.png" alt="PLS-DA" width="600">
</p>
In **PLS-DA**, `Y` encodes class membership (e.g., one-hot/dummy coding).  
**N-PLS-DA** extends PLS-DA to **multi-way (tensor) `X`**, extracting components that respect the tensor modes.


3) A workflow diagram
<p align="center">
  <img src="https://github.com/alejanner/TensorPLS/blob/main/man/figures/overviewTensorPLS.png" alt="Pipeline overview" width="600">
</p>
# 4) Tutorial: How to use TensorPLS?

Here we show a minimal workflow using TensorPLS on example data, this example consider a very large dataset with 136 Subjects, 21285 features (Genes) and 5 time points. So we need to 
expect for those type of dataset a significant computational time in Steps like imputation (if needed). 
In this example we have an external metadata file with information about subjects that user can use to filter the dataset (cohort_filter parameter). I will also include some examples from the GCTOF dataset to highlight specific differences when needed, while the complete pipeline is illustrated using the Gene Expression dataset.

```r
library(TensorPLS)

# Load example data
raw_path    <- system.file("extdata", "GeneExpressionDataProcessed.csv", package = "TensorPLS")
cohort_path <- system.file("extdata", "CohortData.csv", package = "TensorPLS")

# Build tensor (Subjects × Features × Time). 
## Parameters used in the example

X <- prepare_omics(
  data          = raw_path,              # Input CSV with omics data
  id_col        = "Individual.Id",       # Subject identifier column
  time_col      = "Time.to.IA",          # Time-point column
  transpose     = "always",              # Force transposition (raw file is wide-by-feature)
  coercion_mode = "force_numeric",       # Coerce IDs and Time to numeric, features to numeric
  cohort        = cohort_path,           # Metadata file with subject annotations
  cohort_id_col = "Group.Id",            # Column in metadata matching subject IDs
  cohort_filter = "Model.or.Validation=='Model'" # Keep only the "Model" subset
)
#> [1]   136 21285     5
```

## Imputation with Tucker decomposition

Now that we have created the tensor `X`, we need to address **missing data**.  
TensorPLS uses a **Tucker decomposition** combined with **ALS (Alternating Least Squares)** for imputation.

 On large datasets this step can take hours (2–3h for the Gene Expression dataset).  
On smaller datasets (e.g. metabolomics, 136 × 514 × 5), imputation may only take a few minutes.

---

## Choosing the number of components

Before imputation, we must decide how many components to use for each tensor mode:

- **Mode 1** = Subjects  
- **Mode 2** = Features  
- **Mode 3** = Time  

This decision is made using **Pareto (elbow) plots** and **heatmaps**.

---

### Example: Gene Expression dataset

Pareto plot:  
<p align="center">
  <img src="https://github.com/alejanner/TensorPLS/blob/main/man/figures/paretoGE.png" alt="Elbow Gene Expression Dataset" width="600">
</p>

We see a good tradeoff at **11 components**.  
From the heatmap below, we select **(4, 4, 3)**.  

<p align="center">
  <img src="https://github.com/alejanner/TensorPLS/blob/main/man/figures/heatmapComponents.png" alt="Heatmap Gene Expression" width="600">
</p>

---

### At this step, I also provide a plot from the GCTOF dataset (used for testing), where a more distinct elbow is observed at 10 components.

Pareto plot:  
<p align="center">
  <img src="https://github.com/alejanner/TensorPLS/blob/main/man/figures/paretoGCTOFX.png" alt="Elbow GCTOF dataset" width="600">
</p>



## Running the imputation

So now we are ready to perform the **imputation step**!  
We need to pass to the imputation function the number of components to be used for the Tucker decomposition — a vital step in the imputation method.

```r
# Gene Expression dataset
fullarrayGeneExpression <- ImputemethodPackage(
  X        = X,
  fac      = c(4, 4, 3),
  conver   = 1e-07,
  max.iter = 1000
)
```
