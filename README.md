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
In **PLS-DA**, `Y` encodes class membership (e.g., one-hot/dummy coding).  
**N-PLS-DA** extends PLS-DA to **multi-way (tensor) `X`**, extracting components that respect the tensor modes.

> Intuition: PLS searches for directions in `X` that best **co-var** with `Y`, thus focusing on **predictive structure** relevant to the outcome.
