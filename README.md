# SCFA: Semi-confirmatory factor analysis and subnetwork extraction tool

## Installation
You can install the developer version `scfa` package directly from GitHub using the `devtools` package:

```r
# Install devtools if you haven't already
install.packages("devtools")
# Install SCFA from GitHub
devtools::install_github("xavienzo/SCFA")
```

Major versions will be uploaded to CRAN. Updates will be annouced here when available.

## Usage

```r
# A simulated dataset named `sim` with 100 variables and 100 observations
data(sim)
data <- sim
matrix <- cor(data)
```

### Subnetwork extraction
```r
results <- dense(matrix)
# Number of nodes in each dense subnetwork
results$CID
# Index of reordered nodes
results$Clist
# Reordered correlation matrix
results$W_dense
```
### Parameter tuning
```r
# specify a vector of cutout thresholds and a vector of lambdas for grid search
prctile_vec = seq(94, 99, by = 0.5)
lam_vec = seq(0.4, 0.8, length.out = 5)
# grid search
param <- param_tuning_sigmau(matrix, data, prctile_vec, lam_vec)
# use optimal parameters to extract subnetworks
results <- dense(W_original = matrix, threshold = param$cut_out, lambda = param$lambda_out)
# Number of nodes in each dense subnetwork
results$CID
# Index of reordered nodes
results$Clist
# Reordered correlation matrix
results$W_dense
```

### Visualization
```r
# Heatmap of the original correlation structure
plotMatrix(matrix)
# Heatmap of the reordered matrix showing subnetwork structures
plotMatrix(results$W_dense)
```
<div style="display: flex; justify-content: space-between; align-items: center;">
  <div style="text-align: center; width: 48%;">
    <p><strong>Original</strong></p>
    <img src="https://github.com/user-attachments/assets/5da11b49-108e-4965-993d-75f83688281a" alt="sim" style="width: 48%;"/>
  </div>
  <div style="text-align: center; width: 48%;">
    <p><strong>After subnetwork extraction</strong></p>
    <img src="https://github.com/user-attachments/assets/7f8a0133-ba30-4ea5-9cba-03ad104cec84" alt="sim_dense" style="width: 48%;"/>
  </div>
</div>

### SCFA
```r
# perform SCFA
fa <- scfa(data, results$CID, results$Clist)
# factor loadings
fa$loading
# factor scores
fa$factorscore
```

## License
This package is licensed under the MIT License. See the LICENSE file for more details.

## Contact
For questions or comments, please contact `ypan@som.umaryland.edu`

## Reference

1. Wu Q, Huang X, Culbreth AJ, Waltz JA, Hong LE, Chen S. Extracting brain disease-related connectome subgraphs by adaptive dense subgraph discovery. Biometrics. 2022 Dec;78(4):1566-1578. doi: 10.1111/biom.13537. Epub 2021 Aug 22. PMID: 34374075; PMCID: PMC10396394.
2. Shuo Chen, Yuan Zhang, Qiong Wu, Chuan Bi, Peter Kochunov, L Elliot Hong, Identifying covariate-related subnetworks for whole-brain connectome analysis, Biostatistics, 2023;, kxad007, https://doi.org/10.1093/biostatistics/kxad007
3. Yang Y, Ma T, Bi C, Chen S. Semi-confirmatory factor analysis for high-dimensional data with interconnected community structures. arXiv [statME]. 2024. http://arxiv.org/abs/2401.00624.
