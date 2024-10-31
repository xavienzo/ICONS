# ICONS: Integrative analysis of COvariance matrix and Network Structure

## Installation
You can install the development version `ICONS` package directly from GitHub using the `devtools` package:

```r
# Install devtools if you haven't already
install.packages("devtools")
# Install SCFA from GitHub
devtools::install_github("xavienzo/ICONS")
```

Major versions will be uploaded to CRAN. Updates will be annouced here when available.

## Usage

### Example dataset
```r
# A simulated dataset named `sim` with 100 variables and 100 observations
data(sim)
matrix <- cor(sim)
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
param <- param_tuning_sigmau(matrix, sim, prctile_vec, lam_vec)

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
# Plot the original correlation matrix without diagonal values
matrix_wodiag <- matrix - diag(diag(matrix)) #Remove the diagonal elements

# To save the plot as a file, uncomment the lines
plotMatrix(matrix_wodiag, 
           # save.image = T, 
           # filepath = "simulation_orig.png", 
           # format = "png",
           cex.axis = 1.3, cex.lab = 1.3)

# Plot the reordered correlation matrix showing network structures
# To save the plot as a file, uncomment the lines
plotMatrix(results$W_dense, 
           # save.image = T, 
           # filepath = "figure/simulation_dense.png", 
           # format = "png",
           cex.axis = 1.3, cex.lab = 1.3)
```
<div style="display: flex; justify-content: space-between;">
  <div style="text-align: center; width: 40%;">
    <p><strong>Original</strong></p>
    <img src="https://github.com/user-attachments/assets/4536fe79-8d64-4619-98a6-b2fa3fb4495e" alt="sim" style="width: 40%;"/>
  </div>
  <div style="text-align: center; width: 40%;">
    <p><strong>After subnetwork extraction</strong></p>
    <img src="https://github.com/user-attachments/assets/be9b0f12-ac41-4ea0-b153-3d9066d9291e" alt="sim_dense" style="width: 40%;"/>
  </div>
</div>

### SCFA (Semi-confirmatory factor analysis)
```r
# perform SCFA
fa <- scfa(sim, results$CID, results$Clist)
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
