# ICONS

<!-- badges: start -->
[![R-CMD-check](https://github.com/xavienzo/ICONS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/xavienzo/ICONS/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**I**ntegrative analysis of **CO**variance matrices and **N**etwork **S**tructure.

High-dimensional biomedical data — genomics, proteomics, metabolomics,
neuroimaging — tend to have covariance matrices with *interconnected community
structure*: groups of mutually correlated features, the groups themselves
correlated, plus features that belong to no group. ICONS finds that structure
and uses it to fit a factor model.

Confirmatory factor analysis needs you to know in advance which variables load
on which factor, and does not scale past a few hundred variables. ICONS solves
both problems: the communities are learned from the data, and every parameter
has a closed-form maximum likelihood estimator, so there is no optimiser to
converge.

## Installation

```r
# install.packages("devtools")
devtools::install_github("xavienzo/ICONS")
```

## Usage

```r
library(ICONS)
data(sim)                       # 100 observations, 200 variables

W <- cor(sim)
```

### 1. Detect communities

```r
part <- icons_detect(W, threshold = 0.6)
part
#> <ICONS partition>
#>   variables   : 200
#>   communities : 9
#>   singletons  : 1
#>   threshold   : 0.6   lambda: 0.5
#>   sizes       : 20 20 20 15 15 15 15 15 64
```

Give the threshold as a quantile of the edge weights instead, if that is easier
to reason about:

```r
icons_detect(W, probs = 0.95)
```

### 2. Tune, if you want the data to pick the settings

```r
tuned <- icons_tune(W, sim, probs = seq(0.90, 0.99, 0.01),
                    lambda = c(0.4, 0.5, 0.6, 0.7, 0.8))
tuned
plot(tuned)                     # the whole criterion surface

part <- icons_detect(W, threshold = tuned$threshold, lambda = tuned$lambda)
```

The threshold matters far more than `lambda`. A flat criterion surface means
the choice is not well identified — worth knowing before you report it.

### 3. Choose the number of factors

```r
nf <- n_factors(sim, part)
nf$suggested
plot(nf)
```

### 4. Fit the factor model

```r
fit <- scfa(sim, part)
fit
#> <Semi-confirmatory factor analysis>
#>   observations : 100
#>   variables    : 200 (199 modelled, 1 singleton)
#>   factors      : 9
#>   Sigma_u      : mle
#>   frobenius criterion: 16.25  (relative 0.3365)

fit$scores                      # n x K factor scores
fit$Sigma_f                     # K x K factor covariance
fit$a                           # error variance per community

summary(fit)                    # estimates with Wald intervals
confint(fit)                    # exact, from Theorem 3
```

Standard extractors work as expected: `coef()`, `vcov()`, `confint()`,
`fitted()`, `residuals()`, `predict()`, `nobs()`, plus `factor_loadings()` and
`sigma_u()`.

### 5. Look at it

```r
plot_matrix(W, main = "Original")
plot_matrix(reorder_matrix(W, part), partition = part, main = "Reordered")
```

## What changed in 0.2.0

0.2.0 is a rewrite. Three things are worth calling out.

**The estimator is now the published one.** ICONS 0.1.x estimated the factor
covariance as `cov(F_hat)`. Theorem 3 of Yang et al. (2024) shows that
`cov(F_hat) = Sigma_f + diag(a_kk / p_k)` exactly, so that overstates every
factor variance by `a_kk / p_k` — 12.5% for a community of 8 with
`a = 0.5, b = 0.5`, and worse for small or noisy communities. 0.2.0 uses the
closed-form UMVUEs from their equation (4), and adds the exact standard errors
the paper derives, which 0.1.x did not expose at all.

**Nothing forms a `p` by `p` matrix.** The estimators depend on the sample
covariance only through `sum(S_kk')` and `tr(S_kk)`, and both fall out of a
single pass over the data. `scfa()` is `O(np + nK^2)` rather than `O(np^2)`;
the fit criteria use the `n` by `n` Gram matrix. At `p = 50000` a fit takes a
few seconds, where 0.1.x would have needed a 20 GB covariance matrix.

**Several real bugs are fixed** — including `greedy_peeling()` returning a node
list with a duplicate, `dense()` leaving one variable unassigned, and
`param_tuning_sigmau()` aborting the entire grid search whenever one cell
produced a single community. See [NEWS.md](NEWS.md).

Speedups against 0.1.9, at p = 1000–4000: roughly 9–14x for detection,
180–800x for `scfa()`, and 300x for the factor-count path.

Old function names still work and warn once:

| 0.1.x | 0.2.0 |
| --- | --- |
| `dense()` | `icons_detect()` |
| `param_tuning_sigmau()` | `icons_tune()` |
| `k.elbow()` | `n_factors()` |
| `plotMatrix()` | `plot_matrix()` |
| `get_membership()` | `as_membership()` |
| `get_index()` | `block_index()` |
| `get_vectorform()` | `half_vec()` |

## References

1. Yang, Y., Ma, T., Bi, C., & Chen, S. (2024). Semi-confirmatory factor
   analysis for high-dimensional data with interconnected community structures.
   *arXiv:2401.00624*.
2. Yang, Y., Chen, C., & Chen, S. (2024). Covariance matrix estimation for
   high-throughput biomedical data with interconnected communities.
   *The American Statistician*, 78(4), 401–411.
3. Chen, S., Zhang, Y., Wu, Q., Bi, C., Kochunov, P., & Hong, L. E. (2024).
   Identifying covariate-related subnetworks for whole-brain connectome
   analysis. *Biostatistics*, 25(2), 541–558.
4. Wu, Q., Huang, X., Culbreth, A. J., Waltz, J. A., Hong, L. E., & Chen, S.
   (2022). Extracting brain disease-related connectome subgraphs by adaptive
   dense subgraph discovery. *Biometrics*, 78(4), 1566–1578.

## License

MIT. See [LICENSE](LICENSE).

## Contact

`ypan@som.umaryland.edu`
