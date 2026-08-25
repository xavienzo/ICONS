# ICONS 0.2.0

A rewrite. The statistical content is now the published estimator rather than
an approximation of it, the hot paths are compiled, and nothing in the fitting
pipeline allocates a `p` by `p` matrix.

## Correctness

* **`scfa()` now returns the estimators from the paper.** ICONS 0.1.x estimated
  the factor covariance by `cov(F_hat)`. Theorem 3 of Yang, Ma, Bi and Chen
  (2024) gives `cov(F_hat) = Sigma_f + diag(a_kk / p_k)` exactly, so that
  estimator overstates every factor variance, badly for small communities.
  `scfa()` now uses the closed-form maximum likelihood estimators of their
  equation (4), which are uniformly minimum-variance unbiased.

* **Factor scores are the OLS/GLS/FGLS solution.** 0.1.x reweighted the scores
  by a per-variable residual variance. Under the model's `Sigma_u` the three
  estimators coincide (Theorem 2) and reweighting moves the estimator off the
  optimum. Reweighting is still available, and is now internally consistent,
  via `sigma_u = "diagonal"`.

* **Fixed `greedy_peeling()` returning an invalid node list.** When exactly one
  node was peeled, `Node_Seq` was set to the whole recording matrix, so `Clist`
  came back with `N + 1` entries, one of them duplicated, and the duplicated
  node was assigned to two communities.

* **Fixed `dense()` leaving a variable unassigned.** The loop condition stopped
  at `ncol(W) - 1`, so `Clist` could omit the last node. `icons_detect()` always
  returns a permutation of `seq_len(ncol(W))`.

* **Fixed `dense()` allocating a `p^2`-length vector.** `orig_Node <- 1:length(W)`
  used `length()` on a matrix. Harmless at `p = 200`, 800 MB at `p = 10000`.

* **Fixed `param_tuning_sigmau()` aborting the whole grid.** Any parameter pair
  that produced a single community raised `invalid 'times' argument` inside
  `scfa()` and killed the search. `icons_tune()` records such cells as `NA` and
  continues.

* **Fixed `plotMatrix()` leaving the graphics device in a two-panel layout.**

* **Fixed unevenly spaced, non-round axis labels on downsampled heat maps.**
  Tick positions were chosen by `pretty()` on the *downsampled* grid and the
  labels then multiplied by the scale factor, so a 4199-variable matrix was
  labelled 1050, 2100, 3149, 4199. `plot_matrix()` now picks round ticks in the
  original index range and maps them onto the drawn grid, and downsamples by a
  whole-number block size so that mapping is exact. Colour-bar ticks likewise
  step up in density rather than leaving the bar half unlabelled.

* Singletons are now tracked explicitly instead of being smuggled in as a
  trailing pseudo-community that `remove.singletons` had to strip.

## Inference

* Exact standard errors from Theorem 3 for every covariance parameter, exposed
  through `vcov()` and `confint()`. Verified by simulation: 95% Wald intervals
  cover at 95%.
* `cov_scores` gives the exact covariance of a factor score vector.

## Speed

Measured on this package's benchmark script, against 0.1.9:

| operation | p = 1000 | p = 4000 |
| --- | --- | --- |
| community detection | ~9x | ~14x |
| `scfa()` | ~180x | ~800x |
| factor-count path (K = 21) | ~300x | -- |
| tuning grid | ~5x | -- |

The `scfa()` ratio keeps growing with `p`, because the old implementation was
quadratic in `p` and the new one is linear.

* Greedy peeling moved to C++ with a compressed sparse row graph and a
  lazy-deletion heap; the adaptive loop stays in C++ and never copies a
  shrinking submatrix.
* `scfa()` is `O(np + nK^2)` instead of `O(np^2)`. The estimators depend on the
  sample covariance only through `sum(S_kk')` and `tr(S_kk)`, both of which
  follow from one pass over the data.
* Fit criteria use the `n` by `n` Gram matrix, `O(n^2 p)`, instead of forming
  `S`.
* `n_factors()` computes the whole path from one pass plus `O(K^3)`, rather
  than refitting at every `k`.
* `scfa()` now runs at `p = 50000` in a few seconds; 0.1.x needed a 20 GB
  covariance matrix to get there.

## API

Function names are now consistent and the results are S3 objects with `print`,
`summary` and `plot` methods. Old names still work and warn once per session.

| 0.1.x | 0.2.0 |
| --- | --- |
| `dense()` | `icons_detect()` |
| `param_tuning_sigmau()` | `icons_tune()` |
| `k.elbow()` | `n_factors()` |
| `plotMatrix()` | `plot_matrix()` |
| `get_membership()` | `as_membership()` |
| `get_index()` | `block_index()` |
| `get_vectorform()` | `half_vec()` |

`scfa()` keeps its name but takes a partition object or a membership vector
instead of the `CID`/`Clist` pair.

New: `reorder_matrix()`, `factor_loadings()`, `sigma_u()`, `scfa_criterion()`,
and `coef()`, `vcov()`, `confint()`, `fitted()`, `residuals()`, `predict()`,
`nobs()` methods for `"scfa"`.

`icons_detect()` gains `probs` for specifying the threshold as a quantile,
`min_size`, and `max_k`. `icons_tune()` reports the whole criterion surface and
plots it.

## Other

* Compiled code via Rcpp; no Fortran or BLAS dependency beyond base R.
* Test suite (`testthat`) checking the compiled peeling against a literal R
  transcription of the algorithm, the estimators against equation (4) computed
  from the full `p` by `p` covariance, the fast criteria against explicit
  `p` by `p` computations, and unbiasedness and interval coverage by simulation.
* `criterion = "frobenius"`, the new default, is `||S - Sigma_hat||_F`.
  `criterion = "legacy"` reproduces the 0.1.x objective exactly, for
  reproducing earlier results.

# ICONS 0.1.9

* `plotMatrix` now supports rectangle matrices.

# ICONS 0.1.8

* Fixed the bug when there is only one community detected and the closed-form estimation runs into errors.

# ICONS 0.1.7

* Substitute the inversion of matrices in the factor score estimation process (the `scfa` function) with closed-form solutions.

# ICONS 0.1.6

* Added `get_membership` utility function.

# ICONS 0.1.5

* Fixed the factor score bug: `F_HAT` to `F_HAT_FINAL`.

# ICONS 0.1.4

* `plotMatrix` now supports `format = "pdf"` export.
* `plotMatrix` now supports manual colorbar settings through `colorbar.range = NULL`.
* Added `remove.singletons = TRUE` argument in `scfa`.

# ICONS 0.1.3

* Now requires SCALED data input.
* Added `k.elbow.R`.
* Added epsilon argument (to fix unsolvable matrix issue) and  more useful return values in the scfa function for flexibility.
* The parameter tuning function now returns the sigmau results for all prespecified parameter combinations.

# ICONS 0.1.2

* Introduced the option of parallel computing in parameter tuning
* Added option of moment estimator vs. MLE estimator of sigma_u in SCFA
* Removed `entropy.R` and merged SCFA with sigma_u estimation functions

# ICONS 0.1.1

* New version of SCFA: now k-1 factor scores are returned, with the singleton set removed.
* New version of entropy estimation: sigma u now is estimated from cov(X-FL^T), instead of covariance matrix subtraction.
* New version of parameter tuning: the objective function now considers a penalty term on diagnal terms.
* Function `get_sigmau` was separated from utility functions and now stored in `entropy.R`

# ICONS 0.1.0.9006

* Changed package name from ICON to ICONS.
* Changed `plotMatrix` function to a brand new version that has prettier, more flexible axis labeling and visualizes matices in its intrinsic order.
* Changed parameter tuning objective function from `sigmau_norm + log(length(CID_temp))` to `sigmau_norm`.

# ICON 0.1.0.9005

* Changed package name from SCFA to ICON.
* `plot.R` now dynamically shows axis labels and mimic matlab heatmap presentation.

# SCFA 0.1.0.9004

* Added `scfa.R`.

# SCFA 0.1.0.9002

* Added a `NEWS.md` file to track changes to the package.
