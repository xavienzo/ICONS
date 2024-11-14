# ICONS 0.1.2

# ICONS

# ICONS 0.1.2
* Introduced option of parallel computing in parameter tuning
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
