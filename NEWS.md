# ICONS 0.1.6
* Added `get_membership` utility function.

# ICONS 0.1.5
* Fixed the factor score bug: F_HAT to F_HAT_FINAL.

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
* Introduced the option of parallel computing in parameter tuning🎉
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
