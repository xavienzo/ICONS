## Test environments

* local macOS 26.6, R 4.3.0 (aarch64)
* GitHub Actions: macOS-latest (release), windows-latest (release),
  ubuntu-latest (devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 2 notes

* "New submission" — this is the first CRAN submission of ICONS.
* "unable to verify current time" — local network restriction, not reproducible
  on the CI runners.

## Notes for the reviewer

* The package contains compiled code (Rcpp only; no BLAS, LAPACK or Fortran
  linkage).
* Parallelism in `icons_tune()` is opt-in and defaults to a single core
  (`getOption("ICONS.ncores", 1L)`); examples and tests do not use more than
  one.
* `fitted.scfa()` can allocate a p-by-p matrix and therefore refuses to do so
  above `max_p = 5000` unless the user raises the limit.
* The long-running simulation test (interval coverage) is guarded with
  `skip_on_cran()`.
