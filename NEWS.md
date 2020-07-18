# crandep 0.0.2 (2020-07-18)

## New functions

- `dupp()` and `Supp()`: density and survival functions, respectively, of the discrete power law (above a threshold).
- `mcmc_upp()`: fitting the discrete power law (above a threshold) to data using Markov chain Monte Carlo (MCMC).
- `dmix()`, `Smix()`: density and survival functions, respectively, of a discrete extreme value mixture distribution.
- `mcmc_mix()`: fitting the discrete extreme value mixture distribution to data using MCMC.

## Minor changes

- Additional argument in `get_dep_all()` and `get_dep_df()`: `scrape = TRUE` is the same as previous version, while `scrape = FALSE` means `tools::CRAN_package_db()` (thanks to Dirk Eddelbuettel (#1)) will be used instead. Note that changing this argument should still give the same result; the main difference is the time taken.