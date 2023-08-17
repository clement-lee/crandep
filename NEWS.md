# crandep 0.3.2 (2023-08-17)

## Functions
- `dupp()` and `mcmc_upp()` are replaced by `dpol()` and `mcmc_pol()`, respectively, to facilitate a generalisation of the discrete power law, namely the Zipf-polylog distribution.
- `dmix()` and `mcmc_mix()` are replaced by `dmix2()` and `mcmc_mix2()`, respectively, with a new parametrisation, for the 2-component mixture distribution.
- `dmix3()` and `mcmc_mix3()` are added for the 3-component mixture distribution.

## Dependencies
- RcppGSL is no longer needed as a LinkingTo dependency.

## Vignettes and README
- The vignettes and README are updated according to the above changes.
- Also, the previous pipe operator "%>%" is replaced by the native one "|>" throughout the vignettes and README.

# crandep 0.3.1 (2022-06-02)

## Functions

- `get_dep()`, `html_text_vec()` and `get_dep_all_packages()` now return an error message if Internet resources are not available.

# crandep 0.3.0 (2021-06-10)

## Functions

- The previous functionality of `get_dep()` is replaced by that of `get_dep_df()`, while `get_dep_df()` is soft deprecated. This means the former is the single function for obtaining dependencies in a non-igraph object.

- The argument `type` in `get_dep()` and `get_graph_all_packages()` now allows `Enhances` and `Reverse enhances` as the value. These two kind of dependencies are also included in the data frame obtained using `get_dep_all_packages()`.

# crandep 0.2.0 (2021-05-10)

## Functions

- Examples in internal functions `html_text_vec()`, `get_dep_str()` and `get_dep_vec()` removed to minimise the errors due to no internet connection and/or timeout.

- Multiple dependencies are now allowed in the `type` argument in `get_graph_all_packages()`.

- Arguments `give_log` in `dupp()` & `dmix()` changed to `log` without changing the functionality. For uniformity, `Supp()` & `Smix()` are also given the additional argument `log`.

- The ordering of arguments in the funcions `*mix()` is made consistent.

## Vignettes

- For the vignette on dependencies of all CRAN packages, community detection is added.

- Replace https://cran.r-project.org/web/packages/available_packages_by_name.html by https://cran.r-project.org, in the dependency network vignette, to prevent NOTE on possibly invalid URL.

- In the vignette on modelling the number of reverse dependencies, a section on fitting extreme value mixture distribution is added.

## Data

- Added is a citation network of the CHI conference papers, that can serve as a comparison to the CRAN dependency network, in terms of network summaries and characteristics, such as degree distribution.

- Replace https://cran.r-project.org/web/packages/available_packages_by_name.html by https://cran.r-project.org, in the manual of `cran_dependencies`, to prevent NOTE on possibly invalid URL.






# crandep 0.1.0 (2020-08-10)

## New functions

- `get_dep_all_packages()` and `get_graph_all_packages()`: The former is for a data frame of all dependencies of all CRAN packages, while the latter is for the graph of one type of depenedencies of all CRAN packages.

- `get_dep()` replaces `get_dep_all()`, with the same functionality. `get_dep()` gets soft deprecated.

## Minor changes

- The argument `type` in `get_dep()` and `types` in `get_dep_df()` allows input more flexibility. For reverse dependencies, either space or underscore is accepted for separating the words e.g. `type = "reverse suggests"` or `type = "reverse_suggests"`.

- The argument `types` in `get_dep_df()` allows (as before) a character vector of dependency words. Also allowed now is `types = "all"` which means all of the four dependencies (depends, suggests, imports, linking to) and their reverse counterparts.

- Previously, there were issues with string manipulation for some packages without scraping. This is because, when using `tools::CRAN_package_db()`, there might be no space between the package name and the left parenthesis for the version. This is not an issue if `scrape = TRUE` as there is always a space on the CRAN page.

- In the output of `get_dep_df()` and `get_dep_all_packages()`, any `LinkingTo` and `Reverse linking to` dependencies will become "linking to" ("linking_to" previously) in the variable `type`, with "FALSE" and "TRUE" in the variable `reverse`, respectively. This is also updated in the data `cran_dependencies`.

## Vignettes
- The sections on obtaining dependencies of all CRAN packages is now moved to a new vignette. In this vignette, we also provide interactive visualisation of the network of `Depends` of all packages.

- The degree modelling vignette is now for `Imports` network, not `Depends` network. In addition to discrete power law, a discrete extreme value mixture distribution is also used to model the same data set.





# crandep 0.0.2 (2020-07-18)

## New functions

- `dupp()` and `Supp()`: density and survival functions, respectively, of the discrete power law (above a threshold).
- `mcmc_upp()`: fitting the discrete power law (above a threshold) to data using Markov chain Monte Carlo (MCMC).
- `dmix()`, `Smix()`: density and survival functions, respectively, of a discrete extreme value mixture distribution.
- `mcmc_mix()`: fitting the discrete extreme value mixture distribution to data using MCMC.

## Minor changes

- Additional argument in `get_dep_all()` and `get_dep_df()`: `scrape = TRUE` is the same as previous version, while `scrape = FALSE` means `tools::CRAN_package_db()` (thanks to Dirk Eddelbuettel (#1)) will be used instead. Note that changing this argument should still give the same result; the main difference is the time taken.





