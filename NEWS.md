# crandep 0.0.3 (2020-08-)

## New functions

- `get_dep_all_packages()` and `get_graph_all_packages()`: The former is for a data frame of all dependencies of all CRAN packages, while the latter is for the graph of one type of depenedencies of all CRAN packages.

- `get_dep()` replaces `get_dep_all()`, with the same functionality.

## Minor changes

- The argument `type` in `get_dep()` and `types` in `get_dep_df()` allows input more flexibility. For reverse dependencies, either space or underscore is accepted for separating the words e.g. `type = "reverse suggests"` or `type = "reverse_suggests"`.

- The argument `types` in `get_dep_df()` allows (as before) a character vector of dependency words. Also allowed now is `types = "all"` which means all of the four dependencies (depends, suggests, imports, linking to) and their reverse counterparts.

- Previously, there were issues with string manipulation for some packages without scraping. This is because, when using `tools::CRAN_package_db()`, there might be no space between the package name and the left parenthesis for the version. This is not an issue if `scrape = TRUE` as there is always a space on the CRAN page.

- In the output of `get_dep_df()` and `get_dep_all_packages()`, any `LinkingTo` and `Reverse linking to` dependencies will become "linking to" ("linking_to" previously) in the variable `type`, with "FALSE" and "TRUE" in the variable `reverse`, respectively. This is also updated in the data `cran_dependencies`.