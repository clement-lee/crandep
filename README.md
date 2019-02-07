<!-- README.md is generated from README.Rmd. Please edit that file -->
rackage
=======

The goal of rackage is to provide functions for analysing R packages.

Installation
------------

You can install rackage from github with:

    # install.packages("devtools")
    devtools::install_github("clement-lee/rackage")

Example
-------

This is a basic example which shows you how to obtain one kind of
dependencies (Imports) for one package (dplyr):

    library(rackage)

    ## Get the CRAN URL for package dplyr
    foo <- cran_url("dplyr")
    foo
    #> https://CRAN.R-project.org/package=dplyr

    ## Get the html text of the page
    bar <- html_text_vec(foo)

    ## Get the string of imported packages
    qux <- get_dep_str(bar, "Imports")
    qux
    #> [1] "assertthat (≥ 0.2.0), bindrcpp (≥ 0.2.0.9000), glue (≥ 1.1.1), magrittr (≥ 1.5), methods, pkgconfig (≥ 2.0.1), R6 (≥ 2.2.2), Rcpp (≥ 0.12.19), rlang (≥ 0.3.0), tibble (≥ 1.4.2), tidyselect (≥ 0.2.3), utils"

    ## Split the above string into a vector
    quz <- list_of_dep(qux)
    quz
    #>  [1] "assertthat" "bindrcpp"   "glue"       "magrittr"   "methods"   
    #>  [6] "pkgconfig"  "R6"         "Rcpp"       "rlang"      "tibble"    
    #> [11] "tidyselect" "utils"
