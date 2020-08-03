
<!-- README.md is generated from README.Rmd. Please edit that file -->

# crandep

<!-- badges: start -->

<!-- badges: end -->

The goal of crandep is to provide functions for analysing the
dependencies of CRAN packages using social network analysis.

## Installation

You can install crandep from github with:

``` r
# install.packages("devtools")
devtools::install_github("clement-lee/crandep")
```

``` r
library(crandep)
library(dplyr)
library(ggplot2)
library(igraph)
```

## Overview

The functions and example dataset can be divided into the following
categories:

1.  For obtaining data frames of package dependencies, use `get_dep()`,
    `get_dep_df()`, `get_dep_all_packages()`.
2.  For obtaining igraph objects of package dependencies, use
    `get_graph_all_packages()` and `df_to_graph()`.
3.  For modelling the number of dependencies, use `*upp()` and `*mix()`.
4.  There is also an example data set `cran_dependencies`.

## One kind of dependencies

To obtain the information about various kinds of dependencies of a
package, we can use the function `get_dep()` which takes the package
name and the type of dependencies as the first and second arguments,
respectively. Currently, the second argument accepts `Depends`,
`Imports`, `LinkingTo`, `Suggests`, `Reverse_depends`,
`Reverse_imports`, `Reverse_linking_to`, and `Reverse_suggests`, or any
variations in their letter cases, or if the underscore "\_" is replaced
by a space.

``` r
get_dep("dplyr", "Imports")
#>  [1] "ellipsis"   "generics"   "glue"      
#>  [4] "lifecycle"  "magrittr"   "methods"   
#>  [7] "R6"         "rlang"      "tibble"    
#> [10] "tidyselect" "utils"      "vctrs"
get_dep("MASS", "depends")
#> [1] "grDevices" "graphics"  "stats"    
#> [4] "utils"
```

We only consider the 4 most common types of dependencies in R packages,
namely `Imports`, `Depends`, `Suggests` and `LinkingTo`, and their
reverse counterparts. For more information on different types of
dependencies, see [the official
guidelines](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Package-Dependencies)
and <http://r-pkgs.had.co.nz/description.html>.

## Multiple kind of dependencies

As the information all dependencies of one package are on the same page
on CRAN, to avoid scraping the same multiple times, we can use
`get_dep_df()` instead of `get_dep()`. The output will be a data frame
instead of a character vector.

``` r
get_dep_df("dplyr", c("imports", "LinkingTo"))
#>     from         to    type reverse
#> 1  dplyr   ellipsis imports   FALSE
#> 2  dplyr   generics imports   FALSE
#> 3  dplyr       glue imports   FALSE
#> 4  dplyr  lifecycle imports   FALSE
#> 5  dplyr   magrittr imports   FALSE
#> 6  dplyr    methods imports   FALSE
#> 7  dplyr         R6 imports   FALSE
#> 8  dplyr      rlang imports   FALSE
#> 9  dplyr     tibble imports   FALSE
#> 10 dplyr tidyselect imports   FALSE
#> 11 dplyr      utils imports   FALSE
#> 12 dplyr      vctrs imports   FALSE
```

The column `type` is the type of the dependency converted to lower case.
Also, `LinkingTo` is now converted to `linking to` for consistency. For
the four reverse dependencies, the substring `"reverse_"` will not be
shown in `type`; instead the `reverse` column will be `TRUE`. This can
be illustrated by the following:

``` r
get_dep("abc", "depends")
#> [1] "abc.data" "nnet"     "quantreg" "MASS"    
#> [5] "locfit"
get_dep("abc", "reverse_depends")
#> [1] "abctools" "EasyABC"
get_dep_df("abc", c("depends", "reverse_depends"))
#>   from       to    type reverse
#> 1  abc abc.data depends   FALSE
#> 2  abc     nnet depends   FALSE
#> 3  abc quantreg depends   FALSE
#> 4  abc     MASS depends   FALSE
#> 5  abc   locfit depends   FALSE
#> 6  abc abctools depends    TRUE
#> 7  abc  EasyABC depends    TRUE
```

Theoretically, for each forward dependency

    #>   from to type reverse
    #> 1    A  B    c   FALSE

there should be an equivalent reverse dependency

    #>   from to type reverse
    #> 1    B  A    c    TRUE

Aligning the `type` in the forward dependency and the reverse dependency
enables this to be checked easily.

To obtain all 8 types of dependencies, we can use `"all"` in the second
argument, instead of typing a character vector of all 8 words:

``` r
df0.abc <- get_dep_df("abc", "all")
df0.abc
#>    from         to     type reverse
#> 1   abc   abc.data  depends   FALSE
#> 2   abc       nnet  depends   FALSE
#> 3   abc   quantreg  depends   FALSE
#> 4   abc       MASS  depends   FALSE
#> 5   abc     locfit  depends   FALSE
#> 9   abc   abctools  depends    TRUE
#> 10  abc    EasyABC  depends    TRUE
#> 11  abc ecolottery  imports    TRUE
#> 12  abc       ouxy  imports    TRUE
#> 14  abc      coala suggests    TRUE
df0.rstan <- get_dep_df("rstan", "all")
dplyr::count(df0.rstan, type, reverse) # all 8 types
#>         type reverse  n
#> 1    depends   FALSE  2
#> 2    depends    TRUE 21
#> 3    imports   FALSE 10
#> 4    imports    TRUE 63
#> 5 linking to   FALSE  5
#> 6 linking to    TRUE 50
#> 7   suggests   FALSE 12
#> 8   suggests    TRUE 18
```

As of 2020-08-03, the packages that have all 8 types of dependencies are
bigmemory, gRbase, rstan, sf, stochvol, xts.

## Building and visualising a dependency network

To build a dependency network, we have to obtain the dependencies for
multiple packages. For illustration, we choose the [core packages of the
tidyverse](https://www.tidyverse.org/packages/), and find out what each
package `Imports`. We put all the dependencies into one data frame, in
which the package in the `from` column imports the package in the `to`
column. This is essentially the edge list of the dependency network.

``` r
df0.imports <- rbind(
    get_dep_df("ggplot2", "Imports"),
    get_dep_df("dplyr", "Imports"),
    get_dep_df("tidyr", "Imports"),
    get_dep_df("readr", "Imports"),
    get_dep_df("purrr", "Imports"),
    get_dep_df("tibble", "Imports"),
    get_dep_df("stringr", "Imports"),
    get_dep_df("forcats", "Imports")
)
head(df0.imports)
#>      from        to    type reverse
#> 1 ggplot2    digest imports   FALSE
#> 2 ggplot2      glue imports   FALSE
#> 3 ggplot2 grDevices imports   FALSE
#> 4 ggplot2      grid imports   FALSE
#> 5 ggplot2    gtable imports   FALSE
#> 6 ggplot2   isoband imports   FALSE
tail(df0.imports)
#>       from       to    type reverse
#> 59 stringr magrittr imports   FALSE
#> 60 stringr  stringi imports   FALSE
#> 61 forcats ellipsis imports   FALSE
#> 62 forcats magrittr imports   FALSE
#> 63 forcats    rlang imports   FALSE
#> 64 forcats   tibble imports   FALSE
```

## All types of dependencies, in a data frame

The example dataset `cran_dependencies` contains all dependencies as of
2020-05-09.

``` r
data(cran_dependencies)
cran_dependencies
#> # A tibble: 211,381 x 4
#>    from  to             type     reverse
#>    <chr> <chr>          <chr>    <lgl>  
#>  1 A3    xtable         depends  FALSE  
#>  2 A3    pbapply        depends  FALSE  
#>  3 A3    randomForest   suggests FALSE  
#>  4 A3    e1071          suggests FALSE  
#>  5 aaSEA DT             imports  FALSE  
#>  6 aaSEA networkD3      imports  FALSE  
#>  7 aaSEA shiny          imports  FALSE  
#>  8 aaSEA shinydashboard imports  FALSE  
#>  9 aaSEA magrittr       imports  FALSE  
#> 10 aaSEA Bios2cor       imports  FALSE  
#> # … with 211,371 more rows
dplyr::count(cran_dependencies, type, reverse)
#> # A tibble: 8 x 3
#>   type       reverse     n
#>   <chr>      <lgl>   <int>
#> 1 depends    FALSE   11123
#> 2 depends    TRUE     9672
#> 3 imports    FALSE   57617
#> 4 imports    TRUE    51913
#> 5 linking to FALSE    3433
#> 6 linking to TRUE     3721
#> 7 suggests   FALSE   35018
#> 8 suggests   TRUE    38884
```

This is essentially a snapshot of CRAN. We can obtain all the current
dependencies using `get_dep_all_packages()`, which requires no
arguments:

``` r
df0.cran <- get_dep_all_packages()
head(df0.cran)
#>    from             to    type reverse
#> 2 aaSEA             DT imports   FALSE
#> 3 aaSEA      networkD3 imports   FALSE
#> 4 aaSEA          shiny imports   FALSE
#> 5 aaSEA shinydashboard imports   FALSE
#> 6 aaSEA       magrittr imports   FALSE
#> 7 aaSEA       Bios2cor imports   FALSE
dplyr::count(df0.cran, type, reverse) # numbers in general larger than above
#>         type reverse     n
#> 1    depends   FALSE 11184
#> 2    depends    TRUE  9606
#> 3    imports   FALSE 60653
#> 4    imports    TRUE 54264
#> 5 linking to   FALSE  3653
#> 6 linking to    TRUE  3936
#> 7   suggests   FALSE 37265
#> 8   suggests    TRUE 41018
```

## Network of one type of dependencies, as an igraph object

We can build dependency network using `get_graph_all_packages()`.
Furthermore, we can verify that the forward and reverse dependency
networks are (almost) the same, by looking at their size (number of
edges) and order (number of nodes).

``` r
g0.depends <- get_graph_all_packages(type = "depends")
g0.rev_depends <- get_graph_all_packages(type = "reverse depends")
g0.depends
#> IGRAPH b77a083 DN-- 4802 7999 -- 
#> + attr: name (v/c), type (e/c), reverse
#> | (e/l)
#> + edges from b77a083 (vertex names):
#>  [1] A3      ->xtable   A3      ->pbapply 
#>  [3] abc     ->abc.data abc     ->nnet    
#>  [5] abc     ->quantreg abc     ->MASS    
#>  [7] abc     ->locfit   abcdeFBA->Rglpk   
#>  [9] abcdeFBA->rgl      abcdeFBA->corrplot
#> [11] abcdeFBA->lattice  ABCp2   ->MASS    
#> [13] abctools->abc      abctools->abind   
#> + ... omitted several edges
g0.rev_depends
#> IGRAPH c40edf9 DN-- 4802 7999 -- 
#> + attr: name (v/c), type (e/c), reverse
#> | (e/l)
#> + edges from c40edf9 (vertex names):
#>  [1] abc     ->abctools   abc     ->EasyABC   
#>  [3] abc.data->abc        abd     ->tigerstats
#>  [5] abind   ->abctools   abind   ->BCBCSF    
#>  [7] abind   ->CPMCGLM    abind   ->depth     
#>  [9] abind   ->dgmb       abind   ->dynamo    
#> [11] abind   ->fractaldim abind   ->funLBM    
#> [13] abind   ->informR    abind   ->interplot 
#> + ... omitted several edges
```

The dependency words accepted by the argument `type` is the same as in
`get_dep()` and `get_dep_df()`. The two networks’ size and order should
be very close if not identical to each other. Because of the dependency
direction, their edge lists should be the same but with the column names
`from` and `to` swapped.

For verification, the exact same graphs can be obtained by filtering the
data frame for the required dependency and applying `df_to_graph()`:

``` r
g1.depends <- df0.cran %>%
    dplyr::filter(type == "depends" & !reverse) %>%
    df_to_graph(nodelist = dplyr::rename(df0.cran, name = from))
g1.rev_depends <- df0.cran %>%
    dplyr::filter(type == "depends" & reverse) %>%
    df_to_graph(nodelist = dplyr::rename(df0.cran, name = from))
g1.depends # same as g0.depends
#> IGRAPH effe470 DN-- 4802 7999 -- 
#> + attr: name (v/c), type (e/c), reverse
#> | (e/l)
#> + edges from effe470 (vertex names):
#>  [1] A3      ->xtable   A3      ->pbapply 
#>  [3] abc     ->abc.data abc     ->nnet    
#>  [5] abc     ->quantreg abc     ->MASS    
#>  [7] abc     ->locfit   abcdeFBA->Rglpk   
#>  [9] abcdeFBA->rgl      abcdeFBA->corrplot
#> [11] abcdeFBA->lattice  ABCp2   ->MASS    
#> [13] abctools->abc      abctools->abind   
#> + ... omitted several edges
g1.rev_depends # same as g0.rev_depends
#> IGRAPH 7d0fbfb DN-- 4802 7999 -- 
#> + attr: name (v/c), type (e/c), reverse
#> | (e/l)
#> + edges from 7d0fbfb (vertex names):
#>  [1] abc     ->abctools   abc     ->EasyABC   
#>  [3] abc.data->abc        abd     ->tigerstats
#>  [5] abind   ->abctools   abind   ->BCBCSF    
#>  [7] abind   ->CPMCGLM    abind   ->depth     
#>  [9] abind   ->dgmb       abind   ->dynamo    
#> [11] abind   ->fractaldim abind   ->funLBM    
#> [13] abind   ->informR    abind   ->interplot 
#> + ... omitted several edges
```
