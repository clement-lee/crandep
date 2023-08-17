
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
    `get_dep_all_packages()`.
2.  For obtaining igraph objects of package dependencies, use
    `get_graph_all_packages()` and `df_to_graph()`.
3.  For modelling the number of dependencies, use `*pol()` and
    `*mix2()`.
4.  There is also an example data set `cran_dependencies`.

## One or multiple types of dependencies

To obtain the information about various kinds of dependencies of a
package, we can use the function `get_dep()` which takes the package
name and the type of dependencies as the first and second arguments,
respectively. Currently, the second argument accepts a character vector
of one or more of the following words: `Depends`, `Imports`,
`LinkingTo`, `Suggests`, `Enhances`, `Reverse_depends`,
`Reverse_imports`, `Reverse_linking_to`, `Reverse_suggests`, and
`Reverse_enhances`, or any variations in their letter cases, or if the
underscore "\_" is replaced by a space.

``` r
get_dep("dplyr", "Imports")
#>     from         to    type reverse
#> 1  dplyr        cli imports   FALSE
#> 2  dplyr   generics imports   FALSE
#> 3  dplyr       glue imports   FALSE
#> 4  dplyr  lifecycle imports   FALSE
#> 5  dplyr   magrittr imports   FALSE
#> 6  dplyr    methods imports   FALSE
#> 7  dplyr     pillar imports   FALSE
#> 8  dplyr         R6 imports   FALSE
#> 9  dplyr      rlang imports   FALSE
#> 10 dplyr     tibble imports   FALSE
#> 11 dplyr tidyselect imports   FALSE
#> 12 dplyr      utils imports   FALSE
#> 13 dplyr      vctrs imports   FALSE
get_dep("MASS", c("depends", "suggests"))
#>   from        to     type reverse
#> 1 MASS grDevices  depends   FALSE
#> 2 MASS  graphics  depends   FALSE
#> 3 MASS     stats  depends   FALSE
#> 4 MASS     utils  depends   FALSE
#> 5 MASS   lattice suggests   FALSE
#> 6 MASS      nlme suggests   FALSE
#> 7 MASS      nnet suggests   FALSE
#> 8 MASS  survival suggests   FALSE
```

For more information on different types of dependencies, see [the
official
guidelines](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Package-Dependencies)
and <https://r-pkgs.org/description.html>.

In the output, the column `type` is the type of the dependency converted
to lower case. Also, `LinkingTo` is now converted to `linking to` for
consistency.

``` r
get_dep("xts", "LinkingTo")
#>   from  to       type reverse
#> 1  xts zoo linking to   FALSE
get_dep("xts", "linking to")
#>   from  to       type reverse
#> 1  xts zoo linking to   FALSE
```

For the reverse dependencies, the substring `"reverse_"` will not be
shown in `type`; instead the `reverse` column will be `TRUE`. This can
be illustrated by the following:

``` r
get_dep("abc", c("depends", "reverse_depends"))
#>   from       to    type reverse
#> 1  abc abc.data depends   FALSE
#> 2  abc     nnet depends   FALSE
#> 3  abc quantreg depends   FALSE
#> 4  abc     MASS depends   FALSE
#> 5  abc   locfit depends   FALSE
#> 6  abc abctools depends    TRUE
#> 7  abc  EasyABC depends    TRUE
get_dep("xts", c("linking to", "reverse linking to"))
#>   from       to       type reverse
#> 1  xts      zoo linking to   FALSE
#> 2  xts ichimoku linking to    TRUE
#> 3  xts  RcppXts linking to    TRUE
#> 4  xts      TTR linking to    TRUE
```

Theoretically, for each forward dependency

    #>   from to type reverse
    #> 1    A  B    c   FALSE

there should be an equivalent reverse dependency

    #>   from to type reverse
    #> 1    B  A    c    TRUE

Aligning the `type` in the forward dependency and the reverse dependency
enables this to be checked easily.

To obtain all types of dependencies, we can use `"all"` in the second
argument, instead of typing a character vector of all words:

``` r
df0.abc <- get_dep("abc", "all")
df0.abc
#>    from         to     type reverse
#> 1   abc   abc.data  depends   FALSE
#> 2   abc       nnet  depends   FALSE
#> 3   abc   quantreg  depends   FALSE
#> 4   abc       MASS  depends   FALSE
#> 5   abc     locfit  depends   FALSE
#> 10  abc   abctools  depends    TRUE
#> 11  abc    EasyABC  depends    TRUE
#> 12  abc ecolottery  imports    TRUE
#> 14  abc      coala suggests    TRUE
df0.rstan <- get_dep("rstan", "all") # too many rows to display
dplyr::count(df0.rstan, type, reverse) # hence the summary using count()
#>         type reverse   n
#> 1    depends   FALSE   2
#> 2    depends    TRUE  24
#> 3   enhances    TRUE   2
#> 4    imports   FALSE   8
#> 5    imports    TRUE 127
#> 6 linking to   FALSE   5
#> 7 linking to    TRUE 113
#> 8   suggests   FALSE  13
#> 9   suggests    TRUE  31
```

As of 2023-08-17, there are 0 packages that have all 10 types of
dependencies, and 6 packages that have 9 types of dependencies: Matrix,
bigmemory, miceadds, rstan, sf, xts.

## Building and visualising a dependency network

To build a dependency network, we have to obtain the dependencies for
multiple packages. For illustration, we choose the [core packages of the
tidyverse](https://www.tidyverse.org/packages/), and find out what each
package `Imports`. We put all the dependencies into one data frame, in
which the package in the `from` column imports the package in the `to`
column. This is essentially the edge list of the dependency network.

``` r
df0.imports <- rbind(
  get_dep("ggplot2", "Imports"),
  get_dep("dplyr", "Imports"),
  get_dep("tidyr", "Imports"),
  get_dep("readr", "Imports"),
  get_dep("purrr", "Imports"),
  get_dep("tibble", "Imports"),
  get_dep("stringr", "Imports"),
  get_dep("forcats", "Imports")
)
head(df0.imports)
#>      from        to    type reverse
#> 1 ggplot2       cli imports   FALSE
#> 2 ggplot2      glue imports   FALSE
#> 3 ggplot2 grDevices imports   FALSE
#> 4 ggplot2      grid imports   FALSE
#> 5 ggplot2    gtable imports   FALSE
#> 6 ggplot2   isoband imports   FALSE
tail(df0.imports)
#>       from        to    type reverse
#> 73 forcats       cli imports   FALSE
#> 74 forcats      glue imports   FALSE
#> 75 forcats lifecycle imports   FALSE
#> 76 forcats  magrittr imports   FALSE
#> 77 forcats     rlang imports   FALSE
#> 78 forcats    tibble imports   FALSE
```

## All types of dependencies, in a data frame

The example dataset `cran_dependencies` contains all dependencies as of
2020-05-09.

``` r
data(cran_dependencies)
cran_dependencies
#> # A tibble: 211,381 × 4
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
#> # ℹ 211,371 more rows
dplyr::count(cran_dependencies, type, reverse)
#> # A tibble: 8 × 3
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
#>       from         to    type reverse
#> 3 AATtools   magrittr imports   FALSE
#> 4 AATtools      dplyr imports   FALSE
#> 5 AATtools doParallel imports   FALSE
#> 6 AATtools    foreach imports   FALSE
#> 7   ABACUS    ggplot2 imports   FALSE
#> 8   ABACUS      shiny imports   FALSE
dplyr::count(df0.cran, type, reverse) # numbers in general larger than above
#>          type reverse     n
#> 1     depends   FALSE 10560
#> 2     depends    TRUE  9157
#> 3    enhances   FALSE   564
#> 4    enhances    TRUE   583
#> 5     imports   FALSE 94243
#> 6     imports    TRUE 86459
#> 7  linking to   FALSE  5470
#> 8  linking to    TRUE  5859
#> 9    suggests   FALSE 59859
#> 10   suggests    TRUE 66261
```

## Network of one type of dependencies, as an igraph object

We can build dependency network using `get_graph_all_packages()`.
Furthermore, we can verify that the forward and reverse dependency
networks are (almost) the same, by looking at their size (number of
edges) and order (number of nodes).

``` r
g0.depends <- get_graph_all_packages(type = "depends")
g0.depends
#> IGRAPH 4bb852e DN-- 4631 7539 -- 
#> + attr: name (v/c)
#> + edges from 4bb852e (vertex names):
#>  [1] A3          ->xtable     A3          ->pbapply    abc         ->abc.data  
#>  [4] abc         ->nnet       abc         ->quantreg   abc         ->MASS      
#>  [7] abc         ->locfit     ABCp2       ->MASS       abctools    ->abc       
#> [10] abctools    ->abind      abctools    ->plyr       abctools    ->Hmisc     
#> [13] abd         ->nlme       abd         ->lattice    abd         ->mosaic    
#> [16] abodOutlier ->cluster    abundant    ->glasso     Ac3net      ->data.table
#> [19] acc         ->mhsmm      accelmissing->mice       accelmissing->pscl      
#> [22] accessrmd   ->ggplot2    accrual     ->tcltk2     accrualPlot ->lubridate 
#> + ... omitted several edges
```

We could obtain essentially the same graph, but with the direction of
the edges reversed, by specifying `type = "reverse depends"`:

``` r
# Not run
g0.rev_depends <- get_graph_all_packages(type = "reverse depends")
g0.rev_depends
```

The dependency words accepted by the argument `type` is the same as in
`get_dep()`. The two networks’ size and order should be very close if
not identical to each other. Because of the dependency direction, their
edge lists should be the same but with the column names `from` and `to`
swapped.

For verification, the exact same graphs can be obtained by filtering the
data frame for the required dependency and applying `df_to_graph()`:

``` r
g1.depends <- df0.cran |>
  dplyr::filter(type == "depends" & !reverse) |>
  df_to_graph(nodelist = dplyr::rename(df0.cran, name = from))
g1.depends # same as g0.depends
#> IGRAPH 39d30f3 DN-- 4631 7539 -- 
#> + attr: name (v/c), type (e/c), reverse (e/l)
#> + edges from 39d30f3 (vertex names):
#>  [1] A3          ->xtable     A3          ->pbapply    abc         ->abc.data  
#>  [4] abc         ->nnet       abc         ->quantreg   abc         ->MASS      
#>  [7] abc         ->locfit     ABCp2       ->MASS       abctools    ->abc       
#> [10] abctools    ->abind      abctools    ->plyr       abctools    ->Hmisc     
#> [13] abd         ->nlme       abd         ->lattice    abd         ->mosaic    
#> [16] abodOutlier ->cluster    abundant    ->glasso     Ac3net      ->data.table
#> [19] acc         ->mhsmm      accelmissing->mice       accelmissing->pscl      
#> [22] accessrmd   ->ggplot2    accrual     ->tcltk2     accrualPlot ->lubridate 
#> + ... omitted several edges
```
