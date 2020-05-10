---
title: "Introduction to crandep"
date: "2020-05-10"
output: md_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



The goal of crandep is to provide functions for analysing the dependencies of CRAN packages using social network analysis.

## Installation

You can install crandep from github with:


```r
# install.packages("devtools")
devtools::install_github("clement-lee/crandep")
```




```r
library(crandep)
library(dplyr)
library(ggplot2)
library(igraph)
```





## One kind of dependencies 
To obtain the information about various kinds of dependencies of a package, we can use the function `get_dep_all()` which takes the package name and the type of dependencies as the first and second arguments, respectively. Currently, the second argument accepts `Depends`, `Imports`, `LinkingTo`, `Suggests`, `Reverse_depends`, `Reverse_imports`, `Reverse_linking_to`, and `Reverse_suggests`, or any variations in their letter cases.


```r
get_dep_all("dplyr", "Imports")
#>  [1] "ellipsis"   "assertthat" "glue"       "magrittr"   "methods"   
#>  [6] "pkgconfig"  "R6"         "Rcpp"       "rlang"      "tibble"    
#> [11] "tidyselect" "utils"
get_dep_all("MASS", "depends")
#> [1] "grDevices" "graphics"  "stats"     "utils"
get_dep_all("MASS", "dePends")
#> [1] "grDevices" "graphics"  "stats"     "utils"
```

`Imports` and `Depends` are the most common types of dependencies in `R` packages, but there are other types such as `Suggests`. For more information on different types of dependencies, see [the official guidelines](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Package-Dependencies) and [http://r-pkgs.had.co.nz/description.html](http://r-pkgs.had.co.nz/description.html).





## Multiple kind of dependencies
As the information all dependencies of one package are on the same page on CRAN, to avoid scraping the same multiple times, we can use `get_dep_df()` instead of  `get_dep_all()`. The output will be a data frame instead of a character vector.


```r
get_dep_df("dplyr", c("imports", "LinkingTo"))
#>     from         to       type reverse
#> 1  dplyr   ellipsis    imports   FALSE
#> 2  dplyr assertthat    imports   FALSE
#> 3  dplyr       glue    imports   FALSE
#> 4  dplyr   magrittr    imports   FALSE
#> 5  dplyr    methods    imports   FALSE
#> 6  dplyr  pkgconfig    imports   FALSE
#> 7  dplyr         R6    imports   FALSE
#> 8  dplyr       Rcpp    imports   FALSE
#> 9  dplyr      rlang    imports   FALSE
#> 10 dplyr     tibble    imports   FALSE
#> 11 dplyr tidyselect    imports   FALSE
#> 12 dplyr      utils    imports   FALSE
#> 13 dplyr         BH linking_to   FALSE
#> 14 dplyr      plogr linking_to   FALSE
#> 15 dplyr       Rcpp linking_to   FALSE
```
The column `type` is the type of the dependency converted to lower case. Also, `LinkingTo` is now converted to `linking_to` for consistency. For the four reverse dependencies, the substring `"reverse_"` will not be shown in `type`; instead the `reverse` column will be `TRUE`. This can be illustrated by the following:


```r
get_dep_all("abc", "depends")
#> [1] "abc.data" "nnet"     "quantreg" "MASS"     "locfit"
get_dep_all("abc", "reverse_depends")
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

```
#>   from to type reverse
#> 1    A  B    c   FALSE
```
there should be an equivalent reverse dependency

```
#>   from to type reverse
#> 1    B  A    c    TRUE
```
Aligning the `type` in the forward dependency and the reverse dependency enables this to be checked easily.





## Building and visualising a dependency network
To build a dependency network, we have to obtain the dependencies for multiple packages. For illustration, we choose the [core packages of the tidyverse](https://www.tidyverse.org/packages/), and find out what each package `Imports`. We put all the dependencies into one data frame, in which the package in the `from` column imports the package in the `to` column. This is essentially the edge list of the dependency network.


```r
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
#> 61 stringr magrittr imports   FALSE
#> 62 stringr  stringi imports   FALSE
#> 63 forcats ellipsis imports   FALSE
#> 64 forcats magrittr imports   FALSE
#> 65 forcats    rlang imports   FALSE
#> 66 forcats   tibble imports   FALSE
```

