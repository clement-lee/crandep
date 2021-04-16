#' Dependencies of CRAN packages
#'
#' A dataset containing the dependencies of various types (Imports, Depends, Suggests, LinkingTo, and their reverse counterparts) of more than 14600 packages available on CRAN as of 2020-05-09.
#'
#' @format A data frame with 211408 rows and 4 variables:
#' \describe{
#'   \item{from}{the name of the package that introduced the dependencies}
#'   \item{to}{the name of the package that the dependency is directed towards}
#'   \item{type}{the type of dependency, which can take the follow values (all in lowercase): "depends", "imports", "linking to", "suggests"}
#'   \item{reverse}{a boolean representing whether the dependency is a reverse one (TRUE) or a forward one (FALSE)}
#' }
#' @source The CRAN pages of all the packages available on \url{https://cran.r-project.org}
"cran_dependencies"

#' Citation network of CHI papers
#'
#' A dataset containing the citations of conference papers of the ACM Conference on Human Factors in Computing Systems (CHI) from 1981 to 2019, obtained from the ACM digital library. The resulting citation network can be compared to the dependencies network of CRAN packages, in terms of network-related characteristics, such as degree distribution and sparsity.
#'
#' @format A data from with21951 rows and 4 variables:
#' \describe{
#'   \item{from}{the unique identifier (in the digital library) of the paper that cites other papers}
#'   \item{to}{the unique identifier of the paper that is being cited}
#'   \item{year_from}{the publication year of the citing paper}
#'   \item{year_to}{the publication year of the cited paper}
#' }
#' @source \url{https://dl.acm.org/conference/chi}
"chi_citations"
