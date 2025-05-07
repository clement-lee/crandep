#' Construct the giant component of the network from two data frames
#'
#' @param edgelist A data frame with (at least) two columns: from and to
#' @param nodelist NULL, or a data frame with (at least) one column: name, that contains the nodes to include
#' @param gc Boolean, if 'TRUE' (default) then the giant component is extracted, if 'FALSE' then the whole graph is returned
#' @importFrom dplyr semi_join
#' @importFrom igraph graph_from_data_frame decompose.graph gorder
#' @return An igraph object & a connected graph if gc is 'TRUE'
#' @examples
#' from <- c("1", "2", "4")
#' to <- c("2", "3", "5")
#' edges <- data.frame(from = from, to = to, stringsAsFactors = FALSE)
#' nodes <- data.frame(name = c("1", "2", "3", "4", "5"), stringsAsFactors = FALSE)
#' df_to_graph(edges, nodes)
#' @export
df_to_graph <- function(edgelist, nodelist = NULL, gc = TRUE) {
  if (is.null(nodelist)) {
    g <- igraph::graph_from_data_frame(edgelist)
  } else {
    nodelist <- unique(nodelist[, "name", drop = FALSE])
    df <- dplyr::semi_join(edgelist, nodelist, c("to" = "name")) # semi join as some nodes may have become obsolete
    df <- dplyr::semi_join(df, nodelist, c("from" = "name"))
    g <- igraph::graph_from_data_frame(df, vertices = nodelist$name)
  }
  if (gc) {
    l <- igraph::decompose.graph(g)
    v <- sapply(l, igraph::gorder, simplify = TRUE)
    g <- l[[which.max(v)]]
  }
  g
}

#' Graph of dependencies of all CRAN packages
#'
#' \code{get_graph_all_packages} returns an igraph object representing the network of one or more types of dependencies of all CRAN packages.
#' @param type A character vector that contains one or more of the following dependency words: "Depends", "Imports", "LinkingTo", "Suggests", "Enhances", up to letter case and space replaced by underscore. Alternatively, if 'type = "all"', all five dependencies will be obtained; if 'type = "strong"', "Depends", "Imports" & "LinkingTo" will be obtained.
#' @param gc Boolean, if 'TRUE' (default) then the giant component is extracted, if 'FALSE' then the whole graph is returned
#' @param reverse Boolean, whether forward (FALSE, default) or reverse (TRUE) dependencies are requested.
#' @return An igraph object & a connected graph if gc is 'TRUE'
#' @importFrom dplyr inner_join
#' @examples
#' \dontrun{
#' g0.cran.depends <- get_graph_all_packages("depends")
#' g1.cran.imports <- get_graph_all_packages("imports", reverse = TRUE)
#' }
#' @seealso \code{\link{get_dep_all_packages}} for the dependencies of all CRAN packages in a data frame, and \code{\link{df_to_graph}} for constructing the giant component of the network from two data frames
#' @export
get_graph_all_packages <- function(type, gc = TRUE, reverse = FALSE) {
  ## change params to align with others
  type <- check_dep_word(type)
  l0 <- get_dep_all_packages()
  df1 <- data.frame(type = type, reverse = reverse)
  df1$type <- conditional_change(tolower(df1$type), "linkingto", "linking to")
  df2 <- dplyr::inner_join(l0$dependencies, df1, c("type", "reverse")) # edgelist
  df2 <- unique(df2[, c("from", "to")])
  df_to_graph(df2, l0$packages, gc)
}
