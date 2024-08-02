#' Return a sorted vector of nodes id
#'
#' @param g An igraph object of a DAG
#' @param random Boolean, whether the order of selected nodes is randomised in the process
#' @return A data frame with two columns: "id" is the names of nodes in g, and "id_num" is the topological ordering
#' @examples
#' df0 <- data.frame(from = c("a", "b"), to = c("b", "c"), stringsAsFactors = FALSE)
#' g0 <- igraph::graph_from_data_frame(df0, directed = TRUE)
#' topo_sort_kahn(g0)
#' @importFrom igraph as_data_frame
#' @importFrom igraph V
#' @importFrom igraph as_adjacency_matrix
#' @importFrom igraph is_dag
#' @export
topo_sort_kahn <- function(g, random = FALSE) {
  if (!is_dag(g)) {
    stop("g has to be a DAG")
  }
  e0 <- igraph::as_data_frame(g)
  names(e0) <- c("citing", "cited")
  v <- names(igraph::V(g))
  v1 <- sort(unique(c(e0$citing, e0$cited)))
  l <- setdiff(v, v1)
  s0 <- sort(setdiff(e0$citing, e0$cited))
  while (length(s0) > 0L) {
    if (random) {
      n0 <- sample(s0, 1L)
      s0 <- s0[s0 != n0]
    } else {
      n0 <- s0[1L]
      s0 <- s0[-1L]
    }
    l <- c(l, n0)
    ## outgoing edges of n0
    i0 <- e0$citing == n0
    e1 <- e0[i0, , drop = FALSE]
    e0 <- e0[!i0, , drop = FALSE]
    if (nrow(e1) != 0L) {
      e2 <- setdiff(e1$cited, e0$cited)
      if (random) {
        e2 <- sample(e2, length(e2))
      }
      s0 <- c(s0, e2)
    }
  }
  if (nrow(e0) > 0L) {
    stop("topo_sort_kahn: graph has at least 1 cycle")
  }
  o <- match(l, v)
  a <- as_adjacency_matrix(g)[o, o]
  data.frame(
    id = l,
    id_num = seq_along(l),
    stringsAsFactors = FALSE
  )
}

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
    df <- dplyr::semi_join(edgelist, nodelist, c("to" = "name")) # semi join as some nodes may have become obsolete
    g <- igraph::graph_from_data_frame(df)
  }
  if (gc) {
    l <- igraph::decompose.graph(g)
    n <- length(l)
    v <- rep(as.integer(NA), n)
    for (i in seq(n)) {
      v[i] <- igraph::gorder(l[[i]])
    }
    g <- l[[which.max(v)]]
  }
  g
}

#' Graph of dependencies of all CRAN packages
#'
#' \code{get_graph_all_packages} returns an igraph object representing the network of one or more types of dependencies of all CRAN packages.
#' @param type A character vector that contains one or more of the following dependency words: "Depends", "Imports", "LinkingTo", "Suggests", "Enhances", up to letter case and space replaced by underscore. Alternatively, if 'types = "all"', all five dependencies will be obtained.
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
  df0 <- get_dep_all_packages()
  df1 <- data.frame(type = type, reverse = reverse)
  df1$type <- conditional_change(tolower(df1$type), "linkingto", "linking to")
  df2 <- dplyr::inner_join(df0, df1, c("type", "reverse")) # edgelist
  df2 <- unique(df2[, c("from", "to")])
  df3 <- data.frame(name = unique(df0$from))
  df_to_graph(df2, df3, gc)
}
