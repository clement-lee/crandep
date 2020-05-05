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
