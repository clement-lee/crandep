#' Obtain the URL on CRAN using the package name
#'
#' @param name String, name of the package
#' @importFrom glue glue
#' @return The URL for the package on CRAN
#' @export
cran_url <- function(name) {
    glue::glue("https://CRAN.R-project.org/package={name}") # canonical form
}

#' Scrap the page (of a package) as a text vector
#'
#' @param url An URL
#' @importFrom xml2 read_html
#' @importFrom rvest html_text
#' @importFrom stringr str_split
#' @return The html text of the page as a string vector
#' @export
html_text_vec <- function(url) {
    as.vector(stringr::str_split(rvest::html_text(xml2::read_html(url)), "\n", simplify = TRUE))
}

#' Find string corresponding to "Imports", "Depends" etc.
#'
#' @param v A vector of strings
#' @param x One of the following dependency words: "Depends", "Imports", "LinkingTo", "Suggests", "Reverse_depends", "Reverse_imports", "Reverse_linking_to", "Reverse_suggests"
#' @importFrom stringr str_sub str_detect str_replace_all
#' @importFrom glue glue
#' @return A scalar, concatenated string of dependencies
get_dep_str <- function(v, x) {
    x <- stringr::str_replace_all(x, "_", "\u00a0") # \u00a0 is the actual underscore used in pages, while the usual one is \u005F
    x <- glue::glue("{x}:")
    s <- seq_along(v)
    i <- which(stringr::str_sub(v, 1L, nchar(x)) == x)
    if (length(i) == 0L) {
        y <- as.character(NA)
    } else {
        j <- which(s > i & stringr::str_detect(v, ":"))[1L] # look for next line with ":"
        y <- paste(v[(i + 1L):(j - 1L)], collapse = " ") # concatenate all lines in between
    }
    y
}

#' Split a string to a list of dependencies
#'
#' @param x A scalar string, preferably an output of get_dep_str()
#' @importFrom stringr str_split str_locate str_sub
#' @return A string vector of dependencies
#' @export
list_of_dep <- function(x) {
    if (is.na(x)) {
        y <- as.character(NA)
    } else {
        u <- str_split(x, ", ")[[1L]] # the strings of dependencies
        v <- str_locate(u, " ")[, 1L] # the indices of 1st space
        w <- str_sub(u, 1L, ifelse(is.na(v), -1L, v - 1L)) # the substrings
        if (length(w) == 1L && w == "R") {
            y <- as.character(NA)
        } else {
            y <- w[w != "R"]
        }
    }
    y
}

#' Obtain the data frame of one kind of dependencies
#'
#' @param df A data frame with at least two columns: name - a string vector, and <type> - a list of string vectors
#' @param type The dependency type desired
#' @importFrom rlang !! enquo .data
#' @importFrom dplyr select filter rename
#' @importFrom tidyr unnest
#' @return A data frame with two columns: from & to
#' @export
unnest_dep <- function(df, type) {
    type_enquo <- rlang::enquo(type)
    dplyr::rename(dplyr::filter(tidyr::unnest(dplyr::select(df, .data$name, !!type_enquo)), !is.na(!!type_enquo)), from = .data$name, to = !!type_enquo)
}


#' Construct the giant component of the network from two data frames
#'
#' @param edgelist A data frame with (at least) two columns: from and to, preferably an output of unnest_dep()
#' @param nodelist A data frame with (at least) one column: name, that contains the nodes to include
#' @importFrom dplyr semi_join
#' @importFrom igraph graph_from_data_frame decompose.graph V
#' @importFrom purrr map_int map
#' @return An igraph object & a connected graph
df_to_graph <- function(edgelist, nodelist) {
    l <- igraph::decompose.graph(igraph::graph_from_data_frame(dplyr::semi_join(edgelist, nodelist, c("to" = "name")))) # semi join as some nodes (packages) may have become obsolete
    l[[which.max(purrr::map_int(purrr::map(l, igraph::V), length))]]
}
