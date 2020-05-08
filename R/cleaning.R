#' Obtain the URL on CRAN using the package name
#'
#' @param name String, name of the package
#' @return A string of the URL for the package on CRAN
#' @examples
#' \dontrun{cran_url("dplyr")}
#' @keywords internal
cran_url <- function(name) {
    paste0("https://CRAN.R-project.org/package=", name) # canonical form
}

#' Scrap the page (of a package) as a text vector
#'
#' @param url An URL
#' @importFrom xml2 read_html
#' @importFrom rvest html_text
#' @importFrom stringr str_split
#' @return A string vector of the html text of the page according to the url
#' @examples
#' \dontrun{
#' str.dplyr <- "https://cran.r-project.org/web/packages/dplyr/index.html" # CRAN page for dplyr
#' html.dplyr <- html_text_vec(str.dplyr)
#' html.abc <- html_text_vec(cran_url("abc")) # the page for abc on CRAN
#' }
#' @keywords internal
html_text_vec <- function(url) {
    as.vector(stringr::str_split(rvest::html_text(xml2::read_html(url)), "\n", simplify = TRUE))
}

#' Find string corresponding to "Imports", "Depends" etc.
#'
#' @param v A vector of strings
#' @param x One of the following dependency words: "Depends", "Imports", "LinkingTo", "Suggests", "Reverse_depends", "Reverse_imports", "Reverse_linking_to", "Reverse_suggests"
#' @importFrom stringr str_sub str_detect str_replace_all str_to_title
#' @importFrom glue glue
#' @return A string of the concatenation of the dependencies
#' @examples
#' \dontrun{
#' get_dep_str(html_text_vec(cran_url("MASS")), "Depends")
#' get_dep_str(html_text_vec(cran_url("dplyr")), "Imports")
#' get_dep_str(html_text_vec(cran_url("Rcpp")), "Reverse_depends")
#' get_dep_str(html_text_vec(cran_url("lme4")), "LinkingTo")
#' }
#' @keywords internal
get_dep_str <- function(v, x) {
    words <- c("Depends", "Imports", "LinkingTo", "Suggests", "Reverse\u00a0depends", "Reverse\u00a0imports", "Reverse\u00a0linking\u00a0to", "Reverse\u00a0suggests")
    x <- stringr::str_to_title(x)
    x <- stringr::str_replace_all(x, "_", "\u00a0") # \u00a0 is the actual underscore used in pages, while the usual one is \u005F
    if (x == "Linkingto") {
        x <- "LinkingTo"
    }
    if (!(x %in% words)) {
        stop(as.character(glue::glue("x has to be one of the following: {paste(stringr::str_replace_all(words, '\u00a0', '_'), collapse = ', ')}."))) 
    }
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
#' @param x A scalar string, possibly an output of get_dep_str()
#' @importFrom stringr str_split str_locate str_sub
#' @return A string vector of dependencies
#' @examples
#' \dontrun{
#' string <- get_dep_str(html_text_vec(cran_url("MASS")), "Depends") # the packages MASS depends on
#' get_dep_vec(string) # R (<version>) will be removed
#' }
#' @keywords internal
get_dep_vec <- function(x) {
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

#' Obtain one type of dependencies of a package directly. This is essentially the same as chaining get
#'
#' @param name String, name of the package
#' @param type One of the following dependency words: "Depends", "Imports", "LinkingTo", "Suggests", "Reverse_depends", "Reverse_imports", "Reverse_linking_to", "Reverse_suggests"
#' @return A string vector of dependencies
#' @examples
#' get_dep_all("dplyr", "Imports")
#' get_dep_all("MASS", "Depends")
#' @export
get_dep_all <- function(name, type) {
    get_dep_vec(get_dep_str(html_text_vec(cran_url(name)), type))
}

#' Obtain the data frame of one kind of dependencies
#'
#' @param df A data frame with at least two columns: name - a string vector, and <type> - a list of string vectors
#' @param type The dependency type desired
#' @importFrom rlang !! enquo .data
#' @importFrom dplyr select filter rename
#' @importFrom tidyr unnest
#' @return A data frame with two columns: from & to
#' @examples
#' if (requireNamespace("tibble", quietly = TRUE)) {
#'     name.vec <- c("dplyr", "MASS", "Rcpp")
#'     depends.vec <- rep(NA, length(name.vec))
#'     imports.vec <- rep(NA, length(name.vec))
#'     for (i in seq_along(name.vec)) {
#'         depends.vec[i] <- list(get_dep_all(name.vec[i], "Depends"))
#'         imports.vec[i] <- list(get_dep_all(name.vec[i], "Imports"))
#'     }
#'     df0 <- tibble::tibble(name = name.vec, depends = depends.vec, imports = imports.vec)
#'     unnest_dep(df0, depends)
#'     unnest_dep(df0, imports)
#' }
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
#' @examples
#' from <- c("1", "2", "4")
#' to <- c("2", "3", "5")
#' edges <- data.frame(from = from, to = to, stringsAsFactors = FALSE)
#' nodes <- data.frame(name = c("1", "2", "3", "4", "5"), stringsAsFactors = FALSE)
#' df_to_graph(edges, nodes)
#' @export
df_to_graph <- function(edgelist, nodelist) {
    l <- igraph::decompose.graph(igraph::graph_from_data_frame(dplyr::semi_join(edgelist, nodelist, c("to" = "name")))) # semi join as some nodes (packages) may have become obsolete
    l[[which.max(purrr::map_int(purrr::map(l, igraph::V), length))]]
}
