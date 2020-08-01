#' Check and convert dependency word(s)
#'
#' @param x A character vector of dependency words
#' @importFrom stringr str_to_title str_replace_all
#' @return A character vector of modified dependency words
#' @keywords internal
check_dep_word <- function(x) {
    types <- c("Depends", "Imports", "LinkingTo", "Suggests", "Reverse depends", "Reverse imports", "Reverse linking to", "Reverse suggests")
    if (length(x) == 1L && stringr::str_to_title(x) == "All") {
        x <- types
    } else {
        x <- stringr::str_replace_all(x, " ", "_")
        x <- stringr::str_replace_all(x, "\u00a0", "_")
        x <- stringr::str_to_title(x)
        x <- ifelse(x == "Linkingto", "LinkingTo", x)
        x <- ifelse(x == "Linking_to", "LinkingTo", x)
        x <- stringr::str_replace_all(x, "_", " ")
    }
    if (!all(x %in% types)) {
        s <- paste(types, collapse = ", ")
        stop(as.character(paste0("check_dep_word: each dependency word has to be one of the following: ", s, ", up to letter case and space replaced by underscore.")))
    }
    x
}

#' Obtain the URL on CRAN using the package name
#'
#' @param name String, name of the package
#' @return A string of the URL for the package on CRAN
#' @examples
#' crandep:::cran_url("dplyr")
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
#' url.abc <- crandep:::cran_url("abc") # the page for abc on CRAN
#' html.abc <- crandep:::html_text_vec(url.abc)
#' @keywords internal
html_text_vec <- function(url) {
    as.vector(stringr::str_split(rvest::html_text(xml2::read_html(url)), "\n", simplify = TRUE))
}

#' Find string corresponding to "Imports", "Depends" etc.
#'
#' @param v A vector of strings
#' @param x One of the following dependency words: "Depends", "Imports", "LinkingTo", "Suggests", "Reverse_depends", "Reverse_imports", "Reverse_linking_to", "Reverse_suggests"
#' @importFrom stringr str_detect str_replace_all
#' @return A string of the concatenation of the dependencies
#' @examples
#' url.mass <- crandep:::cran_url("MASS")
#' html.mass <- crandep:::html_text_vec(url.mass)
#' crandep:::get_dep_str(html.mass, "Depends")
#' @keywords internal
get_dep_str <- function(v, x) {
    x <- check_dep_word(x)
    x <- stringr::str_replace_all(x, " ", "\u00a0") # \u00a0 is the actual underscore used in pages, while the usual one is \u005F
    x <- paste0(x, ":")
    s <- seq_along(v)
    i <- which(substr(v, 1L, nchar(x)) == x)
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
#' @importFrom stringr str_split str_locate str_replace_all
#' @return A string vector of dependencies
#' @examples
#' url.mass <- crandep:::cran_url("MASS")
#' html.mass <- crandep:::html_text_vec(url.mass)
#' str.mass <- crandep:::get_dep_str(html.mass, "Depends") # the packages MASS depends on
#' crandep:::get_dep_vec(str.mass) # R (<version>) will be removed
#' @keywords internal
get_dep_vec <- function(x) {
    if (is.na(x)) {
        y <- as.character(NA)
    } else {
        x <- stringr::str_replace_all(x, ",", ", ")
        x <- stringr::str_replace_all(x, "\\(", " \\(")
        x <- stringr::str_replace_all(x, "\n", " ")
        x <- stringr::str_replace_all(x, "  ", " ")
        x <- stringr::str_replace_all(x, "  ", " ")
        u <- stringr::str_split(x, ", ")[[1L]] # the strings of dependencies
        v <- stringr::str_locate(u, " ")[, 1L] # the indices of 1st space
        w <- substr(u, 1L, ifelse(is.na(v), nchar(u), v - 1L)) # the substrings
        if (length(w) == 1L && w == "R") {
            y <- as.character(NA)
        } else {
            y <- w[w != "R"]
        }
    }
    y
}
