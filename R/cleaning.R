#' Conditionally change a string
#'
#' @param x A character vector
#' @param from A character vector of words to change from
#' @param to A string to change to
#' @return A string
#' @keywords internal
conditional_change <- function(x, from, to) {
  ifelse(x %in% from, to, x)
}

#' Check and convert dependency word(s)
#'
#' @param x A character vector of dependency words
#' @importFrom stringr str_to_title str_replace_all
#' @return A character vector of modified dependency words
#' @keywords internal
check_dep_word <- function(x) {
  types <- c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances")
  if (length(x) == 1L && stringr::str_to_title(x) == "All") {
    x <- types
  } else if (length(x) == 1L && stringr::str_to_title(x) == "Strong") {
    x <- c("Depends", "Imports", "LinkingTo")
  } else if (any(substr(tolower(x), 1L, 7L) == "reverse")) {
    stop("check_dep_word: some dependency types contain 'Reverse ' up to letter case. Remove these substrings and use the argument 'reverse' in the function that calls check_dep_word.")
  } else {
    x <- stringr::str_replace_all(x, " ", "_")
    x <- stringr::str_replace_all(x, "\u00a0", "_")
    x <- stringr::str_to_title(x)
    x <- conditional_change(x, c("Linkingto", "Linking_to"), "LinkingTo")
    x <- stringr::str_replace_all(x, "_", " ")
  }
  if (!all(x %in% types)) {
    s <- paste(types, collapse = ", ")
    stop(as.character(paste0("check_dep_word: each dependency word has to be one of the following: ", s, ", up to letter case and space replaced by underscore.")))
  }
  x
}

#' Split a string to a list of dependencies
#'
#' @param x A scalar string, possibly an output of get_dep_str()
#' @importFrom stringr str_split str_locate str_replace_all
#' @return A string vector of dependencies
#' @keywords internal
get_dep_vec <- function(x) {
  if (is.na(x)) {
    y <- as.character(NA)
  } else {
    u <- stringr::str_replace_all(x, ",$", "") # the strings of dependencies
    u <- stringr::str_replace_all(u, ",", ", ")
    u <- stringr::str_replace_all(u, "\\(", " \\(")
    u <- stringr::str_replace_all(u, "\n", " ")
    u <- stringr::str_replace_all(u, "  ", " ")
    u <- stringr::str_replace_all(u, "  ", " ")
    u <- stringr::str_split(u, ", ")
    v <- stringr::str_locate(u[[1L]], " ")[, 1L] # the indices of 1st space
    w <- substr(u[[1L]], 1L, ifelse(is.na(v), nchar(u[[1L]]), v - 1L)) # the substrings
    if (length(w) == 1L && w == "R") {
      y <- as.character(NA)
    } else {
      y <- w[w != "R"]
    }
  }
  y
}
