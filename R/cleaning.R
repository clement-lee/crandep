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
#' @importFrom stringr str_sub str_detect str_replace_all str_to_title
#' @return A string of the concatenation of the dependencies
#' @examples
#' url.mass <- crandep:::cran_url("MASS")
#' html.mass <- crandep:::html_text_vec(url.mass)
#' crandep:::get_dep_str(html.mass, "Depends")
#' @keywords internal
get_dep_str <- function(v, x) {
    words <- c("Depends", "Imports", "LinkingTo", "Suggests", "Reverse\u00a0depends", "Reverse\u00a0imports", "Reverse\u00a0linking\u00a0to", "Reverse\u00a0suggests")
    x <- stringr::str_to_title(x)
    x <- stringr::str_replace_all(x, "_", "\u00a0") # \u00a0 is the actual underscore used in pages, while the usual one is \u005F
    if (x == "Linkingto") {
        x <- "LinkingTo"
    }
    if (!(x %in% words)) {
        s <- paste(stringr::str_replace_all(words, '\u00a0', '_'), collapse = ', ')
        stop(as.character(paste0("x has to be one of the following: ", s, "."))) 
    }
    x <- paste0(x, ":")
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
#' url.mass <- crandep:::cran_url("MASS")
#' html.mass <- crandep:::html_text_vec(url.mass)
#' str.mass <- crandep:::get_dep_str(html.mass, "Depends") # the packages MASS depends on
#' crandep:::get_dep_vec(str.mass) # R (<version>) will be removed
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

#' Obtain dependencies of all CRAN packages
#'
#' @importFrom tools CRAN_package_db
#' @return A data frame of dependencies of all CRAN packages
#' @examples
#' df.cran <- crandep:::get_dep_all_packages()
#' @keywords internal
get_dep_all_packages <- function() {
    df0 <- as.data.frame(tools::CRAN_package_db())
    v0 <- v1 <- names(df0)
    v1[v0 == "Reverse depends"] <- "Reverse_depends"
    v1[v0 == "Reverse imports"] <- "Reverse_imports"
    v1[v0 == "Reverse linking to"] <- "Reverse_linking_to"
    v1[v0 == "Reverse suggests"] <- "Reverse_suggests"
    v1[v0 == "Reverse enhances"] <- "Reverse_enhances"
    names(df0) <- v1
    df0
}

#' Obtain one type of dependencies of a package directly
#'
#' @param name String, name of the package
#' @param type One of the following dependency words: "Depends", "Imports", "LinkingTo", "Suggests", "Reverse_depends", "Reverse_imports", "Reverse_linking_to", "Reverse_suggests"
#' @param scrape Boolean. If 'TRUE' (default), the page of the package will be scraped. If 'FALSE', tools::CRAN_package_db() will be used.
#' @importFrom stringr str_replace_all
#' @return A string vector of dependencies
#' @examples
#' get_dep_all("dplyr", "Imports")
#' get_dep_all("MASS", "Depends")
#' get_dep_all("MASS", "Depends", FALSE) # same result as above line
#' @export
get_dep_all <- function(name, type, scrape = TRUE) {
    if (scrape) {
        str0 <- get_dep_str(html_text_vec(cran_url(name)), type)
    } else {
        df0 <- get_dep_all_packages()
        str0 <- df0[df0$Package == name, type][1L] # some packages have multiple rows
        str0 <- stringr::str_replace_all(str0, "\n", " ")
    }
    unique(get_dep_vec(str0))
}

#' Obtain the data frame of multiple kinds of dependencies
#' @param name String, name of the package
#' @param type A character vector that contains one or more of the following dependency words: "Depends", "Imports", "LinkingTo", "Suggests", "Reverse_depends", "Reverse_imports", "Reverse_linking_to", "Reverse_suggests"
#' @param scrape Boolean. If 'TRUE' (default), the page of the package will be scraped. If 'FALSE', tools::CRAN_package_db() will be used.
#' @importFrom tools CRAN_package_db
#' @importFrom stringr str_to_lower str_detect str_sub
#' @return A data frame of dependencies
#' @examples
#' get_dep_df("dplyr", c("Imports", "Depends"))
#' get_dep_df("MASS", c("Suggests", "Depends", "Imports"))
#' get_dep_df("MASS", c("Suggests", "Depends", "Imports"), FALSE) # same result as above line
#' @export
get_dep_df <- function(name, type, scrape = TRUE) {
    l0 <- list()
    if (scrape) {
        html0 <- html_text_vec(cran_url(name))
        for (i in seq_along(type)) {
            typei <- type[i]
            v0 <- get_dep_vec(get_dep_str(html0, typei))
            l0[[i]] <- data.frame(from = name, to = v0, type = typei, stringsAsFactors = FALSE)
        }
    } else {
        df1 <- get_dep_all_packages()
        for (i in seq_along(type)) {
            typei <- type[i]
            str0 <- df1[df1$Package == name, typei]
            str0 <- stringr::str_replace_all(str0[1L], "\n", " ") # some packages have multiple rows
            v0 <- get_dep_vec(str0)
            l0[[i]] <- data.frame(from = name, to = v0, type = typei, stringsAsFactors = FALSE)
        }
    }
    df0 <- do.call(rbind, l0)
    df0 <- df0[!is.na(df0$to),]
    df0$type <- stringr::str_to_lower(df0$type)
    df0$type <- ifelse(df0$type == "linkingto", "linking_to", df0$type)
    df0$reverse <- stringr::str_detect(df0$type, "reverse_")
    df0$type <- ifelse(df0$reverse, stringr::str_sub(df0$type, 9L, -1L), df0$type)
    unique(df0) # there are duplicates
}

#' Construct the giant component of the network from two data frames
#'
#' @param edgelist A data frame with (at least) two columns: from and to
#' @param nodelist NULL, or a data frame with (at least) one column: name, that contains the nodes to include
#' @param gc Boolean, if 'TRUE' (default) then the giant component is extracted, if 'FALSE' then the whole graph is returned
#' @importFrom dplyr semi_join
#' @importFrom igraph graph_from_data_frame decompose.graph gorder
#' @return An igraph object & a connected graph
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
