#' Multiple types of dependencies
#'
#' \code{get_dep} returns a data frame of multiple types of dependencies of a package
#' @param name String, name of the package
#' @param type A character vector that contains one or more of the following dependency words: "Depends", "Imports", "LinkingTo", "Suggests", "Enhances", up to letter case and space replaced by underscore. Alternatively, if 'type = "all"', all five dependencies will be obtained.
#' @param reverse Boolean, whether forward (FALSE, default) or reverse (TRUE) dependencies are requested.
#' @importFrom tools package_dependencies
#' @importFrom utils available.packages findCRANmirror
#' @return A data frame of dependencies
#' @examples
#' get_dep("dplyr", c("Imports", "Depends"))
#' get_dep("MASS", c("Suggests", "Depends", "Imports"), TRUE)
#' @seealso \code{\link{get_dep_all_packages}} for the dependencies of all CRAN packages, and \code{\link{get_graph_all_packages}} for obtaining directly a network of dependencies as an igraph object
#' @export
get_dep <- function(name, type, reverse = FALSE) {
  type <- check_dep_word(type)
  pdb <- utils::available.packages(repos = utils::findCRANmirror("web"))
  if (length(type) == 1L) {
    v0 <- try(tools::package_dependencies(name, db = pdb, which = type, reverse = reverse), silent = TRUE)
    if (inherits(v0, "try-error")) {
      stop("get_dep() uses tools::package_dependencies() which fails. Check Internet connection.")
    } else {
      if (length(v0[[1L]]) == 0L) {
        v0[[1L]] <- as.character(NA)
      }
      df0 <- data.frame(from = name, to = v0[[1L]], type = type, reverse = reverse, stringsAsFactors = FALSE)
    }
  } else {
    l0 <- list()
    for (i in seq_along(type)) {
      typei <- type[i]
      l0[[i]] <- get_dep(name, typei, reverse)
    }
    df0 <- do.call(rbind, l0)
  }
  df0 <- df0[!is.na(df0$to),]
  df0$type <- conditional_change(tolower(df0$type), "linkingto", "linking to")
  unique(df0) # there might be duplicates
}

#' Reshape the data frame of dependencies
#' 
#' @param x A character vector of dependencies, each element of which corresponds to an individual package
#' @param names A character vector of package names of the same length as x 
#' @return A data frame of dependencies
#' @keywords internal
reshape_dep <- function(x, names) {
  y <- lapply(x, get_dep_vec)
  df0 <- data.frame(
    from = rep(names, sapply(y, length)),
    to = unlist(y),
    stringsAsFactors = FALSE
  )
}

#' Dependencies of all CRAN packages
#'
#' \code{get_dep_all_packages} returns the data frame of dependencies of all packages currently available on CRAN.
#'
#' @importFrom tools CRAN_package_db
#' @importFrom dplyr bind_rows
#' @return A data frame of dependencies of all CRAN packages
#' @examples
#' \dontrun{
#' df.cran <- get_dep_all_packages()
#' }
#' @seealso \code{\link{get_dep}} for multiple types of dependencies, and \code{\link{get_graph_all_packages}} for obtaining directly a network of dependencies as an igraph object
#' @export
get_dep_all_packages <- function() {
  db0 <- try(tools::CRAN_package_db(), silent = TRUE)
  if (inherits(db0, "try-error")) {
    stop("get_dep_all_packages() uses tools::CRAN_package_db() which fails. Check Internet connection.")
  } else {
    df0 <- as.data.frame(db0, stringsAsFactors = FALSE)
    pkgnames <- df0$Package
    df1 <-
      dplyr::bind_rows(
               `FALSE` = dplyr::bind_rows(
                                  imports = reshape_dep(df0$Imports, pkgnames),
                                  depends = reshape_dep(df0$Depends, pkgnames),
                                  `linking to` = reshape_dep(df0$LinkingTo, pkgnames),
                                  suggests = reshape_dep(df0$Suggests, pkgnames),
                                  enhances = reshape_dep(df0$Enhances, pkgnames),
                                  .id = "type"
                                ),
               `TRUE` = dplyr::bind_rows(
                                 imports = reshape_dep(df0$`Reverse imports`, pkgnames),
                                 depends = reshape_dep(df0$`Reverse depends`, pkgnames),
                                 `linking to` = reshape_dep(df0$`Reverse linking to`, pkgnames),
                                 suggests = reshape_dep(df0$`Reverse suggests`, pkgnames),
                                 enhances = reshape_dep(df0$`Reverse enhances`, pkgnames),
                                 .id = "type"
                               ),
               .id = "reverse"
             )
    df2 <- data.frame(
      from = df1$from,
      to = df1$to,
      type = df1$type,
      reverse = as.logical(df1$reverse),
      stringsAsFactors = FALSE
    )
    df2 <- df2[!is.na(df2$to),]
    df2 <- df2[df2$to != "",] # there are rows with "" - may need to fix from source get_dep_vec()
    return(unique(df2))
  }
}
