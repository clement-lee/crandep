This is a resubmission of version 0.1.0, which also fixes the ERROR in last released version's CRAN status - see `checking examples ... NOTE / ERROR` below.





## checking CRAN incoming feasibility ... NOTE

- The invalid file URIs in NEWS.md are now removed.

- Regarding the NOTE below (between the two lines of equal signs), the URL concerned is for the list of all packages instead of 1 individual package, and therefore the canonical form is not possible. Also, in the first file where this URL is found, no changes have been made since the previous version.

=========================================
Found the following (possibly) invalid URLs:
  URL: https://cran.r-project.org/web/packages/available_packages_by_name.html
    From: man/cran_dependencies.Rd
          inst/doc/introduction.html
    Status: 200
    Message: OK
    CRAN URL not in canonical form
  The canonical URL of the CRAN page for a package is 
    https://CRAN.R-project.org/package=pkgname
==========================================





## checking examples ... NOTE / ERROR

- The functions with example times > 10s are `get_dep()` or `get_dep_all()`, `get_graph_all_packages()` and `get_dep_all_packages()`, all of which import `tools::CRAN_package_db()`, and therefore take time and require internet connection. The fixes are as below:

- For `get_dep()` and `get_dep_all()`, the lines of examples ` which resulted in errors (due to failure to establish server connection) previously are now removed, thus removing the aforementioned error.

- For `get_graph_all_packages()` and `get_dep_all_packages()`, the examples are now wrapped by \dontrun{}.