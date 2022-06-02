This is a first submission of version 0.3.1. This is mainly in response to Brian Ripley's email (dated 2022-05-22 and cc'd to cran@r-project.org). The functions html_text_vec(), get_dep() and get_dep_all_packages() have now included error messages that compile with the policy "Packages which use Internet resources should fail gracefully with an informative message if the resource is not available or has changed (and not give a check warning nor error)."

There are no ERRORs or WARNINGs, and possible NOTEs are explained below.

## checking installed package size ... NOTE
    installed size is  7.2Mb
    sub-directories of 1Mb or more:
      data   1.2Mb
      doc    2.3Mb
      libs   3.6Mb
