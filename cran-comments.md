This is the first submission of version 0.3.10. The total check time is expected to be below 10min. The changes are partially in response to the email dated 2024-08-10 regarding the CRAN policy: 'Packages which use Internet resources should fail gracefully with an informative message if the resource is not available or has changed (and not give a check warning nor error).' Functions now which use Internet resources will return NULL while printing a message instead of stopping and returning an error message.

There are no ERRORs or WARNINGs, and possible NOTEs are explained below.

## checking installed package size ... NOTE
    installed size is 21.1Mb
    sub-directories of 1Mb or more:
      data   1.2Mb
      doc    2.2Mb
      libs  17.5Mb
